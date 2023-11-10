use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, May 2019
#
# part 2 of a multi-part analysis that tries to detect chimeric RNA sequences
# here we are going through the mapping results of start and end sequences
# and extract pairs together with distances, sequences and qualities
# an optional GFF file can be specified to only store results that overlap with
# annotated features

my $location_file = '';
my $sample = 'SampleName';

&GetOptions(
    'sample=s' => \$sample,
    'locations=s' => \$location_file,
    );

my $start_file = shift;
my $end_file = $start_file;
$end_file =~ s/\.start/.end/;

# if a location file has been specified 
# read in the areas covered by the features of interest:
my %loc = ();
if ($location_file) {
    open (IN, $location_file)
        or die "Can't find location file $location_file: $!";

#NZ_CP008706.1   RefSeq  gene    92960   93224   .       +       .       ID=sRNA1
#NZ_CP008706.1   RefSeq  gene    3810843 3810999 .       +       .       ID=sRNA2

    while (<IN>) {
        chomp;
        my @h = split /\t/, $_;
        my $start = $h[3];
        my $end = $h[4];
        foreach ($start..$end) {
            $loc{$h[0]}{$_}++;
        }
    }
    close IN;
}


# read in all mappings from the start of the originally unmapped reads
# and store position and strand of the location
my %start = ();
open (IN, "samtools view $start_file |")
    or die;

#NS500339:321:HLWL7BGX9:2:21312:15678:15383      16      gi|798817478|ref|NZ_CP008706.1| 7  
while (<IN>) {
    chomp;
    my ($id, $flag, $chr, $pos, $mapq, @rest) = split /\t/, $_;
    if (defined $start{$id}) {
        # keep only the first position in case of multi-mappers
        warn "$id mapped more than once\n";
        next;
    }
    my $strand = '+';
    if ($flag & 16) {
        $strand = '-';
    }
    $start{$id} = "$chr\t$pos\t$strand\t$mapq";
}

close IN;

# now we're reading in the mappings of the ends of the originally unmapped reads
# and keep all those where the according start part has been mapped as well:
open (IN, "samtools view $end_file |")
    or die;

my %out = ();
while (<IN>) {
    chomp;
    my ($id, $flag, $chr, $pos, $mapq, @rest) = split /\t/, $_;
    my $strand = '+';
    if ($flag & 16) {
        $strand = '-';
    }
    if (defined $start{$id}) {
        my ($chr2, $pos2, $strand2, $mapq2) = split /\t/, $start{$id};
        my $dist = $pos - $pos2;
	if ($chr2 ne $chr) {
            $dist = 'NA';
        }
        $out{$id} = "$sample\t$chr2\t$pos2\t$strand2\t$chr\t$pos\t$strand\t$dist\t$mapq2\t$mapq";
    }
}
close IN;

my $num = scalar keys %out;
print STDERR "Found $num mapped pairs for $sample\n";

print "ReadID\tSample\tstart chr\tstart location\tstart strand\tFeature\tend location\tend chr\tend strand\tdistance\tstart MAPQ\tend MAPQ\tstart sequence\tend sequence\n";

# retrieve the 20bp sequences from the start/end fastq files and add them to the collected output:
my $seq_start = $start_file;
$seq_start =~ s/\.bam//;

unless (open (IN, $seq_start)) {
    warn "Can't read sequences from $seq_start: $!";
} else {

    while (<IN>) {
        
        # @NS500339:321:HLWL7BGX9:4:13403:3564:7725 1:N:0:ATTACTCG+AGAGGATA ef:0;if:0 36:19 0:54426:1
        s/^\@//;
        s/ .+//;
        my $seq = <IN>;
        my $head2 = <IN>;
        my $qual = <IN>;
        chomp;
        if (defined $out{$_}) {
            chomp $seq;
            $out{$_} .= "\t$seq";
        }
    }
    close IN;
}

my $seq_end = $end_file;
$seq_end =~ s/\.bam//;
unless (open (IN, $seq_end)) {
    warn "Can't read sequences from $seq_end: $!";
} else {

    while (<IN>) {
        s/^\@//;
        s/ .+//;
        my $seq = <IN>;
        my $head2 = <IN>;
        my $qual = <IN>;
        chomp;
        if (defined $out{$_}) {
            chomp $seq;
            $out{$_} .= "\t$seq";
        }
    }
    close IN;
}

# print out all of the information of mapped start/end pairs
# (or only those that overlap with annotated features):
my $pass = 0;
foreach my $id (keys %out) {
    if (%loc) {
        my @h = split /\t/, $out{$id};
        my $found = 0;
        my $chr1 = $h[1];       
        my $chr2 = $h[4];       
        if (defined $loc{$chr1}{$h[2]}      
            or
            defined $loc{$chr2}{$h[5]}) {
            $found = 1;
        }
        next unless ($found);
        $pass++;
    }
    print "$id\t$out{$id}\n";
}

if (%loc) {
    print STDERR "$sample: $pass pairs overlap with locations in $location_file\n";
}
