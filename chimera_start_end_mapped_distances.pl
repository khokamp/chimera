use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, May 2019
#
# part 2 of a multi-part analysis that tries to detect chimeric RNA sequences
# here we are preparing input files by splitting reads that could not be mapped 
# against the reference genome into a start and an end library:

my $location_file = '';

&GetOptions(
    'locations=s' => \$location_file,
    );

my $start_file = shift;
my $end_file = $start_file;
$end_file =~ s/\.start/.end/;

# read in all mappings from the start of the originally unmapped reads
# and store position and strand of the location
my %start = ();
open (IN, "samtools view $start_file |")
    or die;

#NS500339:321:HLWL7BGX9:2:21312:15678:15383      16      gi|798817478|ref|NZ_CP008706.1| 7  
while (<IN>) {
    chomp;
    my ($id, $flag, $chr, $pos, @rest) = split /\t/, $_;
    if (defined $start{$id}) {
	warn "$id mapped more than once\n";
	next;
    }
    my $strand = '+';
    if ($flag == 16) {
	$strand = '-';
    }
    $start{$id} = "$pos\t$strand";
}

close IN;

# if a location file has been specified 
# read in the areas covered by the features of interest:
my %loc = ();
if ($location_file) {
    open (IN, $location_file)
	or die;

#NZ_CP008706.1   RefSeq  gene    92960   93224   .       +       .       ID=sRNA1
#NZ_CP008706.1   RefSeq  gene    3810843 3810999 .       +       .       ID=sRNA2

    while (<IN>) {
	chomp;
	my @h = split /\t/, $_;
	my $start = $h[3];
	my $end = $h[4];
	foreach ($start..$end) {
	    $loc{$_}++;
	}
    }
    close IN;
}

# now we're reading in the mappings of the ends of the originally unmapped reads
# and keep all those where the according start part has been mapped as well:
open (IN, "samtools view $end_file |")
    or die;

my %out = ();
while (<IN>) {
    chomp;
    my ($id, $flag, $chr, $pos, @rest) = split /\t/, $_;
    my $strand = '+';
    if ($flag == 16) {
        $strand = '-';
    }
    if (defined $start{$id}) {
	my ($pos2, $strand2) = split /\t/, $start{$id};
	my $dist = $pos - $pos2;
	$out{$id} = "$pos2\t$strand2\t$pos\t$strand\t$dist";
    }
}
close IN;

my $num = scalar keys %out;
print STDERR "Found $num mapped pairs\n";

print "ReadID\tstart location\tstart strand\tend location\tend strand\tdistance\tstart sequence\tend sequence\n";

# retrieve the 20bp sequences from the start/end fastq files and add them to the collected output:
my $seq_start = $start_file;
$seq_start =~ s/\.bam//;
open (IN, $seq_start)
    or die;
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

my $seq_end = $end_file;
$seq_end =~ s/\.bam//;
open (IN, $seq_end)
    or die;
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

# print out all of the information of mapped start/end pairs
# (or only those that overlap with annotated features):
my $pass = 0;
foreach my $id (keys %out) {
    if (%loc) {
	my @h = split /\t/, $out{$id};
	my $found = 0;
	if (defined $loc{$h[0]}	    
	    or
	    defined $loc{$h[2]}) {
	    $found = 1;
	}
	next unless ($found);
	$pass++;
    }
    print "$id\t$out{$id}\n";
}

if (%loc) {
    print STDERR "$pass pairs overlap with locations in $location_file\n";
}

