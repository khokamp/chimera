use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, May 2019                                                                                                                                                                                                 
#                                                                                                                                                                                                                                                  
# part 1 of a multi-part analysis that tries to detect chimeric RNA sequences                                                                                                                                                                      
# here we are preparing input files by splitting reads that could not be mapped                                                                                                                                                                    
# against the reference genome into a start and an end library:                                                                                                                                                                                    

# length of sequence string to use from each end                                                                                                                                                                                                   
my $len = 20;

# if dealing with Fastq instead of Fasta file, set this to '1':                                                                                                                                                                                    
# or use the -fastq flag when running the script                                                                                                                                                                                                   
my $fastq_file = '';

&GetOptions(
    'length=i' => \$len,
    'fastq' => \$fastq_file,
    );

foreach my $file (@ARGV) {

    # all files in bzip format, uncompress on the fly, adjust if needed                                                                                                                                                                            
    if ($file =~ /\.bz2$/) {
        open (IN, "bunzip2 -c $file |")
            or die;
    } elsif ($file =~ /\.gz$/) {
        open (IN, "gunzip -c $file |")
            or die;
    } else {
        open (IN, $file)
            or die;
    }

    # name output files after input file                                                                                                                                                                                                           
    # attache suffixes '.start' and '.end':                                                                                                                                                                                                        
    my $out1 = $file;
    $out1 =~ s/(.+)\..+/$1/;
    my $out2 = $out1;
    $out1 .= '.start';
    $out2 .= '.end';
    open (OUT1, ">$out1")
        or die;
    open (OUT2, ">$out2")
        or die;

    while (<IN>) {
        my $head = $_;
        $head =~ s/^>/\@/;
        my $seq = <IN>;
        chomp $seq;
        my $head2 = '+';

        # generate dummy quality values                                                                                                                                                                                                            
        my $qual = 'I' x length($seq);

        if ($fastq_file) {
            $head2 = <IN>;
            chomp $head2;
            $qual = <IN>;
            chomp $qual;
        }

        my $start = substr $seq, 0, $len;
        my $end = substr $seq, -$len;
        my $qual_start = substr $qual, 0, $len;
        my $qual_end = substr $qual, -$len;
        print OUT1 "$head$start\n$head2\n$qual_start\n";
        print OUT2 "$head$end\n$head2\n$qual_end\n";
    }
    close IN;
    close OUT1;
    close OUT2;
}
