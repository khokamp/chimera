use warnings;
use strict;

# Karsten Hokamp, Trinity College Dublin, May 2019
#
# part 1 of a multi-part analysis that tries to detect chimeric RNA sequences
# here we are preparing input files by splitting reads that could not be mapped 
# against the reference genome into a start and an end library:

# if dealing with Fastq instead of Fasta file, set this to '1':
my $fastq_file = 0;

foreach my $file (@ARGV) {

    # all files in bzip format, uncompress on the fly, adjust if needed
    if ($file =~ /\.bz2$/) {
	open (IN, "bunzip2 -c $file |")
	    or die;
    } else {
	die;
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

	my $start = substr $seq, 0, 20;
	my $end = substr $seq, -20;
	my $qual_start = substr $qual, 0, 20;
	my $qual_end = substr $qual, -20;
	print OUT1 "$head$start\n$head2\n$qual_start\n";	
	print OUT2 "$head$end\n$head2\n$qual_end\n";
    }
    close IN;
    close OUT1;
    close OUT2;
}
