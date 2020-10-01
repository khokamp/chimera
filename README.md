# chimera
Perl scripts for analysis of unmapped RNA-seq reads with the aim to find chimeras, consisting of part sRNA part mRNA.

## part 1:
generate_start_and_end_from_shortreads.pl
here we are preparing input files by splitting reads that could not be mapped 
against the reference genome into a start and an end library:

## part 2:
chimera_start_end_mapped_distances.pl
here we are preparing input files by splitting reads that could not be mapped 
against the reference genome into a start and an end library:

## part 3:
chimera_annotation.pl
here we focus on a subset of annotated features (sRNAs) and create output files for each of them
the output contains all reads where either the start or end maps to the according sRNA
and the feature that is either overlapped by the other start/end or lies nearby its mapping position
