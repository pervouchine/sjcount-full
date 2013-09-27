Copyright 2013 Dmitri Pervouchine (dp@crg.eu), Lab Roderic Guigo
Bioinformatics and Genomics Group @ Centre for Genomic Regulation
Parc de Recerca Biomedica: Dr. Aiguader, 88, 08003 Barcelona

This file is a part of the 'sjcount' package. 
'sjcount' package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

'sjcount' package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License 
along with 'sjcount' package.  If not, see <http://www.gnu.org/licenses/>.

============================================================================

DESCRIPTION

sjcount is a utility for fast SJ quantification. It is the annotation-agnostic 
vesrion of bam2ssj (see https://github.com/pervouchine/bam2ssj)

============================================================================

INSTALLATION

To install, type 'make all'

Prerequisites:
	You need to install samtools

	Get it by svn:
	svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools
	enter the directory and type 'make all'

	Samtools package Version: 0.1.18-dev (r982:313) is compartible (and likely so are oder versions)

NOTE that

1.	You need to to update the SAMDIR varibale in the makefile

2.	Some users get error messages even when compiling sjcount with a correct SAMDIR path, something like

	/centos6/samtools-9.3.2013/libbam.a(bgzf.o): In function `mt_destroy`:
	/centos6/samtools-9.3.2013/bgzf.c:458: undefined reference to `pthread_join`
	/centos6/samtools-9.3.2013/libbam.a(bgzf.o): In function `bgzf_mt`:
	/centos6/samtools-9.3.2013/bgzf.c:445: undefined reference to `pthread_create`

This error has to do with big zip libraries, not with samtools

============================================================================

USAGE

./sjcount -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file] [-maxlen max_intron_length] [-minlen min_intron_length] [-margin length] [-read1 0|1] [-read2 0|1] [-nbins number_of_bins] [-binsize bin_size] [-lim number_of_lines] [-quiet]

Input:   a (sorted) BAM file

Options:

	-maxlen upper limit on intron length, 0 = no limit (default=0)
	-minlen lower limit on intron length, 0 = no limit (default=0)
	-margin length, minimum number of flanking nucleotides in the read in order to support SJ or EB, (default=0)
	-read1 0/1, reverse complement read1 no/yes (default=0)
	-read2 0/1, reverse complement read2 no/yes (default=0)
	-binsize size of the overhang bin, (default=+INFTY)
	-nbins number of overhang bins, (default=1)
	-lim nreads stop after nreads, (default=no limit)
	-quiet, suppress verbose output

Output:

	-ssj: Splice junction counts, tab-delimited  (default = stdout)
	Columns are: chr, begin, end, offset, count (+), count (-)
	-ssc: Splice boundary counts, tab-delimited  (default = none)
	Columns are: chr, position, count (+), count (-)


============================================================================

DETAILS

	Each 'N' in the CIGAR string defines a splice junction. The lengths of the continuous portions of the upstream 
	and downstream alignment (M, match only!) have to be at least <margin> nucleotides. The break has to be between
	<min_intron_length> and <max_intron_length>. Next, each 'N' in the CIGAR has a continuous portion upstream from 
	the start of the alignment or previous N (i.e., only M/I/D). It is called 'offset'. Offset is used to bin reads 
	by the formula {bin = (int) offset/binsize}. The max value of 'bin' is enforced to be nbins - 1. By default, 
	nbins = 1, i.e., all SJ counts fall into one bin.
 








