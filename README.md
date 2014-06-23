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
	You need samtools, but it will be installed automatically by the makefile
	Samtools package Version: 0.1.18

NOTE that

1.	If uyou have samtools already installed, you can simply update the SAMTOOLS_DIR varibale in the makefile

2.	Some users get error messages even when compiling sjcount with a correct SAMDIR path, something like

	/centos6/samtools-9.3.2013/libbam.a(bgzf.o): In function `mt_destroy`:
	/centos6/samtools-9.3.2013/bgzf.c:458: undefined reference to `pthread_join`
	/centos6/samtools-9.3.2013/libbam.a(bgzf.o): In function `bgzf_mt`:
	/centos6/samtools-9.3.2013/bgzf.c:445: undefined reference to `pthread_create`

This error has to do with big zip libraries, not with samtools. 

============================================================================

USAGE

 ./sjcount_v3 -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file] [-maxlen max_intron_length] [-minlen min_intron_length] [-read1 0|1] [-read2 0|1] [-nbins number_of_bins] [-binsize bin_size] [-lim number_of_lines] [-quiet]

Input:   a (sorted) BAM file

Options:

	-read1 0/1, reverse complement read1 no/yes (default=1)
	-read2 0/1, reverse complement read2 no/yes (default=0)
	-nbins number of bins for offsets, (default=1)
	-lim nreads stop after nreads, (default=no limit)
	-maxnh, the max value of the NH tag, (default=none)
	-continuous, no mismatches when overlapping splice boundaries
	-unstranded, force strand=0
	-quiet, suppress verbose output

Output:	-ssj: Splice Junction counts, tab-delimited  (default=stdout)
	Columns are: chr, begin, end, strand, offset, count
	-ssc: Splice boundary counts, tab-delimited  (default=none)
	Columns are: chr, position, position, strand, offset, count

============================================================================

DETAILS
	See documentation in latex/sjcount_v3.pdf
