LATEXDIR=latex/
GCC=g++
SAMTOOLS_DIR=samtools-0.1.18/

.PHONY: all clean test

all: sjcount-deprecated sjcount

clean ::
	rm -f -r sjcount-deprecated sjcount progressbar.o 

${SAMTOOLS_DIR}libbam.a:
	wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
	tar -xf samtools-0.1.18.tar.bz2
	rm -f samtools-0.1.18.tar.bz2
	make -C samtools-0.1.18 all
	# If SAMTOOLS is already installed, you might want to update SAMTOOLS_DIR path without installing a fresh copy

progressbar.o :	progressbar.c progressbar.h
	$(GCC) -c progressbar.c 

sjcount-deprecated : sjcount-deprecated.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount-deprecated.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount-deprecated

sjcount : sjcount.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount

${LATEXDIR}sjcount.pdf : ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex

######################################################################################################################

TESTDIR=test/

test :: ${TESTDIR}control.ssj ${TESTDIR}control.ssc
	#===> tests passed successfully <===#

PARAMS=-lim 100000 -nbins 50 -read1 0 -read2 0

${TESTDIR}test.bam : 
	wget http://gasv.googlecode.com/files/Example.bam -O ${TESTDIR}test.bam

${TESTDIR}test.ssj ${TESTDIR}test.ssc : ${TESTDIR}test.bam sjcount
	sjcount -bam ${TESTDIR}test.bam -ssj ${TESTDIR}test.ssj -ssc ${TESTDIR}test.ssc ${PARAMS}

${TESTDIR}control.ssj : ${TESTDIR}test.bam ${TESTDIR}test.ssj ${TESTDIR}sam2sj.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sj.pl ${PARAMS} |sort -k1,1 -k2,3n > ${TESTDIR}control.ssj
	cmp ${TESTDIR}test.ssj ${TESTDIR}control.ssj

${TESTDIR}control.ssc : ${TESTDIR}test.bam ${TESTDIR}test.ssc ${TESTDIR}control.ssj ${TESTDIR}sam2sb.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sb.pl -ssj ${TESTDIR}control.ssj ${PARAMS} |sort -k1,1 -k2,3n > ${TESTDIR}control.ssc
	cmp ${TESTDIR}test.ssc ${TESTDIR}control.ssc

clean ::
	rm -f ${TESTDIR}test.bam ${TESTDIR}test.ssj ${TESTDIR}test.ssc ${TESTDIR}control.ssj ${TESTDIR}control.ssc



