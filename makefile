LATEXDIR=latex/
GCC=g++
SAMTOOLS_DIR=samtools-0.1.18/

.PHONY: all clean test

all:  sjcount_v3 sjcount_v1 sjcount_v2

clean ::
	rm -f -r sjcount_v3 progressbar.o sjcount_v1 sjcount_v2

${SAMTOOLS_DIR}libbam.a:
	wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2
	tar -xf samtools-0.1.18.tar.bz2
	rm -f samtools-0.1.18.tar.bz2
	make -C samtools-0.1.18 all
	# If SAMTOOLS is already installed, you might want to update SAMTOOLS_DIR path without installing a fresh copy

progressbar.o : progressbar.c progressbar.h
	$(GCC) -c progressbar.c 

sjcount_v1 : sjcount_v1.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount_v1.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount_v1

sjcount_v2 : sjcount_v2.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount_v2.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount_v2

sjcount_v3 : sjcount_v3.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount_v3.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount_v3


${LATEXDIR}sjcount_v3.pdf : ${LATEXDIR}sjcount_v3.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount_v3.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount_v3.tex

######################################################################################################################

TESTDIR=test/

PARAMS=-lim 1000000 -nbins 50 -read1 0 -read2 0

${TESTDIR}test.bam : 
	wget genome.crg.es/~dmitri/export/sjcount/test.bam -O ${TESTDIR}test.bam

${TESTDIR}test.ssj ${TESTDIR}test.ssc : ${TESTDIR}test.bam sjcount
	sjcount_v3 -bam ${TESTDIR}test.bam ${PARAMS} -ssj ${TESTDIR}test.ssj -ssc ${TESTDIR}test.ssc

${TESTDIR}control.ssj : ${TESTDIR}test.bam  ${TESTDIR}sam2sj3.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sj3.pl ${PARAMS} | sort > ${TESTDIR}control.ssj

${TESTDIR}control.ssc : ${TESTDIR}test.bam ${TESTDIR}sam2sb3.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  |  perl ${TESTDIR}sam2sb3.pl -ssj ${TESTDIR}test.ssj ${PARAMS} | sort > ${TESTDIR}control.ssc


test :: ${TESTDIR}test.ssj ${TESTDIR}control.ssj ${TESTDIR}test.ssc ${TESTDIR}control.ssc
	sort ${TESTDIR}test.ssj | cmp ${TESTDIR}control.ssj
	sort ${TESTDIR}test.ssc | cmp ${TESTDIR}control.ssc
	#===> tests v3 passed successfully <===#

clean ::
	rm -f ${TESTDIR}test.ssj ${TESTDIR}control.ssj
