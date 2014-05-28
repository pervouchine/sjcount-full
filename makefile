LATEXDIR=latex/
GCC=g++
SAMTOOLS_DIR=samtools-0.1.18/

.PHONY: all clean test

all: sjcount-deprecated sjcount sjcount3

clean ::
	rm -f -r sjcount-deprecated sjcount progressbar.o sjcount3

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

sjcount3 : sjcount3.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount3.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount3


${LATEXDIR}sjcount.pdf : ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex

######################################################################################################################

TESTDIR=test/

test :: ${TESTDIR}control.ssj ${TESTDIR}control.ssc
	sort ${TESTDIR}test.ssj | cmp ${TESTDIR}control.ssj
	sort ${TESTDIR}test.ssc | cmp ${TESTDIR}control.ssc
	#===> tests v2 passed successfully <===#

PARAMS=-lim 1000000 -nbins 50 -read1 0 -read2 0

${TESTDIR}test.bam : 
	wget genome.crg.es/~dmitri/export/sjcount/test.bam -O ${TESTDIR}test.bam

${TESTDIR}test.ssj ${TESTDIR}test.ssc : ${TESTDIR}test.bam sjcount
	sjcount -bam ${TESTDIR}test.bam -ssj ${TESTDIR}test.ssj -ssc ${TESTDIR}test.ssc ${PARAMS}

${TESTDIR}control.ssj : ${TESTDIR}test.bam ${TESTDIR}test.ssj ${TESTDIR}sam2sj.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sj.pl ${PARAMS} | sort > ${TESTDIR}control.ssj

${TESTDIR}control.ssc : ${TESTDIR}test.bam ${TESTDIR}test.ssc ${TESTDIR}control.ssj ${TESTDIR}sam2sb.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sb.pl -ssj ${TESTDIR}control.ssj ${PARAMS} | sort > ${TESTDIR}control.ssc

clean ::
	rm -f ${TESTDIR}test.bam ${TESTDIR}test.ssj ${TESTDIR}test.ssc ${TESTDIR}control.ssj ${TESTDIR}control.ssc

######################################################################################################################

${TESTDIR}test.ssj3 ${TESTDIR}test.ssc3 : ${TESTDIR}test.bam sjcount3
	sjcount3 -bam ${TESTDIR}test.bam ${PARAMS} -ssj ${TESTDIR}test.ssj3 -ssc ${TESTDIR}test.ssc3

${TESTDIR}control.ssj3 : ${TESTDIR}test.bam  ${TESTDIR}sam2sj3.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  | perl ${TESTDIR}sam2sj3.pl ${PARAMS} | sort > ${TESTDIR}control.ssj3

${TESTDIR}control.ssc3 : ${TESTDIR}test.bam ${TESTDIR}sam2sb3.pl
	${SAMTOOLS_DIR}samtools view ${TESTDIR}test.bam  |  perl ${TESTDIR}sam2sb3.pl -ssj ${TESTDIR}test.ssj3 ${PARAMS} | sort > ${TESTDIR}control.ssc3


test3 :: ${TESTDIR}test.ssj3 ${TESTDIR}control.ssj3 ${TESTDIR}test.ssc3 ${TESTDIR}control.ssc3
	sort ${TESTDIR}test.ssj3 | cmp ${TESTDIR}control.ssj3
	sort ${TESTDIR}test.ssc3 | cmp ${TESTDIR}control.ssc3
	#===> tests v3 passed successfully <===#


clean ::
	rm -f ${TESTDIR}test.ssj3 ${TESTDIR}control.ssj3


