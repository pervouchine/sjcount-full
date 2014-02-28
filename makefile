LATEXDIR=latex/
GCC=g++
SAMTOOLS_DIR=samtools-0.1.18/

.PHONY: all clean

all: sjcount sjcount2

clean ::
	rm -f -r ${SAMTOOLS_DIR} sjcount sjcount2 progressbar.o 

${SAMTOOLS_DIR}libbam.a:
	wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download
	tar -xf samtools-0.1.18.tar.bz2
	rm -f samtools-0.1.18.tar.bz2
	make -C samtools-0.1.18 all
	# If SAMTOOLS is already installed, you might want to update SAMTOOLS_DIR path without installing a fresh copy

progressbar.o :	progressbar.c progressbar.h
	$(GCC) -c progressbar.c 

sjcount : sjcount.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount

sjcount2 : sjcount2.c progressbar.o ${SAMTOOLS_DIR}libbam.a
	$(GCC) -I ${SAMTOOLS_DIR} sjcount2.c progressbar.o ${SAMTOOLS_DIR}libbam.a -lz -o sjcount2

${LATEXDIR}sjcount.pdf : ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex


