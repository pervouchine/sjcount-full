SAMDIR=~/software/samtools/
LATEXDIR=latex/
GCC=g++

.PHONY: all

all: sjcount sjcount2 ${LATEXDIR}sjcount.pdf

$(SAMDIR)libbam.a:
	# You need to install samtools
	# Get it by svn:
	# svn co https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools
	# enter the dir and type 'make all'
	# don't forget to update the SAMDIR varibale in this makefile
	exit 1	

progressbar.o:	progressbar.c progressbar.h
	$(GCC) -c progressbar.c 

sjcount : sjcount.c progressbar.o $(SAMDIR)libbam.a
	$(GCC) -I $(SAMDIR) sjcount.c progressbar.o $(SAMDIR)libbam.a -lz -o sjcount

sjcount2 : sjcount2.c progressbar.o $(SAMDIR)libbam.a
	$(GCC) -I $(SAMDIR) sjcount2.c progressbar.o $(SAMDIR)libbam.a -lz -o sjcount2

${LATEXDIR}sjcount.pdf : ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex
	pdflatex -output-directory=${LATEXDIR} ${LATEXDIR}sjcount.tex


clean:
	rm -f -r progressbar.o sjcount sjcount2

