LIBS_GSL = -lgsl -lgslcblas

LIBS = $(LIBS_GSL)
INCL = -I./

CXX = g++
CXXFLAGS += -O2 -Wall -c
DEBUG = -g

PDFLATEX = pdflatex

default: all

.SUFFIXES: .o .cpp

all: test

clean:
	rm -f *~
	rm -f *.{o,x}
	rm -f *.{aux,log,out,dat,pdf}

test: \
	test.o
	$(CXX) \
	test.o \
	-o test.x $(LIBS)

look: look.pdf


%.o: %.cpp
	$(CXX) $(DEBUG) $(CXXFLAGS) $(INCL) $<

%.pdf: %.tex
	$(PDFLATEX) $<
