all: RDAnalysis

RDAnalysis: RDAnalysis.cc MyParser.h
	g++ `root-config --libs` -lMinuit `root-config --cflags` --std=c++11 -g \
		RDAnalysis.cc -o RDAnalysis
