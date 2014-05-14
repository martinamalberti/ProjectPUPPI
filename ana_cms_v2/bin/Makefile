CC        = g++
CCFLAGS   = -Wall
SOURCES   =
ROOTFLAGS = `root-config --cflags --glibs` -lTMVA


all:
	$(CC) $(CCFLAGS) -o JetTreeAnalysis $(ROOTFLAGS) JetTreeAnalysis.cc JetTree.h

JetTreeAnalysis:
	$(CC) $(CCFLAGS) -o JetTreeAnalysis $(ROOTFLAGS) JetTreeAnalysis.cc JetTree.cc

clean:
	rm JetTreeAnalysis

