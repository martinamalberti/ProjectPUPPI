CC        = g++
CCFLAGS   = -Wall
SOURCES   =
ROOTFLAGS = `root-config --cflags --glibs` -lTMVA


all:
	$(CC) $(CCFLAGS) -o JetTreeAnalysis $(ROOTFLAGS) JetTreeAnalysis.cc JetTree.h
	$(CC) $(CCFLAGS) -o TestJetTreeAnalysis $(ROOTFLAGS) TestJetTreeAnalysis.cc JetTreeAnalyzer.cc

JetTreeAnalysis:
	$(CC) $(CCFLAGS) -o JetTreeAnalysis $(ROOTFLAGS) JetTreeAnalysis.cc JetTree.cc

TestJetTreeAnalysis:
	$(CC) $(CCFLAGS) -o TestJetTreeAnalysis $(ROOTFLAGS) TestJetTreeAnalysis.cc JetTreeAnalyzer.cc

clean:
	rm JetTreeAnalysis TestJetTreeAnalysis

