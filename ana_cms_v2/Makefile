# ----------------------------------------
# inputs 

#FASTJETLOCATION=/afs/cern.ch/sw/lcg/external/fastjet/3.0.6/x86_64-slc5-gcc46-opt/bin/
FASTJETLOCATION=/afs/cern.ch/user/m/malberti/fastjet-3.0.6/
FASTJETINC :=`$(FASTJETLOCATION)/fastjet-config --cxxflags --plugins`
FASTJETLIB :=`$(FASTJETLOCATION)/fastjet-config --libs --plugins`
#BACONLIB   :=-L/Users/Phil/Physics/JWJ/git_v2/ProjectPUPPI/ana_cms/Bacon/BaconAna/DataFormatsOffline/ -lBacon  -Wl,-rpath,/Users/Phil/Physics/JWJ/git_v2/ProjectPUPPI/ana_cms/Bacon/BaconAna/DataFormatsOffline/ -lm -ldl
BACONLIB   :=-LBacon/BaconAna/DataFormatsOffline/ -lBacon  -Wl,-rpath,Bacon/BaconAna/DataFormatsOffline/ -lm -ldl
ROOTINC := `root-config --glibs --cflags` -lTMVA

CC = g++

# ----------------------------------------


all: puppiSimpleAnalyzer

## FOR MAVERICKS (CLANG)

puppiSimpleAnalyzer: puppiSimpleAnalyzer
	$(CC) -o puppiSimpleAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) $(BACONLIB) -lvectorDict -lJetCleanser puppiSimpleAnalyzer.cpp puppiCleanContainer.cc NoTrees.cc -I.

## FOR MAC OSX <= 10.8
## puppiSimpleAnalyzer: puppiSimpleAnalyzer
##	$(CC) -o puppiSimpleAnalyzer $(ROOTINC) $(FASTJETINC) $(FASTJETLIB) $(BACONLIB) -lvectorDict -lJetCleanser puppiSimpleAnalyzer.cpp puppiCleanContainer.cc puppiCleanContainer.hh NoTrees.cc NoTrees.hh -I.

clean:
	rm puppiSimpleAnalyzer
