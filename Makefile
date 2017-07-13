CC 		= g++
CFLAGS		= -c -Wall -Wextra -pedantic -O2 -std=gnu++14
LFLAGS		= -Wall -fopenmp
BOOST		= -I $(HOME)/Boost/include/
SEQAN		= -I $(HOME)/seqan-library-2.0.0/include
VIENNA		= -I $(HOME)/Vienna/include
VIENNALIB	= -L $(HOME)/Vienna/lib
BOOSTLIB	= -L $(HOME)/Boost/lib
CXX		= /usr/bin/gcc
LIBS		= -lboost_filesystem -lboost_system -lboost_iostreams
LIBSPOLYFUD	= -lboost_filesystem -lboost_system -lboost_iostreams -lRNA
LIBSTEST	= -lboost_unit_test_framework -lboost_filesystem -lboost_system -lboost_iostreams
TPATH		= bin/
SRCPATH		= src/

OBJS		= $(TPATH)hgvsParser.o \
	lib/polarUtility.o \
	$(TPATH)seqStruct.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)polyFud.o \

.PHONY : all
all : $(TPATH)uTests $(TPATH)perf_polyFud $(TPATH)perf_polyFudBpp

lib/polarUtility.o : src/polarUtility.hpp src/polarUtility.cpp 
	@echo "[Compile] polar utility functions"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) src/polarUtility.cpp -o lib/polarUtility.o

#unit tests
$(TPATH)uTests : $(TPATH)uTests.o $(OBJS) 
	@echo "[Link] uTests"
	@$(CC) $^ $(BOOSTLIB) $(LFLAGS) $(LIBSTEST) -o bin/uTests

$(TPATH)uTests.o : unit_tests/uTests.cpp
	@echo "[Compile] uTests"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) unit_tests/uTests.cpp -o $(TPATH)uTests.o

$(TPATH)hgvsParser.o : src/hgvsParser.hpp src/hgvsParser.cpp
	@echo "[Compile] HGVS Parser"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) src/hgvsParser.cpp -o $(TPATH)hgvsParser.o

$(TPATH)readKnownPolyA.o : perf_testing/readKnownPolyA.hpp perf_testing/readKnownPolyA.cpp 
	@echo "[Compile] readKnownPolyA"
	@$(CC) $(BOOST) $(SEQAN) $(LIBPATH) $(CFLAGS) $(LIBS) perf_testing/readKnownPolyA.cpp -o $(TPATH)readKnownPolyA.o

$(TPATH)seqStruct.o : src/seqStruct.hpp src/seqStruct.cpp src/hgvsParser.hpp
	@echo "[Compile] SeqStruct"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) src/seqStruct.cpp -o $(TPATH)seqStruct.o

$(TPATH)utr3Finder.o : src/utr3Finder.hpp src/utr3Finder.cpp src/hgvsParser.hpp src/seqStruct.hpp
	@echo "[Compile] UTR3Finder"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) src/utr3Finder.cpp -o $(TPATH)utr3Finder.o

$(TPATH)polyFud.o : src/polyFud.hpp src/polyFud.cpp src/utr3Finder.hpp
	@echo "[Compile] PolyFud"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) src/polyFud.cpp -o $(TPATH)polyFud.o

$(TPATH)polyFudBpp.o : src/polyFudBpp.hpp src/polyFudBpp.cpp src/polyFud.hpp src/utr3Finder.hpp
	@echo "[Compile] PolyFud-BPP"
	@$(CC) $(BOOST) $(SEQAN) $(VIENNA) $(BOOSTLIB) $(CFLAGS) $(LIBSPOLYFUD) src/polyFudBpp.cpp -o $(TPATH)polyFudBpp.o

$(TPATH)createTPset.o : perf_testing/createTPset.hpp perf_testing/createTPset.cpp
	@echo "[Compile] TP data set"@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) perf_testing/createTPset.cpp -o $(TPATH)createTPset.o
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) perf_testing/createTPset.cpp -o $(TPATH)createTPset.o

$(TPATH)createTNset.o : perf_testing/createTNset.hpp perf_testing/createTNset.cpp
	@echo "[Compile] TN data set"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) perf_testing/createTNset.cpp -o $(TPATH)createTNset.o

$(TPATH)perf_polyFud : $(TPATH)perf_polyFud.o \
	lib/polarUtility.o \
	$(TPATH)createTPset.o \
	$(TPATH)createTNset.o \
	$(TPATH)hgvsParser.o \
	$(TPATH)polyFud.o \
	$(TPATH)refGeneParser.o \
	$(TPATH)readKnownPolyA.o \
	$(TPATH)seqStruct.o \
	$(TPATH)utr3Finder.o 
	@echo "[Link] perf_polyFud"
	@$(CC) $(BOOST) $(SEQAN) $^ $(BOOSTLIB) $(BOOSTLIBS) $(LFLAGS) $(LIBS) -o $(TPATH)perf_polyFud

$(TPATH)perf_polyFud.o : perf_testing/perf_polyFud.hpp perf_testing/perf_polyFud.cpp src/refGeneParser.hpp src/polarUtility.hpp
	@echo "[Compile] perf_polyFud"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) perf_testing/perf_polyFud.cpp -o $(TPATH)perf_polyFud.o

$(TPATH)refGeneParser.o : src/refGeneParser.hpp src/refGeneParser.cpp
	@echo "[Compile] RefGeneParser"
	@$(CC) $(BOOST) $(SEQAN) $(BOOSTLIB) $(CFLAGS) $(LIBS) src/refGeneParser.cpp -o $(TPATH)refGeneParser.o

$(TPATH)perf_polyFudBpp : $(TPATH)perf_polyFudBpp.o \
	$(TPATH)hgvsParser.o \
	lib/polarUtility.o \
	$(TPATH)refGeneParser.o \
	$(TPATH)utr3Finder.o \
	$(TPATH)polyFud.o \
	$(TPATH)seqStruct.o \
	$(TPATH)polyFudBpp.o
	@echo "[Link] perf_polyFud-BPP"
	@$(CC) $(BOOST) $(SEQAN) $^ $(BOOSTLIB) $(VIENNA) $(VIENNALIB) $(LFLAGS) $(LIBSPOLYFUD) -o $(TPATH)perf_polyFudBpp

$(TPATH)perf_polyFudBpp.o : perf_testing/perf_polyFudBpp.hpp perf_testing/perf_polyFudBpp.cpp src/polarUtility.hpp
	@echo "[Compile] perf_polyFud-BPP"
	@$(CC) $(BOOST) $(BOOSTLIB) $(SEQAN) $(VIENNA) $(VIENNALIB) $(CFLAGS) $(LIBSPOLYFUD) perf_testing/perf_polyFudBpp.cpp -o $(TPATH)perf_polyFudBpp.o

.PHONY : clean
clean :
	@echo "[Delete] object and binary files"
	@rm $(TPATH)perf_polyFud $(TPATH)uTests $(TPATH)perf_polyFudBpp $(TPATH)*.o 
