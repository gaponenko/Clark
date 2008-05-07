# First get the general variable definitions.

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

DEPENDENCIES  = $(shell ls *.C)
DEPENDENCIES2 = $(shell ls *.h)

# includes and libs
INCLUDE = -I/usr/include -I/twist/local/include

CXX = g++

# Flags for the compiler
CXXFLAGS = -g -O2 -m32 
CXXFLAGS += $(ROOTCFLAGS) -DEXE

# Flags for the linker
LDFLAGS = -g -O2 -m32

LINKER = g++

AR = ar  -rc
RANLIB = ranlib


LIBS		+= $(ROOTLIBS) -L/twist/local/lib -llog4cpp -L/usr/lib -lgsl -lgslcblas -lm -lboost_regex

OBJ = CommandLine.o ConfigFile.o DetectorGeo.o EventClass.o HistogramFactory.o TreeClass.o Users.o FuncLib.o EventLib.o DefaultConfig.o Clark.o
EXE = Clark

all : link



Clark.o : Clark.cpp $(DEPENDENCIES2)
	@echo "  CC  $<"
	@$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} $<

DefaultConfig.o : DefaultConfig.cpp
	@echo "  CC  $<"
	@$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} $<

%.o : %.C $(DEPENDENCIES) $(DEPENDENCIES2)
	@echo "  CC  $<"
	@$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} $<

%.o : %.cpp %.h
	@echo "  CC  $<"
	@$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} $<

link: $(OBJ)
	@echo "  LD  ${EXE}"
	@${LINKER} -o $(EXE) ${LDFLAGS} $(OBJ) $(OBJ2) ${LIBS}


clean:
	@echo "  RM  ${EXE} ${OBJ} core.*"
	@rm -rf ${EXE} ${OBJ} core.*

# install:
# 	@mkdir -p ${BINPREFIX}
# 	@chmod 755 ${EXE}
# 	@cp -f ${EXE} ${BINPREFIX}
# 	@echo installed ${EXE} to ${BINPREFIX}
	

