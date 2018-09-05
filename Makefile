# First get the general variable definitions.

OPTFLAGS      := -O2

ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLDFLAGS    := $(shell root-config --ldflags)
ROOTLIBS      := -L$(shell root-config --libdir) -lGenVector $(shell root-config --libs)

# includes and libs
# On detsim we use products from the UPS areas
ifeq ($(UPS_DIR),)
        #mu2esoft gone: INCLUDE = -I/twist/local/include -I/twist/datap/muminus/mu2esoft/RooUnfold/trunk
        #mu2esoft gone: LIBS += $(ROOTLIBS) -L/twist/local/lib -L/twist/datap/muminus/mu2esoft/RooUnfold/trunk -lRooUnfold -llog4cpp -lgsl -lgslcblas -lm -lboost_regex
	INCLUDE = -I/twist/local/include
	LIBS += $(ROOTLIBS) -L/twist/local/lib -llog4cpp -lgsl -lgslcblas -lm -lboost_regex
else
	INCLUDE =  $(shell log4cpp-config --cflags) $(shell gsl-config --cflags) -I$(BOOST_INC)
	LIBS += $(ROOTLIBS) $(shell log4cpp-config --libs) $(shell gsl-config --libs) -L$(BOOST_LIB) -lboost_regex
endif

# Flags for the compiler
CXX = g++
CXXFLAGS := -g $(OPTFLAGS) $(ROOTCFLAGS) -Wall -Wno-parentheses -Wno-sign-compare -std=c++0x

# Flags for the linker
LINKER = g++
#mu2esoft gone: LDFLAGS := -g $(OPTFLAGS) $(ROOTLDFLAGS) -Xlinker -rpath=/twist/local/lib -Xlinker -rpath=$(shell root-config --libdir)  -Xlinker -rpath=/twist/datap/muminus/mu2esoft/RooUnfold/trunk
LDFLAGS := -g $(OPTFLAGS) $(ROOTLDFLAGS) -Xlinker -rpath=/twist/local/lib -Xlinker -rpath=$(shell root-config --libdir)

#
HEADERS := $(shell ls *.h)
SRCS := $(shell ls *.cpp)
OBJS := $(SRCS:%.cpp=%.o)
EXE := Clark

all : $(EXE)

# Users.C includes all other .C files, so we must compile and link just this object.
MODULES  := $(shell ls *.C)
Users.o : $(MODULES) $(HEADERS)
	$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} Users.C

# Rule for the "normal" C++ sources, with excessive dependencies
%.o : %.cpp $(HEADERS)
	$(CXX) -o $@ $(CXXFLAGS) -c ${INCLUDE} $<

$(EXE): $(OBJS) Users.o
	${LINKER} -o $(EXE) ${LDFLAGS} $^ ${LIBS}

clean:
	rm -rf ${EXE} ${OBJS} Users.o core.*

# install:
# 	@mkdir -p ${BINPREFIX}
# 	@chmod 755 ${EXE}
# 	@cp -f ${EXE} ${BINPREFIX}
# 	@echo installed ${EXE} to ${BINPREFIX}
