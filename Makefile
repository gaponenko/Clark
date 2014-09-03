# First get the general variable definitions.

ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLDFLAGS    := $(shell root-config --ldflags)
ROOTLIBS      := -L$(shell root-config --libdir) -lGenVector $(shell root-config --libs)

# includes and libs
# On detsim we use products from the UPS areas
ifeq ($(UPS_DIR),)
	INCLUDE = -I/home/ahillair/local/include -I/home/ahillair/TWIST/Software/2014_02_02_Unfolding/RooUnfold-1.1.1
	LIBS += $(ROOTLIBS) -L/home/ahillair/local/lib -L/home/ahillair/TWIST/Software/2014_02_02_Unfolding/RooUnfold-1.1.1 -lRooUnfold -llog4cpp -lgsl -lgslcblas -lm -lboost_regex
else
	INCLUDE = -I$(LOG4CPP_DIR)/include -I$(GSL_DIR)/include -I$(BOOST_INC)
	LIBS += $(ROOTLIBS) -L$(LOG4CPP_DIR)/lib  -L$(GSL_DIR)/lib -L$(BOOST_LIB) -llog4cpp -lgsl -lgslcblas -lm -lboost_regex
endif

# Flags for the compiler
CXX = g++
CXXFLAGS := -g -O2 $(ROOTCFLAGS) -Wall -Wno-parentheses -Wno-sign-compare

# Flags for the linker
LINKER = g++
LDFLAGS := -g -O2 $(ROOTLDFLAGS) -Xlinker -rpath=/home/ahillair/local/lib -Xlinker -rpath=$(shell root-config --libdir)  -Xlinker -rpath=/home/ahillair/TWIST/Software/2014_02_02_Unfolding/RooUnfold-1.1.1

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
