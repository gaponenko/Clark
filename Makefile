# First get the general variable definitions.

ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLDFLAGS    := $(shell root-config --ldflags)
ROOTLIBS      := $(shell root-config --libs)

CXX = g++

# Flags for the compiler
CXXFLAGS = -g -O2 $(ROOTCFLAGS)

# Flags for the linker
LDFLAGS = -g -O2 $(ROOTLDFLAGS)

LINKER = g++

# includes and libs
INCLUDE = -I/usr/include -I/twist/local/include
LIBS += $(ROOTLIBS) -L/twist/local/lib -llog4cpp -L/usr/lib -lgsl -lgslcblas -lm -lboost_regex

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
