EXECUTABLE := gmres
TEST 	   := test
CU_FILES   := gmres.cu test.cu
CU_DEPS    :=
CC_FILES   := main.cpp
LOGS       := logs

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')
OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O3 -Wall -g -std=c++11
HOSTNAME=$(shell hostname)

LIBS       :=
FRAMEWORKS :=

#NVCCFLAGS=-O3 -m64 --gpu-architecture compute_35
#LIBS += GL glut cudart
INC=-Icusplibrary

ifneq ($(wildcard /opt/cuda-8.0/.*),)
# Latedays
NVCCFLAGS=-O3 -m64 --gpu-architecture compute_35
LDFLAGS=-L/opt/cuda-8.0/lib64/ -lcudart
else
# GHC or AWS
NVCCFLAGS=-O3 -m64 --gpu-architecture compute_30
LDFLAGS=-L/usr/local/cuda/lib64/ -lcudart
endif

LDLIBS  := $(addprefix -l, $(LIBS))
LDFRAMEWORKS := $(addprefix -framework , $(FRAMEWORKS))

NVCC=nvcc
OBJS=gmres.o mtxvec.o qr.o utils.o
TESTOBJS=test.o gmres.o mtxvec.o qr.o utils.o

.PHONY: clean

default: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(CC_FILES) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

test: $(TESTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

%.o: %.cu
	$(NVCC) $(INC) $< $(NVCCFLAGS) -c -o $@

clean:
	rm -rf *~ $(EXECUTABLE) $(TEST) $(OBJS) $(TESTOBJS) $(LOGS)
