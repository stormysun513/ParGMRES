EXECUTABLE := test
CU_FILES   := test.cu
CU_DEPS    :=
CC_FILES   := main.cpp
LOGS       := logs

###########################################################

ARCH=$(shell uname | sed -e 's/-.*//g')
OBJDIR=objs
CXX=g++ -m64
CXXFLAGS=-O3 -Wall -g
HOSTNAME=$(shell hostname)

LIBS       :=
FRAMEWORKS :=

NVCCFLAGS=-O3 -m64 --gpu-architecture compute_35
LIBS += GL glut cudart
INC=-Icusplibrary

ifneq ($(wildcard /opt/cuda-8.0/.*),)
# Latedays
LDFLAGS=-L/opt/cuda-8.0/lib64/ -lcudart
else
# GHC
LDFLAGS=-L/usr/local/cuda/lib64/ -lcudart
endif

LDLIBS  := $(addprefix -l, $(LIBS))
LDFRAMEWORKS := $(addprefix -framework , $(FRAMEWORKS))

NVCC=nvcc
OBJS=test.o

.PHONY: clean

default: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS) $(LDFRAMEWORKS)

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

%.o: %.cu
	$(NVCC) $(INC) $< $(NVCCFLAGS) -c -o $@

clean:
	rm -rf *~ $(EXECUTABLE) $(OBJS) $(LOGS) 
