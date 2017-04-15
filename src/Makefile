CXX=g++ -m64
CXXFLAGS=-std=c++11 -Icommon -Iobjs/ -O2 -Wall
# CXXFLAGS=-I../common -Iobjs/ -O2 -Wall -mavx2 -mfma

APP_NAME=pgmres
OBJDIR=objs
CFILE=sin


# default: $(APP_NAME) $(APP_NAME).da $(CFILE).s
default: $(APP_NAME)

.PHONY: dirs clean

dirs:
	/bin/mkdir -p $(OBJDIR)/

clean:
	/bin/rm -rf $(OBJDIR) *~ $(APP_NAME) $(APP_NAME).da

# OBJS=$(OBJDIR)/main.o $(OBJDIR)/sin.o
OBJS=$(OBJDIR)/main.o $(OBJDIR)/utils.o $(OBJDIR)/loadmtx.o

$(APP_NAME): dirs $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm -lpthread

# $(APP_NAME).da: $(APP_NAME)
# 	objdump -d $(APP_NAME) > $(APP_NAME).da

# $(CFILE).s: $(CFILE).cpp
# 	$(CXX) $(CXXFLAGS) -S $(CFILE).cpp

$(OBJDIR)/%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

# $(OBJDIR)/%.o: $(COMMONDIR)/%.cpp
# 	$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/main.o: main.hpp loadmtx.hpp Makefile
$(OBJDIR)/utils.o: utils.hpp
$(OBJDIR)/loadmtx.o: loadmtx.hpp

test: test.cpp loadmtx.cpp
	g++ -std=c++11 -o test test.cpp loadmtx.cpp

