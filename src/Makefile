CXX=g++ -m64
CXXFLAGS=-g -std=c++11 -Icommon -Iobjs/ -IEigen/ -O0 -Wall

APP_NAME=pgmres
OBJDIR=objs
SRCS=utils.cpp loadmtx.cpp mtxvec.cpp main.cpp
OBJS=$(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

.PHONY: dirs clean

default: $(APP_NAME)

dirs:
	@/bin/mkdir -p $(OBJDIR)/

clean:
	/bin/rm -rf $(OBJDIR) *~ $(APP_NAME)

$(APP_NAME): dirs $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm -lpthread

$(OBJDIR)/%.o: %.cpp $(DEPS)
	$(CXX) -MMD $(CXXFLAGS) -c -o $@ $<

-include $(OBJDIR)/*.d
