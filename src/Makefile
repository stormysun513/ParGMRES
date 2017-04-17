CXX=g++ -m64
CXXFLAGS=-std=c++11 -Iobjs/ -IEigen/ -O3 -Wall

APP_NAME=pgmres
OBJDIR=objs
SRCS=utils.cpp loadmtx.cpp mtxvec.cpp main.cpp
OBJS=$(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

.PHONY: dirs clean

default: dirs $(APP_NAME)

dirs:
	@/bin/mkdir -p $(OBJDIR)

clean:
	/bin/rm -rf $(OBJDIR) *~ $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm -lpthread

$(OBJDIR)/%.o: %.cpp $(DEPS)
	$(CXX) -MMD $(CXXFLAGS) -c -o $@ $<

-include $(OBJDIR)/*.d
