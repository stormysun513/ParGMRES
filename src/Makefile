CXX=g++ -m64
CXXFLAGS=-std=c++11 -Iobjs/ -IEigen/ -O3 -Wall -fopenmp

APP=pgmres
TEST=svdtest
OBJDIR=objs

APP_MAIN=main.cpp
TEST_MAIN=svdtest.cpp
SRCS=utils.cpp loadmtx.cpp mtxvec.cpp eigen.cpp
OBJS=$(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

.PHONY: dirs clean

default: dirs $(APP)

dirs:
	@/bin/mkdir -p $(OBJDIR)

clean:
	/bin/rm -rf $(OBJDIR) *~ $(APP) $(TEST)

test: $(TEST)

$(TEST): $(TEST_MAIN) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(TEST_MAIN) $(OBJS)

$(APP): $(APP_MAIN) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(APP_MAIN) $(OBJS)

$(OBJDIR)/%.o: %.cpp $(DEPS)
	$(CXX) -MMD $(CXXFLAGS) -c -o $@ $<

-include $(OBJDIR)/*.d
