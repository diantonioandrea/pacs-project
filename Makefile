.PHONY: test clean distclean
CXXFLAGS = -Wall -pedantic -std=c++20 -I./include -O3

# Further optimization.
# CXXFLAGS += -DNDEBUG

# Verbosity.
CXXFLAGS += -DVERBOSE

# Parallel computing using OpenMP.
ifneq ($(OpenMP),) # $(OpenMP) set to /path/to/libomp.
ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -Xclang
endif
CXXFLAGS += -fopenmp
CPPFLAGS += -I$(OpenMP)/include
LDFLAGS += -L$(OpenMP)/lib
LDLIBS += -lomp
else
CXXFLAGS += -Wno-unknown-pragmas
endif

OBJECTS = $(subst .cpp,.o,$(shell ls ./src))
HEADERS = ./include/*.hpp

TEST_EXECS = $(subst .cpp,.out,$(shell ls ./test))
TEST_OBJECTS = $(subst .cpp,.o,$(shell ls ./test))

# Test.
test: $(TEST_EXECS)
	@echo "Done!"

$(TEST_EXECS): %.out: %.o $(OBJECTS)
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $^ to $@"; else echo "Linking $^ to $@ with the following flags: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(TEST_OBJECTS): %.o: test/%.cpp
	@echo "Compiling $< using $(CXX) with the following flags: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJECTS): %.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with the following flags: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Clean.
clean:
	@echo "Cleaning the repo."
	@$(RM) ./*.o
	@$(RM) ./*.poly

distclean: clean
	@$(RM) ./*.out
