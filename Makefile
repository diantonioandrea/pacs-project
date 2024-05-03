.PHONY: tests clean distclean
CXXFLAGS = -Wall -pedantic -std=c++20 -I./include -O3

# Further optimization.
# CXXFLAGS += -DNDEBUG

# Parallel computing using OpenMP.
ifneq ($(OpenMP),) # $(OpenMP) set to /path/to/libomp.
ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -Xclang
endif
CXXFLAGS += -fopenmp -I$(OpenMP)/include
LDFLAGS += -L$(OpenMP)/lib
LDLIBS += -lomp
endif

# Headers, recompiling purposes.
HEADERS = ./include/*.hpp

# Tests.
TEST_SOURCES = $(shell find tests -name "*.cpp")
TEST_OBJECTS = $(TEST_SOURCES:tests/%.cpp=%.o)
TEST_EXECS = $(TEST_SOURCES:tests/%.cpp=%.out)

tests: $(TEST_EXECS)

$(TEST_EXECS): $(TEST_OBJECTS)
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $^ to $@"; else echo "Linking $< to $@ with the following flags: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $< -o $@

$(TEST_OBJECTS): $(TEST_SOURCES) $(HEADERS)
	@echo "Compiling $< using $(CXX) with the following flags: $(CXXFLAGS)"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean.
clean:
	$(RM) $(TEST_OBJECTS)

distclean: clean
	$(RM) $(TEST_EXECS)
