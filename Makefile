.PHONY: test run clean distclean
CXXFLAGS = -Wall -pedantic -std=c++20 -I./include -O3

# Further optimization.
# CXXFLAGS += -DNDEBUG

# Verbosity.
# CPPFLAGS += -DVERBOSE

# Dynamic compression.
CPPFLAGS += -DDYNAMIC_SPARSE

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

TEST_RUN = $(subst .cpp,,$(shell ls ./test))
TEST_EXECS = $(subst .cpp,.out,$(shell ls ./test))
TEST_OBJECTS = $(subst .cpp,.o,$(shell ls ./test))

OUTPUT_DIR = ./output

# Test.
test: $(TEST_EXECS)
	@echo "Done!"

run: $(TEST_RUN)
	@echo "Done!"

$(TEST_RUN): $(TEST_EXECS) $(OUTPUT_DIR)
	@echo "Executing ./$@ and redirecting the output to $(OUTPUT_DIR)/$@.txt"
	@./$@.out > $(OUTPUT_DIR)/$@.txt

$(TEST_EXECS): %.out: %.o $(OBJECTS)
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $^ to $@"; else echo "Linking $^ to $@ with the following flags: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(TEST_OBJECTS): %.o: test/%.cpp
	@echo "Compiling $< using $(CXX) with the following flags: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OBJECTS): %.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with the following flags: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(OUTPUT_DIR):
	@mkdir -p $(OUTPUT_DIR)

# Clean.
clean:
	@echo "Cleaning the repo."
	@$(RM) ./*.o
	@$(RM) ./*.poly
	@$(RM) -r $(OUTPUT_DIR)

distclean: clean
	@$(RM) ./*.out
