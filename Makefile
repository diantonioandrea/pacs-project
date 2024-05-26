.PHONY: test example testrun clean distclean
CXXFLAGS = -Wall -Wno-sign-compare -pedantic -std=c++20 -I./include -O3 -fPIC

# Verbosity.
# CPPFLAGS += -DNVERBOSE

# Further optimization.
# CPPFLAGS += -DNDEBUG

# Parallel computing using OpenMP.
# ifneq ($(OpenMP),) # $(OpenMP) set to /path/to/libomp.
# ifeq ($(shell uname),Darwin) # Apple's clang.
# CXXFLAGS += -Xclang
# endif
# CXXFLAGS += -fopenmp
# CPPFLAGS += -I$(OpenMP)/include
# LDFLAGS += -L$(OpenMP)/lib
# LDLIBS += -lomp
# else
# CXXFLAGS += -Wno-unknown-pragmas
# endif

# Parallel computing using STL and modules.
ifneq ($(mkTbbLib),)
LDFLAGS += -L$(mkTbbLib)
LDLIBS += -ltbb
endif

# Files.
OBJECTS = $(subst src/,objects/,$(subst .cpp,.o,$(shell find src -name "*.cpp")))
HEADERS = ./include/*.hpp

EXAMPLE_EXECS = $(subst example/,executables/,$(subst .cpp,.out,$(shell find example -name "*.cpp")))
EXAMPLE_OBJECTS = $(subst example/,objects/,$(subst .cpp,.o,$(shell find example -name "*.cpp")))

TEST_RUN = $(subst .cpp,,$(shell ls ./test))
TEST_EXECS = $(subst test/,executables/,$(subst .cpp,.out,$(shell find test -name "*.cpp")))
TEST_OBJECTS = $(subst test/,objects/,$(subst .cpp,.o,$(shell find test -name "*.cpp")))

# Directories.
OUTPUT_DIR = ./output
OBJECT_DIR = ./objects
EXEC_DIR = ./executables

# All.
all: test example $(OUTPUT_DIR)

# Test.
test: $(OBJECT_DIR) $(EXEC_DIR) $(TEST_EXECS)
	@echo "Compiled tests!"

testrun: $(OBJECT_DIR) $(EXEC_DIR) $(OUTPUT_DIR) $(TEST_RUN) 
	@echo "Run tests!"

$(TEST_RUN): $(TEST_EXECS)
	@echo "Executing $(EXEC_DIR)/$@ and redirecting the output to $(OUTPUT_DIR)/$@.txt"
	@$(EXEC_DIR)/$@.out > $(OUTPUT_DIR)/$@.txt

$(TEST_EXECS): executables/%.out: objects/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst objects/,,$<) and base objects to $@"; else echo "Linking $(subst objects/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(TEST_OBJECTS): objects/%.o: test/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Examples.
example: $(OBJECT_DIR) $(EXEC_DIR) $(EXAMPLE_EXECS)
	@echo "Compiled examples!"

$(EXAMPLE_EXECS): executables/%.out: objects/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst objects/,,$<) and base objects to $@"; else echo "Linking $(subst objects/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(EXAMPLE_OBJECTS): objects/%.o: example/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Objects.
$(OBJECTS): ./objects/%.o: src/%.cpp $(HEADERS)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Directories.
$(OUTPUT_DIR):
	@mkdir -p $(OUTPUT_DIR)

$(EXEC_DIR):
	@mkdir -p $(EXEC_DIR)

$(OBJECT_DIR):
	@mkdir -p $(OBJECT_DIR)

# Clean.
clean:
	@echo "Cleaning the repo."
	@$(RM) -r $(OBJECT_DIR)
	@$(RM) -r $(OUTPUT_DIR)

distclean: clean
	@$(RM) -r $(EXEC_DIR)
