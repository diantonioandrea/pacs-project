.PHONY: all lib install tests examples domains testrun clean distclean
CXXFLAGS = -Wall -Wno-sign-compare -pedantic -std=c++20 -march=native -fPIC -I./include -O2 -fno-unsafe-math-optimizations -fno-fast-math

ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -ffp-model=precise
endif


# FLAGS.

# Disables solution output.
CPPFLAGS += -DNSOLUTIONS

# Disables verbosity.
# CPPFLAGS += -DNVERBOSE

# Disables debugging. Enhances performance.
# CPPFLAGS += -DNDEBUG


# PARALLEL COMPUTING.

# # Parallel computing using STL with modules. Disabled by default.
# ifneq ($(mkTbbLib),)
# CPPFLAGS += -DPARALLEL
# CXXFLAGS += -I$(mkTbbInc)
# LDFLAGS += -L$(mkTbbLib)
# LDLIBS += -ltbb
# endif

# Parallel computing using OpenMP. Enabled by default.
ifneq ($(OpenMP),) # $(OpenMP) set to /path/to/libomp.
ifeq ($(shell uname),Darwin) # Apple's clang.
CXXFLAGS += -Xclang
LDLIBS += -lomp
else
LDLIBS += -lgomp
endif
CXXFLAGS += -fopenmp
CPPFLAGS += -I$(OpenMP)/include
LDFLAGS += -L$(OpenMP)/lib
else
ifneq ($(mkPrefix),) # Parallel computing using OpenMP with modules.
CXXFLAGS += -fopenmp
LDFLAGS += -L$(mkToolchainPrefix)/lib
LDLIBS += -lgomp
else
CXXFLAGS += -Wno-unknown-pragmas
endif
endif


# Directories.
OUTPUT_DIR = ./output
OBJECT_DIR = ./objects
LIB_DIR = ./lib
EXEC_DIR = ./executables

# Library.
LIBRARY_NAME = PacsHPDG
LIBRARY = $(LIB_DIR)/lib$(LIBRARY_NAME).a
INCLUDE_DESTINATION = $(HOME)/include
LIB_DESTINATION = $(HOME)/lib

# Files.
OBJECTS = $(subst src/,$(OBJECT_DIR)/,$(subst .cpp,.o,$(shell find src -name "*.cpp")))

HEADERS = ./include/*.hpp # Recompilation purposes.
HEADERS += ./include/PacsHPDG/*.hpp
HEADERS += ./include/PacsHPDG/Algebra/*.hpp
HEADERS += ./include/PacsHPDG/Algebra/Methods/*.hpp
HEADERS += ./include/PacsHPDG/Geometry/*.hpp
HEADERS += ./include/PacsHPDG/Fem/*.hpp
HEADERS += ./include/PacsHPDG/Laplacian/*.hpp
HEADERS += ./include/PacsHPDG/Statistics/*.hpp

EXAMPLE_EXECS = $(subst examples/,$(EXEC_DIR)/,$(subst .cpp,.out,$(shell find examples -name "*.cpp")))
EXAMPLE_OBJECTS = $(subst examples/,$(OBJECT_DIR)/,$(subst .cpp,.o,$(shell find examples -name "*.cpp")))
HEADERS += ./examples/*.hpp

DOMAIN_EXECS = $(subst domains/,$(EXEC_DIR)/,$(subst .cpp,.out,$(shell find domains -name "*.cpp")))
DOMAIN_OBJECTS = $(subst domains/,$(OBJECT_DIR)/,$(subst .cpp,.o,$(shell find domains -name "*.cpp")))

TEST_RUN = $(subst .cpp,,$(shell ls ./test))
TEST_EXECS = $(subst test/,$(EXEC_DIR)/,$(subst .cpp,.out,$(shell find test -name "*.cpp")))
TEST_OBJECTS = $(subst test/,$(OBJECT_DIR)/,$(subst .cpp,.o,$(shell find test -name "*.cpp")))

# All.
all: tests examples domains

# Library.
lib: $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	@echo "Archiving the library to $(LIBRARY)"
	@ar rcs $(LIBRARY) $(OBJECTS)

ifeq ($(shell find . -wholename $(LIBRARY)), $(LIBRARY)) # Available only after compilation.
install: # Manual spacing for consistency between platforms.
	@echo "Installing the library"
	@echo "    - Includes under $(INCLUDE_DESTINATION)"
	@echo "    - Lib under $(LIB_DESTINATION)"
	@echo "Remember: -I$(INCLUDE_DESTINATION) -L$(LIB_DESTINATION) -l$(LIBRARY_NAME)"

	@mkdir -p $(INCLUDE_DESTINATION)
	@mkdir -p $(LIB_DESTINATION)

	@cp -r ./include/* $(INCLUDE_DESTINATION)/
	@cp $(LIBRARY) $(LIB_DESTINATION)/
endif

# Test.
tests: $(EXEC_DIR) $(OUTPUT_DIR) $(TEST_EXECS)
	@echo "Compiled tests!"

testrun: $(EXEC_DIR) $(OUTPUT_DIR) $(TEST_RUN) 
	@echo "Run tests!"

$(TEST_RUN): $(TEST_EXECS)
	@echo "Executing $(EXEC_DIR)/$@ and redirecting the output to $(OUTPUT_DIR)/$@.txt"
	@$(EXEC_DIR)/$@.out > $(OUTPUT_DIR)/$@.txt

$(TEST_EXECS): executables/%.out: $(OBJECT_DIR)/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@"; else echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(TEST_OBJECTS): $(OBJECT_DIR)/%.o: test/%.cpp $(HEADERS) $(OBJECT_DIR)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Examples.
examples: $(OBJECT_DIR) $(EXEC_DIR) $(OUTPUT_DIR) $(EXAMPLE_EXECS)
	@echo "Compiled examples!"

$(EXAMPLE_EXECS): executables/%.out: $(OBJECT_DIR)/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@"; else echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(EXAMPLE_OBJECTS): $(OBJECT_DIR)/%.o: examples/%.cpp $(HEADERS) $(OBJECT_DIR)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Domains.
domains: $(OBJECT_DIR) $(EXEC_DIR) $(OUTPUT_DIR) $(DOMAIN_EXECS)
	@echo "Compiled domains!"

$(DOMAIN_EXECS): executables/%.out: $(OBJECT_DIR)/%.o $(OBJECTS) 
	@if [ "$(LDFLAGS) $(LDLIBS)" = " " ]; then echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@"; else echo "Linking $(subst $(OBJECT_DIR)/,,$<) and base objects to $@ with: $(LDFLAGS) $(LDLIBS)"; fi
	@$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(DOMAIN_OBJECTS): $(OBJECT_DIR)/%.o: domains/%.cpp $(HEADERS) $(OBJECT_DIR)
	@echo "Compiling $< using $(CXX) with: $(CXXFLAGS) $(CPPFLAGS)"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Objects.
$(OBJECTS): ./$(OBJECT_DIR)/%.o: src/%.cpp $(HEADERS) $(OBJECT_DIR)
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
	@$(RM) -r $(LIB_DIR)
	@$(RM) -r $(EXEC_DIR)
