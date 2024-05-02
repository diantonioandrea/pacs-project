.PHONY: tests clean distclean
CXXFLAGS = -Wall -pedantic -std=c++20 -I./include -O3

ifneq ($(PACS_ROOT),)
LDFLAGS += -L$(PACS_ROOT)/lib
LDLIBS += -lpacs
endif

# Further optimization.
# CXXFLAGS += -DNDEBUG

# Parallel computing.
ifneq ($(mkTbbLib),)
CXXFLAGS += -DPARALLEL_PACS
LDFLAGS += -L$(mkTbbLib)
LDLIBS += -ltbb
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
