CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libbifurcation.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_bifurcation.c
TEST_BIN = build/test_bifurcation
.PHONY: all builddir clean test
all: builddir $(LIB)
builddir:; @mkdir -p build
$(LIB): $(OBJS); $(AR) rcs $@ $^
build/%.o: src/%.c; $(CC) $(CFLAGS) -c $< -o $@
test: all; $(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lbifurcation $(LDFLAGS) && ./$(TEST_BIN)
clean:; rm -rf build
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
examples: all $(EXAMPLE_BINS); @echo "All examples built."
build/%: examples/%.c $(LIB); $(CC) $(CFLAGS) -o $@ $< -Lbuild -lbifurcation $(LDFLAGS)
demo: examples; @for e in $(EXAMPLE_BINS); do echo "=== $$(basename $$e) ===" && ./$$e; done
