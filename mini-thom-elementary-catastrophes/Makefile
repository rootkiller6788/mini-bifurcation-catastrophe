CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libtec.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(patsubst tests/%.c,build/test_%,$(TEST_SRCS))
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/example_%,$(EXAMPLE_SRCS))

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

build/test_%: tests/%.c $(LIB)
	$(CC) $(CFLAGS) $< -Lbuild -ltec $(LDFLAGS) -o $@

build/example_%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) $< -Lbuild -ltec $(LDFLAGS) -o $@

test: $(LIB) $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do $$t || exit 1; done
	@echo "=== All tests completed ==="

examples: $(LIB) $(EXAMPLE_BINS)
	@echo "=== Running Examples ==="
	@for e in $(EXAMPLE_BINS); do $$e; done

clean:
	rm -rf build
