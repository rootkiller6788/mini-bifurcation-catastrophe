CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libpdepattern.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_pde.c
TEST_BIN = build/test_pde

all: builddir $(LIB)
builddir:
	@mkdir -p build
$(LIB): $(OBJS)
	$(AR) rcs $@ $^
build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
test: builddir $(LIB)
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lpdepattern $(LDFLAGS)
	./$(TEST_BIN)
examples: builddir $(LIB)
	for f in examples/ex*.c; do \
	  $(CC) $(CFLAGS) -o build/$$(basename $$f .c) $$f -Lbuild -lpdepattern $(LDFLAGS); \
	  ./build/$$(basename $$f .c); \
	done
clean:
	rm -rf build
