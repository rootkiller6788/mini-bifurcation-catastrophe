CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libhomoclinic_bifurcation.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
.PHONY: all builddir clean test
all: builddir $(LIB)
builddir:;@mkdir -p build
$(LIB): $(OBJS);$(AR) rcs $@ $^
build/%.o: src/%.c | builddir;$(CC) $(CFLAGS) -c $< -o $@
test: $(LIB);$(CC) $(CFLAGS) tests/test_homoclinic.c -Lbuild -lhomoclinic_bifurcation $(LDFLAGS) -o build/test.exe;./build/test.exe
clean:;rm -rf build
