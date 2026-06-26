CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
AR = ar
LIB = build/libnormalform.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_nft.c
TEST_BIN = build/test_nft
.PHONY: all builddir clean test
all: builddir $(LIB)
builddir:; mkdir -p build
$(LIB): $(OBJS); $(AR) rcs $@ $^
build/%.o: src/%.c; $(CC) $(CFLAGS) -c $< -o $@
test: all; $(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lnormalform $(LDFLAGS) && ./$(TEST_BIN)
clean:; rm -rf build
