CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
AR = ar rcs
LIB = libbcd.a
SRCS = $(wildcard src/*.c)
OBJS = $(SRCS:.c=.o)

.PHONY: all clean test

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB)
	@for t in tests/test_*.c; do \
		echo "--- $$t ---"; \
		$(CC) $(CFLAGS) $$t -L. -lbcd -lm -o $$t.out && ./$$t.out; \
		rm -f $$t.out; \
	done

clean:
	rm -f src/*.o $(LIB) tests/*.out
