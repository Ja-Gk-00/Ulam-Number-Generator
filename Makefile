CC     = gcc
CFLAGS = -O3 -std=c11 -march=native -fPIC -Wall -Wextra
SRC    = c_src/ulam.c
HDR    = c_src/ulam.h

UNAME := $(shell uname -s 2>/dev/null || echo Windows)

ifeq ($(UNAME),Darwin)
  LIB     = ulam/libulam.dylib
  LDFLAGS = -shared
else ifeq ($(UNAME),Windows)
  LIB     = ulam/ulam.dll
  LDFLAGS = -shared
  CFLAGS  := $(filter-out -fPIC,$(CFLAGS))
else
  LIB     = ulam/libulam.so
  LDFLAGS = -shared -lm
endif

.PHONY: all clean test format format-check

all: $(LIB)

$(LIB): $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<
	@echo "Built: $@"

clean:
	rm -f ulam/ulam.dll ulam/libulam.so ulam/libulam.dylib

test: $(LIB)
	python -m pytest tests/ -v

format:
	clang-format -i $(SRC) $(HDR)

format-check:
	clang-format --dry-run --Werror $(SRC) $(HDR)
