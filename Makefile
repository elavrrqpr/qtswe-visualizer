CC      ?= gcc
CFLAGS  ?= -O3 -Wall -Wextra -std=gnu11
OPENMP  ?= 1
ifeq ($(OPENMP),1)
CFLAGS  += -fopenmp
endif
LDLIBS  += -lm
EXE     :=
ifeq ($(OS),Windows_NT)
CFLAGS  += -D__USE_MINGW_ANSI_STDIO=1
LDLIBS  += -lws2_32
EXE     := .exe
endif

qtswe$(EXE): main.c qtswe.c qtswe.h
	$(CC) $(CFLAGS) -o $@ main.c qtswe.c $(LDLIBS)

clean:
	rm -f qtswe qtswe.exe
.PHONY: clean
