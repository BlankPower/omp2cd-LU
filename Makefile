include ../../Makefile.def

SHELL=/bin/bash

TARGET=lu.exe
SOURCES= lu.c c_print_results.c c_timers.c wtime.c

.PHONY: all darts omp clean-darts clean

all: omp

darts:
	$(OMP2CD) $(OMP2CDFLAGS) $(SOURCES) -- $(CLANGFLAGS)
	@clang-format -i -style="{BasedOnStyle: WebKit, ColumnLimit: 100}" ./*.output.darts.*

omp: $(SOURCES:.c=.o)
	$(CC) -O3 $(SOURCES:.c=.o) $(LDFLAGS_OPENMP) -o $(TARGET)

%.o: %.c
	$(CC) -c $< $(CFLAGS_OPENMP) -o $@
	
clean-darts:
	rm -f ./*.output.darts.*

clean:
	rm -f *~ *.o $(TARGET) *.dat.0 *.logf.* *.view.dat.* *.view.sz.*
