IDIR =../include
CC=g++

#Use these flags if you don't want to use optimization profiles
CFLAGS=-I$(IDIR) -Wall -O3 -fopenmp

#CFLAGS=-I$(IDIR) -O3 -fprofile-generate=./profiles
#CFLAGS=-I$(IDIR) -O3 -fprofile-use=./profiles

ODIR=.
LDIR =.

LIBS=-lm -lgsl -lgslcblas -lstdc++

_DEPS = xy_class.h init_functions.h params.h pa_class.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ = main.o xy_class.o init_functions.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

.PHONY: all

all : main auto pa_sim

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
wolff_main: wolff_main.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

main: main.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

pa_sim: pa_main.o xy_class.o init_functions.o pa_class.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

square: square.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

ribbon: ribbon.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#auto: correlation.o xy_class.o init_functions.o
#auto: correlation.o xy_class.o init_functions.o
#	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

auto: correlation.o xy_class.o init_functions.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
