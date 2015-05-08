CC=gcc
CFLAGS=-I.
DEPS = ../include/hellomake.h
OBJ = hellomake.o hellofunc.o 

%.o: %.c $(DEPS)
	echo $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hellomake: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)
