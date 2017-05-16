CC=mpicc
LD = mpicc
CFLAGS= -Ofast -Wall -std=c99
LIBS=-lmpi -lm

OBJS = wave.o
EXECUTABLE = wave
all: $(EXECUTABLE)
$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $(EXECUTABLE)
wave.o: wave.c
	$(CC) $(CFLAGS) -c wave.c
clean: 
	-rm -f wave.o
	-rm -f wave
