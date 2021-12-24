CC = mpicc 
UNDFAST= -funroll-all-loops -ffast-math -fprefetch-loop-arrays
CFLAGS = -I$(UNDMOL13DIR)/include $(UNDCFLAGS) $(UNDFAST) 
LDFLAGS = -L$(UNDMOL13DIR)/lib -lundmol -lundmolmpi -lm

EXECUTABLE = tqcorr_parallel.exe

SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c, %.o, $(SRCS))

$(EXECUTABLE): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -g -o $(EXECUTABLE) 
	mv $(EXECUTABLE) $(UNDMOL13DIR)/bin
	rm -f *.o

%.o: %.c
	$(CC) $(CFLAGS) -c -g -o $@ $< 

clean:
	rm -f *.o
