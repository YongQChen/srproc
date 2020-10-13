
CC = g++
CFLAGS = -Wall -O2
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: srproc

srproc: $(OBJFILES)
	$(CC) -o srproc $(OBJFILES)

%.o: %.cpp
	$(COMPILE) -o $@ $<

#postrun: $(srproc)
all: 
	rm *.o;
