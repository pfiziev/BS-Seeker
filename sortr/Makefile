CC = gcc
CFLAGS = -O3 -Wall 
CFLAGS += -D_FILE_OFFSET_BITS=64
CFLAGS += -D_LARGEFILE_SOURCE
CFLAGS += -fno-exceptions
CFLAGS += -finline-functions
CFLAGS += -funroll-loops

all: clean sortr
sortr:
	$(CC) $(CFLAGS) sortr.c utils.c -o sortr
	
clean:  
	rm -f sortr 
