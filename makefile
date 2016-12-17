CC=gcc 
CFLAGS=-O3 -c -Wall 
SOURCES=kdtree.c pfield.c main.c 
OBJECTS=$(SOURCES:.c=.o) 
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ -lm
	rm  $(OBJECTS)
.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm  $(EXECUTABLE)
