# Ben Postlethwaite 
# Makefile version 3
#
# Set variables and compiler for make file.

CC=g++  
CFLAGS=-c -Wall
PROJECT=driven
OBJECTS=$(PROJECT)Main.o $(PROJECT)Defs.o  

$(PROJECT): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $(PROJECT)

$(PROJECT)Main.o: $(PROJECT)Headers.h

$(PROJECT)Defs.o: $(PROJECT)Headers.h

clean:
	-rm -rf $(OBJECTS)

bigclean:
	-rm -rf $(OBJECTS) *.dat $(PROJECT)