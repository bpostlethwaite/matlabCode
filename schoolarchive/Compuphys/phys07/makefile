# Ben Postlethwaite 
# Makefile version 3
#
# Set variables and compiler for make file.

CC=g++  
CFLAGS= -g -Wall
PROJECT=propane
OBJECTS=$(PROJECT)Main.o $(PROJECT)Defs.o  

$(PROJECT): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o $(PROJECT)

$(PROJECT)Main.o: $(PROJECT)Main.cpp $(PROJECT)Headers.h
	$(CC) $(CFLAGS) -c  $(PROJECT)Main.cpp 

$(PROJECT)Defs.o: $(PROJECT)Defs.cpp $(PROJECT)Headers.h
	$(CC) $(CFLAGS) -c $(PROJECT)Defs.cpp

clean:
	-rm -rf $(OBJECTS)

bigclean:
	-rm -rf $(OBJECTS) *.dat $(PROJECT)