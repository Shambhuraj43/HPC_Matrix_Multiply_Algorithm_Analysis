## Compiler, tools and options
CC = gcc
LINK = gcc
OPT = -O1

CCFLAGS = $(OPT) -Wall
LDFLAGS = $(OPT)

## PAPI setup
PAPI_ROOT = /opt/cray/pe/papi/6.0.0.7
#PAPI_LIB = $(PAPI_ROOT)/lib/libpapi.a

## Libraries
LIBS = -L $(PAPI_ROOT)/lib -lpapi
INC = -I $(PAPI_ROOT)/include

## Files
OBJECTS = main.o
TARGET = main

## Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $<

## Build Rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core

