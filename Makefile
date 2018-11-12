# Makefile for compiling the HADESthreshscan

VERSION = v1
TRBNETDIR = /home/hadaq/git/trbnettools
ROOTDIR = /cvmfs/hades.gsi.de/install/root-5.34.34
BOOSTDIR = /usr

# Use this compiler
CC = g++

# Includes
TRBNETINCDIR = -I$(TRBNETDIR)/include
ROOTINCDIR = -I$(ROOTDIR)/include
BOOSTINCDIR = -I$(BOOSTDIR)/include

INCLUDEDIRS = $(TRBNETINCDIR) $(ROOTINCDIR) $(BOOSTINCDIR)

# Libraries
TRBNETLIBDIR = -L$(TRBNETDIR)/trbnetd
TRBNETLIB = -ltrbnet
ROOTLIBDIR = -L$(ROOTDIR)/lib
ROOTLIB = -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -lGui
BOOSTLIBDIR =  -L$(BOOSTDIR)/lib
BOOSTLIB =  -lboost_program_options -lboost_filesystem -lboost_system

LIBS = $(TRBNETLIB) $(ROOTLIB) $(BOOSTLIB)
LIBDIRS = $(TRBNETLIBDIR) $(ROOTLIBDIR) $(BOOSTLIBDIR)

# Options
OPT = -c -Wall -Wextra -pedantic -O3 -std=c++11

# Make rules
HADESthreshscan_$(VERSION): HADESthreshscan_$(VERSION).o
	$(CC) -o HADESthreshscan_$(VERSION) HADESthreshscan_$(VERSION).o $(LIBDIRS) $(LIBS)
	$(shell echo 'export LD_LIBRARY_PATH=$(ROOTDIR)/lib:$(TRBNETDIR)/lib:$(BOOSTDIR)/lib:"$$LD_LIBRARY_PATH"' > setLD)
clean:
	/bin/rm -f *.o HADESthreshscan_$(VERSION) setLD
.C.o:  $*.C
	$(CC) $*.C $(INCLUDEDIRS) $(OPT)
