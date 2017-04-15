.PHONY: all clean

CCMPI = mpicxx
CCFLAGS = -O3 -m64  -Wall
CCFLAGS += -fopenmp -std=c++11
C_PROG = phyNGSC
D_PROG = phyNGSD
OBJSC = bit_stream.o huffman.o tasks.o
OBJSD = bit_stream.o huffman.o tasks.o

.cpp.o:
	$(CCMPI) $(CCFLAGS) -c $< -o $@

$(C_PROG): phyNGSC.o $(OBJSC)
	$(CCMPI) $(CCFLAGS) -o $(C_PROG) phyNGSC.o $(OBJSC)

$(D_PROG): phyNGSD.o $(OBJSD)
	$(CCMPI) $(CCFLAGS) -o $(D_PROG) phyNGSD.o $(OBJSD)

all: $(C_PROG) $(D_PROG)

clean:
	-rm *.o
	-rm $(C_PROG) phyNGSD 