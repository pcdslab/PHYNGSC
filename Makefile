CCMPI = mpicxx
CCFLAGS = -O3 -m64  -Wall
CCFLAGS += -fopenmp -std=c++11
PROG = phyNGSC
OBJS = bit_stream.o huffman.o tasks.o

.cpp.o:
	$(CCMPI) $(CCFLAGS) -c $< -o $@

$(PROG): phyNGSC.o $(OBJS)
	$(CCMPI) $(CCFLAGS) -o $(PROG) phyNGSC.o $(OBJS)

clean:
	-rm *.o
	-rm $(PROG)