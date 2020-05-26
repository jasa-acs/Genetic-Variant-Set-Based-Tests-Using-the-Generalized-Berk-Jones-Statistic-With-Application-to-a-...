CXX=g++ -std=c++11
INCLUDES = -I/PATH/TO/BOOST
CXXFLAGS = -g -Wall $(INCLUDES)

all: ebb_crossprob_cor

# 
# It doesn't automatically know the name of the binary
#

ebb_crossprob_cor: ebb_crossprob_cor.o 
	$(CXX) ebb_crossprob_cor.o -o ebb_crossprob_cor 

#
# It automatically knows to make the .o from the .c, that's why we only need $^
#

ebb_crossprob_cor.o: ebb_crossprob_cor.cpp 
	$(CXX) $(CXXFLAGS) -c $^ 

clean:
	rm -f *.o ebb_crossprob_cor

