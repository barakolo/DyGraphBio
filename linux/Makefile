DEPS = node.hpp node.cpp edge.hpp edge.cpp parser.hpp parser.cpp dygraph.hpp dygraph.cpp propgraph.hpp propgraph.cpp
CFLAGS = -I Eigen
LDFLAGS = 
MFLAGS = -std=c++11 -s -O3 -DNDEBUG
CC = g++ 

clean:
	rm -f test tests main

tests: clean
	$(CC) $(MFLAGS) $(CFLAGS) $(DEPS) ./tests.cpp -o tests $(LDFLAGS)

