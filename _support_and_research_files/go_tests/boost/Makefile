CXX=g++
CXXFLAGS=-fPIC -I ~/git/boost.math/
LDFLAGS=-static-libgcc -static-libstdc++ -shared

LIB=boost_hyp1f1

OBJS=instantiated.o

$(LIB): $(OBJS)
	$(CXX) $(LDFLAGS) -o lib$(LIB).so $?
