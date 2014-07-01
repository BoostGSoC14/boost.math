CC=gcc
FC=gfortran

LIB=cephes_scipy

CFLAGS=-fPIC -I/usr/include/python2.7
LDFLAGS=-static-libgfortran -static-libgcc -static-libstdc++ -shared

C_FILES=specfun_wrappers.o sf_error.o
F_FILES=specfun.o

$(LIB) : $(C_FILES) $(F_FILES)
	$(FC) $(LDFLAGS) -o lib$(LIB).so $?
