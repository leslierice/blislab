CC        = icc
OPTS      = -std=c99 -O3 

#MKL_SERIAL= -lmkl_sequential -lmkl_core -lmkl_intel_lp64
#MKL_PARALLEL= -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64
#LIB       = ${MKL_SERIAL} -lpthread
#LIB       = ${MKL_PARALLEL} -lpthread -openmp
LIB       = -mkl=parallel -openmp

#LIBDIR   += -L${MKLROOT}/lib/intel64 -I${MKLROOT}/include

all: test_dgemm_mkl.x


%.x: %.o
	$(CC) $(OPTS) $(LIBDIR) $< -o $@ $(LIB)

#.SUFFIXES: .c .cpp .cu .o

%.o : %.c
	$(CC) $(OPTS) $(INC) -c $< $(LIB)

clean:
	rm -f *.x *~ *.o 

##.PHONY: clean

