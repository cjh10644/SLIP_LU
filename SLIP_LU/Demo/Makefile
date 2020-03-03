default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../Include -I../Source -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include

# LDFLAGS = -L../../lib

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd
CS = ../Lib/libslip.a $(LDLIBS)

all: lib example example2 example3 example4 example5 lu test
	- ./example
	- ./example2
	- ./example3
	- ./example4
	- ./example5
	- ./SLIPLU
	- ./test

runtest: lib test
	- ./test

runtest1: lib test1
	- ./test1

runlu: lib lu
	- ./SLIPLU

compare1: lib lu test test1
	- ./SLIPLU f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./test f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./test1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt

compare2: lib lu test test1
	- ./SLIPLU f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt
	- ./test f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt
	- ./test1 f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt

compare3: lib lu test test1
	- ./SLIPLU f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt
	- ./test f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt
	- ./test1 f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt

runlub: lib lu
	- ./SLIPLU f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./SLIPLU f ../ExampleMats/prob159_mat.txt ../ExampleMats/prob159_v.txt
	- ./SLIPLU f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/10teams.mat ../ExampleMats/BasisLIB_ALL/RHS/10teams.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.mat ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/co9.mat ../ExampleMats/BasisLIB_ALL/RHS/co9.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.mat ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/d18512.mat ../ExampleMats/BasisLIB_ALL/RHS/d18512.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/d15112.mat ../ExampleMats/BasisLIB_ALL/RHS/d15112.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.mat ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/dbic1.mat ../ExampleMats/BasisLIB_ALL/RHS/dbic1.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/de080285.mat ../ExampleMats/BasisLIB_ALL/RHS/de080285.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/fome12.mat ../ExampleMats/BasisLIB_ALL/RHS/fome12.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/model11.mat ../ExampleMats/BasisLIB_ALL/RHS/model11.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/model6.mat ../ExampleMats/BasisLIB_ALL/RHS/model6.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.mat ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.mat ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.v.rhs
	- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/rail4284.mat ../ExampleMats/BasisLIB_ALL/RHS/rail4284.v.rhs
	
runtestb: lib test
	- ./test f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./test f ../ExampleMats/prob159_mat.txt ../ExampleMats/prob159_v.txt
	- ./test f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/10teams.mat ../ExampleMats/BasisLIB_ALL/RHS/10teams.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.mat ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/co9.mat ../ExampleMats/BasisLIB_ALL/RHS/co9.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.mat ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/d18512.mat ../ExampleMats/BasisLIB_ALL/RHS/d18512.v.rhs
	- ./test  f ../ExampleMats/BasisLIB_ALL/RHS/d15112.mat ../ExampleMats/BasisLIB_ALL/RHS/d15112.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.mat ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/dbic1.mat ../ExampleMats/BasisLIB_ALL/RHS/dbic1.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/de080285.mat ../ExampleMats/BasisLIB_ALL/RHS/de080285.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/fome12.mat ../ExampleMats/BasisLIB_ALL/RHS/fome12.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/model11.mat ../ExampleMats/BasisLIB_ALL/RHS/model11.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/model6.mat ../ExampleMats/BasisLIB_ALL/RHS/model6.v.rhs
	- ./test f ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.mat ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.v.rhs
	- ./test  f ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.mat ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.v.rhs
	- ./test  f ../ExampleMats/BasisLIB_ALL/RHS/rail4284.mat ../ExampleMats/BasisLIB_ALL/RHS/rail4284.v.rhs
	- ./test  f ../ExampleMats/BasisLIB_ALL/RHS/rl5934.mat ../ExampleMats/BasisLIB_ALL/RHS/rl5934.v.rhs- ./SLIPLU f ../ExampleMats/BasisLIB_ALL/RHS/rl5934.mat ../ExampleMats/BasisLIB_ALL/RHS/rl5934.v.rhs

runtest1b: lib test1
	- ./test1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./test1 f ../ExampleMats/prob159_mat.txt ../ExampleMats/prob159_v.txt
	- ./test1 f ../ExampleMats/NSR8K_mat.txt ../ExampleMats/NSR8K_v.txt
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/10teams.mat ../ExampleMats/BasisLIB_ALL/RHS/10teams.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.mat ../ExampleMats/BasisLIB_ALL/RHS/30_70_4.5_0.5_100.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/co9.mat ../ExampleMats/BasisLIB_ALL/RHS/co9.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.mat ../ExampleMats/BasisLIB_ALL/RHS/cont1_l.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/d18512.mat ../ExampleMats/BasisLIB_ALL/RHS/d18512.v.rhs
	- ./test1  f ../ExampleMats/BasisLIB_ALL/RHS/d15112.mat ../ExampleMats/BasisLIB_ALL/RHS/d15112.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.mat ../ExampleMats/BasisLIB_ALL/RHS/d2q06c.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/dbic1.mat ../ExampleMats/BasisLIB_ALL/RHS/dbic1.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/de080285.mat ../ExampleMats/BasisLIB_ALL/RHS/de080285.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/fome12.mat ../ExampleMats/BasisLIB_ALL/RHS/fome12.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/model11.mat ../ExampleMats/BasisLIB_ALL/RHS/model11.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/model6.mat ../ExampleMats/BasisLIB_ALL/RHS/model6.v.rhs
	- ./test1 f ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.mat ../ExampleMats/BasisLIB_ALL/RHS/stormg2_1000.v.rhs
	- ./test1  f ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.mat ../ExampleMats/BasisLIB_ALL/RHS/scfxm1-2r-64.v.rhs
	- ./test1  f ../ExampleMats/BasisLIB_ALL/RHS/rail4284.mat ../ExampleMats/BasisLIB_ALL/RHS/rail4284.v.rhs
	- ./test1  f ../ExampleMats/BasisLIB_ALL/RHS/rl5934.mat ../ExampleMats/BasisLIB_ALL/RHS/rl5934.v.rhs


lib:
	( cd ../Lib ; $(MAKE) )

example: lib example.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example example.c $(CS)

example2: lib example2.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example2 example2.c demos.c $(CS)

example3: lib example3.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example3 example3.c  $(CS)

example4: lib example4.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example4 example4.c  $(CS)

example5: lib example5.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example5 example5.c demos.c $(CS)

lu: lib SLIPLU.c demos.h demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SLIPLU SLIPLU.c demos.c $(CS)

test: lib test.c demos.h demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o test test.c demos.c $(CS)

test1: lib test1.c demos.h demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o test1 test1.c demos.c $(CS)
clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example example2 example3 example4 example5 SLIPLU test *.a *.dSYM *.obj *.dll
