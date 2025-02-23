
#######################################################################
#                                                                     #
#             Makefile for ORTHPOL                                    #
#                                                                     #
#  Usage:                                                             #
#                                                                     #
#  Type  make                  to compile and run all test programs   #
#  Type  make testN.out.local  to compile and run the test program    #
#                              number N                               #
#  Type  make clean            to remove all but the generated output #
#                                                                     #
#######################################################################

#----------------------------------------------------------------------
# This section specifies the name of the linker and compilation options
# for various systems. All but one are commented out.
#----------------------------------------------------------------------

FC = gfortran
FCFLAGS = -Wall -fbacktrace #Do not use -O2; it may introduce errors
FLFLAGS =

#-------------------------------------------------------------------
# This section specifies which object files are needed for each test
# program.
#-------------------------------------------------------------------

TEST1 = test1.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o
TEST2 = test2.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o
TEST3 = test3.o r1mach.o sti.o lancz.o
TEST4 = test4.o r1mach.o mcdis.o qgp.o sti.o lancz.o recur.o gauss.o
TEST5 = test5.o r1mach.o d1mach.o recur.o drecur.o mcdis.o qgp.o sti.o lancz.o gauss.o
TEST6 = test6.o r1mach.o d1mach.o mcdis.o dmcdis.o qgp.o dqgp.o sti.o dsti.o lancz.o dlancz.o recur.o drecur.o gauss.o dgauss.o
TEST7 = test7.o r1mach.o d1mach.o mcdis.o dmcdis.o qgp.o dqgp.o sti.o dsti.o lancz.o dlancz.o
TEST8 = test8.o r1mach.o recur.o mccheb.o qgp.o cheb.o
TEST9 = test9.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o chri.o dchri.o
TEST10= test10.o r1mach.o d1mach.o recur.o drecur.o gauss.o dgauss.o chri.o dchri.o
TEST11= test11.o r1mach.o d1mach.o recur.o drecur.o knum.o nu0jac.o dknum.o gchri.o dgchri.o cheb.o dcheb.o chri.o dchri.o

#--------------------------------------------------------------
# This is the default target. It will cause all the tests to be
# generated and run, saving their outputs in files.
#--------------------------------------------------------------

.phony : all
all :   test1.out.local test2.out.local test3.out.local test4.out.local test5.out.local test6.out.local test7.out.local test8.out.local test9.out.local test10.out.local test11.out.local

#----------------------------------------------------------------------
# This section shows how to run each test program. The output for testN
# is put in testN.out.local.
#----------------------------------------------------------------------

test1.out.local  : test1 
	./test1  > $@
test2.out.local  : test2 
	./test2  > $@
test3.out.local  : test3 
	./test3  > $@
test4.out.local  : test4 
	./test4  > $@
test5.out.local  : test5
	./test5  > $@
test6.out.local  : test6
	./test6  > $@
test7.out.local  : test7
	./test7  > $@
test8.out.local  : test8
	./test8  > $@
test9.out.local  : test9
	./test9  > $@
test10.out.local : test10
	./test10 > $@
test11.out.local : test11
	./test11 > $@

#--------------------------------------------------
# This section shows how to link each test program.
#--------------------------------------------------

test1 : $(TEST1)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test2 : $(TEST2)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test3 : $(TEST3)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test4 : $(TEST4)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test5 : $(TEST5)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test6 : $(TEST6)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test7 : $(TEST7)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test8 : $(TEST8)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test9 : $(TEST9)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test10: $(TEST10)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
test11: $(TEST11)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#----------------------------------------------
# This section shows how to clean up afterward.
#----------------------------------------------

.phony : clean
clean :
	rm -f test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 *.o

