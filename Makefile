
##########
## Set Variables
##########

CXX= g++
CXXFLAGS = -fPIC -D NDEBUG
#CXXFLAGS = -fsched-spec-load -fPIC -D NDEBUG
#For debugging mode use the following instead:
#CXXFLAGS = -g -fsched-spec-load -fPIC -D NDEBUG

ifeq (${LIBFLAGS},-shared)
	LIBRARY = libRNA.so
	HYBRID_LIBRARY = libHybridRNA.so
else
	LIBRARY = RNA_Library.dll
	HYBRID_LIBRARY = HybridRNA_Library.dll
endif



LINK = ${CXX} ${CXXFLAGS} -Wall -o $@
#For debugging mode use the following instead:
#LINK=${CXX} ${CXXFLAGS} -g -Wall -o $@


all:
	@echo 'Use `make MonteCarlo` to MonteCarlo'


COMMON_OBJ_FILES= MonteCarlo.o MonteCarloClass.o  
RNAstructureZero_FILES= ../RNAstructureZero/RNA_class/RNA.o ../RNAstructureZero/src/algorithm.o ../RNAstructureZero/src/alltrace.o ../RNAstructureZero/src/arrayclass.o ../RNAstructureZero/src/dotarray.o ../RNAstructureZero/src/draw.o ../RNAstructureZero/src/forceclass.o ../RNAstructureZero/src/MaxExpect.o ../RNAstructureZero/src/MaxExpectStack.o ../RNAstructureZero/src/outputconstraints.o ../RNAstructureZero/src/pfunction.o ../RNAstructureZero/src/probknot.o ../RNAstructureZero/src/random.o ../RNAstructureZero/src/rna_library.o ../RNAstructureZero/src/stackclass.o ../RNAstructureZero/src/stackstruct.o ../RNAstructureZero/src/stochastic.o ../RNAstructureZero/src/structure.o ../RNAstructureZero/RNA_class/thermodynamics.o ../RNAstructureZero/src/TProgressDialog.o

MonteCarlo: ${RNAstructureZero_FILES} ${COMMON_OBJ_FILES} 
	${LINK} ${RNAstructureZero_FILES} ${COMMON_OBJ_FILES} 

# Sequential Object Files (SEQ_OBJ_FILES)
MonteCarlo.o: MonteCarlo.cpp MonteCarloClass.h  ../RNAstructureZero/RNA_class/RNA.h ../RNAstructureZero/src/algorithm.h ../RNAstructureZero/src/alltrace.h ../RNAstructureZero/src/arrayclass.h ../RNAstructureZero/src/dotarray.h ../RNAstructureZero/src/draw.h ../RNAstructureZero/src/forceclass.h ../RNAstructureZero/src/MaxExpect.h ../RNAstructureZero/src/MaxExpectStack.h ../RNAstructureZero/src/outputconstraints.h ../RNAstructureZero/src/pfunction.h ../RNAstructureZero/src/probknot.h ../RNAstructureZero/src/random.h ../RNAstructureZero/src/rna_library.h ../RNAstructureZero/src/stackclass.h ../RNAstructureZero/src/stackstruct.h ../RNAstructureZero/src/stochastic.h ../RNAstructureZero/src/structure.h ../RNAstructureZero/RNA_class/thermodynamics.h ../RNAstructureZero/src/TProgressDialog.h

# Common Object Files (COMMON_OBJ_FILES)
MonteCarloClass.o: MonteCarloClass.cpp MonteCarloClass.h

#RNAstructureZero files (RNAstructureZero_FILES)
-g ../RNAstructureZero/RNA_class/RNA.o: ../RNAstructureZero/RNA_class/RNA.cpp ../RNAstructureZero/RNA_class/RNA.h
../RNAstructureZero/src/algorithm.o: ../RNAstructureZero/src/algorithm.cpp ../RNAstructureZero/src/algorithm.h
../RNAstructureZero/src/alltrace.o: ../RNAstructureZero/src/alltrace.cpp ../RNAstructureZero/src/alltrace.h
../RNAstructureZero/src/arrayclass.o: ../RNAstructureZero/src/arrayclass.cpp ../RNAstructureZero/src/arrayclass.h
../RNAstructureZero/src/dotarray.o: ../RNAstructureZero/src/dotarray.cpp ../RNAstructureZero/src/dotarray.h
../RNAstructureZero/src/draw.o: ../RNAstructureZero/src/draw.cpp ../RNAstructureZero/src/draw.h
../RNAstructureZero/src/forceclass.o: ../RNAstructureZero/src/forceclass.cpp ../RNAstructureZero/src/forceclass.h
../RNAstructureZero/src/MaxExpect.o: ../RNAstructureZero/src/MaxExpect.cpp ../RNAstructureZero/src/MaxExpect.h
../RNAstructureZero/src/MaxExpectStack.o: ../RNAstructureZero/src/MaxExpectStack.cpp ../RNAstructureZero/src/MaxExpectStack.h
../RNAstructureZero/src/outputconstraints.o: ../RNAstructureZero/src/outputconstraints.cpp ../RNAstructureZero/src/outputconstraints.h
../RNAstructureZero/src/pfunction.o: ../RNAstructureZero/src/pfunction.cpp ../RNAstructureZero/src/pfunction.h
../RNAstructureZero/src/probknot.o: ../RNAstructureZero/src/probknot.cpp ../RNAstructureZero/src/probknot.h
../RNAstructureZero/src/random.o: ../RNAstructureZero/src/random.cpp ../RNAstructureZero/src/random.h
../RNAstructureZero/src/rna_library.o: ../RNAstructureZero/src/rna_library.cpp ../RNAstructureZero/src/rna_library.h
../RNAstructureZero/src/stackclass.o: ../RNAstructureZero/src/stackclass.cpp ../RNAstructureZero/src/stackclass.h
../RNAstructureZero/src/stackstruct.o: ../RNAstructureZero/src/stackstruct.cpp ../RNAstructureZero/src/stackstruct.h
../RNAstructureZero/src/stochastic.o: ../RNAstructureZero/src/stochastic.cpp ../RNAstructureZero/src/stochastic.h
../RNAstructureZero/src/structure.o: ../RNAstructureZero/src/structure.cpp ../RNAstructureZero/src/structure.h
../RNAstructureZero/RNA_class/thermodynamics.o: ../RNAstructureZero/RNA_class/thermodynamics.cpp ../RNAstructureZero/RNA_class/thermodynamics.h
../RNAstructureZero/src/TProgressDialog.o: ../RNAstructureZero/src/TProgressDialog.cpp ../RNAstructureZero/src/TProgressDialog.h


##########
## Object and Executable Cleanup
##########

# Remove object files
clean:
	rm -f ../RNAstructureZero/RNA_class/*.o
	rm -f ../RNAstructureZero/src/*.o
	rm -f *.o

# Remove object files and executables
realclean: clean
	rm MonteCarlo
