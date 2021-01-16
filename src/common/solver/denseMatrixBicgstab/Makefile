CC = g++
CFLAG = -c -O3
CCFLAG =-lm -o  
OBJECTS =denseMatrixBicgstab.o denseMatrixILUdcmp.o denseMatrixLinearEquationSolver.o \
denseMatrixOperation.o MatrixOperation.o 

all : $(OBJECTS)
	$(CC) check.cpp $(CCFLAG) check.out $(OBJECTS)
	$(CC) check2.cpp $(CCFLAG) check2.out $(OBJECTS)
	$(CC) check3.cpp $(CCFLAG) check3.out $(OBJECTS)

denseMatrixBicgstab.o : denseMatrixBicgstab.cpp denseMatrixBicgstab.h
	$(CC) denseMatrixBicgstab.cpp $(CFLAG)

denseMatrixILUdcmp.o : denseMatrixILUdcmp.cpp denseMatrixILUdcmp.h
	$(CC) denseMatrixILUdcmp.cpp $(CFLAG)

denseMatrixLinearEquationSolver.o : denseMatrixLinearEquationSolver.cpp denseMatrixLinearEquationSolver.h
	$(CC) denseMatrixLinearEquationSolver.cpp $(CFLAG)

denseMatrixOperation.o : denseMatrixOperation.cpp denseMatrixOperation.h
	$(CC) denseMatrixOperation.cpp $(CFLAG)

MatrixOperation.o : MatrixOperation.cpp MatrixOperation.h
	$(CC) MatrixOperation.cpp $(CFLAG)

clean : 
	rm *.o *.out
