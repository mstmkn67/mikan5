CC = g++
CFLAG = -c -O3
CCFLAG =-lm -o 
OBJECT=Vector2d.o Vector3d.o
OBJECTS2D=AprxDoubleSquareLatticeField2d.o AprxVector2dSquareLatticeField2d.o
OBJECTS3D=AprxDoubleSquareLatticeField3d.o AprxVector3dSquareLatticeField3d.o

all : $(OBJECT) $(OBJECTS2D) $(OBJECTS3D) 
	$(CC) check2d.cpp $(CCFLAG) check2d.out $(OBJECT) $(OBJECTS2D)
	$(CC) check3d.cpp $(CCFLAG) check3d.out $(OBJECT) $(OBJECTS3D)

Vector2d.o : ../Vector2d.cpp ../Vector2d.h
	$(CC) ../Vector2d.cpp $(CFLAG)

Vector3d.o : ../Vector3d.cpp ../Vector3d.h
	$(CC) ../Vector3d.cpp $(CFLAG)

AprxDoubleSquareLatticeField2d.o : AprxDoubleSquareLatticeField2d.cpp AprxDoubleSquareLatticeField2d.h
	$(CC) AprxDoubleSquareLatticeField2d.cpp $(CFLAG)

AprxVector2dSquareLatticeField2d.o : AprxVector2dSquareLatticeField2d.cpp AprxVector2dSquareLatticeField2d.h
	$(CC) AprxVector2dSquareLatticeField2d.cpp $(CFLAG)

AprxDoubleSquareLatticeField3d.o : AprxDoubleSquareLatticeField3d.cpp AprxDoubleSquareLatticeField3d.h
	$(CC) AprxDoubleSquareLatticeField3d.cpp $(CFLAG)

AprxVector3dSquareLatticeField3d.o : AprxVector3dSquareLatticeField3d.cpp AprxVector3dSquareLatticeField3d.h
	$(CC) AprxVector3dSquareLatticeField3d.cpp $(CFLAG)

clean : 
	rm *.o *.out
