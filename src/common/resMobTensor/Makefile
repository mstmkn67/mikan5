CC = g++  
CFLAG = -c -O3
CCFLAG = -lm -llapack 
OBJECT1=Vector3d.o EulerAngle.o ResistanceMobilityTensor.o AxisymmetricResistanceTensor.o \
LapackFunctions.o
OBJECT2=MobilityTensorOfBrownianSimulation.o random.o rng.o twister.o 

all : $(OBJECT1) $(OBJECT2)
	$(CC) checkAxisymmetricTensor.cpp -o checkAxis.out $(OBJECT1) $(CCFLAG)
	$(CC) checkBrownianMotion.cpp -o checkBrown.out $(OBJECT1) $(OBJECT2) $(CCFLAG)
	$(CC) checkZeroEigenValue.cpp -o checkZero.out $(OBJECT1) $(CCFLAG)
	
Vector3d.o : ../Vector3d.cpp ../Vector3d.h
	$(CC) ../Vector3d.cpp $(CFLAG)

EulerAngle.o : ../EulerAngle.cpp ../EulerAngle.h
	$(CC) ../EulerAngle.cpp $(CFLAG)

AxisymmetricResistanceTensor.o : AxisymmetricResistanceTensor.cpp AxisymmetricResistanceTensor.h
	$(CC) AxisymmetricResistanceTensor.cpp $(CFLAG)

ResistanceMobilityTensor.o : ResistanceMobilityTensor.cpp ResistanceMobilityTensor.h
	$(CC) ResistanceMobilityTensor.cpp $(CFLAG)

MobilityTensorOfBrownianSimulation.o : MobilityTensorOfBrownianSimulation.cpp MobilityTensorOfBrownianSimulation.h
	$(CC) MobilityTensorOfBrownianSimulation.cpp $(CFLAG)

LapackFunctions.o : ../solver/LapackFunctions.cpp ../solver/LapackFunctions.h
	$(CC) ../solver/LapackFunctions.cpp $(CFLAG)

random.o : ../random/random.cpp ../random/random.h
	$(CC) ../random/random.cpp $(CFLAG)

rng.o : ../random/rng.cpp ../random/rng.h
	$(CC) ../random/rng.cpp $(CFLAG)

twister.o : ../random/twister.cpp ../random/twister.h
	$(CC) ../random/twister.cpp $(CFLAG)

clean : 
	rm *.o *.out *.dat
