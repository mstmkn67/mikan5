CC = g++ 
CFLAG = -c -O3 -I $(PF_FILES)/include 
CCFLAG = -o particle_simulator5 -L $(PF_FILES)/lib/linux_64  -llapack -lblas -lplatform -fopenmp
OBJECT=Vector3d.o Vector2d.o EulerAngle.o Timer.o \
	AprxDoubleSquareLatticeField2d.o AprxDoubleSquareLatticeField3d.o \
	AprxVector2dSquareLatticeField2d.o AprxVector3dSquareLatticeField3d.o \
	gourmain.o cmdlineopts.o random.o rng.o twister.o \
	ResistanceMobilityTensor.o MobilityTensorOfBrownianSimulation.o \
	LapackFunctions.o UdfFlow.o ConstantForce.o Particle.o ParticleSimulator.o \
	EllipticRotationField.o SawtoothField.o DynamicFlow.o \
	ElectricMagneticTorque.o\
	particleSimulatorMain.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Vector3d.o : ../common/Vector3d.cpp ../common/Vector3d.h
	$(CC) ../common/Vector3d.cpp $(CFLAG)

Vector2d.o : ../common/Vector2d.cpp ../common/Vector2d.h
	$(CC) ../common/Vector2d.cpp $(CFLAG)

EulerAngle.o : ../common/EulerAngle.cpp ../common/EulerAngle.h
	$(CC) ../common/EulerAngle.cpp $(CFLAG)

Timer.o : ../common/Timer.cpp ../common/Timer.h
	$(CC) ../common/Timer.cpp $(CFLAG)

AprxDoubleSquareLatticeField2d.o : ../common/aprxField/AprxDoubleSquareLatticeField2d.cpp  ../common/aprxField/AprxDoubleSquareLatticeField2d.h
	$(CC)  ../common/aprxField/AprxDoubleSquareLatticeField2d.cpp $(CFLAG)

AprxDoubleSquareLatticeField3d.o : ../common/aprxField/AprxDoubleSquareLatticeField3d.cpp  ../common/aprxField/AprxDoubleSquareLatticeField3d.h
	$(CC)  ../common/aprxField/AprxDoubleSquareLatticeField3d.cpp $(CFLAG)

AprxVector2dSquareLatticeField2d.o : ../common/aprxField/AprxVector2dSquareLatticeField2d.cpp  ../common/aprxField/AprxVector2dSquareLatticeField2d.h
	$(CC)  ../common/aprxField/AprxVector2dSquareLatticeField2d.cpp $(CFLAG)

AprxVector3dSquareLatticeField3d.o : ../common/aprxField/AprxVector3dSquareLatticeField3d.cpp  ../common/aprxField/AprxVector3dSquareLatticeField3d.h
	$(CC)  ../common/aprxField/AprxVector3dSquareLatticeField3d.cpp $(CFLAG)

gourmain.o : ../common/udf/gourmain.cpp ../common/udf/gourmain.h
	$(CC) ../common/udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : ../common/udf/cmdlineopts.cpp ../common/udf/cmdlineopts.h
	$(CC) ../common/udf/cmdlineopts.cpp $(CFLAG)

random.o : ../common/random/random.cpp ../common/random/random.h
	$(CC) ../common/random/random.cpp $(CFLAG)

rng.o : ../common/random/rng.cpp ../common/random/rng.h
	$(CC) ../common/random/rng.cpp $(CFLAG)

twister.o : ../common/random/twister.cpp ../common/random/twister.h
	$(CC) ../common/random/twister.cpp $(CFLAG)

ResistanceMobilityTensor.o : ../common/resMobTensor/ResistanceMobilityTensor.cpp ../common/resMobTensor/ResistanceMobilityTensor.h
	$(CC) ../common/resMobTensor/ResistanceMobilityTensor.cpp $(CFLAG)

MobilityTensorOfBrownianSimulation.o : ../common/resMobTensor/MobilityTensorOfBrownianSimulation.cpp ../common/resMobTensor/MobilityTensorOfBrownianSimulation.h
	$(CC) ../common/resMobTensor/MobilityTensorOfBrownianSimulation.cpp $(CFLAG)

LapackFunctions.o : ../common/solver/LapackFunctions.cpp ../common/solver/LapackFunctions.h
	$(CC) ../common/solver/LapackFunctions.cpp $(CFLAG)

UdfFlow.o : ../common/velocityField/UdfFlow.cpp ../common/velocityField/UdfFlow.h
	$(CC) ../common/velocityField/UdfFlow.cpp $(CFLAG)

DynamicFlow.o : ../common/velocityField/DynamicFlow.cpp ../common/velocityField/DynamicFlow.h
	$(CC) ../common/velocityField/DynamicFlow.cpp $(CFLAG)

ConstantForce.o : ConstantForce.cpp ConstantForce.h
	$(CC) ConstantForce.cpp $(CFLAG)

particleSimulatorMain.o : particleSimulatorMain.cpp
	$(CC) particleSimulatorMain.cpp $(CFLAG)

Particle.o : Particle.cpp Particle.h
	$(CC) Particle.cpp $(CFLAG)

ParticleSimulator.o : ParticleSimulator.cpp ParticleSimulator.h
	$(CC) ParticleSimulator.cpp $(CFLAG)

EllipticRotationField.o : ../common/electricMagneticField/EllipticRotationField.h  ../common/electricMagneticField/EllipticRotationField.cpp
	$(CC)  ../common/electricMagneticField/EllipticRotationField.cpp $(CFLAG)

SawtoothField.o : ../common/electricMagneticField/SawtoothField.h  ../common/electricMagneticField/SawtoothField.cpp
	$(CC)  ../common/electricMagneticField/SawtoothField.cpp $(CFLAG)

ElectricMagneticTorque.o : ElectricMagneticTorque.h ElectricMagneticTorque.cpp
	$(CC) ElectricMagneticTorque.cpp $(CFLAG)

clean : 
	rm -f *.o particle_simulator5

install :
	strip particle_simulator5
	cp particle_simulator5 ../../bin
