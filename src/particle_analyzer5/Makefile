CC = g++ 
#CC = g++ -mno-cygwin -DMINGW -Wno-deprecated
#CFLAG = -c -O3 -I $(PF_FILES)/include
#CCFLAG =-lm -o particle_analyzer5.exe  -llapack -lblas -L$(PF_FILES)/lib/win64  -lplatform 
CFLAG = -c -O3 -I $(PF_FILES)/include 
CCFLAG = -o particle_analyzer5  -L $(PF_FILES)/lib/linux_64  -llapack -lblas -lplatform -fopenmp
#CCFLAG =-lm -o particle_analyzer5.exe  -L$(PF_FILES)/lib/win32  -lplatform -llapack_win -lblas_win -lg2c 

OBJECT=Vector3d.o EulerAngle.o ResistanceMobilityTensor.o AxisymmetricResistanceTensor.o \
	ParticleModel.o BeadsModel.o BoundaryElementModel.o \
	SuperpositionModel.o Timer.o  StokesGreenFunc.o BeadsInteraction.o \
	LapackFunctions.o gourmain.o cmdlineopts.o \
	ParticleAnalyzer.o particleAnalyzerMain.o

all :$(OBJECT)
	$(CC) $(OBJECT) $(CCFLAG)

Vector3d.o : ../common/Vector3d.cpp ../common/Vector3d.h
	$(CC) ../common/Vector3d.cpp $(CFLAG)

EulerAngle.o : ../common/EulerAngle.cpp ../common/EulerAngle.h
	$(CC) ../common/EulerAngle.cpp $(CFLAG)

ResistanceMobilityTensor.o : ../common/resMobTensor/ResistanceMobilityTensor.cpp ../common/resMobTensor/ResistanceMobilityTensor.h
	$(CC) ../common/resMobTensor/ResistanceMobilityTensor.cpp $(CFLAG)

AxisymmetricResistanceTensor.o : ../common/resMobTensor/AxisymmetricResistanceTensor.cpp ../common/resMobTensor/AxisymmetricResistanceTensor.h
	$(CC) ../common/resMobTensor/AxisymmetricResistanceTensor.cpp $(CFLAG)

ParticleModel.o : ParticleModel.cpp ParticleModel.h
	$(CC) ParticleModel.cpp $(CFLAG)

BoundaryElementModel.o : BoundaryElementModel.cpp BoundaryElementModel.h
	$(CC) BoundaryElementModel.cpp $(CFLAG)

BeadsModel.o : BeadsModel.cpp BeadsModel.h
	$(CC) BeadsModel.cpp $(CFLAG)

SuperpositionModel.o : SuperpositionModel.cpp SuperpositionModel.h
	$(CC) SuperpositionModel.cpp $(CFLAG)

Timer.o : ../common/Timer.cpp ../common/Timer.h
	$(CC) ../common/Timer.cpp $(CFLAG)

StokesGreenFunc.o : ../common/GreenFunc/StokesGreenFunc.cpp ../common/GreenFunc/StokesGreenFunc.h
	$(CC) ../common/GreenFunc/StokesGreenFunc.cpp $(CFLAG)

BeadsInteraction.o : ../common/GreenFunc/BeadsInteraction.cpp ../common/GreenFunc/BeadsInteraction.h
	$(CC) ../common/GreenFunc/BeadsInteraction.cpp $(CFLAG)

LapackFunctions.o : ../common/solver/LapackFunctions.cpp ../common/solver/LapackFunctions.h
	$(CC) ../common/solver/LapackFunctions.cpp $(CFLAG)

gourmain.o : ../common/udf/gourmain.cpp ../common/udf/gourmain.h
	$(CC) ../common/udf/gourmain.cpp $(CFLAG)

cmdlineopts.o : ../common/udf/cmdlineopts.cpp ../common/udf/cmdlineopts.h
	$(CC) ../common/udf/cmdlineopts.cpp $(CFLAG)

ParticleAnalyzer.o : ParticleAnalyzer.cpp ParticleAnalyzer.h
	$(CC) ParticleAnalyzer.cpp $(CFLAG)

particleAnalyzerMain.o : particleAnalyzerMain.cpp
	$(CC) particleAnalyzerMain.cpp $(CFLAG)

clean : 
	rm -f *.o particle_analyzer5

install :
	strip particle_analyzer5
	cp particle_analyzer5 ../../bin
