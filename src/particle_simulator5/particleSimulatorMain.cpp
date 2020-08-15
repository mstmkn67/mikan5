//                  Masato Makino

#include "ParticleSimulator.h"
#include "../common/udf/gourmain.h"


void udfHeaderCheck()
{
	string version("V5.0"),engine("MIKAN");
	cout << "**************************************************************" << endl;
	cout << "             Particle Simulator   " << version << "           " << endl;
	cout << "                                        Masato Makino         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: particle_simulator5 -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	ParticleSimulator* sim=new ParticleSimulator(in,out);
	sim->run();
	delete sim;
	return 0;
}
