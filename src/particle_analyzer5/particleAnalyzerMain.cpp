
#include "ParticleAnalyzer.h"
#include "../common/udf/gourmain.h"
#include <iostream>
#include <fstream>
using namespace std;

void header(){
	string version("V5.0");
	cout << "****************************************************************" << endl;
	cout << "            Particle Analyzer "     << version << "             " << endl;
	cout << "                                        Masato MAKINO           " << endl;
	cout << "****************************************************************" << endl;
}

void error_massage(){
	cout << "usage: particle_analyzer5 -I input_udf [-O output_udf] " << endl;
}

int gourmain(UDFManager* in,UDFManager* out)
{
	header();
	ParticleAnalyzer* pa=new ParticleAnalyzer(in,out);
	pa->update();
	delete pa;
	return 0;
}

