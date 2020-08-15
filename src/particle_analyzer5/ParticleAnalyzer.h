//ParticleAnalyzer
//ó±éqâêÕÉNÉâÉX

#ifndef _PARTICLE_ANALYZER_H_
#define _PARTICLE_ANALYZER_H_

#include "udfmanager.h"
#include "BeadsModel.h"
#include "BoundaryElementModel.h"
#include "SuperpositionModel.h"

#include <string>
#include <vector>
using namespace std;

class ParticleAnalyzer{
public:
	ParticleAnalyzer(UDFManager* in_udf,UDFManager* out_udf);
	virtual ~ParticleAnalyzer();
	virtual void update();
protected:
	virtual void analyze_superposition();
	virtual void analyze_BEM();
	virtual void analyze_beads();

	virtual void output_common(double time,Vector3d& translate,Tensor3x3& rotate,ResistanceTensor& rt);
private:
	UDFManager* in_udf;UDFManager* out_udf;
	Location loc_type;
};

#endif // _PARTICLE_ANALYZER_H_
