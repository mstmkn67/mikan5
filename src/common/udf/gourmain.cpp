#include "gourmain.h"
#include <omp.h>

int main(int argc,char* argv[]){
	CmdLineOpts opts(argc,argv,"I:O:");
	string in_udf_file_name=opts.arg('I');
	string out_udf_file_name=opts.arg('O');
	string num_thread=opts.arg('N');
	if(in_udf_file_name.size()==0){
		error_massage();
		return 1;
	}
	if(out_udf_file_name.size()==0)out_udf_file_name=in_udf_file_name;
	if(num_thread.size()==0){
		omp_set_num_threads(1);
	}else{
		int n=atoi(num_thread.c_str());
		omp_set_num_threads(n);
	}
	try {
		UDFManager* inudf=new UDFManager(in_udf_file_name);
		UDFManager* outudf=new UDFManager(out_udf_file_name,in_udf_file_name);
		gourmain(inudf,outudf);
		delete outudf;
		delete inudf;
	}
	catch(UDFStreamException &e) {
		cerr << e.GetErrorMessage() << endl;
		return 4;
	}
	catch (UDFObjectException &e) {
		cerr << e.GetErrorMessage() << endl;
		return 5;
	}
	catch (Location::LocationException &e) {
		cerr << e.what() << endl;
		return 6;
	}
	catch (UDFManager::UDFManagerException &e) {
		cerr << e.what() << endl;
		return 7;
	}
	catch (PFException &e) {
		cerr << e.GetErrorMessage() << endl;
		return 3;
	}
	catch (...) {
		cerr << "Undefined exception" << endl;
		return 99;
	}

	return 0;
}
