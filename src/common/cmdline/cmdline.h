//==========================================================================
//	cmdline.h
//
//	class CmdLine 
//
//		get input_udf_name, output_udf_name, def_udf_name and
//		control_file_name from command line args.
//
//	PASTA-2.6	Nov.05, 2001	J.Takimoto, DOI project
#ifndef _CMDLINE_H_
#define _CMDLINE_H_
#ifdef _WIN32
#pragma warning (disable : 4786)	// VC++ でデバッガーに関する不要の警告を押さえる
#endif

#include <string>
#include <strstream>
using namespace std;

class CmdLine {
public:				// data mambers are public
	string input_udf;
	string output_udf;
	string def_udf;		// output def udf
	string log_file;
	string summ_udf;	// summary udf
	string ctrl_file;	// control file

public:

	// the following two methods return 0 if success, non-zero otherwise.

	// get various file names from command line args
	int readArgs(int argc, const char* const argv[]);

	// check whether we have necessary info or not
	int check();

	// error messages
	string message() { return errmessage.str(); }

private:
	ostrstream	errmessage;
};

#endif // _CMDLINE_H_
