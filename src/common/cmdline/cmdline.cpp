//===========================================================================
//
//	cmdline.cpp
//
// PASTA-2.6	Nov.05, 2001	J.Takimoto, Doi-project
//===========================================================================

#include "cmdline.h"
#include "cmdlineopts.h"

//===== readArgs() ==============================
//	get various file names from comand line args
//
//	returns 0 if success, non-zero otherwise.

int
CmdLine::readArgs(int argc, const char* const argv[])
{
	CmdLineOpts opts(argc,argv,"I:O:D:L:S:M:");
	
	if( ! opts.good() ) {
		errmessage << opts.message() << endl;
		return 1;
	}

	input_udf  = opts.arg('I');
	output_udf = opts.arg('O');
	def_udf    = opts.arg('D');
	log_file   = opts.arg('L');
	summ_udf   = opts.arg('S');
	ctrl_file  = opts.arg('M');

	return 0;
}

//===== check() ==========================
//	check whether necessary info are there
//

int
CmdLine::check()
{
	int errcount = 0;

	if( input_udf == "" ) {
		errmessage << "Input UDF not specified." << endl;
		++errcount;
	}
	if( output_udf == "" ) {
		errmessage << "Output UDF not specified." << endl;
		++errcount;
	}
#if 0
	if( def_udf == "" ) {
		errmessage << "UDF def_file not specified." << endl;
		++errcount;
	}
#endif
	return errcount;
}
