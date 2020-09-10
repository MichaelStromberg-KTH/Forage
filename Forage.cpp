// Forage.cpp : Defines the entry point for the console application.
//

#include "Ace.h"

void ShowHelp(char* execFilename);

int main(int argc, char* argv[]) {

	char aceDir[FILENAME_BUFFER_SIZE+1];
	char aceFile[FILENAME_BUFFER_SIZE+1];

	aceDir[0]  = 0;
	aceFile[0] = 0;

	printf("Forage: SNP Discovery Software (c) 2004 Michael Stromberg\n=========================================================\n");
	
	if(argc == 1) ShowHelp(argv[0]);

	// evaluate all of the switches
	for(int i=1;i<argc;i++) {

		// help
		if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"-?") || !strcmp(argv[i],"\\h") || !strcmp(argv[i],"\\?")) ShowHelp(argv[0]);

		// verbose
		if(!strcmp(argv[i],"-v")) CRuntimeParameters::IsVerbose = true;
		
		// define directory
		if(!strcmp(argv[i],"-dir")) {	
			
			int arglen = (int)strlen(argv[i+1]);
			if(arglen > FILENAME_BUFFER_SIZE) {
				printf("ERROR: Maximum phrap cluster pathname is %u characters. Current length: %u\n",FILENAME_BUFFER_SIZE,arglen);
				exit(1);
			}

			memcpy(aceDir,argv[i+1],arglen);
			
			if(aceDir[arglen-1] == '\\') aceDir[arglen-1] = 0;
				else aceDir[arglen] = 0;
			
			i++;
		}

		// define file
		if(!strcmp(argv[i],"-ace")) {

			int arglen = (int)strlen(argv[i+1]);
			if(arglen > FILENAME_BUFFER_SIZE) {
				printf("ERROR: Maximum ace pathname is %u characters. Current length: %u\n",FILENAME_BUFFER_SIZE,arglen);
				exit(1);
			}

			memcpy(aceFile,argv[i+1],arglen);
			aceFile[arglen] = 0;
			i++;
		}

		// define base filtration threshold
		if(!strcmp(argv[i],"-basefilterthreshold")) {
					
			CRuntimeParameters::EnableBaseFiltration = true;
			int threshold = atoi(argv[i+1]);
			i++;

			if((threshold < 1) || (threshold > 255)) {
				printf("ERROR: Invalid base filter threshold selected [1-255]. Default: 20\n");
				exit(1);
			}

			CRuntimeParameters::BaseFilterThreshold = threshold;
		}

		// define minimum sequence threshold
		if(!strcmp(argv[i],"-minsequencethreshold")) {
					
			int threshold = atoi(argv[i+1]);
			i++;

			if(threshold < 2) {
				printf("ERROR: Invalid minimum sequence threshold selected [2+]. Default: 2\n");
				exit(1);
			}

			CRuntimeParameters::MinSequenceThreshold = threshold;
		}

		// enable XML output
		if(!strcmp(argv[i],"-xml")) {
			
			CRuntimeParameters::EnableXmlOutput = true;

			int arglen = (int)strlen(argv[i+1]);
			if(arglen > FILENAME_BUFFER_SIZE) {
				printf("ERROR: Maximum XML pathname is %u characters. Current length: %u\n",FILENAME_BUFFER_SIZE,arglen);
				exit(1);
			}

			memcpy(CRuntimeParameters::XmlFilename,argv[i+1],arglen);
			i++;
		}

		// specify the polymorphism rate
		if(!strcmp(argv[i],"-ppoly")) {
			
			double ppoly = atof(argv[i+1]);

			if((ppoly <= 0) || (ppoly > 0.05)) {
				printf("ERROR: Invalid rate of polymorphism specified [0-0.05]. Default: 0.003 (1 in 333 bases)\n");
				exit(1);
			}

			CRuntimeParameters::P_POLY = ppoly;
			i++;
		}
	}

	//
	// test to see if the ace directory and filename are correctly specified
	//

	if(strlen(aceDir) == 0) {
		printf("ERROR: Phrap cluster directory not specified. Specify with the -dir (directory) flag.\n");
		exit(1);
	}

	if(strlen(aceFile) == 0) {
		printf("ERROR: ace filename not specified. Specify with the -ace (ace filename) flag.\n");
		exit(1);
	}

	// try to read the acefile
	char buffer[522];
	sprintf(buffer,"%s%cedit_dir%c%s",aceDir,OS_SLASH,OS_SLASH,aceFile);

	ifstream ace(buffer, ios::in | ios::binary);

	// if our file is missing, bomb
	if(ace.fail()) {
		printf("ERROR: Could not open the ace file: %s\n",buffer);
		exit(1);
	}

	ace.close();

	//
	// execute Forage
	//
	CBayesianUtils::PreCalculate(CRuntimeParameters::P_POLY);
	CMeanStdDev::Load();
	CDualCodebook::LoadCodebooks();

	CAce cgap(aceDir,aceFile);
	cgap.LoadPhdFiles(aceDir);
	cgap.AlignSequences();
	cgap.FindParalogs();
	cgap.FindSNPs();

	return 0;
}

void ShowHelp(char* execFilename) {

	printf("Usage: %s -dir (phrap cluster directory) -ace (ace filename)\n\t\t  [-xml (xml filename)] [-ppoly (#)] [-basefilterthreshold (#)]\n\t\t  [-minsequencethreshold (#)] [-verbose] [-h]\n\n",execFilename);

	printf("Required parameters\n===================\n");
	printf("   -dir  specifies the phrap cluster directory. e.g. the directory that\n\t contains the edit_dir and phd_dir subdirectories.\n\n");
	printf("   -ace  specifies the ace filename. e.g. CLONE.fasta.screen.ace\n\n");
	
	printf("Optional parameters\n===================\n");
	printf("   -basefilterthreshold  enables base filtration which filters out all bases\n\t\t\t having a base quality less than the specified amount.\n\t\t\t Default: 20. Valid: [1-255]\n\n");
	printf("   -minsequencethreshold specifies the minimum number of sequences featured in\n\t\t\t an alignment slice. Default: 2. Valid [2+]\n\n");
	printf("   -ppoly  specifies a new a priori polymorphism rate. Default: 0.003\n\t   (1 in 333 bases). Valid (0.05,1]\n\n");
	printf("   -xml    saves results to the specified filename in Forage's XML format\n\n");
	printf("   -v      enables verbose output\n\n");
	exit(1);
}
