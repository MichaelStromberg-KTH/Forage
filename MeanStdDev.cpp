#include "MeanStdDev.h"
#include <fstream>
#include <cstring>

using namespace std;

// declare the static variables
double CMeanStdDev::Mean[]   = { 0,0,0,0,0 };
double CMeanStdDev::StdDev[] = { 0,0,0,0,0 };

// loads the mean & standard deviations
void CMeanStdDev::Load(void) {

	ifstream in("MEANSTDDEV.DAT", ios::in | ios::binary);

	// if our file is missing, bomb
	if(in.fail()) {
		printf("ERROR: Could not open the mean and standard deviation file.\n");
		exit(1);
	}

	// retrieve the number of vector dimensions
	unsigned char dim;
	in.read((char*)&dim,sizeof(dim));

	// read in the data points
	char bigEndianDouble[8];

	// read in the mean data points
	for(int j=0;j<dim;j++) {
		in.read(bigEndianDouble,sizeof(double));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
                memcpy(&CMeanStdDev::Mean[j],bigEndianDouble,sizeof(double));
#else
		CMeanStdDev::Mean[j] = GetDouble(bigEndianDouble);
#endif
	}

	// read in the standard deviation data points
	for(int j=0;j<dim;j++) {
		in.read(bigEndianDouble,sizeof(double));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
		memcpy(&CMeanStdDev::StdDev[j],bigEndianDouble,sizeof(double));	
#else
		CMeanStdDev::StdDev[j] = GetDouble(bigEndianDouble);
#endif		
	}

	in.close();
}

// converts a character array to a little-endian double
double CMeanStdDev::GetDouble(char* buffer) {

	char c[8];
	
	c[0] = buffer[7];
	c[1] = buffer[6];
	c[2] = buffer[5];
	c[3] = buffer[4];
	c[4] = buffer[3];
	c[5] = buffer[2];
	c[6] = buffer[1];
	c[7] = buffer[0];

	return *(double*)c;
}
