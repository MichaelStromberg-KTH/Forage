#include "DualCodebook.h"

// declare the static variables
unsigned char CDualCodebook::m_numDimensions = 4;
const char* CDualCodebook::m_classDescriptors = "ACGT";

double** CDualCodebook::m_obd_points  = NULL;
double** CDualCodebook::m_lvq3_points = NULL;

unsigned char* CDualCodebook::m_obd_classifications  = NULL;
unsigned char* CDualCodebook::m_lvq3_classifications = NULL;

unsigned int CDualCodebook::m_obd_numEntries  = 0;
unsigned int CDualCodebook::m_lvq3_numEntries = 0;

char** CDualCodebook::m_obd_names  = NULL;
char** CDualCodebook::m_lvq3_names = NULL;

void CDualCodebook::LoadCodebooks(void) {	

	// load the OBD codebook from file
	LoadObdCodebook("OBD_CODEBOOK.DAT");

	// load the LVQ3 codebook from file
	LoadLvq3Codebook("LVQ3_CODEBOOK.DAT");
}

void CDualCodebook::Dispose(void) {

	// clean up the names and the points (2D)
	for(unsigned int i=0;i<m_obd_numEntries;i++) {
		delete [] m_obd_names[i];	
		delete [] m_obd_points[i];

	}

	for(unsigned int i=0;i<m_lvq3_numEntries;i++) {
		delete [] m_lvq3_names[i];	
		delete [] m_lvq3_points[i];
	}

	// clean up the names and the points (1D)
	delete [] m_obd_names;
	delete [] m_obd_points;
	delete [] m_lvq3_names;
	delete [] m_lvq3_points;

	// clean up the classifications
	delete [] m_obd_classifications;
	delete [] m_lvq3_classifications;
}

// loads the OBD codebook from file
void CDualCodebook::LoadObdCodebook(const char* filename) {

	ifstream in(filename, ios::in | ios::binary);

	// if our file is missing, bomb
	if(in.fail()) {
		printf("ERROR: Could not open the codebook.\n");
		exit(1);
	}

	// retrieve the number of vector entries
	char bigEndianUnsignedInt[4];
	in.read(bigEndianUnsignedInt,sizeof(unsigned int));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
	unsigned int entries;
	memcpy(&entries,bigEndianUnsignedInt,sizeof(unsigned int));
#else
	unsigned int entries  = GetUnsignedInt(bigEndianUnsignedInt);
#endif
	m_obd_numEntries = entries;

	// retrieve the number of vector dimensions
	unsigned char dim;
	in.read((char*)&dim,sizeof(dim));

	// initialize our data points
	m_obd_points = new double*[entries];
	for(unsigned int i=0;i<entries;i++) m_obd_points[i] = new double[dim];

	// initialize the classifications
	m_obd_classifications = new unsigned char[entries];

	// initialize the names
	m_obd_names = new char*[entries];

	// read in the data points
	char bigEndianDouble[8];

	for(unsigned int i=0;i<entries;i++) {

		// read in the data points
		for(int j=0;j<dim;j++) {
			in.read(bigEndianDouble,sizeof(double));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
			memcpy(&m_obd_points[i][j],bigEndianDouble,sizeof(double));
#else
			m_obd_points[i][j] = GetDouble(bigEndianDouble);
#endif
		}

		// read in the classification
		in.read((char*)&m_obd_classifications[i],sizeof(char));
		
		// read in name length
		char nameLength;
		in.read((char*)&nameLength,sizeof(nameLength));

		// read in the string
		char* utfString = new char[nameLength*2];
		m_obd_names[i] = new char[nameLength+1];
		in.read(utfString,nameLength*2);
		GetAnsiString(utfString,m_obd_names[i],nameLength);

		// clean up temporary variables
		delete [] utfString;
	}

	in.close();
}

// loads the LVQ3 codebook from file
void CDualCodebook::LoadLvq3Codebook(const char* filename) {

	ifstream in(filename, ios::in | ios::binary);

	// if our file is missing, bomb
	if(in.fail()) {
		printf("ERROR: Could not open the codebook.\n");
		exit(1);
	}

	// retrieve the number of vector entries
	char bigEndianUnsignedInt[4];
	in.read(bigEndianUnsignedInt,sizeof(unsigned int));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
	unsigned int entries;
	memcpy(&entries,bigEndianUnsignedInt,sizeof(unsigned int));
#else
	unsigned int entries  = GetUnsignedInt(bigEndianUnsignedInt);
#endif
	m_lvq3_numEntries = entries;

	// retrieve the number of vector dimensions
	unsigned char dim;
	in.read((char*)&dim,sizeof(dim));

	// initialize our data points
	m_lvq3_points = new double*[entries];
	for(unsigned int i=0;i<entries;i++) m_lvq3_points[i] = new double[dim];

	// initialize the classifications
	m_lvq3_classifications = new unsigned char[entries];

	// initialize the names
	m_lvq3_names = new char*[entries];

	// read in the data points
	char bigEndianDouble[8];

	for(unsigned int i=0;i<entries;i++) {

		// read in the data points
		for(int j=0;j<dim;j++) {
			in.read(bigEndianDouble,sizeof(double));
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
			memcpy(&m_lvq3_points[i][j],bigEndianDouble,sizeof(double));
#else
			m_lvq3_points[i][j] = GetDouble(bigEndianDouble);
#endif
		}

		// read in the classification
		in.read((char*)&m_lvq3_classifications[i],sizeof(char));
	
		// read in name length
		char nameLength;
		in.read((char*)&nameLength,sizeof(nameLength));

		// read in the string
		char* utfString = new char[nameLength*2];
		m_lvq3_names[i] = new char[nameLength+1];
		in.read(utfString,nameLength*2);
		GetAnsiString(utfString,m_lvq3_names[i],nameLength);

		// clean up temporary variables
		delete [] utfString;
	}

	in.close();
}

// converts a character array to a little-endian double
double CDualCodebook::GetDouble(char* buffer) {

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

// converts a character array to a little-endian unsigned integer
unsigned int CDualCodebook::GetUnsignedInt(char* buffer) {

	char c[4];
	
	c[0] = buffer[3];
	c[1] = buffer[2];
	c[2] = buffer[1];
	c[3] = buffer[0];
	
	return *(unsigned int*)c;
}

// converts a UTF string to an ANSI string
void CDualCodebook::GetAnsiString(char* buffer, char* c, int len) {
	for(int i=0;i<len;i++) c[i] = buffer[2*i+1];
	c[len] = 0x0;
}

// checks the classification of the specified point and returns true if quorum is reached
bool CDualCodebook::IsSNP(double* points) {
	
	unsigned char obdWinnerClassification = UCHAR_MAX;
	double obdWinnerDifference            = DBL_MAX;

	// Go through all OBD code vectors
	for(unsigned int i=0;i<m_obd_numEntries;i++) {

		double difference = 0.0;
		double tempDiff   = 0.0;

		for(unsigned int j=0;j<m_numDimensions;j++) {
			tempDiff = m_obd_points[i][j] - points[j];
			difference += tempDiff * tempDiff;
		}
		
		if(difference < obdWinnerDifference) {
			obdWinnerClassification = m_obd_classifications[i];
			obdWinnerDifference     = difference;
		}
	}

	unsigned char lvq3WinnerClassification = UCHAR_MAX;
	double lvq3WinnerDifference            = DBL_MAX;

	// Go through all LVQ3 code vectors
	for(unsigned int i=0;i<m_lvq3_numEntries;i++) {

		double difference = 0.0;
		double tempDiff   = 0.0;

		for(unsigned int j=0;j<m_numDimensions;j++) {
			tempDiff = m_lvq3_points[i][j] - points[j];
			difference += tempDiff * tempDiff;
		}
		
		if(difference < lvq3WinnerDifference) {
			lvq3WinnerClassification = m_lvq3_classifications[i];
			lvq3WinnerDifference     = difference;
		}
	}

	// if both networks classify this as a SNP, return true
	if((lvq3WinnerClassification == 1) && (obdWinnerClassification == 1)) return true;

	return false;
}
