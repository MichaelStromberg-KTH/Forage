#include <fstream>
#include <limits.h>
#include <float.h>
#include <cstring>
#include "RuntimeParameters.h"
using namespace std;

class CDualCodebook {
public:
	static void LoadCodebooks(void);
	static void Dispose(void);
	// checks the classification of the specified point and returns true if quorum is reached
	static bool IsSNP(double* points);
private:
	// OBD codebook data points
	static double** m_obd_points;
	// LVQ3 codebook data points
	static double** m_lvq3_points;
	// OBD vector classifications
	static unsigned char* m_obd_classifications;
	// LVQ3 vector classifications
	static unsigned char* m_lvq3_classifications;
	// number of entries in the OBD codebook
	static unsigned int m_obd_numEntries;
	// number of entries in the LVQ3 codebook
	static unsigned int m_lvq3_numEntries;
	// descriptors for each class present in the codebooks
	static const char* m_classDescriptors;
	// identification strings for each OBD vector
	static char** m_obd_names;
	// identification strings for each LVQ3 vector
	static char** m_lvq3_names;
	// the number of dimensions present in each codebook vector
	static unsigned char m_numDimensions;
	// loads the OBD codebook from file
	static void LoadObdCodebook(const char* filename);
	// loads the LVQ3 codebook from file
	static void LoadLvq3Codebook(const char* filename);
	// converts a character array to a little-endian double
	static double GetDouble(char* buffer);
	// converts a character array to a little-endian unsigned integer
	static unsigned int GetUnsignedInt(char* buffer);
	// converts a UTF string to an ANSI string
	static void GetAnsiString(char* buffer, char* c, int len);
};
