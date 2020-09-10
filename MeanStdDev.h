#include "RuntimeParameters.h"

class CMeanStdDev {
public:
	// specifies the mean for each dimension
	static double Mean[];
	// specifies the standard deviation for each dimension
	static double StdDev[];
	// loads the mean & standard deviations
	static void Load(void);
private:
	// converts a character array to a little-endian double
	static double GetDouble(char* buffer);
};
