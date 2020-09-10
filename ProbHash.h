#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class CProbHash
{
public:
	CProbHash(void);
	~CProbHash(void);
	// Places the given value into the hash
	void Add(unsigned char variation, double value);
	// Retrieves the given value from the hash
	double Get(unsigned char variation);
	// returns the hash value of the string
	double m_hash[16];
};
