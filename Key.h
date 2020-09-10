#ifndef KEY_H_SEEN
#define KEY_H_SEEN

#include <stdio.h>
#include <stdlib.h>

#define mix(a,b,c) { \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8);  \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12); \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5);  \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

class CKey {
public:
	CKey();
	CKey(double s, unsigned char v);
	CKey(unsigned int a, unsigned int c, unsigned int g, unsigned int t);
	~CKey(void);
	// The quantity of 'A' alleles in the key
	unsigned int m_a;
	// The quantity of 'C' alleles in the key
	unsigned int m_c;
	// The quantity of 'G' alleles in the key
	unsigned int m_g;
	// The quantity of 'T' alleles in the key
	unsigned int m_t;
	// the key variation derived from the counts
	unsigned char m_Variation;
	// the posterior sum represented by the key
	double m_Sum;
	// creates a new key by incrementing the given allele
	CKey MakeNewKey(unsigned char base);
	// calculates a new posterior sum
	void CalculateSum(unsigned char base, unsigned int depth, double priorVar, double newPriorVar, double oldSum, double Perror);
	// get the hashcode for the key
	unsigned int HashCode(void);
	// get the hashcode for a given key
	static unsigned int HashCode(unsigned int depth_a, unsigned int depth_c, unsigned int depth_g, unsigned int depth_t);
	// compares this key to another key k
	//bool IsEqual(CKey* k);
	// operator used in STL sort
	bool operator<(const CKey& key) const;
	// used to make life easier
	//void operator=(const CKey& key);
	// used to make life easier
	void operator=(const CKey* key);
private:
	// the saved hash
	unsigned int m_Hash;
};

#endif /* !KEY_H_SEEN */
