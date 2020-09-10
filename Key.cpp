#include "Key.h"

CKey::CKey()
: m_a(0)
, m_c(0)
, m_g(0)
, m_t(0)
, m_Variation(0)
, m_Sum(0)
, m_Hash(0)
{}

CKey::CKey(double s, unsigned char v)
: m_a(0)
, m_c(0)
, m_g(0)
, m_t(0)
, m_Variation(v)
, m_Sum(s)
, m_Hash(0)
{}

// TODO: make sure that in the Bayesian algorithm, that we never encounter a null variation
CKey::CKey(unsigned int a, unsigned int c, unsigned int g, unsigned int t)
: m_a(a)
, m_c(c)
, m_g(g)
, m_t(t)
, m_Variation(0)
, m_Sum(0)
, m_Hash(0)
{}

CKey::~CKey(void) {}

// creates a new key by incrementing the given allele
CKey CKey::MakeNewKey(unsigned char base) {
	CKey tempKey(m_a,m_c,m_g,m_t);

	switch(base) {
		case 0:
			tempKey.m_a++;
			break;
		case 1:
			tempKey.m_c++;
			break;
		case 2:
			tempKey.m_g++;
			break;
		case 3:
			tempKey.m_t++;
			break;
		default:
			printf("ERROR: invalid base found: %d\n",base);
			exit(1);
	}

	unsigned char position = 0;

	// m_Variation is reset in the constructor
	if(tempKey.m_a > 0) tempKey.m_Variation += 1;
	if(tempKey.m_c > 0) tempKey.m_Variation += 2;
	if(tempKey.m_g > 0) tempKey.m_Variation += 4;
	if(tempKey.m_t > 0) tempKey.m_Variation += 8;

	return tempKey;
}

// calculates a new posterior sum
void CKey::CalculateSum(unsigned char base, unsigned int depth, double priorVar, double newPriorVar, double oldSum, double Perror) {
	unsigned int count;

	switch(base) {
		case 0:
			count = m_a;
			break;
		case 1:
			count = m_c;
			break;
		case 2:
			count = m_g;
			break;
		case 3:
			count = m_t;
			break;
		default:
			printf("ERROR: invalid base found: %d\n",base);
			exit(1);
	}

	// calculate prior modification factor
	// - take into account change in variation
	// - take into account base count change
	depth++;
	m_Sum += oldSum * Perror * newPriorVar * count / (priorVar * depth);
}

// get the hashcode for the key
unsigned int CKey::HashCode(void) {

	unsigned int a = 0x9e3779b9 + m_a;
	unsigned int b = 0x9e3779b9 + m_c;
	unsigned int c = m_g;

	mix(a,b,c);

	a += m_t;
	b += m_a;
	c += m_c;

	mix(a,b,c);

	//printf("hash: %11u a: %3d c: %3d g: %3d t: %3d\n",c,m_a,m_c,m_g,m_t);

	return c;
}

// get the hashcode for the key
unsigned int CKey::HashCode(unsigned int depth_a, unsigned int depth_c, unsigned int depth_g, unsigned int depth_t) {

	unsigned int a = 0x9e3779b9 + depth_a;
	unsigned int b = 0x9e3779b9 + depth_c;
	unsigned int c = depth_g;

	mix(a,b,c);

	a += depth_t;
	b += depth_a;
	c += depth_c;

	mix(a,b,c);

	return c;
}

// compares this key to another key k
/*
bool CKey::IsEqual(CKey* k) {

	if(m_Sum != k->m_Sum) return false;

	if(m_a != k->m_a) return false;
	if(m_c != k->m_c) return false;
	if(m_g != k->m_g) return false;
	if(m_t != k->m_t) return false;

	return true;
}
*/
bool CKey::operator<(const CKey& key) const {
    return m_Sum > key.m_Sum;
}
/*
void CKey::operator=(const CKey& key) {
    
	m_a         = key.m_a;
	m_c         = key.m_c;
	m_g         = key.m_g;
	m_t         = key.m_t;
	m_Sum       = key.m_Sum;
	m_Variation = key.m_Variation;
}
*/
void CKey::operator=(const CKey* key) {
    
	m_a         = key->m_a;
	m_c         = key->m_c;
	m_g         = key->m_g;
	m_t         = key->m_t;
	m_Sum       = key->m_Sum;
	m_Variation = key->m_Variation;
}
