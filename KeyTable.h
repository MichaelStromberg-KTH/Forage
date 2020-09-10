#include "Key.h"
#include <algorithm>
#include <cstring>
using namespace std;

#define CAPACITY 257

class CKeyTable {
public:
	CKeyTable(void);
	~CKeyTable(void);
	// size of the key hash table
	unsigned int m_Size;
	unsigned int m_Collisions;
	unsigned int m_Same;
	// key array used in hash table
	CKey m_Keys[CAPACITY];
	unsigned int m_KeyHashes[CAPACITY];
	// retrieves the key having the given hash
	CKey* Get(unsigned int hash);
	// returns true if the table contains the key
	bool Contains(unsigned int hash);
	// puts a key into the hash table
	void Put(CKey* key);
	// erases all keys from the hash table
	void Clear(void);
	// returns the keys in the hash table
	void GetKeys(CKey* keys);
	// Retains the 50 highest probabilities and adds them to the specified key table
	void SortAndAdd(CKeyTable* kt, unsigned int depth);
private:
	// returns the index given a key with a certain hash
	unsigned int IndexFor(unsigned int hash);
};
