#include "KeyTable.h"

CKeyTable::CKeyTable(void)
: m_Size(0)
, m_Collisions(0)
, m_Same(0)
{
	memset(m_KeyHashes,0,sizeof(unsigned int)*CAPACITY);
}

CKeyTable::~CKeyTable(void) {}

// returns the index given a key with a certain hash
unsigned int CKeyTable::IndexFor(unsigned int hash) {

	//printf("Hash: %d\n",hash);
	return hash % CAPACITY;
	
	//hash += ~(hash << 9);
	//hash ^=  (hash >> 14);
	//hash +=  (hash << 4);
	//hash ^=  (hash >> 10);

	//return hash & (CAPACITY - 1);
}

// retrieves the key having the given hash
CKey* CKeyTable::Get(unsigned int hash) {

	// if we have no elements, return null
	if(m_Size == 0) return NULL;

	int i = IndexFor(hash);
	
	while(m_KeyHashes[i] > 0) {
		if(hash == m_KeyHashes[i]) return &m_Keys[i];
		if(++i == CAPACITY) i = 0;
	}
	
	return NULL;
}

// returns true if the table contains the key
bool CKeyTable::Contains(unsigned int hash) {

	// if we have no elements, return false
	if(m_Size == 0) return false;

	int i = IndexFor(hash);
	
	while(m_KeyHashes[i] > 0) {
		if(hash == m_KeyHashes[i]) return true;
		if(++i == CAPACITY) i = 0;
	}
	
	return false;
}

// puts a key into the hash table
void CKeyTable::Put(CKey* key) {

	// check how many entries we already have
	if(m_Size == CAPACITY) {
		printf("ERROR: Attempted to add more keys to key table than capacity allows.\n");
		exit(1);
	}

	// save the key's hash code and find the correlated hash table index
	unsigned int hash = key->HashCode();
	int i = IndexFor(hash);

	// if the index is already occupied, find an empty pocket or update old entry
	while(m_KeyHashes[i] > 0) {
		m_Collisions++;

		// check to see if this is the same key
		if(m_KeyHashes[i] == hash) {

			// since hash codes are a,c,g,t-based, only the sum and variation need to be updated
			m_Keys[i].m_Sum       = key->m_Sum;
			m_Keys[i].m_Variation = key->m_Variation;
			m_Same++;

			return;
		}

		// wrap around if we've reached capacity
		if(++i == CAPACITY) i = 0;
	}

	// found an empty pocket, update the hash array
	m_KeyHashes[i] = hash;

	// clone the key
	m_Keys[i].m_a         = key->m_a;
	m_Keys[i].m_c         = key->m_c;
	m_Keys[i].m_g         = key->m_g;
	m_Keys[i].m_t         = key->m_t;
	m_Keys[i].m_Sum       = key->m_Sum;
	m_Keys[i].m_Variation = key->m_Variation;
	
	// increment the hash table size
	m_Size++;

	return;
}

// erases all keys from the hash table
void CKeyTable::Clear(void) {
	
	// to clear the hash table, all we have to do is reset 
	// the hash array and the counter variables
	memset(m_KeyHashes,0,sizeof(unsigned int)*CAPACITY);
	m_Collisions = 0;
	m_Size = 0;
	m_Same = 0;
}

// returns the keys in the hash table
void CKeyTable::GetKeys(CKey* keys) {

	// if we have no elements, return
	if(m_Size == 0) return;

	for(unsigned int i=0,numEntries=0;i<CAPACITY;i++) 
		if(m_KeyHashes[i] > 0) {
			keys[numEntries].m_a         = m_Keys[i].m_a;
			keys[numEntries].m_c         = m_Keys[i].m_c;
			keys[numEntries].m_g         = m_Keys[i].m_g;
			keys[numEntries].m_t         = m_Keys[i].m_t;
			keys[numEntries].m_Sum       = m_Keys[i].m_Sum;
			keys[numEntries++].m_Variation = m_Keys[i].m_Variation;
		}
}

// Retains the 50 highest probabilities and adds them to the specified key table
void CKeyTable::SortAndAdd(CKeyTable* kt, unsigned int depth) {

	CKey* newKeys              = kt->m_Keys;
	unsigned int* newKeyHashes = kt->m_KeyHashes;

	unsigned int newLength = CAPACITY;
	unsigned int newSize   = kt->m_Size;

	// if the source Key Table is empty, all we have to do is clear this one.
	if(newSize == 0) {
		Clear();
		return;
	}

	// copy the valid (non-null) keys into a new array so that we can sort them
	CKey* newValidKeys = new CKey[newSize];

	unsigned int numEntries = 0;
	for(unsigned int i=0;i<CAPACITY;i++)
		if(newKeyHashes[i] > 0) newValidKeys[numEntries++] = newKeys[i];

	sort(newValidKeys,newValidKeys+numEntries);

	// erase the contents from the host Key Table
	Clear();

	// select the keys with the highest 50 probabilities
	numEntries = newSize;
	if(numEntries > 50) numEntries = 50;

	for(unsigned int i=0,movedEntries=0;i<newSize;i++) {
		Put(&newValidKeys[i]);
		if(++movedEntries == numEntries) break;
	}

	// dispose of the new valid entries
	delete [] newValidKeys;

	//
	// keep non-polymorphic combinations
	//

	unsigned int depthInt = depth + 1;
	unsigned int hashCode;

	//printf("depthInt: %u.\n",depthInt);

	// A monomorphic case
	hashCode = CKey::HashCode(depthInt,0,0,0); 
	if(!Get(hashCode)) Put(kt->Get(hashCode));

	// C monomorphic case
	hashCode = CKey::HashCode(0,depthInt,0,0); 
	if(!Get(hashCode)) Put(kt->Get(hashCode));

	// G monomorphic case
	hashCode = CKey::HashCode(0,0,depthInt,0); 
	if(!Get(hashCode)) Put(kt->Get(hashCode));

	// T monomorphic case
	hashCode = CKey::HashCode(0,0,0,depthInt); 
	if(!Get(hashCode)) Put(kt->Get(hashCode));

	// erase the contents from the source Key Table
	kt->Clear();
}
