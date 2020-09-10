#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "RuntimeParameters.h"

class CSequence {
public:
	CSequence();
	~CSequence(void);
	// the sequence name
	char* m_Name;
	// if true, the sequence needs to be complemented
	bool m_IsReverseComplement;
	// the sequence offset in reference to the consensus
	int m_Offset;
	// the sequence filename
	char* m_Filename;
	// the bases in the sequence
	char* m_Bases;
	// the sequence qualities
	char* m_Qualities;
	// True when the sequence has been aligned with the contig
	bool m_IsAligned;
	// position where sequence clip starts
	int m_ClipStart;
	// position where sequence clip ends
	int m_ClipEnd;
	// number of bases in the sequence
	unsigned int m_NumBases;
	// Aligns the sequence with respect to the contig
	void Align(unsigned int alignBases);
	// Aligns the sequence (old version) with respect to the contig
	void AlignOldVersion(unsigned int alignBases);
};
