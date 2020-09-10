#ifndef ACE_H_SEEN
#define ACE_H_SEEN

#include "Contig.h"
#include "Sequence.h"
#include "BayesianAlgorithm.h"
#include "SynchronizedArrays.h"
//#include "RuntimeParameters.h"
#include "MeanStdDev.h"
#include "DualCodebook.h"
#include <fstream>
#include <math.h>
#include <cstring>
using namespace std;

#define BUFFER_SIZE		2048
#define SEQ_BUFFER_SIZE	16384
#define FILENAME_SIZE	512
#define FLANK_RADIUS	32

class CAce {
public:
	CAce(char* dirPrefix, char* ace);
	~CAce(void);
	// Loads the specified contig from the ace file
	void BuildContigs(char* dirPrefix, char* ace);
	// Loads the specified contig from the old version ace file
	void BuildOldAceContigs(char* dirPrefix, char* ace);
	// Loads the sequences found in the ace file
	void LoadPhdFiles(char* dirPrefix);
	// Reverses a sequence (for reverse complements)
	void ReverseSequence(char* str, unsigned int length);
	// Aligns all of the sequences in the ace file
	void AlignSequences(void);
	// Dumps all of the sequences in the ace file
	void DumpSequences(void);
	// Dumps all of the sequences in the ace file at the specified position
	void DumpPosition(unsigned int contigPos);
	// Scours the cluster for SNPs
	void FindSNPs(void);
	// Scours the cluster for paralogs
	void FindParalogs(void);
	// the contigs read from the ace file
	CContig* m_Contigs;
	// temporary container for sequences read in ace file
	CSequence** m_ContigSequences;
	// masking vector for the temporary contig sequence container
	bool** m_ContigSeqsMask;
	// the number of contigs present in the ace file
	unsigned int m_NumContigs;
	// the number of sequences in each contig
	unsigned int* m_NumContigSeqs;
	// the number of sequences in each contig that are masked
	unsigned int* m_NumContigSeqsMasked;
	// is true if an old version ace file is encountered
	bool m_IsOldVersion;
	// is true if the cluster is aligned
	bool m_IsAligned;
	// Retrieves a IUPAC ambiguity code based on the variation string
	char GetIupacAmbiguityCode(unsigned char variation);
};

#endif /* !ACE_H_SEEN */
