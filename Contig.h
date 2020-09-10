#ifndef NULL
#define NULL 0
#endif

class CContig {
public:
	CContig();
	~CContig(void);
	// The contig name
	char* m_Name;
	// The consensus base sequence
	char* m_Consensus;
	// the consensus quality values
	char* m_ConsensusQuality;
	// The consensus length
	unsigned int m_ConsensusLength;
};
