#include "Contig.h"

CContig::CContig()
: m_Consensus(NULL)
, m_ConsensusQuality(NULL)
, m_Name(NULL)
, m_ConsensusLength(0)
{}

CContig::~CContig(void) {
	//printf("Contig destructor.\n");

	// delete the name
	if(m_Name) delete [] m_Name;

	// delete the consensus
	if(m_Consensus) delete [] m_Consensus;

	// delete the consensus quality
	if(m_ConsensusQuality) delete [] m_ConsensusQuality;
}
