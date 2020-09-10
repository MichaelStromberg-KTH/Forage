#include "Sequence.h"

CSequence::CSequence()
: m_Name(NULL)
, m_IsReverseComplement(false)
, m_Offset(0)
, m_Filename(NULL)
, m_Bases(NULL)
, m_Qualities(NULL)
, m_IsAligned(false)
, m_ClipStart(0)
, m_ClipEnd(0)
, m_NumBases(0)
{}

CSequence::~CSequence(void) {
	//printf("Sequence destructor: %s.\n",m_Name);

	// delete the name
	if(m_Name) delete [] m_Name;

	// delete the filename
	if(m_Filename) delete [] m_Filename;

	// delete the bases
	if(m_Bases) delete [] m_Bases;

	// delete the qualities
	if(m_Qualities) delete [] m_Qualities;
}

// Aligns the sequence (old version) with respect to the contig
void CSequence::AlignOldVersion(unsigned int alignBases) {

	// TODO: make sure the qualities have been read
	if(m_IsAligned) {
		printf("ERROR: Sequence is already aligned.\n");
		exit(1);
	}

	//printf("Offset: %d\n",m_Offset);

	// adjust the offset
	int offset = m_Offset;

	// the consensus sequence starts at position 1
	unsigned int consensusLength = alignBases + 1;

	// how many bases should we copy from the original
	unsigned int numBasesToCopy = m_NumBases;

	int startPosition    = 0;
	int newStartPosition = offset;
	
	int endPosition      = numBasesToCopy;
	int newEndPosition   = numBasesToCopy + offset;

	// if sequence starts before the consensus, adjust the beginning of the sequence
	if(newStartPosition < 0) {
		newStartPosition = 1;
		newEndPosition  += offset;
		numBasesToCopy   = newEndPosition - newStartPosition;
	}

	// if sequence is longer than consensus, adjust the end of the sequence
	if(newEndPosition > (signed)consensusLength) {
		newEndPosition = consensusLength;
		numBasesToCopy   = newEndPosition - newStartPosition;
	}

	// align the bases
	char* pStartPos    = m_Bases + startPosition;
	char* pNewStartPos = m_Bases + newStartPosition;

	memmove(pNewStartPos,pStartPos,numBasesToCopy);

	// align the qualities
	pStartPos    = m_Qualities + startPosition;
	pNewStartPos = m_Qualities + newStartPosition;
	
	memmove(pNewStartPos,pStartPos,numBasesToCopy);

	unsigned int numBaseTrimBegin = newStartPosition;
	unsigned int numBaseTrimEnd   = consensusLength - newEndPosition;

	if(numBaseTrimBegin > 0) {
		memset(m_Bases,0,sizeof(char)*numBaseTrimBegin);
		memset(m_Qualities,0,sizeof(char)*numBaseTrimBegin);
	}
	
	if(numBaseTrimEnd > 0) {
		memset(m_Bases+newEndPosition,0,sizeof(char)*numBaseTrimEnd);
		memset(m_Qualities+newEndPosition,0,sizeof(char)*numBaseTrimEnd);
	}

	// perform base filtration
	if(CRuntimeParameters::EnableBaseFiltration) {
		unsigned int threshold = CRuntimeParameters::BaseFilterThreshold;
		for(unsigned int i=0;i<alignBases;i++) 
			if(m_Qualities[i] < (signed)threshold) {
				m_Bases[i]     = 0;
				m_Qualities[i] = 0;
			}
	}

	// set the aligned flag
	m_IsAligned = true;

	// set the new length
	m_NumBases = alignBases;
}

// Aligns the sequence with respect to the contig
void CSequence::Align(unsigned int alignBases) {

	// the consensus sequence starts at position 1
	unsigned int consensusLength = alignBases + 1;

	// how many bases should we copy from the original
	unsigned int numBasesToCopy = m_NumBases;

	// all contigs begin at position 1
	int contigStartPos = 1;

	// this is where the sequence starts on the contig scale
	contigStartPos += m_Offset;

	//printf("Contig start position: %u\n",contigStartPos);

	// this is where the sequence starts on the sequence scale
	int sequenceStartPos = 1;

	// if the starting position is off the scale, adjust
	if(contigStartPos < 1) {
		contigStartPos    = 1;
		sequenceStartPos -= m_Offset;
	}

	unsigned int sequencePosDiff = numBasesToCopy - sequenceStartPos;
	unsigned int contigPosDiff   = consensusLength - contigStartPos;

	if(sequencePosDiff < contigPosDiff) numBasesToCopy = sequencePosDiff;
		else numBasesToCopy = contigPosDiff;

	//printf("consensus length: %u, sequence length: %u, bases being copied: %u.\n",alignBases,m_NumBases,numBasesToCopy);

	//for(unsigned int i=0;i<consensusLength;i++) printf("%c",m_Bases[i]);
	//printf("\n\n");

	// align the bases
	char* pSeqStartPos    = m_Bases + sequenceStartPos;
	char* pContigStartPos = m_Bases + contigStartPos;
	
	//printf("* first base: %c, second base: %c, third base: %c\n",m_Bases[0],m_Bases[1],m_Bases[2]);
	memmove(pContigStartPos,pSeqStartPos,numBasesToCopy);

	// align the qualities
	pSeqStartPos    = m_Qualities + sequenceStartPos;
	pContigStartPos = m_Qualities + contigStartPos;
	
	memmove(pContigStartPos,pSeqStartPos,numBasesToCopy);

	//for(unsigned int i=0;i<consensusLength;i++) printf("%c",m_Bases[i]);
	//printf("\n\n");

	unsigned int endPos = contigStartPos+numBasesToCopy;

	//printf("start: %u, end: %u, consensus length: %u\n",contigStartPos,contigStartPos+numBasesToCopy,consensusLength);

	unsigned int numBaseTrimBegin = contigStartPos;
	unsigned int numBaseTrimEnd   = consensusLength - endPos;

	if(numBaseTrimBegin > 0) {
		memset(m_Bases,0,sizeof(char)*numBaseTrimBegin);
		memset(m_Qualities,0,sizeof(char)*numBaseTrimBegin);
	}
	
	if(numBaseTrimEnd > 0) {	
		memset(m_Bases+endPos,0,sizeof(char)*numBaseTrimEnd);
		memset(m_Qualities+endPos,0,sizeof(char)*numBaseTrimEnd);
	}

	// perform base filtration
	if(CRuntimeParameters::EnableBaseFiltration) {
		unsigned int threshold = CRuntimeParameters::BaseFilterThreshold;
		for(unsigned int i=0;i<alignBases;i++) 
			if(m_Qualities[i] < (signed)threshold) {
				m_Bases[i]     = 0;
				m_Qualities[i] = 0;
			}
	}

	// set the aligned flag
	m_IsAligned = true;

	// set the new length
	m_NumBases = alignBases;
}
