#include "Ace.h"
#include <algorithm>

CAce::CAce(char* dirPrefix, char* ace)
: m_Contigs(NULL)
, m_ContigSequences(NULL)
, m_ContigSeqsMask(NULL)
, m_NumContigs(0)
, m_NumContigSeqs(NULL)
, m_NumContigSeqsMasked(NULL)
, m_IsOldVersion(false)
, m_IsAligned(false)
{
	BuildContigs(dirPrefix,ace);
}

CAce::~CAce(void) {
	//printf("Ace destructor.\n");

	// delete our contigs
	if(m_Contigs) delete [] m_Contigs;

	// delete our temporary contig sequence container and masking vector
	for(unsigned int i=0;i<m_NumContigs;i++) {
		if(m_ContigSequences[i]) delete [] m_ContigSequences[i];
		if(m_ContigSeqsMask[i])  delete [] m_ContigSeqsMask[i];
	}
	
	if(m_ContigSequences) delete [] m_ContigSequences;
	if(m_ContigSeqsMask)  delete [] m_ContigSeqsMask;
	
	// delete the contig sequences count vector
	if(m_NumContigSeqs) delete [] m_NumContigSeqs;

	// delete the contig masked sequences count vector
	if(m_NumContigSeqsMasked) delete [] m_NumContigSeqsMasked;
}

// Loads the specified contig from the ace file
void CAce::BuildContigs(char* dirPrefix, char* ace) {

	char buffer[BUFFER_SIZE];
	char sequenceBuffer[SEQ_BUFFER_SIZE];
	char absoluteFilename[FILENAME_SIZE];

	char* consensusPosition;
	char* sequencePosition;
	char* token;

	unsigned int positionsRead;

	//printf("Loading %s.\n",ace);

	sprintf(absoluteFilename,"%s%cedit_dir%c%s",dirPrefix,OS_SLASH,OS_SLASH,ace);

	//printf("Absolute filename: %s\n",absoluteFilename);

	ifstream in(absoluteFilename);

	// if our file is missing, bomb
	if(in.fail()) {
		printf("ERROR: Could not open the ace file: %s\n",absoluteFilename);
		exit(1);
	}

	// find the first non-empty line
	do {
		in.getline(buffer,BUFFER_SIZE);
	} while(strlen(buffer) == 0);

	token = strtok(buffer," ");

	if(strncmp(token,"AS",2) == 0) {

		// extract the # of contigs
		token = strtok(NULL," ");
		m_NumContigs = atoi(token);

		//printf("contigs: %u\n",m_NumContigs);

		// initialize main variables
		m_Contigs = new CContig[m_NumContigs];

		m_NumContigSeqs       = new unsigned int[m_NumContigs];
		m_NumContigSeqsMasked = new unsigned int[m_NumContigs];
		m_ContigSequences     = new CSequence*[m_NumContigs];
		m_ContigSeqsMask      = new bool*[m_NumContigs];

		for(unsigned int i=0;i<m_NumContigs;i++) {
			m_ContigSequences[i]     = NULL;
			m_ContigSeqsMask[i]      = NULL;
			m_NumContigSeqsMasked[i] = 0;
		}

		// investigate all of our contigs
		for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {

			// find the next contig tag
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strncmp(buffer,"CO",2) != 0);

			// extract the contig name
			token = strtok(buffer," ");
			token = strtok(NULL," ");

			unsigned int nameLen = (unsigned int)strlen(token);
			m_Contigs[currentContig].m_Name = new char[nameLen+1];
			strncpy(m_Contigs[currentContig].m_Name,token,nameLen);
			m_Contigs[currentContig].m_Name[nameLen] = 0;
			//printf("Contig name: %s.\n",m_Contigs[currentContig].m_Name);

			// extract the consensus length
			token = strtok(NULL," ");
			m_Contigs[currentContig].m_ConsensusLength = atoi(token);
			//printf("Consensus length: %d.\n",m_Contigs[currentContig].m_ConsensusLength);

			// extract the # of contig sequences
			token = strtok(NULL," ");
			m_NumContigSeqs[currentContig] = atoi(token);
			//printf("Contig sequences: %u.\n",m_NumContigSeqs[currentContig]);

			// read in the entire consensus sequence
			sequencePosition = sequenceBuffer;

			do {
				in.getline(sequencePosition,SEQ_BUFFER_SIZE);
				positionsRead = (unsigned int)strlen(sequencePosition);
				sequencePosition += positionsRead;
			} while(positionsRead != 0);

			// replace all asterisks with hyphens and capitalize
			for(sequencePosition = sequenceBuffer;*sequencePosition;sequencePosition++) {
				if((*sequencePosition > 96) && (*sequencePosition < 123)) *sequencePosition -= 32;
				if(*sequencePosition == 42) *sequencePosition = 45;
			}

			// TODO: is this right?
			m_Contigs[currentContig].m_Consensus        = new char[m_Contigs[currentContig].m_ConsensusLength+2];
			m_Contigs[currentContig].m_ConsensusQuality = new char[m_Contigs[currentContig].m_ConsensusLength+2];

			// TODO: do we really need to clear the entire base and quality sequence
			memset(m_Contigs[currentContig].m_Consensus,0,m_Contigs[currentContig].m_ConsensusLength+2);
			memset(m_Contigs[currentContig].m_ConsensusQuality,0,m_Contigs[currentContig].m_ConsensusLength+2);

			sequencePosition = m_Contigs[currentContig].m_Consensus + 1;
			memcpy(sequencePosition,sequenceBuffer,m_Contigs[currentContig].m_ConsensusLength);
			//printf("Consensus: %s.\n",sequencePosition);

			// find the BQ line
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strlen(buffer) == 0);

			if(strncmp(buffer,"BQ",2) != 0) {
				printf("Not BQ.\n");
				exit(1);
			}

			//printf("BQ buffer: %s.\n",buffer);

			// read in the entire consensus quality sequence
			sequencePosition = sequenceBuffer;

			do {
				in.getline(sequencePosition,SEQ_BUFFER_SIZE);
				positionsRead = (unsigned int)strlen(sequencePosition);
				sequencePosition += positionsRead;
			} while(positionsRead != 0);

			// convert the qualities to values
			sequencePosition  = m_Contigs[currentContig].m_ConsensusQuality + 1;
			consensusPosition = m_Contigs[currentContig].m_Consensus + 1;

			token = strtok(sequenceBuffer," ");
			for(;*consensusPosition;sequencePosition++) {

				// hop over all gaps
				while(*consensusPosition == 45) {
					*sequencePosition = 0;
					consensusPosition++;
					sequencePosition++;
				}

				*sequencePosition = atoi(token);

				// increment the consensus position
				consensusPosition++;

				// stop if this is the end of the consensus
				if(!*consensusPosition) break;

				// extract the next token
				token = strtok(NULL," ");

				if(!token) {
					printf("ERROR: No more tokens. Current consensus: %c\n",*consensusPosition);
					exit(1);
				}
			}

			// find the AF line
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strlen(buffer) == 0);

			if(strncmp(buffer,"AF",2) != 0) {
				printf("Not AF.\n");
				exit(1);
			}

			//printf("AF buffer: %s.\n",buffer);

			// create our contig sequences and prepare a masking array
			m_ContigSequences[currentContig] = new CSequence[m_NumContigSeqs[currentContig]];
			m_ContigSeqsMask[currentContig]  = new bool[m_NumContigSeqs[currentContig]];

			for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) m_ContigSeqsMask[currentContig][i] = false;

			// read all of the AF positions
			unsigned int contigSeqCounter = 0;

			do {

				if(strncmp(buffer,"BS",2) == 0) break;

				// skip past the AF token
				token = strtok(buffer," ");

				// extract the sequence name
				token = strtok(NULL," ");
				unsigned int nameLen = (unsigned int)strlen(token);

				m_ContigSequences[currentContig][contigSeqCounter].m_Name = new char[nameLen+1];
				strncpy(m_ContigSequences[currentContig][contigSeqCounter].m_Name,token,nameLen);
				m_ContigSequences[currentContig][contigSeqCounter].m_Name[nameLen] = 0;

				//printf("Sequence name: %s\n",m_ContigSequences[currentContig][contigSeqCounter].m_Name);

				// extract the complement direction
				token = strtok(NULL," ");
				if(*token == 'C') m_ContigSequences[currentContig][contigSeqCounter].m_IsReverseComplement = true;

				//if(m_ContigSequences[currentContig][contigSeqCounter].m_IsReverseComplement) printf("REVERSE.\n");

				// extract offset
				token = strtok(NULL," ");
				m_ContigSequences[currentContig][contigSeqCounter++].m_Offset = atoi(token);

				//printf("Offset: %d\n",atoi(token));

				// get new line
				in.getline(buffer,BUFFER_SIZE);

			} while(strlen(buffer) != 0);

			//for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) printf("Name: %s %d\n",m_ContigSequences[currentContig][i].m_Name,m_ContigSequences[currentContig][i].m_Offset);

			// skip over BS
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strlen(buffer) > 0);

			// find next line
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strlen(buffer) == 0);

			//
			// READ SEQUENCES FROM ACE FILE
			//

			int numSeqsRead = 0;

			while(strncmp(buffer,"RD",2) == 0) {

				// skip past the RD token
				token = strtok(buffer," ");

				// check the sequence name
				token = strtok(NULL," ");
				if(strcmp(token,m_ContigSequences[currentContig][numSeqsRead].m_Name) != 0) {
					printf("ERROR: Not the same sequence name. Token: %s (%s)\n",token,m_ContigSequences[currentContig][numSeqsRead].m_Name);
					exit(1);
				}

				// extract the # of sequence bases
				token = strtok(NULL," ");
				unsigned int numSeqBases = atoi(token);
				//printf("# of bases in sequence: %u\n",numSeqBases);

				// pick the larger of the two: consensus length or sequence size
				unsigned int allocationSize;
				if((numSeqBases+1) > (m_Contigs[currentContig].m_ConsensusLength+2)) allocationSize = numSeqBases + 1;
					else allocationSize = m_Contigs[currentContig].m_ConsensusLength + 2;

				// initialize our sequence arrays
				m_ContigSequences[currentContig][numSeqsRead].m_Bases     = new char[allocationSize];
				m_ContigSequences[currentContig][numSeqsRead].m_Qualities = new char[allocationSize];
				m_ContigSequences[currentContig][numSeqsRead].m_NumBases  = numSeqBases;

				unsigned int filenameLength = (unsigned int)strlen(m_ContigSequences[currentContig][numSeqsRead].m_Name)+6;
				m_ContigSequences[currentContig][numSeqsRead].m_Filename = new char[filenameLength+1];
				sprintf(m_ContigSequences[currentContig][numSeqsRead].m_Filename,"%s.phd.1",m_ContigSequences[currentContig][numSeqsRead].m_Name);

				//printf("Filename: %s\n",m_ContigSequences[currentContig][numSeqsRead].m_Filename);

				// read in the entire sequence
				sequencePosition = sequenceBuffer;

				do {
					in.getline(sequencePosition,SEQ_BUFFER_SIZE);
					positionsRead = (unsigned int)strlen(sequencePosition);
					sequencePosition += positionsRead;
				} while(positionsRead != 0);

				// replace all asterisks with hyphens and capitalize
				for(sequencePosition = sequenceBuffer;*sequencePosition;sequencePosition++) {
					if((*sequencePosition > 96) && (*sequencePosition < 123)) *sequencePosition -= 32;
					if(*sequencePosition == 42) *sequencePosition = 45;
				}

				// copy the sequence bases to the sequence object
				// TODO: Maybe we should do this directly in the sequence object?
				//printf("numSeqBases: %u, base length: %u\n",numSeqBases,m_Contigs[currentContig].m_ConsensusLength+2);
				strncpy(m_ContigSequences[currentContig][numSeqsRead].m_Bases,sequenceBuffer,numSeqBases);

				//printf("first base: %c\n",m_ContigSequences[currentContig][numSeqsRead].m_Bases[0]);

				// evaluate the QA tag
				do {
					in.getline(buffer,BUFFER_SIZE);
				} while(strncmp(buffer,"QA",2) != 0);

				// skip past the QA token
				token = strtok(buffer," ");

				// extract the quality clip start
				token = strtok(NULL," ");
				int qualClipStart = atoi(token);

				// extract the quality clip end
				token = strtok(NULL," ");
				int qualClipEnd = atoi(token);

				// extract the alignment clip start
				token = strtok(NULL," ");
				int alignClipStart = atoi(token);

				// extract the alignment clip end
				token = strtok(NULL," ");
				int alignClipEnd = atoi(token);

				// analyze sequence clip starting and ending points
				if(qualClipStart > alignClipStart) m_ContigSequences[currentContig][numSeqsRead].m_ClipStart = qualClipStart;
				else m_ContigSequences[currentContig][numSeqsRead].m_ClipStart = alignClipStart;

				if(qualClipEnd < alignClipEnd) m_ContigSequences[currentContig][numSeqsRead].m_ClipEnd = qualClipEnd;
				else m_ContigSequences[currentContig][numSeqsRead].m_ClipEnd = alignClipEnd;

				// mask the sequence if it has low quality
				if(((qualClipStart == -1) && (qualClipEnd == -1)) || ((alignClipStart == -1) && (alignClipEnd == -1))) {
					//printf("MASKING: %s.\n",m_ContigSequences[currentContig][numSeqsRead].m_Name);
					m_ContigSeqsMask[currentContig][numSeqsRead] = true;
					m_NumContigSeqsMasked[currentContig]++;

				} else {
					// check to see if file exists
					sprintf(absoluteFilename,"%s%cphd_dir%c%s",dirPrefix,OS_SLASH,OS_SLASH,m_ContigSequences[currentContig][numSeqsRead].m_Filename);

					ifstream seq(absoluteFilename);

					// if our sequence is missing, remove it
					if(seq.fail()) {
						printf("- MISSING: %s\n",absoluteFilename);
						m_ContigSeqsMask[currentContig][numSeqsRead] = true;
						m_NumContigSeqsMasked[currentContig]++;
					}
				}

				// increment our sequence counter
				numSeqsRead++;

				// find next sequence
				if(numSeqsRead == m_NumContigSeqs[currentContig]) break;

				do {
					in.getline(buffer,BUFFER_SIZE);
				} while(strncmp(buffer,"RD",2) != 0);
			}

			//printf("Sequences read: %d.\n",numSeqsRead);

		} // end contig for loop

		in.close();

	} else {

		// This is most likely an old style ace file then
		in.close();
		printf("ERROR: AS tag not found. Trying to load as old style ace file.\n");
		BuildOldAceContigs(dirPrefix,ace);
	}
}

// Loads the specified contig from the ace file
void CAce::BuildOldAceContigs(char* dirPrefix, char* ace) {

	char buffer[BUFFER_SIZE];
	char sequenceBuffer[SEQ_BUFFER_SIZE];
	char absoluteFilename[FILENAME_SIZE];

	char* consensusPosition;
	char* sequencePosition;
	char* token;

	unsigned int positionsRead;

	printf("Loading %s.\n",ace);

	sprintf(absoluteFilename,"%s%cedit_dir%c%s",dirPrefix,OS_SLASH,OS_SLASH,ace);

	//printf("Absolute filename: %s\n",absoluteFilename);

	ifstream in(absoluteFilename);

	// if our file is missing, bomb
	if(in.fail()) {
		printf("ERROR: Could not open the old version ace file: %s\n",absoluteFilename);
		exit(1);
	}

	// find the first non-empty line
	// TODO: might have to check for EOF here
	do {
		in.getline(buffer,BUFFER_SIZE);
	} while(strlen(buffer) == 0);

	// check to see if this is a valid old style ace file
	token = strtok(buffer," ");
	if(strncmp(token,"DNA",3) != 0) {
		printf("ERROR: This is neither an ace file nor an old version ace file.\n");
		exit(1);
	}

	// mark this as an old version ace object
	m_IsOldVersion = true;

	//
	// first pass (count contigs)
	//

	unsigned int numContigSeqs = 0;
	unsigned int numDnaTags    = 0;

	while(!in.eof()) {
		in.getline(buffer,BUFFER_SIZE);
		if(strncmp(buffer,"Assembled_from*",15) == 0) numContigSeqs++;
		if(strncmp(buffer,"DNA",3) == 0) numDnaTags++;
	}

	m_NumContigs = numDnaTags - numContigSeqs + 1;
	printf("Contigs: %u.\n",m_NumContigs);

	// initialize main variables
	m_Contigs = new CContig[m_NumContigs];

	m_NumContigSeqs       = new unsigned int[m_NumContigs];
	m_NumContigSeqsMasked = new unsigned int[m_NumContigs];
	m_ContigSequences     = new CSequence*[m_NumContigs];
	m_ContigSeqsMask      = new bool*[m_NumContigs];

	for(unsigned int i=0;i<m_NumContigs;i++) {
		m_ContigSequences[i]     = NULL;
		m_ContigSeqsMask[i]      = NULL;
		m_NumContigSeqsMasked[i] = 0;
	}

	//
	// second pass (count sequences per contig)
	//

	// go back to the beginning of the file
	in.clear();
	in.seekg(0,ios_base::beg);

	for(unsigned char currentContig=0;currentContig<m_NumContigs;currentContig++) {

		// find the next DNA tag (contig)
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strncmp(buffer,"DNA",3) != 0);

		// find next sequence tag
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strncmp(buffer,"Sequence",8) != 0);

		// count the number of Assembled_from* tags
		unsigned int numContigSeqs = 0;

		do {
			in.getline(buffer,BUFFER_SIZE);
			if(strncmp(buffer,"Assembled_from*",15) == 0) numContigSeqs++;
		} while(strncmp(buffer,"Base_segment",12) != 0);

		printf("Contig sequences: %u\n",numContigSeqs);
		
		// initialize sequence related arrays
		m_NumContigSeqs[currentContig]   = numContigSeqs;
		m_ContigSequences[currentContig] = new CSequence[numContigSeqs];
		m_ContigSeqsMask[currentContig]  = new bool[numContigSeqs];

		for(unsigned int i=0;i<numContigSeqs;i++) m_ContigSeqsMask[currentContig][i] = false;

		// skip over the sequences
		unsigned int skippedSeqs = 0;

		do {
			in.getline(buffer,BUFFER_SIZE);
			if(strncmp(buffer,"DNA",3) == 0) skippedSeqs++;
		} while(skippedSeqs < numContigSeqs);
	}

	//
	// third pass (read content)
	//

	// go back to the beginning of the file
	in.clear();
	in.seekg(0,ios_base::beg);

	for(unsigned char currentContig=0;currentContig<m_NumContigs;currentContig++) {

		// find the next DNA tag (contig)
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strncmp(buffer,"DNA",3) != 0);

		// extract the contig name
		token = strtok(buffer," ");
		token = strtok(NULL," ");

		unsigned int nameLen = (unsigned int)strlen(token);
		m_Contigs[currentContig].m_Name = new char[nameLen+1];
		strncpy(m_Contigs[currentContig].m_Name,token,nameLen);
		m_Contigs[currentContig].m_Name[nameLen] = 0;
		//printf("Contig name: %s.\n",m_Contigs[currentContig].m_Name);

		// read in the entire consensus sequence
		sequencePosition = sequenceBuffer;

		do {
			in.getline(sequencePosition,SEQ_BUFFER_SIZE);
			positionsRead = (unsigned int)strlen(sequencePosition);
			sequencePosition += positionsRead;
		} while(positionsRead != 0);

		// replace all asterisks with hyphens and capitalize
		for(sequencePosition = sequenceBuffer;*sequencePosition;sequencePosition++) {
			if((*sequencePosition > 96) && (*sequencePosition < 123)) *sequencePosition -= 32;
			if(*sequencePosition == 42) *sequencePosition = 45;
		}

		// extract the consensus length
		m_Contigs[currentContig].m_ConsensusLength = (unsigned int)strlen(sequenceBuffer);
		//printf("Consensus length: %d.\n",m_Contigs[currentContig].m_ConsensusLength);

		// TODO: is this right?
		m_Contigs[currentContig].m_Consensus        = new char[m_Contigs[currentContig].m_ConsensusLength+2];
		m_Contigs[currentContig].m_ConsensusQuality = new char[m_Contigs[currentContig].m_ConsensusLength+2];

		// TODO: do we really need to clear the entire base and quality sequence
		memset(m_Contigs[currentContig].m_Consensus,0,m_Contigs[currentContig].m_ConsensusLength+2);
		memset(m_Contigs[currentContig].m_ConsensusQuality,0,m_Contigs[currentContig].m_ConsensusLength+2);

		sequencePosition = m_Contigs[currentContig].m_Consensus + 1;
		memcpy(sequencePosition,sequenceBuffer,m_Contigs[currentContig].m_ConsensusLength);
		//printf("Consensus: %s.\n",sequencePosition);

		// find the BaseQuality line
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strlen(buffer) == 0);

		if(strncmp(buffer,"BaseQuality",11) != 0) {
			printf("Not BaseQuality.\n");
			exit(1);
		}

		//printf("BaseQuality buffer: %s.\n",buffer);

		// read in the entire consensus quality sequence
		sequencePosition = sequenceBuffer;

		do {
			in.getline(sequencePosition,SEQ_BUFFER_SIZE);
			positionsRead = (unsigned int)strlen(sequencePosition);
			sequencePosition += positionsRead;
		} while(positionsRead != 0);

		//printf("BaseQuality: %s.\n",sequenceBuffer);

		// convert the qualities to values
		sequencePosition  = m_Contigs[currentContig].m_ConsensusQuality + 1;
		consensusPosition = m_Contigs[currentContig].m_Consensus + 1;

		token = strtok(sequenceBuffer," ");
		for(;*consensusPosition;sequencePosition++) {

			// hop over all gaps
			while(*consensusPosition == 45) {
				*sequencePosition = 0;
				consensusPosition++;
				sequencePosition++;
			}

			*sequencePosition = atoi(token);

			// increment the consensus position
			consensusPosition++;

			// stop if this is the end of the consensus
			if(!*consensusPosition) break;

			// extract the next token
			token = strtok(NULL," ");

			if(!token) {
				printf("ERROR: No more tokens. Current consensus: %c\n",*consensusPosition);
				exit(1);
			}
		}

		// find the Sequence line
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strlen(buffer) == 0);

		if(strncmp(buffer,"Sequence",8) != 0) {
			printf("Not Sequence.\n");
			exit(1);
		}

		//printf("Sequence buffer: %s.\n",buffer);

		// read all of the Assembled_from* positions
		unsigned int contigSeqCounter = 0;

		do {

			if(strncmp(buffer,"Base_segment",12) == 0) break;

			// only evaluate the Assembled_from* tags
			if(strncmp(buffer,"Assembled_from*",15) == 0) {

				// skip past the Assembled_from* token
				token = strtok(buffer," ");

				// extract the sequence name
				token = strtok(NULL," ");
				unsigned int nameLen = (unsigned int)strlen(token);

				m_ContigSequences[currentContig][contigSeqCounter].m_Name = new char[nameLen+1];
				strncpy(m_ContigSequences[currentContig][contigSeqCounter].m_Name,token,nameLen);
				m_ContigSequences[currentContig][contigSeqCounter].m_Name[nameLen] = 0;

				//printf("Sequence name: %s.\n",m_ContigSequences[currentContig][contigSeqCounter].m_Name);

				// extract the complement direction
				sequencePosition = token + nameLen - 5;
				if(strncmp(sequencePosition,".comp",5) == 0) m_ContigSequences[currentContig][contigSeqCounter].m_IsReverseComplement = true;
				//if(m_ContigSequences[currentContig][contigSeqCounter].m_IsReverseComplement) printf("REVERSE.\n");

				// extract offset
				token = strtok(NULL," ");
				m_ContigSequences[currentContig][contigSeqCounter++].m_Offset = atoi(token);

				//printf("Offset: %d\n",atoi(token)-1);
			}

			// get new line
			in.getline(buffer,BUFFER_SIZE);

		} while(strlen(buffer) != 0);

		//for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) printf("Name: %s %d\n",m_ContigSequences[currentContig][i].m_Name,m_ContigSequences[currentContig][i].m_Offset);

		// skip over Base_segment
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strlen(buffer) > 0);

		// find next line
		do {
			in.getline(buffer,BUFFER_SIZE);
		} while(strlen(buffer) == 0);

		//
		// READ SEQUENCES FROM ACE FILE
		//

		int numSeqsRead = 0;

		while(strncmp(buffer,"DNA",3) == 0) {

			// skip past the DNA token
			token = strtok(buffer," ");

			// check the sequence name
			token = strtok(NULL," ");
			if(strcmp(token,m_ContigSequences[currentContig][numSeqsRead].m_Name) != 0) {
				printf("ERROR: Not the same sequence name. Token: %s (%s)\n",token,m_ContigSequences[currentContig][numSeqsRead].m_Name);
				exit(1);
			}

			// set the sequence filename
			int nameLength = (int)strlen(m_ContigSequences[currentContig][numSeqsRead].m_Name);
			strncpy(buffer,m_ContigSequences[currentContig][numSeqsRead].m_Name,BUFFER_SIZE);
			if(m_ContigSequences[currentContig][numSeqsRead].m_IsReverseComplement) nameLength -= 5;
			buffer[nameLength] = 0;

			unsigned int filenameLength = nameLength + 6;
			m_ContigSequences[currentContig][numSeqsRead].m_Filename = new char[filenameLength+1];
			sprintf(m_ContigSequences[currentContig][numSeqsRead].m_Filename,"%s.phd.1",buffer);

			//printf("Name: %s, filename: %s\n",m_ContigSequences[currentContig][numSeqsRead].m_Name,m_ContigSequences[currentContig][numSeqsRead].m_Filename);

			// read in the entire sequence
			sequencePosition = sequenceBuffer;

			do {
				in.getline(sequencePosition,SEQ_BUFFER_SIZE);
				positionsRead = (unsigned int)strlen(sequencePosition);
				sequencePosition += positionsRead;
			} while(positionsRead != 0);

			// replace all asterisks with hyphens and capitalize
			for(sequencePosition = sequenceBuffer;*sequencePosition;sequencePosition++) {
				if((*sequencePosition > 96) && (*sequencePosition < 123)) *sequencePosition -= 32;
				if(*sequencePosition == 42) *sequencePosition = 45;
			}

			// extract the # of sequence bases
			unsigned int numSeqBases = (unsigned int)strlen(sequenceBuffer);
			//printf("# of bases in sequence: %u\n",numSeqBases);

			// pick the larger of the two: consensus length or sequence size
			unsigned int allocationSize;
			if((numSeqBases+1) > (m_Contigs[currentContig].m_ConsensusLength+2)) allocationSize = numSeqBases + 1;
				else allocationSize = m_Contigs[currentContig].m_ConsensusLength + 2;

			// initialize the base array
			m_ContigSequences[currentContig][numSeqsRead].m_Bases     = new char[allocationSize];
			m_ContigSequences[currentContig][numSeqsRead].m_Qualities = new char[allocationSize];
			m_ContigSequences[currentContig][numSeqsRead].m_NumBases  = numSeqBases;

			// TODO: temporary
			memset(m_ContigSequences[currentContig][numSeqsRead].m_Bases,'_',sizeof(char)*m_Contigs[currentContig].m_ConsensusLength+2);
			
			// copy the sequence bases to the sequence object
			strncpy(m_ContigSequences[currentContig][numSeqsRead].m_Bases,sequenceBuffer,numSeqBases);
			//printf("first base: %c, second base: %c, third base: %c\n",m_ContigSequences[currentContig][numSeqsRead].m_Bases[0],m_ContigSequences[currentContig][numSeqsRead].m_Bases[1],m_ContigSequences[currentContig][numSeqsRead].m_Bases[2]);
			//printf("Bases: %s\n",m_ContigSequences[currentContig][numSeqsRead].m_Bases);

			// evaluate the Clipping* tag
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strncmp(buffer,"Clipping*",9) != 0);

			// skip past the Clipping* token
			token = strtok(buffer," ");

			// extract the clip start
			token = strtok(NULL," ");
			m_ContigSequences[currentContig][numSeqsRead].m_ClipStart = atoi(token);

			// extract the clip end
			token = strtok(NULL," ");
			m_ContigSequences[currentContig][numSeqsRead].m_ClipEnd   = atoi(token);

			// check to see if file exists
			sprintf(absoluteFilename,"%s%cphd_dir%c%s",dirPrefix,OS_SLASH,OS_SLASH,m_ContigSequences[currentContig][numSeqsRead].m_Filename);

			ifstream seq(absoluteFilename);

			// if our sequence is missing, remove it
			if(seq.fail()) {
				printf("- MISSING: %s\n",absoluteFilename);
				m_ContigSeqsMask[currentContig][numSeqsRead] = true;
				m_NumContigSeqsMasked[currentContig]++;
			}

			// increment our sequence counter
			numSeqsRead++;

			// find next sequence
			if(numSeqsRead == m_NumContigSeqs[currentContig]) break;

			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strncmp(buffer,"DNA",3) != 0);
		}

		printf("Sequences read: %d.\n",numSeqsRead);

	} // end contig for loop
	
	in.close();
}

// Loads the specified contig from the ace file
void CAce::LoadPhdFiles(char* dirPrefix) {

	char qualityBuffer[BUFFER_SIZE];
	char buffer[BUFFER_SIZE];
	char absoluteFilename[FILENAME_SIZE];

	char* qualityPosition;
	char* token;

	unsigned int seqLength;

	// return if we have no contigs to load
	if(m_NumContigs == 0) {
		printf("ERROR: Cannot load sequences for cluster with no contigs.\n");	
		exit(1);
	}

	// return if we have no sequences to load
	if(!m_ContigSequences) {
		printf("ERROR: Contig sequences not properly initialized.\n");	
		exit(1);
	}

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {

		//printf("Loading sequences for %s.\n",m_Contigs[currentContig].m_Name);

		for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) {

			// if the sequence is masked, skip sequence
			if(m_ContigSeqsMask[currentContig][i]) continue;
				
			//printf("- Loading sequence %u: %s.\n",i+1,m_ContigSequences[currentContig][i].m_Name);

			// construct our absolute filename
			sprintf(absoluteFilename,"%s%cphd_dir%c%s",dirPrefix,OS_SLASH,OS_SLASH,m_ContigSequences[currentContig][i].m_Filename);

			ifstream in(absoluteFilename);

			// if we can't open the file, the prefiltering has failed
			if(in.fail()) {
				printf("ERROR: Could not open %s.\n",absoluteFilename);
				exit(1);
			}
			
			// initialize our variables
			qualityPosition = qualityBuffer;
			seqLength       = 0;

			// find next BEGIN_DNA line
			do {
				in.getline(buffer,BUFFER_SIZE);
			} while(strncmp(buffer,"BEGIN_DNA",9) != 0);

			// keep going until we hit END_DNA
			in.getline(buffer,BUFFER_SIZE);

			do {
			
				// skip over the base
				token = strtok(buffer," ");

				// extract the quality
				token = strtok(NULL," ");
				*qualityPosition = atoi(token);
				
				// increment positions
				qualityPosition++;

				// get next line
				in.getline(buffer,BUFFER_SIZE);

				seqLength++;

			} while(strncmp(buffer,"END_DNA",7) != 0);

			// close the sequence file
			in.close();

			// spit out an error if our buffer size is too small
			if(seqLength > BUFFER_SIZE) {
				printf("ERROR: Sequence being parsed is larger than the allocated buffer size: %u (%u).\n",seqLength,BUFFER_SIZE);
				exit(1);
			}

			// if the sequence is reverse complement, reverse it
			if(m_ContigSequences[currentContig][i].m_IsReverseComplement) ReverseSequence(qualityBuffer,seqLength);

			char* bases           = m_ContigSequences[currentContig][i].m_Bases;
			char* qualities       = m_ContigSequences[currentContig][i].m_Qualities;
			unsigned int numBases = m_ContigSequences[currentContig][i].m_NumBases;

			unsigned int qualityCounter = 0;
			unsigned int baseCounter    = 0;
			
			while((qualityCounter < numBases) && (baseCounter < numBases)) {
				while((bases[baseCounter] == '-') && (baseCounter < numBases)) qualities[baseCounter++] = 0;
				qualities[baseCounter++] = qualityBuffer[qualityCounter++];
			}

			// trim the sequence
			//printf("Sequence length: %u, clip start: %u, clip end: %u\n",storedSeqLength,m_ContigSequences[currentContig][i].m_ClipStart,m_ContigSequences[currentContig][i].m_ClipEnd);
			
			// the trimming seems to work perfectly
			for(int j=0;j<m_ContigSequences[currentContig][i].m_ClipStart-1;j++) {
				m_ContigSequences[currentContig][i].m_Bases[j]     = 0;
				m_ContigSequences[currentContig][i].m_Qualities[j] = 0;
			}

			for(unsigned int j=m_ContigSequences[currentContig][i].m_ClipEnd;j<numBases;j++) {
				m_ContigSequences[currentContig][i].m_Bases[j]     = 0;
				m_ContigSequences[currentContig][i].m_Qualities[j] = 0;
			}
		
			// debugging
			//for(unsigned int j=0;j<20;j++) printf("%c",bases[j]);
			//printf("\n");

		} // end contig seqs for loop
	} // end contig for loop
}

// Reverses a sequence (for reverse complements)
void CAce::ReverseSequence(char* str, unsigned int length) {

	unsigned int end   = length - 1;
	unsigned int start = 0;

	while(start<end) {
		str[start]   ^= str[end];
		str[end]     ^= str[start];
		str[start++] ^= str[end--];
	}
}

// Aligns all of the sequences in the ace file
void CAce::AlignSequences(void) {
	
	// bomb if we're trying to align an aligned cluster
	if(m_IsAligned) {
		printf("ERROR: Cluster is already aligned.\n");
		exit(1);
	}

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {
		
		unsigned int consensusLength = m_Contigs[currentContig].m_ConsensusLength;
		//printf("Aligning %s\n",m_Contigs[currentContig].m_Name);
		
		for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) {

			// if the sequence is masked, skip sequence
			if(m_ContigSeqsMask[currentContig][i]) continue;
		
			// align the sequence
			if(m_IsOldVersion) m_ContigSequences[currentContig][i].AlignOldVersion(consensusLength);
				else m_ContigSequences[currentContig][i].Align(consensusLength);

		} // // end contig seqs for loop
	} // // end contig for loop

	m_IsAligned = true;
}

// Dumps all of the sequences in the ace file
void CAce::DumpSequences(void) {
	
	char dashBuffer[SEQ_BUFFER_SIZE];

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {
		
		unsigned int consensusLength = m_Contigs[currentContig].m_ConsensusLength + 1;
		
		printf("%s\n",m_Contigs[currentContig].m_Name);
				
		memset(dashBuffer,'=',sizeof(char)*consensusLength);
		dashBuffer[consensusLength] = 0;
		
		printf("================%s\n",dashBuffer);

		printf("%16s",m_Contigs[currentContig].m_Name);
		for(unsigned int j=0;j<consensusLength;j++) printf("%c",m_Contigs[currentContig].m_Consensus[j]);
		printf("\n");
		
		printf("================%s\n",dashBuffer);

		for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) {

			// if the sequence is masked, skip sequence
			if(m_ContigSeqsMask[currentContig][i]) continue;
	
			// align the sequence
			printf("%16s",m_ContigSequences[currentContig][i].m_Name);
		
			for(unsigned int j=0;j<consensusLength;j++) printf("%c",m_ContigSequences[currentContig][i].m_Bases[j]);
			printf("\n");
		} // end contig seqs for loop

		printf("\n");
	} // end contig for loop
}

// Dumps all of the sequences in the ace file at the specified position
void CAce::DumpPosition(unsigned int contigPos) {

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {
	
		for(unsigned int i=0;i<m_NumContigSeqs[currentContig];i++) {

			// if the sequence is masked, skip sequence
			if(m_ContigSeqsMask[currentContig][i]) continue;
	
			// display the position
			printf("%16s %c%c [%c] %c%c: %u\n",m_ContigSequences[currentContig][i].m_Name,m_ContigSequences[currentContig][i].m_Bases[contigPos-2],m_ContigSequences[currentContig][i].m_Bases[contigPos-1],m_ContigSequences[currentContig][i].m_Bases[contigPos],m_ContigSequences[currentContig][i].m_Bases[contigPos+1],m_ContigSequences[currentContig][i].m_Bases[contigPos+2],m_ContigSequences[currentContig][i].m_Qualities[contigPos]);
	
		} // end contig seqs for loop

		printf("\n");
	} // end contig for loop
}

// Scours the cluster for SNPs
void CAce::FindSNPs(void) {

	// initialize our variables
	unsigned int alleleCount[4];
	unsigned char alleleBestQuality[4];

	double alleleFreq[4];
	double lowestError[4];

	double sortedAlleleFreq[4];
	double sortedLowestError[4];

	const char* alleleDescription = "ACGT";

	char tmpBase;
	unsigned char tmpQuality;
		
	double dataPoint[5];

	unsigned int numFilteredBases;
	unsigned int numGaps;

	unsigned char* locusColumnBase;
	unsigned char* locusColumnQuality;

	unsigned int minSequenceThreshold = CRuntimeParameters::MinSequenceThreshold;

	CBayesianAlgorithm bayes;

	// open the xml file if requested
	bool useXML = CRuntimeParameters::EnableXmlOutput;
	ofstream xml;

	if(useXML) {
		xml.open(CRuntimeParameters::XmlFilename, ios_base::out | ios_base::trunc);
		
		// if we can't write the file, bomb
		if(xml.fail()) {
			printf("ERROR: Could not write to the XML file: %s\n",CRuntimeParameters::XmlFilename);
			exit(1);
		}

		xml << "<?xml version=\"1.0\"?>\n<ForageAnalysis>\n";
	}

	// bomb if we're trying to find SNPs in an unaligned cluster
	if(!m_IsAligned) {
		printf("ERROR: Attempted to find SNPs before aligning cluster.\n");
		exit(1);
	}

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {
		
		// skip this contig if the contig does not have enough sequences
		unsigned int numSeqs = m_NumContigSeqs[currentContig];
		unsigned int unmaskedNumSeqs = numSeqs - m_NumContigSeqsMasked[currentContig];

		/*
		printf("# seqs: %u, # unmasked: %u, # masked: %u\n",numSeqs,unmaskedNumSeqs,m_NumContigSeqsMasked[currentContig]);

		for(int pp=0;pp<numSeqs;pp++) {
			printf("%2u: %16s ",pp,m_ContigSequences[currentContig][pp].m_Name);
			if(m_ContigSeqsMask[currentContig][pp]) printf("masked\n");
				else printf("\n");
		}
		*/
		if(CRuntimeParameters::IsVerbose) printf("Looking for SNPs in %s (%u sequences)\n",m_Contigs[currentContig].m_Name,unmaskedNumSeqs);

		if(unmaskedNumSeqs < minSequenceThreshold) {
			if(CRuntimeParameters::IsVerbose) printf("INFO: %s contains less than %u sequences.\n",m_Contigs[currentContig].m_Name,CRuntimeParameters::MinSequenceThreshold);
			continue;
		}

		unsigned int consensusLength = m_Contigs[currentContig].m_ConsensusLength;
		
		locusColumnBase    = new unsigned char[unmaskedNumSeqs];
		locusColumnQuality = new unsigned char[unmaskedNumSeqs];
		
		// reset the contig gaps counter
		numGaps = 0;

		// localize our consensus variable
		char* consensus = m_Contigs[currentContig].m_Consensus;

		// check each position in the contig
		for(unsigned int contigPos=1;contigPos<consensusLength;contigPos++) {

			// if the consensus position is a gap, then continue
			if(consensus[contigPos] == '-') {
				numGaps++;
				continue;
			}

			// reset the counters
			memset(alleleCount,0,sizeof(unsigned int)*4);
			memset(alleleBestQuality,0,sizeof(unsigned char)*4);
			memset(locusColumnBase,0,sizeof(unsigned char)*4);
			memset(locusColumnQuality,0,sizeof(unsigned char)*4);
			numFilteredBases = 0;

			// check each sequence in the position
			for(unsigned int contigSeq=0;contigSeq<numSeqs;contigSeq++) {

				// if the sequence is masked, skip sequence
				if(m_ContigSeqsMask[currentContig][contigSeq]) continue;

				tmpBase    = m_ContigSequences[currentContig][contigSeq].m_Bases[contigPos];
				tmpQuality = m_ContigSequences[currentContig][contigSeq].m_Qualities[contigPos];	

				// if the sequence base is zero, skip the base
				if(tmpBase == 0) continue;

				// if the sequence quality is zero, skip the base
				if(tmpQuality == 0) continue;

				// count the occurrences of each allele also find the best quality found for each allele
				if(tmpBase == alleleDescription[0]) {
					
					alleleCount[0]++;
					if(tmpQuality > alleleBestQuality[0]) alleleBestQuality[0] = tmpQuality;
					locusColumnBase[numFilteredBases] = 0;
					//printf("A");

				} else if(tmpBase == alleleDescription[1]) {
					
					alleleCount[1]++;
					if(tmpQuality > alleleBestQuality[1]) alleleBestQuality[1] = tmpQuality;
					locusColumnBase[numFilteredBases] = 1;
					//printf("C");

				} else if(tmpBase == alleleDescription[2]) {
					
					alleleCount[2]++;
					if(tmpQuality > alleleBestQuality[2]) alleleBestQuality[2] = tmpQuality;
					locusColumnBase[numFilteredBases] = 2;
					//printf("G");
					
				} else if(tmpBase == alleleDescription[3]) {
					
					alleleCount[3]++;
					if(tmpQuality > alleleBestQuality[3]) alleleBestQuality[3] = tmpQuality;
					locusColumnBase[numFilteredBases] = 3;
					//printf("T");
					
				} else continue;
				
				locusColumnQuality[numFilteredBases++] = tmpQuality;

			} // end contig seqs for loop

			// skip this contig position if it does not have enough sequences
			if(numFilteredBases < minSequenceThreshold) continue;

			// get allele frequency and lowest error probability
			for(unsigned char k=0;k<4;k++) {
				alleleFreq[k]  = alleleCount[k] / (double)numFilteredBases;			
				lowestError[k] = pow(10.0,(double)alleleBestQuality[k]/-10.0);
			}

			// sort the frequencies and errors
			memcpy(sortedAlleleFreq,alleleFreq,sizeof(double)*4);
			memcpy(sortedLowestError,lowestError,sizeof(double)*4);
			CSynchronizedArrays::Sort(sortedAlleleFreq,sortedLowestError);

			// if the slice is monomorphic, continue
			if(sortedAlleleFreq[2] == 0) continue;

			// use the Bayesian algorithm
			double Psnp = bayes.Analyze(locusColumnBase,locusColumnQuality,numFilteredBases);
			
			//
			// Use the neural network
			//			
			dataPoint[0] = Psnp;
			dataPoint[1] = sortedAlleleFreq[3];
			dataPoint[2] = sortedLowestError[3];
			dataPoint[3] = sortedLowestError[2];
			dataPoint[4] = sortedLowestError[1];
			
			// normalize the data point variance
			for(int k=0;k<5;k++) dataPoint[k] = (dataPoint[k] - CMeanStdDev::Mean[k]) / CMeanStdDev::StdDev[k];			
			
			if(CDualCodebook::IsSNP(dataPoint)) {
				unsigned int B3 = (int)((-10.0 * log10(sortedLowestError[3])) + 0.5);
				unsigned int B2 = (int)((-10.0 * log10(sortedLowestError[2])) + 0.5);
				unsigned int B1 = (int)((-10.0 * log10(sortedLowestError[1])) + 0.5);
				unsigned int B0 = (int)((-10.0 * log10(sortedLowestError[0])) + 0.5);
	
				// 
				// copy the left flank from the consensus
				//

				unsigned int numFlankPositions = 0;
				unsigned int currentFlankPosition = contigPos - 1;

				while(numFlankPositions < FLANK_RADIUS && currentFlankPosition > 0) {
					if(m_Contigs[currentContig].m_Consensus[currentFlankPosition--] == '-') continue;
					numFlankPositions++;
				}

				currentFlankPosition++;
				char* leftFlank = new char[numFlankPositions+1];
				leftFlank[numFlankPositions] = 0;

				for(unsigned int tempPosition=0;tempPosition<numFlankPositions;currentFlankPosition++) {
					if(m_Contigs[currentContig].m_Consensus[currentFlankPosition] == '-') continue;
					leftFlank[tempPosition++] = m_Contigs[currentContig].m_Consensus[currentFlankPosition];
				}

				// 
				// copy the right flank from the consensus
				//

				numFlankPositions = 0;
				currentFlankPosition = contigPos + 1;

				while(numFlankPositions < FLANK_RADIUS && currentFlankPosition < m_Contigs[currentContig].m_ConsensusLength) {
					if(m_Contigs[currentContig].m_Consensus[currentFlankPosition++] == '-') continue;
					numFlankPositions++;
				}

				char* rightFlank = new char[numFlankPositions+1];
				rightFlank[numFlankPositions] = 0;

				currentFlankPosition = contigPos + 1;
				
				for(unsigned int tempPosition=0;tempPosition<numFlankPositions;currentFlankPosition++) {
					if(m_Contigs[currentContig].m_Consensus[currentFlankPosition] == '-') continue;
					rightFlank[tempPosition++] = m_Contigs[currentContig].m_Consensus[currentFlankPosition];
				}

				// 
				// find the observed variation
				//

				char obsVar = '*';

				bool a = false;
				bool c = false;
				bool g = false;
				bool t = false;

				if(alleleCount[0] > 0) a = true;
				if(alleleCount[1] > 0) c = true;
				if(alleleCount[2] > 0) g = true;
				if(alleleCount[3] > 0) t = true;

				if(a && c && !g && !t) obsVar = 'M';
				if(a && !c && g && !t) obsVar = 'R';
				if(a && !c && !g && t) obsVar = 'W';
				if(!a && c && g && !t) obsVar = 'S';
				if(!a && c && !g && t) obsVar = 'Y';
				if(!a && !c && g && t) obsVar = 'K';

				if(a && c && g && !t) obsVar = 'V';
				if(a && c && !g && t) obsVar = 'H';
				if(a && !c && g && t) obsVar = 'D';
				if(!a && c && g && t) obsVar = 'B';

				if(a && c && g && t) obsVar = 'N';

				// 
				// print information
				//

				printf("Forage SNP @ %u (unpadded: %u) P(SNP): %.4f P(VAR): %.4f\n       Freq3: %.2f | B3: %u B2: %u | A: %u C: %u G: %u T: %u\n",contigPos,contigPos-numGaps,Psnp,bayes.m_Pvar,sortedAlleleFreq[3],B3,B2,alleleCount[0],alleleCount[1],alleleCount[2],alleleCount[3]);
				printf("       %32s %c %s\n\n",leftFlank,GetIupacAmbiguityCode(bayes.m_Variation),rightFlank);

				// write to the xml file
				if(useXML) {
					xml << "   <snp contig=\"" << m_Contigs[currentContig].m_Name << "\" locus=\"" << contigPos << "\" unpadded_locus=\"" << (contigPos-numGaps) << "\">\n";
					xml << "      <psnp>" << Psnp << "</psnp>\n";
					xml << "      <pvar>" << bayes.m_Pvar << "</pvar>\n";
					xml << "      <var>" << GetIupacAmbiguityCode(bayes.m_Variation) << "</var>\n";
					xml << "      <obsvar>" << obsVar << "</obsvar>\n";
					xml << "      <depth>" << numFilteredBases << "</depth>\n";
					xml << "      <quality maj=\"" << B3 << "\" min1=\"" << B2 << "\" min2=\"" << B1 << "\" min3=\"" << B0 << "\"></quality>\n";
					xml << "      <leftflank>" << leftFlank << "</leftflank>\n";
					xml << "      <rightflank>" << rightFlank << "</rightflank>\n";
					xml << "      <frequency maj=\"" << sortedAlleleFreq[3] << "\" min1=\"" << sortedAlleleFreq[2] << "\" min2=\"" << sortedAlleleFreq[1] << "\" min3=\"" << sortedAlleleFreq[0] << "\">";
					
					// figure out what sort of secondary allele frequency is exhibited
					if(sortedAlleleFreq[2] >= 0.20) xml << "frequent";
						else if(sortedAlleleFreq[2] >= 0.05) xml << "common";
							else xml << "rare";

					xml << "</frequency>\n";
					xml << "   </snp>\n";
				}

				// clean up
				delete [] leftFlank;
			}

		} // end consensus postion for loop

		// do some housekeeping
		delete [] locusColumnBase;
		delete [] locusColumnQuality;

	} // end contig for loop

	// close the xml file
	if(useXML) {
	
		xml << "</ForageAnalysis>\n";
		xml.close();
	}
}

// Scours the cluster for paralogs
void CAce::FindParalogs(void) {

	// initialize our variables
	char consensusBase, tmpBase;
	unsigned char consensusQuality, tmpQuality;

	double P_error_consensus;
	double P_error_buffer;

	unsigned int* numDiscs;
	unsigned int* comparedLength;
	double* expDiscTotal;

	double pairwisePPOLY = CBayesianUtils::PAIRWISE_P_POLY;
	
	// bomb if we're trying to find paralogs in an unaligned cluster
	if(!m_IsAligned) {
		printf("ERROR: Attempted to find paralogs before aligning cluster.\n");
		exit(1);
	}

	for(unsigned int currentContig=0;currentContig<m_NumContigs;currentContig++) {
		
		if(CRuntimeParameters::IsVerbose) printf("Looking for paralogs in %s\n",m_Contigs[currentContig].m_Name);

		// skip this contig if all of the sequences are already masked
		unsigned int numSeqs = m_NumContigSeqs[currentContig];
		unsigned int unmaskedNumSeqs = numSeqs - m_NumContigSeqsMasked[currentContig];

		if(unmaskedNumSeqs == 0) continue;

		//printf("Looking for paralogs in %s\n",m_Contigs[currentContig].m_Name);

		// initialize contig variables
		numDiscs       = new unsigned int[numSeqs];
		comparedLength = new unsigned int[numSeqs];
		expDiscTotal   = new double[numSeqs];

		memset(numDiscs,0,sizeof(unsigned int)*numSeqs);
		memset(comparedLength,0,sizeof(unsigned int)*numSeqs);
		memset(expDiscTotal,0,sizeof(double)*numSeqs);

		unsigned int consensusLength = m_Contigs[currentContig].m_ConsensusLength;
		
		// localize our consensus variables
		char* consensusBases     = m_Contigs[currentContig].m_Consensus;
		char* consensusQualities = m_Contigs[currentContig].m_ConsensusQuality;

		//
		// gather pair-wise statistics on discrepencies
		//
		for(unsigned int contigPos=1;contigPos<consensusLength;contigPos++) {

			consensusBase    = consensusBases[contigPos];
			consensusQuality = consensusQualities[contigPos];
			
			// if the consensus position is a gap, then continue
			if(consensusBase == '-') continue;

			// continue if the consensus position doesn't meet the minimum quality
			if(consensusQuality < 6) continue;
			
			P_error_consensus = pow(10.0,(double)consensusQuality/-10.0);
			
			// check each sequence in the position
			for(unsigned int contigSeq=0;contigSeq<numSeqs;contigSeq++) {

				// if the sequence is masked, skip sequence
				if(m_ContigSeqsMask[currentContig][contigSeq]) continue;

				tmpBase    = m_ContigSequences[currentContig][contigSeq].m_Bases[contigPos];
				tmpQuality = m_ContigSequences[currentContig][contigSeq].m_Qualities[contigPos];
				
				// continue if the sequence position is a gap, masked, or empty
				if(tmpBase == '-' || tmpBase == 'X' || tmpBase == 0) continue;
				
				// continue if the sequence position doesn't meet the minimum quality level
				if(tmpQuality < 6) continue;
				
				P_error_buffer = pow(10.0,(double)tmpQuality/-10.0);
				
				expDiscTotal[contigSeq] += 1.0 - (P_error_consensus * P_error_buffer) / 3.0 - (1.0 - P_error_consensus) * (1.0 - P_error_buffer);
				
				if(tmpBase != consensusBase) numDiscs[contigSeq]++;
				 
				comparedLength[contigSeq]++;

			} // end contig seqs for loop

		} // end consensus postion for loop

		//
		// calculate the probability that each sequence is native
		//
		double lambda_native, lambda_paralog;
		double ratio;
		double P_native;
		
		for(unsigned int contigSeq=0;contigSeq<numSeqs;contigSeq++) {

			// if the sequence is masked, skip sequence
			if(m_ContigSeqsMask[currentContig][contigSeq]) continue;

			lambda_native  = comparedLength[contigSeq] * pairwisePPOLY + expDiscTotal[contigSeq];
			lambda_paralog = comparedLength[contigSeq] * 0.02  + expDiscTotal[contigSeq];
			
			// TODO: I don't think this should be the case.
			if(lambda_paralog >= 1) {
			
				ratio = exp(lambda_paralog - lambda_native) * pow(lambda_native / lambda_paralog, (double)numDiscs[contigSeq]);		
				P_native = ratio / (1.0 + ratio);
				
				if(P_native < CBayesianUtils::P_NAT_THRESHOLD) {
					if(CRuntimeParameters::IsVerbose) printf("  - %s marked as a paralog: P(NAT): %f\n",m_ContigSequences[currentContig][contigSeq].m_Name,P_native);
					m_ContigSeqsMask[currentContig][contigSeq] = true;
					m_NumContigSeqsMasked[currentContig]++;
				}

			} else {
				// in this case we have a paralog as well
				// this is just to save time from calculating ratio, etc.
				if(CRuntimeParameters::IsVerbose) printf("  - %s marked as a paralog: P(NAT): %f\n",m_ContigSequences[currentContig][contigSeq].m_Name,P_native);
				m_ContigSeqsMask[currentContig][contigSeq] = true;
				m_NumContigSeqsMasked[currentContig]++;
			}

		} // end contig seq for loop
	} // end contig for loop
}

// Retrieves a IUPAC ambiguity code based on the variation string
char CAce::GetIupacAmbiguityCode(unsigned char variation) {

	switch(variation) {
		case 1:
			return 'A';
			break;
		case 2:
			return 'C';
			break;
		case 3:
			return 'G';
			break;
		case 4:
			return 'T';
			break;
		case 5:
			return 'M';
			break;
		case 6:
			return 'R';
			break;
		case 7:
			return 'W';
			break;
		case 8:
			return 'S';
			break;
		case 9:
			return 'Y';
			break;
		case 10:
			return 'K';
			break;
		case 11:
			return 'V';
			break;
		case 12:
			return 'H';
			break;
		case 13:
			return 'D';
			break;
		case 14:
			return 'B';
			break;
		case 15:
			return 'N';
			break;
		default:
			printf("ERROR: unknown variation encountered: %d\n",variation);
			exit(1);
	}
}
