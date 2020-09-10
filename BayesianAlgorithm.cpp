#include "BayesianAlgorithm.h"

CBayesianAlgorithm::CBayesianAlgorithm() {}

CBayesianAlgorithm::~CBayesianAlgorithm() {}

// Returns the Bayesian probability that an alignment position is polymorphic

// Performance data (iterations / s) with 5 sequences 
// ==================================================
// Perl reference      171
// Java reference	  3738 (22x)
// Microsoft C++	  4364
// Intel C++		  5826 (34x)
// Intel C++ Opt2	  6085 (36x)

// Microsoft C++ char 5976
// Intel C++ char     9086 (53x) (2.4x Java)

double CBayesianAlgorithm::Analyze(unsigned char* bases, unsigned char* qualities, unsigned int numSeqs) {

	// If we are using more sequences than has been tested with our hash function
	// exit the program.
	if(numSeqs > 1023) {
		printf("Number of sequences: %u\n",numSeqs);
		printf("The hash function in the Key class has only been tested with 1023 sequences.\n");
		exit(1);
	}

	unsigned int depth = numSeqs;

	// reset our main variables
	m_Psnp      = 0.0;
	m_Pvar      = 0.0;
	m_Variation = 0;
	
	//
	// calculate base error probability
	//
	double** Perrors = new double*[depth];

	for(unsigned int i=0;i<depth;i++) {
		Perrors[i]      = new double[4];

		double error    = pow(10.0,(double)qualities[i]/-10.0);
		double diverror = error / 3.0;

		for(unsigned char j=0;j<4;j++)
			if(bases[i] == j) Perrors[i][j] = 1.0 - error;
				else Perrors[i][j] = diverror;
	}
	
	// assign variation-specific absolute priors	
	CProbHash priorVariation = CBayesianUtils::GetPriorVariations(depth);

	// processed variation probabilities
	CProbHash probVariation;

	CKeyTable keys, newKeys;
	CKey seedKey(1.0,0);
	keys.Put(&seedKey);

	CKey tempKey, newKey;

	unsigned int maxterms = 50;
	double priorVar, newPriorVar;

	unsigned int hashCode;
	double sum;
	CKey* keyHashes = NULL;
	
	// consider each position in the column		
	for(unsigned int depthIteration=0;depthIteration<depth;depthIteration++) {
	
		//printf("Starting depth: %u, keys: %u.\n",depthIteration,keys.m_Size);
		unsigned int numKeyHashes = keys.m_Size;
		
		// clear previous key hashes
		if(keyHashes) delete [] keyHashes;

		keyHashes = new CKey[numKeyHashes];
		keys.GetKeys(keyHashes);
		
		// consider each old variation
		for(unsigned int i=0;i<numKeyHashes;i++) {

			tempKey = keyHashes[i];
			//printf("* Using key %u:%u:%u:%u sum: %f variation: %u.\n",tempKey.m_a,tempKey.m_c,tempKey.m_g,tempKey.m_t,tempKey.m_Sum,tempKey.m_Variation);

			priorVar	 = priorVariation.Get(tempKey.m_Variation);
			sum			 = tempKey.m_Sum;
	
			// consider contribution of each nucleotide at the current level
			for(unsigned char base=0;base<4;base++) {

				// create new variation
				newKey = tempKey.MakeNewKey(base);

				// find out if there is a duplicate
				hashCode = newKey.HashCode();
				if(newKeys.Contains(hashCode)) newKey = newKeys.Get(hashCode);

				// get the prior variation data
				newPriorVar = priorVariation.Get(newKey.m_Variation);

				// calculate the sum of probabilities					
				newKey.CalculateSum(base,depthIteration,priorVar,newPriorVar,sum,Perrors[depthIteration][base]);

				//printf("Adding key %u:%u:%u:%u sum: %f variation: %u to newKeys.\n",newKey.m_a,newKey.m_c,newKey.m_g,newKey.m_t,newKey.m_Sum,newKey.m_Variation);

				// add the new key to the key set
				newKeys.Put(&newKey);
			}
		}
		
		// Reduce number of terms by keeping the top 50 probabilities
		keys.SortAndAdd(&newKeys,depthIteration);
	}

	// do some cleanup
	for(unsigned int i=0;i<depth;i++) delete [] Perrors[i];
	delete [] Perrors;

	delete [] keyHashes;

	// calculate the sum of the probabilities so that we can normalize them
	// and update total posterior probability for each term
	double normSum = 0;
	double tmpDouble;
	unsigned char variation;

	unsigned int numKeys = keys.m_Size;
	CKey* keyArray = new CKey[numKeys];
	keys.GetKeys(keyArray);

	for(unsigned int i=0;i<numKeys;i++) {
		variation 	= keyArray[i].m_Variation;
		sum 		= keyArray[i].m_Sum;

		// weight the variation probability with prior variation
		tmpDouble = probVariation.Get(variation) + sum;
		probVariation.Add(variation,tmpDouble);

		// aggregate variation probability weighted by prior variation
		normSum += sum;
	}

	// clear previous key array
	delete [] keyArray;

	// Normalize probabilities		
	double probMultiplicity[5];
	memset(probMultiplicity,0,sizeof(double)*5);

	double maxProbVariation = 0;
	char maxVariation       = 0;

	//double[] hashProbabilities 	= probVariation.
		
	//	probVariation.hash;
	//String[] hashVariations 	= probVariation.hashLabels;

	// public String[] hashLabels = { "","a","c","g","t","ac","ag","at","cg","ct","gt","acg","act","agt","cgt","acgt" };

	// normalize the monomorphic variations
	for(unsigned int i=1;i<5;i++) {
		tmpDouble = probVariation.m_hash[i] / normSum;
		probMultiplicity[1] += tmpDouble;

		if(tmpDouble > maxProbVariation) {
			maxProbVariation = tmpDouble;
			maxVariation     = i; 
		}
	}

	// normalize the bi-allelic variations
	for(unsigned int i=5;i<11;i++) {
		tmpDouble = probVariation.m_hash[i] / normSum;
		probMultiplicity[2] += tmpDouble;

		if(tmpDouble > maxProbVariation) {
			maxProbVariation = tmpDouble;
			maxVariation     = i; 
		}
	}

	// normalize the tri-allelic variations
	for(unsigned int i=11;i<15;i++) {
		tmpDouble = probVariation.m_hash[i] / normSum;
		probMultiplicity[3] += tmpDouble;

		if(tmpDouble > maxProbVariation) {
			maxProbVariation = tmpDouble;
			maxVariation     = i; 
		}
	}

	// normalize the tetra-allelic variation
	tmpDouble = probVariation.m_hash[15] / normSum;
	probMultiplicity[4] = tmpDouble;

	if(tmpDouble > maxProbVariation) {
		maxProbVariation = tmpDouble;
		maxVariation     = 15; 

	}

	m_Psnp      = 1 - probMultiplicity[1];
	m_Pvar      = maxProbVariation;
	m_Variation = maxVariation;

	return m_Psnp;
}
