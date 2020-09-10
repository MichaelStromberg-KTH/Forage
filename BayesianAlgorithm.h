#include <math.h>
#include "BayesianUtils.h"
#include "Key.h"
#include "KeyTable.h"

class CBayesianAlgorithm {
public:
	CBayesianAlgorithm();
	~CBayesianAlgorithm();
	// Returns the Bayesian probability that an alignment position is polymorphic
	double Analyze(unsigned char* bases, unsigned char* qualities, unsigned int numSeqs);
	// The most probable variation
	unsigned char m_Variation;
	// The Bayesian probability that the variation is correctly classified
	double m_Pvar;
	// The Bayesian probability that an alignment position is polymorphic
	double m_Psnp;
};
