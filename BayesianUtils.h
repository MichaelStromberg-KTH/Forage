#include "ProbHash.h"
#include <math.h>

class CBayesianUtils {
public:
	// Precalculate nucleotide multiplicity
	static void PreCalculate(double ppoly);
	// Returns a probability hash filled with our prior probabilities
	static CProbHash GetPriorVariations(unsigned int depth);
	// SNP probability
	static double P_POLY;
	// pairwise SNP probability
	static double PAIRWISE_P_POLY;
	// the threshold for the probability that a sequence is native
	static double P_NAT_THRESHOLD;
private:
	// nucleotide multiplicity 1 probability
	static double P_MULT1;
	// nucleotide multiplicity 2 probability
	static double P_MULT2;
	// nucleotide multiplicity 3 probability
	static double P_MULT3;
	// nucleotide multiplicity 4 probability
	static double P_MULT4;
};
