#include "BayesianUtils.h"

// declare the static variables
double CBayesianUtils::P_POLY          = 0.0;
double CBayesianUtils::PAIRWISE_P_POLY = 0.0;

double CBayesianUtils::P_MULT1         = 0.0;
double CBayesianUtils::P_MULT2         = 0.0;
double CBayesianUtils::P_MULT3         = 0.0;
double CBayesianUtils::P_MULT4         = 0.0;

double CBayesianUtils::P_NAT_THRESHOLD = 0.75;

// Precalculate nucleotide multiplicity
void CBayesianUtils::PreCalculate(double ppoly) {

	// construct our likelihood estimates for each nucleotide multiplicity
	// multiply each estimate by the polymorphism likelihood - P_prior(# alleles)

	// calculate our probability u
	double phi = 28.0 + 108.0 * ppoly + 12.0 * sqrt(9.0+42.0*ppoly+81.0*(ppoly*ppoly));
	double u   = (pow(phi,2.0/3.0) - 8.0 - 2.0 * pow(phi,1.0/3.0)) / (6.0 * pow(phi,1.0/3.0));

	P_MULT1 = 1.0 - u;
	P_MULT2 = u;
	P_MULT3 = u*u;
	P_MULT4 = u*u*u;

	// initialize the SNP probability variables
	P_POLY          = ppoly;
	PAIRWISE_P_POLY = ppoly / 3.0;
}

// Returns a probability hash filled with our prior probabilities
CProbHash CBayesianUtils::GetPriorVariations(unsigned int depth) {

	// P_prior(S) = P_prior(S|variation type) * P_prior(variation type|# alleles) * P_prior(# alleles)
	// P_prior(AAGAG) = P_prior(AAGAG|AG) * P_prior(AG|2 alleles) * P_prior(2 alleles)

	// base composition prior probabilities - P_prior(S|variation type)		
	double P_BC_MULT1 = 1.0;
	double P_BC_MULT2 = 1.0 / (depth + 1.0);
	double P_BC_MULT3 = 1.0 / (depth + 1.0) / (depth + 2.0);
	double P_BC_MULT4 = 1.0 / (depth + 1.0) / (depth + 2.0) / (depth + 3.0);

	// variation type prior probabilities - P_prior(variation type|# alleles)
	double P_VAR_MULT1 = 0.25;
	//double P_VAR_MULT2 = 1.0 / 6.0;
	double P_VAR_MULT3 = 0.25;
	double P_VAR_MULT4 = 1.0;

	// based on empirical variations found in the CGAP validated SNP database
	double P_VAR_MULT2_M_K = 0.075;
	double P_VAR_MULT2_R_Y = 0.355;
	double P_VAR_MULT2_W   = 0.09;
	double P_VAR_MULT2_S   = 0.05;

	// assign variation-specific absolute priors	
	CProbHash ph;

	ph.Add(0,1.0);

	double tempDouble = P_VAR_MULT1 * P_MULT1;
	ph.Add(1,tempDouble);
	ph.Add(2,tempDouble);
	ph.Add(4,tempDouble);
	ph.Add(8,tempDouble);

	tempDouble = P_BC_MULT2 * P_MULT2;
	ph.Add(3,tempDouble*P_VAR_MULT2_M_K);
	ph.Add(5,tempDouble*P_VAR_MULT2_R_Y);
	ph.Add(9,tempDouble*P_VAR_MULT2_W);
	ph.Add(6,tempDouble*P_VAR_MULT2_S);
	ph.Add(10,tempDouble*P_VAR_MULT2_R_Y);
	ph.Add(12,tempDouble*P_VAR_MULT2_M_K);

	tempDouble = P_BC_MULT3 * P_VAR_MULT3 * P_MULT3;
	ph.Add(7,tempDouble);
	ph.Add(11,tempDouble);
	ph.Add(13,tempDouble);
	ph.Add(14,tempDouble);

	ph.Add(15,P_BC_MULT4 * P_VAR_MULT4 * P_MULT4);

	return ph;
}
