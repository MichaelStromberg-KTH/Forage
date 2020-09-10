#include "SynchronizedArrays.h"

void CSynchronizedArrays::Sort(double* alleleFreq, double* lowestError) {
	for(unsigned int i=0;i<4;i++) {
		for(unsigned int j=i;j>0 && alleleFreq[j-1]>alleleFreq[j];j--) {
			Swap(alleleFreq,j);
			Swap(lowestError,j);
		}
	}
}

void CSynchronizedArrays::Swap(double* x, unsigned int a) {
	double t = x[a];
	x[a]     = x[a-1];
	x[a-1]   = t;
}
