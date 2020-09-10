#include "ProbHash.h"

CProbHash::CProbHash(void) {
	// initialize our hash array
	memset(m_hash,0,sizeof(double)*16);
}

CProbHash::~CProbHash(void) {}

// release vc++ benchmarking
// with string hash: 4456 ms / 22439597
// with switch:      2994 ms / 33396725

// intel c++ benchmarking
// with switch:      1682 ms / 59438218

// release vc++ benchmarking
// with switch/char  1252 ms / 79884966

// intel c++ benchmarking
// with switch:       551 ms / 181556740 (8 x)

// Places the given value into the hash
void CProbHash::Add(unsigned char variation, double value) {

	switch(variation) {
		// empty space entry
		case 0:
			m_hash[0] = value;
			break;
		
		// monomorphic variations
		case 1: // a
			m_hash[1] = value;
			break;

		case 2: // c
			m_hash[2] = value;
			break;

		case 4: // g
			m_hash[3] = value;
			break;

		case 8: // t
			m_hash[4] = value;
			break;

		// bi-allelic variations
		case 3: // ac
			m_hash[5] = value;
			break;

		case 5: // ag
			m_hash[6] = value;
			break;

		case 9: // at
			m_hash[7] = value;
			break;

		case 6: // cg
			m_hash[8] = value;
			break;

		case 10: // ct
			m_hash[9] = value;
			break;

		case 12: // gt
			m_hash[10] = value;
			break;

		// tri-allelic variations
		case 7: // acg
			m_hash[11] = value;
			break;

		case 11: // act
			m_hash[12] = value;
			break;

		case 13: // agt
			m_hash[13] = value;
			break;

		case 14: // cgt
			m_hash[14] = value;
			break;

		// tetra-allelic variations
		case 15: // acgt
			m_hash[15] = value;
			break;

		// default
		default:
			printf("Invalid hash entry adding to probability hash.\n");
			exit(1);
	}
}

// Retrieves the given value from the hash
double CProbHash::Get(unsigned char variation) {

	double value = 0.0;

	switch(variation) {
		// empty space entry
		case 0:
			value = m_hash[0];
			break;
		
		// monomorphic variations
		case 1: // a
			value = m_hash[1];
			break;

		case 2: // c
			value = m_hash[2];
			break;

		case 4: // g
			value = m_hash[3];
			break;

		case 8: // t
			value = m_hash[4];
			break;

		// bi-allelic variations
		case 3: // ac
			value = m_hash[5];
			break;

		case 5: // ag
			value = m_hash[6];
			break;

		case 9: // at
			value = m_hash[7];
			break;

		case 6: // cg
			value = m_hash[8];
			break;

		case 10: // ct
			value = m_hash[9];
			break;

		case 12: // gt
			value = m_hash[10];
			break;

		// tri-allelic variations
		case 7: // acg
			value = m_hash[11];
			break;

		case 11: // act
			value = m_hash[12];
			break;

		case 13: // agt
			value = m_hash[13];
			break;

		case 14: // cgt
			value = m_hash[14];
			break;

		// tetra-allelic variations
		case 15: // acgt
			value = m_hash[15];
			break;

		// default
		default:
			printf("Invalid hash entry adding to probability hash.\n");
			exit(1);
	}

	return value;
}

/*
// Places the given value into the hash
void CProbHash::Add(char* variation, double value) {

	unsigned int len = (unsigned int)strlen(variation);

	switch(len) {
		// empty space entry
		case 0:
			m_hash[0] = value;
			break;
		
		// monomorphic variations
		case 1:
			switch(variation[0]) {
				case 'a':
					m_hash[1] = value;
					break;
				case 'c':
					m_hash[2] = value;
					break;
				case 'g':
					m_hash[3] = value;
					break;
				case 't':
					m_hash[4] = value;
					break;
				default:
					printf("Invalid hash entry a1.\n");
					exit(1);
			}
			break;

		// bi-allelic variations
		case 2:
			switch(variation[0]) {
				case 'a':
					switch(variation[1]) {
						case 'c':
							m_hash[5] = value;
							break;
						case 'g':
							m_hash[6] = value;
							break;
						case 't':
							m_hash[7] = value;
							break;
						default:
							printf("Invalid hash entry a2.\n");
							exit(1);
					}
					break;
				case 'c':
					switch(variation[1]) {
						case 'g':
							m_hash[8] = value;
							break;
						case 't':
							m_hash[9] = value;
							break;
						default:
							printf("Invalid hash entry a3.\n");
							exit(1);
					}
					break;
				case 'g':
					switch(variation[1]) {
						case 't':
							m_hash[10] = value;
							break;
						default:
							printf("Invalid hash entry a4.\n");
							exit(1);
					}
					break;
			}
			break;

		// tri-allelic variations
		case 3:
			switch(variation[0]) {
				case 'a':
					if((variation[1] == 'c') && (variation[2] == 'g')) {
						m_hash[11] = value;
					} else if((variation[1] == 'c') && (variation[2] == 't')) {
						m_hash[12] = value;
					} else if((variation[1] == 'g') && (variation[2] == 't')) {
						m_hash[13] = value;
					} else {
						printf("Invalid hash entry a5.\n");
						exit(1);
					}

					break;
				case 'c':
					if((variation[1] == 'g') && (variation[2] == 't')) {
						m_hash[14] = value;
					} else {
						printf("Invalid hash entry a6.\n");
						exit(1);
					}
					break;
				default:
					printf("Invalid hash entry a7.\n");
					exit(1);
			}
			break;

		// tetra-allelic variations
		case 4:
			m_hash[15] = value;
			break;

		// default
		default:
			printf("Invalid hash entry a8.\n");
			exit(1);
	}
}

// Retrieves the given value from the hash
double CProbHash::Get(char* variation) {

	double value = 0.0;

	unsigned int len = (unsigned int)strlen(variation);

	switch(len) {
		// empty space entry
		case 0:
			value = m_hash[0];
			break;
		
		// monomorphic variations
		case 1:
			switch(variation[0]) {
				case 'a':
					value = m_hash[1];
					break;
				case 'c':
					value = m_hash[2];
					break;
				case 'g':
					value = m_hash[3];
					break;
				case 't':
					value = m_hash[4];
					break;
				default:
					printf("Invalid hash entry 1.\n");
					exit(1);
			}
			break;

		// bi-allelic variations
		case 2:
			switch(variation[0]) {
				case 'a':
					switch(variation[1]) {
						case 'c':
							value = m_hash[5];
							break;
						case 'g':
							value = m_hash[6];
							break;
						case 't':
							value = m_hash[7];
							break;
						default:
							printf("Invalid hash entry 2.\n");
							exit(1);
					}
					break;
				case 'c':
					switch(variation[1]) {
						case 'g':
							value = m_hash[8];
							break;
						case 't':
							value = m_hash[9];
							break;
						default:
							printf("Invalid hash entry 3.\n");
							exit(1);
					}
					break;
				case 'g':
					switch(variation[1]) {
						case 't':
							value = m_hash[10];
							break;
						default:
							printf("Invalid hash entry 4.\n");
							exit(1);
					}
					break;
			}
			break;

		// tri-allelic variations
		case 3:
			switch(variation[0]) {
				case 'a':
					if((variation[1] == 'c') && (variation[2] == 'g')) {
						value = m_hash[11];
					} else if((variation[1] == 'c') && (variation[2] == 't')) {
						value = m_hash[12];
					} else if((variation[1] == 'g') && (variation[2] == 't')) {
						value = m_hash[13];
					} else {
						printf("Invalid hash entry 5: %s\n",variation);
						exit(1);
					}

					break;
				case 'c':
					if((variation[1] == 'g') && (variation[2] == 't')) {
						value = m_hash[14];
					} else {
						printf("Invalid hash entry 6.\n");
						exit(1);
					}
					break;
				default:
					printf("Invalid hash entry 7.\n");
					exit(1);
			}
			break;

		// tetra-allelic variations
		case 4:
			value = m_hash[15];
			break;

		// default
		default:
			printf("Invalid hash entry 8.\n");
			exit(1);
	}

	return value;
}
*/
