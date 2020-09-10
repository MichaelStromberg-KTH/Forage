#include "RuntimeParameters.h"

// declare the static variables
bool CRuntimeParameters::IsVerbose = false;

bool CRuntimeParameters::EnableBaseFiltration = false;
unsigned char CRuntimeParameters::BaseFilterThreshold = 20;

bool CRuntimeParameters::EnableXmlOutput = false;
char CRuntimeParameters::XmlFilename[FILENAME_BUFFER_SIZE+1];

unsigned int CRuntimeParameters::MinSequenceThreshold = 2;

double CRuntimeParameters::P_POLY = 0.003;
