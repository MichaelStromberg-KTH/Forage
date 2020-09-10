#ifndef RUNTIMEPARAMTERS_H_SEEN
#define RUNTIMEPARAMTERS_H_SEEN

#define FILENAME_BUFFER_SIZE 255

#ifdef _WIN32
#define OS_SLASH '\\'
#else
#define OS_SLASH '/'
#endif

class CRuntimeParameters {
public:
	// adjusts the level of output
	static bool IsVerbose;
	// enables base filtration 
	static bool EnableBaseFiltration;
	// enables xml output 
	static bool EnableXmlOutput;
	// specifies the base filtration threshold 
	static unsigned char BaseFilterThreshold;
	// specifies the minimum number of sequences
	static unsigned int MinSequenceThreshold;
	// specifies the XML output filename
	static char XmlFilename[];
	// specifies the polymorphism rate
	static double P_POLY;
};

#endif /* RUNTIMEPARAMTERS_H_SEEN */
