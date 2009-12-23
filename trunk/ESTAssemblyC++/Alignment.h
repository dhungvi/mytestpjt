#ifndef ZY_Alignment_HH
#define ZY_Alignment_HH

#include <limits.h>
#include <string>
#include <vector>
#include "Param.h"
#include "AlignmentAlgorithm.h"

class Alignment {
	int alignmentThreshold; //It's the threshold for alignment. That is, all the alignment with
							// the distance which is bigger than the value will be seen as infinity.

public:
	Alignment();
	int getDistance(std::string s1, std::string s2);
	AlignResult getLocalAlignment(std::string s1, std::string s2);

private:
	int getSimlarityScore(std::string s1, std::string s2);
	AlignmentAlgorithm alignAlgo;

};


#endif
