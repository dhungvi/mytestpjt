#ifndef ZY_AlignmentAlgorithm_HH
#define ZY_AlignmentAlgorithm_HH

#include <string>
#include "IntMatrix.h"

struct AlignResult {
	int score;
	std::string str1;
	std::string str2;
};

/*
 * This class implements Global and local alignment algorithms
 */
class AlignmentAlgorithm {
private:
	IntMatrix scoreMatrix; //A, C, G, T, N
	int gapPenalty;
	inline int encodeBase(char c) {
	    switch (c) {
	       case 'A' : return 0;
	       case 'C' : return 1;
	       case 'G' : return 2;
	       case 'T' : return 3;
	       case 'N' : return 4;
	       case 'P' : return 5;
	   }
	   return -1;
	}


public:
	AlignmentAlgorithm() {};
	AlignmentAlgorithm(int match, int mismatch, int gap);
	inline void setScoringMatrix(const IntMatrix& m) {scoreMatrix = m;}
	inline void setGapPenalty(int gap) {gapPenalty = gap;}
	int getNWScore(const std::string& s1, const std::string& s2);
	AlignResult getNWAlignment(const std::string& s1, const std::string& s2);
	int getSWScore(const std::string& s1, const std::string& s2);
	AlignResult getSWAlignment(const std::string& s1, const std::string& s2);
	std::vector<int> encodeString(const std::string& s1);
};



#endif
