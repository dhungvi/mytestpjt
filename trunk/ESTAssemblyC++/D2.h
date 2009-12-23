#ifndef JK_D2_HH
#define JK_D2_HH

#include <vector>
#include <limits.h>   // Holds a variable INT_MAX
#include <string>
#include "Param.h"



class BestWindowMatches {
public:
	inline BestWindowMatches(const std::vector<int> &leftStart, int numBestLeft, int leftD2, const std::vector<int> &rightStart, int numBestRight, int rightD2) {
		bestLeftStart = leftStart;
		numBestLeftWindows = numBestLeft;
		bestLeftD2 = leftD2;

		bestRightStart = rightStart;
		numBestRightWindows = numBestRight;
		bestRightD2 = rightD2;
	}

	std::vector<int> bestLeftStart;
	int numBestLeftWindows;
	int bestLeftD2;

	std::vector<int> bestRightStart;
	int numBestRightWindows;
	int bestRightD2;
};




/**
* This class implements the D2 algorithm.
*
*/
class D2 {
public:
	D2();

	// inspectors
	int getWindowSize();
	int getd2WordSize();
	int getHeuristicWordSize();

	// computation methods
	BestWindowMatches matchEndWindows(std::string s1, std::string s2);
	bool uv_tv_Heuristic(std::string s1, std::string s2);

private:
	int windowSize;	// the size of window
	int THRESHOLD;	// THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
	// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).


	// d2 parameters
	int d2WordSize;
	int d2NumWords;     // number of different words;
	int d2WordFilter;

	// heuristic parameters
	int heuristicWordSize;
	int heuristicNumWords;
	int heuristicWordFilter;
	int u;
	int uv_skip;
	int t;
	int tv_max;

	// helper functions
	inline int encodeBase(char c) {
		switch (c) {
	case 'A' :
	case 'a' : return 0;
	case 'C' :
	case 'c' : return 1;
	case 'G' :
	case 'g' : return 2;
	case 'T' :
	case 't' : return 3;
	case 'n' :
	case 'N' : return -1;
		}
		return -2;
	}

	// Returns the word starting at base leftCoord o
	inline int encodeWord(std::string s, int leftCoord, int wordSize, int wordFilter) {
		int code = 0;
		for (int i=0; i < wordSize; i++) {
			int c = encodeBase(s[i + leftCoord]);
			if (c < 0)
				return -1*i - 1;
			else
				code = ((code << 2) | c) & wordFilter;
		}
		return code;
	}

	std::vector<int> createWindowHash(std::string s, int leftCoord, int windowSize, int wordSize, int wordFilter, int numWords);
};



#endif
