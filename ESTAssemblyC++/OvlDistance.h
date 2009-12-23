#ifndef ZY_OvlDistance_HH
#define ZY_OvlDistance_HH

#include <vector>
#include <limits.h>   // Holds a variable INT_MAX
#include <string>
#include "Param.h"
#include "D2.h"
#include "Alignment.h"

/**
 * This class includes all the methods related to overlap distance.
 *
 */

class OvlDistance {
	int windowSize;	// the size of window
	int InclusionThreshold;	// use this value to define overlap distance of two inclusion subsequence.
	// this value is used in getOVLDistance for judging inclusion(s1 includes s2, or versa).
	// When there is no error in est, we can set it to be zero;
	// When error occur, if the average overlap length is len, we can set it to be (1-(len-4)/len)*100; here 4 means
	//	allowing 2 different bases in two ests. That is, if the different bases<=2, we assume them to be inclusion.
	// the distance which is bigger than the value will be seen as infinity.
	Alignment alignment;


public:
	D2 d2;
	OvlDistance();
	std::vector<int> getOVLDistance(const std::string& tS1, const std::string& tS2);
	bool checkInclusion(const std::string& s1, const std::string& s2);

private:
	std::vector<int> reducePos(const std::vector<int>& input, int len);
};

#endif
