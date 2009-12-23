#include "Alignment.h"

using namespace std;

Alignment::Alignment() {
	alignmentThreshold = ALIGNMENT_THRESHOLD;
	alignAlgo = AlignmentAlgorithm(1,-1,-2);
}

/*
 * Calculate the distance of two strings.
 * dis = (1 - similarityScore/lengthOfLongerString)*a, actually, in our case, s1 has the same length as s2.
 * Now, we set a=100. So the return value would be [0, 100]
 * @param s1, s2
 * @return int distance.
 */
int Alignment::getDistance(string s1, string s2) {
	int score = getSimlarityScore(s1, s2);
	int retVal = INT_MAX;
	if (score != 0) {
		int length = s1.length();
		retVal = (int) ((1 - (double) score / length) * 100);
	}

	if (retVal > alignmentThreshold) {
		retVal = INT_MAX;
	}
	return retVal;
}

/*
 * Use Needleman-Wunsch algorithm to calculate similarity score of two string.
 * @param s1, s2
 * @return int similarity score(>=0), if the value is less than 0, it's set to be zero.
 */
int Alignment::getSimlarityScore(string s1, string s2) {
	return alignAlgo.getNWScore(s1, s2);
}

/*
 * Use Smith-Waterman algorithm to get local alignment.
 * @param s1, s2
 * @return string[], [0] and [1] are the two aligned sequences, [2] is the pairwise alignment.
 */
AlignResult Alignment::getLocalAlignment(string s1, string s2) {
	return alignAlgo.getSWAlignment(s1, s2);
}

