#ifndef ZY_PARAM_HH
#define ZY_PARAM_HH

#include <string>

//Parameters for D2

// The window size for d2
const int WINDONW_SIZE=100;
//Threshold for d2. if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
//	THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
const int D2_THRESHOLD=9025;
//Use this value to define overlap distance of two inclusion subsequence. When there is no error in est, we can set it to be zero;
//  When error occurs, if the average overlap length is len, we can set it to be (1-(len-8)/len)*100; here 8 means  allowing 4 different
//   bases in two ests. That is, if the different bases <= 4, we assume them to be inclusion.
const int INCLUSION_THRESHOLD=22;
//Bound of word for d2 distance. This is the word length we use in d2.
const int D2_WORD_SIZE=6;
//The threshold for alignment. That is, all the alignment with the distance which is bigger than the value will be seen as infinity.
const int ALIGNMENT_THRESHOLD=40;


//Parameters for D2 Heuristics
const int HEURISTIC_WORD_SIZE=8;
const int U=6;
const int UV_SKIP=16;
const int T=65;
const int TV_MAX=100;

//The number of levels we will do in order to verify a left end. "0" means doing until leaves.
const int NUMOFLEVELS=0;

//The possible longest length of est. It's used in Reconstruction to reduce the length of strings doing local alignment
const int COMPARISON_LENGTH=800;

#endif
