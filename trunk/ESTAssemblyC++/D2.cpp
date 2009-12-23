#include "D2.h"
#include <math.h>
#include "Param.h"

using namespace std;

#define MAX_VAL 40000

D2::D2() {
	this->windowSize = WINDONW_SIZE;
	this->THRESHOLD = D2_THRESHOLD;

	// d2 parameters
	this->d2WordSize = D2_WORD_SIZE;
	d2WordFilter = (1 << (d2WordSize << 1)) - 1; // 2^(2*d2WordSize) - 1
	d2NumWords = 1 << (d2WordSize << 1); // 2^(2*d2WordSize)

	// heuristic parameters
	this->heuristicWordSize = HEURISTIC_WORD_SIZE;
	heuristicWordFilter = (1 << (heuristicWordSize << 1)) - 1;
	heuristicNumWords = 1 << (heuristicWordSize << 1);

	this->u = U;
	this->uv_skip = UV_SKIP;
	this->t = T;
	this->tv_max = TV_MAX;
}

int D2::getWindowSize() {
	return windowSize;
}

int D2::getd2WordSize() {
	return d2WordSize;
}

int D2::getHeuristicWordSize() {
	return heuristicWordSize;
}

vector<int> D2::createWindowHash(string s, int leftCoord, int windowSize, int wordSize, int wordFilter, int numWords) {
	vector<int> H(numWords);

	int currentWordCode = 0;
	int currentWordSize = 0;
	for (int i=0; i < windowSize; i++) {
		int c = encodeBase(s[i + leftCoord]);
		if (c < 0) {
			currentWordCode = 0;
			currentWordSize = 0;
		}
		else {
			currentWordCode = ((currentWordCode << 2) | c) & wordFilter;
			currentWordSize = min(currentWordSize+1, wordSize);
		}
		if (currentWordSize == wordSize) {
			H[currentWordCode]++;
		}
	}
	return H;
}


BestWindowMatches D2::matchEndWindows(string s1, string s2) {
	if (!uv_tv_Heuristic(s1, s2))
		return BestWindowMatches(vector<int>(0), 0, 0, vector<int>(0), 0, 0);

	vector<int> H1_left = createWindowHash(s1, 0, windowSize, d2WordSize, d2WordFilter, d2NumWords);
	vector<int> H1_right = createWindowHash(s1, s1.length() - windowSize, windowSize, d2WordSize, d2WordFilter, d2NumWords);
	vector<int> H2 = createWindowHash(s2, 0, windowSize, d2WordSize, d2WordFilter, d2NumWords);

	int d2_left = 0;
	int d2_right = 0;
	for (int i=0; i < d2NumWords; i++) {
		if (H1_left[i] != H2[i]) {
		  d2_left += (int)pow(H1_left[i] - H2[i], 2);
		}
		if (H1_right[i] != H2[i]) {
		  d2_right += (int)pow(H1_right[i] - H2[i], 2);
		}
	}

	vector<int> bestLeftWindow(s2.length() - windowSize +1);
	bestLeftWindow[0] = 0;
	int numBestLeft = 1;
	int bestLeftScore = d2_left;

	vector<int> bestRightWindow(s2.length() - windowSize + 1);
	bestRightWindow[0] = 0;
	int numBestRight = 1;
	int bestRightScore = d2_right;

	for (int i=0; i < s2.length() - windowSize; i++) {
		int firstWord = encodeWord(s2, i, d2WordSize, d2WordFilter);
		int lastWord = encodeWord(s2, i + windowSize - d2WordSize + 1, d2WordSize, d2WordFilter);

		if (firstWord != lastWord) {
			if (firstWord >= 0 && lastWord >= 0) {
				// This is what the adjustment to d2 should be:
				//d2 = d2 - (int)pow(H1_left[firstWord] - H2[firstWord], 2) + (int)pow(H1_left[firstWord] - (H2[firstWord]-1), 2) -
				//(int)pow(H1_left[lastWord] - H2[lastWord], 2) + (int)pow(H1_left[lastWord] - (H2[lastWord]+1), 2);
				// Hazelhurst proves the following is equivilent, IF I am understanding correctly -- must be checked.
				d2_left += (((H2[lastWord] - H1_left[lastWord]) - (H2[firstWord] - H1_left[firstWord]) + 1) << 1);
				d2_right += (((H2[lastWord] - H1_right[lastWord]) - (H2[firstWord] - H1_right[firstWord]) + 1) << 1);
				H2[firstWord]--;
				H2[lastWord]++;
			}
			else {
				if (firstWord >= 0) {
					d2_left = (int) (d2_left - pow(H1_left[firstWord] - H2[firstWord], 2) + pow(H1_left[firstWord] - (H2[firstWord]-1), 2));
					d2_right = (int) (d2_right - pow(H1_right[firstWord] - H2[firstWord], 2) + pow(H1_right[firstWord] - (H2[firstWord]-1), 2));
					H2[firstWord]--;
				}
				if (lastWord >= 0) {
					d2_left = (int) (d2_left - pow(H1_left[lastWord] - H2[lastWord], 2) + pow(H1_left[lastWord] - (H2[lastWord]+1), 2));
					d2_right = (int) (d2_right - pow(H1_right[lastWord] - H2[lastWord], 2) + pow(H1_right[lastWord] - (H2[lastWord]+1), 2));
					H2[lastWord]++;
				}
			}

		}
		if (d2_left == bestLeftScore) {
			bestLeftWindow[numBestLeft++] = i+1;
		}
		else if (d2_left < bestLeftScore) {
			bestLeftScore = d2_left;
			bestLeftWindow[0] = i+1;
			numBestLeft = 1;
		}

		if (d2_right == bestRightScore) {
			bestRightWindow[numBestRight++] = i+1;
		}
		else if (d2_right < bestRightScore) {
			bestRightScore = d2_right;
			bestRightWindow[0] = i+1;
			numBestRight = 1;
		}
	}

	if (bestLeftScore > THRESHOLD) {
		bestLeftWindow = vector<int>(0);
		numBestLeft = 0;
		bestLeftScore = MAX_VAL;
	}
	if (d2_right > THRESHOLD) {
		bestRightWindow = vector<int>(0);
		numBestRight = 0;
		bestRightScore = MAX_VAL;
	}
	return BestWindowMatches(bestLeftWindow, numBestLeft, bestLeftScore, bestRightWindow, numBestRight, bestRightScore);
}


bool D2::uv_tv_Heuristic(string s1, string s2) {

	// The u/v heuristic
	// Look at every (uv_skip) word on s1 and count the number of instances on s2.
	// Return false if the value is less than u.
	vector<int> H = createWindowHash(s2, 0, s2.length(), heuristicWordSize, heuristicWordFilter, heuristicNumWords);
	int total = 0;
	for (int i=0; total < u && i <= s1.length() - heuristicWordSize; i += uv_skip) {
		int code = encodeWord(s1, i, heuristicWordSize, heuristicWordFilter);
		if (code >= 0) {
			total += H[code];
		}
	}
	if (total < u)
		return false;

	// the t/v heursitc
	// Must find at least t words on s2 that occur within 100 bases of eachother on s1.
	vector<int> arr(tv_max, 0);
	int current_position = 0;
	int current_code = encodeWord(s1, 0, heuristicWordSize, heuristicWordFilter);
	while (current_code < 0 && current_position <=  s1.length() - heuristicWordSize) {
		current_position += -current_code;
		current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);  // Shift over past the N
	}
	if (current_code >= 0) {
		total = H[current_code];
		arr[current_position % tv_max] = total;
	}
	else
		total = 0;

	while (current_position < s1.length() - heuristicWordSize ) {
		if (total >= t)
			return true;

		int next_char = encodeBase(s1[current_position+heuristicWordSize]);
		if (next_char >= 0)  {
			current_code = ((current_code << 2) | next_char) & heuristicWordFilter;
			current_position++;
		}
		else {
			current_position += 2;
			current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);
			while (current_code < 0) {
				current_position += -current_code;
				if (current_position > s1.length() - heuristicWordSize)
					break;
				current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);
			}
		}

		int current_index = current_position % tv_max;
		int score = current_code >= 0 ? H[current_code] : 0;
		total = total - arr[current_index] + score;
		arr[current_index] = score;
	}
	return false;
}



