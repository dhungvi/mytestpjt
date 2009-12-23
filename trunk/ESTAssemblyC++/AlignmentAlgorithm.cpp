#include <math.h>
#include <algorithm>
#include "AlignmentAlgorithm.h"
using namespace std;

AlignmentAlgorithm::AlignmentAlgorithm(int match, int mismatch, int gap) {
	 IntMatrix score(5,5);
	 score.set(0, 0, match);
	 score.set(0, 1, mismatch);
	 score.set(0, 2, mismatch);
	 score.set(0, 3, mismatch);
	 score.set(0, 4, mismatch);
	 score.set(1, 0, mismatch);
	 score.set(1, 1, match);
	 score.set(1, 2, mismatch);
	 score.set(1, 3, mismatch);
	 score.set(1, 4, mismatch);
	 score.set(2, 0, mismatch);
	 score.set(2, 1, mismatch);
	 score.set(2, 2, match);
	 score.set(2, 3, mismatch);
	 score.set(2, 4, mismatch);
	 score.set(3, 0, mismatch);
	 score.set(3, 1, mismatch);
	 score.set(3, 2, mismatch);
	 score.set(3, 3, match);
	 score.set(3, 4, mismatch);
	 score.set(4, 0, mismatch);
	 score.set(4, 1, mismatch);
	 score.set(4, 2, mismatch);
	 score.set(4, 3, mismatch);
	 score.set(4, 4, match);

	 this->setScoringMatrix(score);
	 this->setGapPenalty(gap);

}

int AlignmentAlgorithm::getNWScore(const string& s1, const string& s2) {
	int numOfRows = s1.length();
	int numOfCols = s2.length();

	IntMatrix alignMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 0; i <= numOfCols; i++) {
		alignMatrix.set(0, i, (this->gapPenalty) * i);
	}
	for (int i = 0; i <= numOfRows; i++) {
		alignMatrix.set(i, 0, (this->gapPenalty) * i);
	}

	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodeBase(s1[i - 1]);
			int base2 = encodeBase(s2[j - 1]);
			int max = alignMatrix.get(i - 1, j) + this->gapPenalty;

			// See if the second term is larger (WEST).
			int west = alignMatrix.get(i, j - 1) + this->gapPenalty;
			if (max <= west) {
				max = west;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix.get(i - 1, j - 1)
					+ (this->scoreMatrix).get(base1, base2);
			if (max <= northwest) {
				max = northwest;
			}

			alignMatrix.set(i, j, max);
		}
	}
	return alignMatrix.get(numOfRows, numOfCols);
}

AlignResult AlignmentAlgorithm::getNWAlignment(const std::string& s1,
		const std::string& s2) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();
	IntMatrix trace(numOfRows + 1, numOfCols + 1);
	trace.set(0, 0, 0);

	IntMatrix alignMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix.set(0, i, (this->gapPenalty) * i);
		trace.set(0, i, 2); //west
	}
	for (int i = 1; i <= numOfRows; i++) {
		alignMatrix.set(i, 0, (this->gapPenalty) * i);
		trace.set(i, 0, 1); //north
	}

	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodeBase(s1[i - 1]);
			int base2 = encodeBase(s2[j - 1]);
			int max = alignMatrix.get(i - 1, j) + this->gapPenalty;

			// See if the second term is larger (WEST).
			int west = alignMatrix.get(i, j - 1) + this->gapPenalty;
			if (max <= west) {
				max = west;
				flag = 2;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix.get(i - 1, j - 1)
					+ (this->scoreMatrix).get(base1, base2);
			if (max <= northwest) {
				max = northwest;
				flag = 3;
			}

			alignMatrix.set(i, j, max);
			trace.set(i, j, flag);
		}
	}

	result.score = alignMatrix.get(numOfRows, numOfCols);

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = numOfRows;
	int col = numOfCols;
	while ((row != 0) || (col != 0)) {
		int flag = trace.get(row, col);
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	string::reverse_iterator rit;
	for (rit = tStr1.rbegin(); rit < tStr1.rend(); rit++)
		result.str1.append(1, *rit);
	for (rit = tStr2.rbegin(); rit < tStr2.rend(); rit++)
		result.str2.append(1, *rit);

	return result;
}

int AlignmentAlgorithm::getSWScore(const string& s1, const string& s2) {
	int numOfRows = s1.length();
	int numOfCols = s2.length();

	IntMatrix alignMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 0; i <= numOfCols; i++) {
		alignMatrix.set(0, i, 0);
	}
	for (int i = 0; i <= numOfRows; i++) {
		alignMatrix.set(i, 0, 0);
	}

	int maxScore = 0;
	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodeBase(s1[i - 1]);
			int base2 = encodeBase(s2[j - 1]);
			int innerMax = max(0, alignMatrix.get(i - 1, j) + this->gapPenalty);

			// See if the second term is larger (WEST).
			int west = alignMatrix.get(i, j - 1) + this->gapPenalty;
			if (innerMax <= west) {
				innerMax = west;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix.get(i - 1, j - 1)
					+ (this->scoreMatrix).get(base1, base2);
			if (innerMax <= northwest) {
				innerMax = northwest;
			}

			alignMatrix.set(i, j, innerMax);

			maxScore = max(maxScore, innerMax);
		}
	}
	return maxScore;
}

AlignResult AlignmentAlgorithm::getSWAlignment(const std::string& s1,
		const std::string& s2) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();
	IntMatrix trace(numOfRows + 1, numOfCols + 1);
	trace.set(0, 0, 0);

	IntMatrix alignMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix.set(0, i, 0);
		trace.set(0, i, 0); //start point
	}
	for (int i = 1; i <= numOfRows; i++) {
		alignMatrix.set(i, 0, 0);
		trace.set(i, 0, 0); //start point
	}

	int maxScore = 0;
	int maxRow = 0;
	int maxCol = 0;
	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodeBase(s1[i - 1]);
			int base2 = encodeBase(s2[j - 1]);
			int max = alignMatrix.get(i - 1, j) + this->gapPenalty;

			// See if the second term is larger (WEST).
			int west = alignMatrix.get(i, j - 1) + this->gapPenalty;
			if (max <= west) {
				max = west;
				flag = 2;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix.get(i - 1, j - 1)
					+ (this->scoreMatrix).get(base1, base2);
			if (max <= northwest) {
				max = northwest;
				flag = 3;
			}

			if (max <= 0) {
				alignMatrix.set(i, j, 0);
				trace.set(i, j, 0); //start point
			} else {
				alignMatrix.set(i, j, max);
				trace.set(i, j, flag);
			}
			if (max > maxScore) {
				maxScore = max;
				maxRow = i;
				maxCol = j;
			}
		}
	}

	result.score = alignMatrix.get(maxRow, maxCol);

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = maxRow;
	int col = maxCol;
	int flag = trace.get(row, col);
	while (flag != 0) {
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
		flag = trace.get(row, col);
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	string::reverse_iterator rit;
	for (rit = tStr1.rbegin(); rit < tStr1.rend(); rit++)
		result.str1.append(1, *rit);
	for (rit = tStr2.rbegin(); rit < tStr2.rend(); rit++)
		result.str2.append(1, *rit);

	return result;
}
