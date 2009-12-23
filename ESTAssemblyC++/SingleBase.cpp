#include "SingleBase.h"
#include <iostream>
using namespace std;

SingleBase::SingleBase(char c) {
	curBase = '\0';
	curBaseNum = 0;
	numA = 0;
	numC = 0;
	numG = 0;
	numT = 0;
	numDash = 0;
	changeBaseNum(c);
	changeCurBase();
}

SingleBase::SingleBase(char c1, char c2) {
	curBase = '\0';
	curBaseNum = 0;
	numA = 0;
	numC = 0;
	numG = 0;
	numT = 0;
	numDash = 0;
	changeBaseNum(c1);
	changeBaseNum(c2);
	changeCurBase();
}

void SingleBase::addOneBase(char base) {
	changeBaseNum(base);
	changeCurBase();
}

char SingleBase::getCurBase() {
	if (curBase == '-') {
		return 'P';
	} else {
		return curBase;
	}
}

void SingleBase::changeBaseNum(char base) {
	switch (base) {
	case 'A':
		numA++;
		break;
	case 'T':
		numT++;
		break;
	case 'C':
		numC++;
		break;
	case 'G':
		numG++;
		break;
	case '-':
		numDash++;
		break;
	}
}

void SingleBase::changeCurBase() {
	int maxNum = numDash;
	char base = '-';

	if (maxNum < numA) {
		base = 'A';
		maxNum = numA;
	}
	if (maxNum < numC) {
		base = 'C';
		maxNum = numC;
	}
	if (maxNum < numG) {
		base = 'G';
		maxNum = numG;
	}
	if (maxNum < numT) {
		base = 'T';
		maxNum = numT;
	}

	curBase = base;
	curBaseNum = maxNum;
}


