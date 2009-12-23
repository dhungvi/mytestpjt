#include "SixTuple.h"

SixTuple::SixTuple() {
	this->curNode = 0;
	this->leftNode = 0;
	this->lOvlLen = 0;
	this->lDis = 0;
	this->rightNode = 0;
	this->rOvlLen = 0;
	this->rDis = 0;
	this->isNull = true;
}

SixTuple::SixTuple(int c, int i1, int i2, int i3, int i4, int i5, int i6, bool b) {
	this->curNode = c;
	this->leftNode = i1;
	this->lOvlLen = i2;
	this->lDis = i3;
	this->rightNode = i4;
	this->rOvlLen = i5;
	this->rDis = i6;
	this->isNull = b;
}

SixTuple::SixTuple(int c) {
	this->curNode = c;
	this->isNull = true;
	this->leftNode = 0;
	this->lOvlLen = 0;
	this->lDis = 0;
	this->rightNode = 0;
	this->rOvlLen = 0;
	this->rDis = 0;
}

void SixTuple::setSixTuple(int i1, int i2, int i3, int i4, int i5, int i6, bool b) {
	this->leftNode = i1;
	this->lOvlLen = i2;
	this->lDis = i3;
	this->rightNode = i4;
	this->rOvlLen = i5;
	this->rDis = i6;
	this->isNull = b;
}
