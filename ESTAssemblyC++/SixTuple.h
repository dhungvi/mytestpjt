#ifndef ZY_SixTuple_HH
#define ZY_SixTuple_HH

class SixTuple {
public:
	int curNode;	//index of the current node
	int leftNode;	//index of node on the left
	int lOvlLen;	//the overlap length with + or -, now it's -
	int lDis; 		//the distance with + or -, now it's -
	int rightNode;	//the index of node on the right
	int rOvlLen;	//the overlap length with + or -, now it's +
	int rDis;		//the distance with + or -, now it's +
	bool isNull;

	SixTuple();

	SixTuple(int c, int i1, int i2, int i3, int i4, int i5, int i6, bool b);

	SixTuple(int c);

	void setSixTuple(int i1, int i2, int i3, int i4, int i5, int i6, bool b);
};


#endif
