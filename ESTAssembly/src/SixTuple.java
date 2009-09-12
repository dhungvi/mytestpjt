
public class SixTuple {
	int leftNode;	//index of node on the left
	int lOvlLen;	//the overlap length with + or -, now it's -
	int lDis; 		//the distance with + or -, now it's -
	int rightNode;	//the index of node on the right
	int rOvlLen;	//the overlap length with + or -, now it's +
	int rDis;		//the distance with + or -, now it's +
	
	public SixTuple() {
	}
	
	public SixTuple(int i1, int i2, int i3, int i4, int i5, int i6) {
		leftNode = i1;
		lOvlLen = i2;
		lDis = i3;
		rightNode = i4;
		rOvlLen = i5;
		rDis = i6;
	}
}
