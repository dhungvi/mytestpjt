#ifndef ZY_SingleBase_HH
#define ZY_SingleBase_HH

/*
 * record one base in the consensus
 */
class SingleBase {
public:
	char curBase;
	int curBaseNum;

	int numA;
	int numC;
	int numG;
	int numT;
	int numDash;

	SingleBase(char c);

	SingleBase(char c1, char c2);

	void addOneBase(char base);

	char getCurBase();

private:
	void changeBaseNum(char base);
	void changeCurBase();
};

#endif
