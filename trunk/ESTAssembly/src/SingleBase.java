public class SingleBase {
	char curBase;
	int curBaseNum;
	
	int numA;
	int numC;
	int numG;
	int numT;
	int numDash;
	
	public SingleBase(char c) {
		changeBaseNum(c);
		changeCurBase();
	}

	public SingleBase(char c1, char c2) {
		changeBaseNum(c1);
		changeBaseNum(c2);
		changeCurBase();
	}
	
	public void addOneBase(char base) {
		changeBaseNum(base);
		changeCurBase();
	}
	
	public char getCurBase() {
		if (curBase == '-') {
			return 'P';
		} else {
			return curBase;
		}
	}
	
	private void changeBaseNum(char base) {
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
	
	private void changeCurBase() {
		int max = numDash;
		char base = '-';
		
		if (max < numA) {
			base = 'A';
			max = numA;
		}
		if (max < numC) {
			base = 'C';
			max = numC;
		}
		if (max < numG) {
			base = 'G';
			max = numG;
		}
		if (max < numT) {
			base = 'T';
			max = numT;
		}

		curBase = base;
		curBaseNum = max;
	}
}
