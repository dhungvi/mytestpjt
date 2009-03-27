import java.util.ArrayList;

import neobio.alignment.BasicScoringScheme;
import neobio.alignment.IncompatibleScoringSchemeException;
import neobio.alignment.NeedlemanWunsch;

public class NewD2 extends D2{
	private int alignmentThreshold = 40;
	
	
	public NewD2() {
		
	}
	
	/**
	 * Overload the method of parent class
	 * Get the overlap distance of the two strings. 
	 * This function is different from the same-name method in its parent class in that 
	 * 	1) It tries to find the position with the minimal d2 distance instead of that with 
	 * 		the same value as the edge;
	 * 	2) It uses alignment to get the similarity value of two substrings instead of using
	 * 		d2. And it sets the similarity value to be the overlap distance of s1 and s2.
	 * The function finds the position in s2 which has the smallest d2 distance between s1's first window and 
	 * s2 or between s1's last window and s2. If it finds the position, it returns the overlap length and the
	 * overlap distance of the two strings. If not, it returns INT_MAX. If it finds two positions for both the 
	 * first and last window in s1, it chooses the one with the smaller distance.
	 * Specifically, if the function finds the position, it returns three kinds of values.
	 * 		If s2 is to the left of s1, the values(overlap length and distance) are negative integer;
	 * 		If s2 is to the right of s1, the values(overlap length and distance) are positive integer;
	 * 		If s2 in included in s1, the overlap distance is INT_MAX.
	 * 
	 * @param s1 String the first string, s2 String the second string, d2Dis int d2 distance.
	 * @return the first element is the overlap length, the second is the overlap distance.
	 * If s2 is to the left of s1, the length are negative, the distance is zero or negative.
	 * If s2 is to the right of s1, the length are positive, the distance is zero or positive.
	 * If no overlap is found, the distance is INT_MAX.
	 * If s2 is included in s1, the distance is INT_MAX.
	 * If s1 is included in s2, the distance is INT_MAX.
	 */
	protected int[] getOVLDistance(String s1, String s2) {
		int[] returnValues = new int[2];
		String leftWindow = s1.substring(0, windowSize);
		String rightWindow = s1.substring(s1.length()-windowSize);
		Object[] leftPos = getWindowPos(leftWindow, s2);
		Object[] rightPos = getWindowPos(rightWindow, s2);
		int disLeft = INT_MAX;
		int disRight = INT_MAX;
		int ovlDis = INT_MAX;
 		int lenOverlap = 0;
 		int lLenOverlap = 0;
 		int rLenOverlap = 0;
  		
		// if all leftPos[i] are -1, disLeft will be kept to be INT_MAX.
 		for (int i=0; i<leftPos.length; i++) {
 			int lPos = ((Integer)(leftPos[i])).intValue();
 			int tLenOverlap = s2.length() - lPos; 
	 		int tmpDis = INT_MAX;
			if (tLenOverlap > s1.length()) {	//if s1 is included in s2
				tmpDis = this.getDistance(s1.substring(0, s1.length()), s2.substring(lPos));
			} else {
				tmpDis = this.getDistance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if ((tmpDis < disLeft) && (tLenOverlap > lLenOverlap)){ // do we need to use two conditions or just one?
				disLeft = tmpDis;
				lLenOverlap = tLenOverlap;
			}
 		}
		
		// if all rightPos[i] are -1, disRight will be kept to be INT_MAX.
 		for (int i=0; i<rightPos.length; i++) {
 			int rPos = ((Integer)(rightPos[i])).intValue();
			int tLenOverlap = rPos + windowSize; 
			int lenInS1 = s1.length()-tLenOverlap;

	 		int tmpDis = INT_MAX;
			if (lenInS1 < 0) {	//if s1 is included in s2
				tmpDis = this.getDistance(s1.substring(0), s2.substring(0, tLenOverlap));
			} else {
				tmpDis = this.getDistance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
			}
			if ((tmpDis < disRight) && (tLenOverlap > rLenOverlap)) {// do we need to use two conditions or just one?
				disRight = tmpDis;
				rLenOverlap = tLenOverlap;
			}
 		}

		
		// compare disLeft and disRight, select the one with smaller value.
		if (disLeft < disRight) {
			ovlDis = -1 * disLeft;	//minus represents that s2 is to the left of s1
			lenOverlap = -1 * lLenOverlap;
		} else {
			ovlDis = disRight;	//s2 is to the right of s1
			lenOverlap = rLenOverlap;
		} 
			
		/*if s2 is included in s1, we will ignore s2 by assigning INT_MAX to the overlap distance.
		 * Here, we only judge two situations: s1: s_1 ... e_1; s2: s'_1 ... e'_1.
		 * 		s_1 = s'_1, or e_1 = e'_1.
		 * If s'_1<s_1 and e'_1>e_1, their ovl distance would be INT_MAX according to above program.
		 */
		if ((Math.abs(lenOverlap) == s2.length()) && 
				(s2.length() <= s1.length()) && 
				((disRight == InclusionThreshold)||	(Math.abs(disLeft) == InclusionThreshold))) {
			lenOverlap = 0;
			ovlDis = INT_MAX;
		}
		
		/*if s1 is included in s2, we will also ignore s2 by assigning INT_MAX to the overlap distance.
		 * 
		 * Note: if ovlWindowSize > length of s1, disRight and disLeft would not be equal to zero.
		 */
		if ((s1.length() <= s2.length()) && 
				(disRight == disLeft) &&
				(disRight == InclusionThreshold)) {
			lenOverlap = 0;
			ovlDis = INT_MAX;
		}
		
		returnValues[0] = lenOverlap;
		returnValues[1] = ovlDis;
		return returnValues;
	}

	/**
	 * Overload the method in parent class.
	 * This function is different from the same-name method in its parent class in that It tries to find the position 
	 * with the minimal d2 distance instead of that with the same value as the edge.
	 * 
	 * Get the window position in string 's2' such that the d2 distance between the two input string is minimal. 
	 * If multiple positions are found, the function will return all of them.
	 * If no position is found, the function will return an empty array.
	 * 
	 * 'w' has the length of windowSize, it only has one window which is itself.
	 * 's2' is a string. The function finds starting position of the window on
	 * 's2' which makes the distance of 'w' and 's2' is 'd2Dis'.
	 * 
	 * @param w String the first string, s2 String the second string, d2Dis int d2 distance.
	 * @return an Object array with two elements.
	 */
	protected Object[] getWindowPos(String w, String s2) {
		ArrayList<Integer> aPos = new ArrayList<Integer> ();
		int minDis = INT_MAX;
		
		int l = s2.length() - w.length() + 1;
		int[] disArray = new int[l];	//store value of d2Dis for all the windows in s2.
		
		for (int i=0; i<l; i++) {
			int dis = super.getD2Distance(w, s2.substring(i, i+w.length()));
			if (dis < minDis) {
				minDis = dis;
			} 
			disArray[i] = dis;
		}

		// get positions with the value of minDis
		if (minDis < INT_MAX) {
			for (int j=0; j<disArray.length; j++) {
				if (minDis == disArray[j]) {
					aPos.add(Integer.valueOf(j));
				}
			}
		}
	
		Object[] positions = aPos.toArray();
		return positions;
	}
	
	/*
	 * calculate the distance of two strings.
	 * dis = (1 - similarityScore/lengthOfLongerString)*a, actually, in our case, s1 has the same length as s2. 
	 * Now, we set a=100. So the return value would be [0, 100]
	 * @param s1, s2
	 * @return int distance.
	 */
	public int getDistance(String s1, String s2) {
		int score = getSimlarityScore(s1, s2);
/*		if (s1.length() > s2.length()) {
			length = s1.length();
		} else {
			length = s2.length();
		}
	*/	
		int retVal = INT_MAX;
		if (score != 0) {
			//int length = s1.length() + s2.length() - score;
			int length = s1.length();
			retVal = (int)((1 - (double)score/length) * 100);
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
	private int getSimlarityScore(String s1, String s2) {
		int match = 1;
		int mismatch = -1;
		int gap = -1;
		NeedlemanWunsch algorithm = new NeedlemanWunsch();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		
		int score = INT_MAX;
		try {
			score = algorithm.getScore();
			//System.out.println(algorithm.getPairwiseAlignment());
		} catch (IncompatibleScoringSchemeException e) {
			// TODO Auto-generated catch block
			System.err.println(e.getMessage());
			System.exit(1);
		}
		if (score < 0) {
			score = 0;
		}
		return score;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		NewD2 d2 = new NewD2();
		String s1 = "CAATCGGGGAACTA";
		String s2 = "CAATCGGGGATACT";
		System.out.println(d2.getSimlarityScore(s1,s2));
		//System.out.println(d2.getDistance(s1,s2));
		int tmp[] = new int[2];
		tmp = d2.getOVLDistance(s1, s2);
		System.out.println(tmp[1]);

	}

}
