import java.util.ArrayList;



/**
 * This class implements the D2 algorithm. Feb 15, 2009.
 * In this version, the author does not use any skills to make D2 faster.
 * It is a naive implementation of D2. That is, compare all the windows between
 * s1 and s2, and then get the minimum value of sed.
 * 
 * Note: we need change the method 'initWords' to generate more words if boundOfWord changes.
 */

public class D2 {
	private final int INT_MAX = 2147483647;
	private int windowSize = 10;	// the size of window
	//private int ovlWindowSize = 0;// the window size for overlap distance, it is assigned 
								// in 'getOVLDistance' with the value of overlap length.
	private final int boundOfWord = 5; 	// the upper bound and the lower bound have the same value
	private int THRESHOLD = 36;	// THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
						// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
	private final int THRESHOLD_OVL = 0;	// THRESHOLD for overlap distance
	private final int InclusionThreshold = 0;	// use this value to define overlap distance of two inclusion subsequence.
						// this value is used in getOVLDistance for judging inclusion(s1 includes s2, or versa).
	private char[] alphabet = new char[]{'A', 'T', 'C', 'G'};	// alphabet of the two compared strings
	private int[][] v2Array;//	the array is used to store all the frequency 
							//	of all the possible words for all the windows in s2 (the second string).
	private String[] words;	// store all the possible words.
	
	public D2() {
		initWords();	// initialize the variable 'words'
	}
	
	public int getWindowSize() {
		return windowSize;
	}
	
	/**
	 * Get the d2 distance of the two strings. 
	 * 
	 * @param s1 String the first string, s2 String the second string.
	 * @return >=0 the distance; -1 length of s1 or s2 is less than the window size.
	 */
	public int getD2Distance(String s1, String s2) {
		if ((s1.length() < windowSize) || (s2.length() < windowSize)) {
			System.out.println("The length of input string is less than the window size!");
			return -1;
		}
		
		getVectorsForS2(s2);	// put frequency of all the words on all the windows in string s2 into 'v2Array'
		
		int l = s1.length() - windowSize + 1;
		int[] v1 = new int[words.length];	// store frequency of words on one window of s1
		int minSed = INT_MAX;	//minSed is initialized to the maximum value of int
		
		for (int i=0; i<l; i++) {
			v1 = getFreqInWindow(s1.substring(i, i+windowSize));
			// calculate sed for the windows in s1 and all the windows in s2 and get minimum value of sed.
			for (int j=0; j<v2Array.length; j++) {
				int sed = 0;
				// calculate sed for two windows in s1 and s2 respectively
				for (int k=0; k<words.length; k++) {
					sed = sed + (v1[k] - v2Array[j][k])* (v1[k] - v2Array[j][k]);
				}
				
				if (sed < minSed) {
					minSed = sed;
				}
			}
		}
		
		// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
		if (minSed > THRESHOLD) {
			minSed = INT_MAX;
		}
		return minSed;
	}

	/**
	 * Get the overlap distance of the two strings. 
	 * The function finds the position in s2 which has the given d2 distance between s1's first window and 
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
	public int[] getOVLDistance(String s1, String s2, int d2Dis) {
		int[] returnValues = new int[2];
		String leftWindow = s1.substring(0, windowSize);
		String rightWindow = s1.substring(s1.length()-windowSize);
		//int[] leftPos = getWindowPos(leftWindow, s2, d2Dis);
		//int[] rightPos = getWindowPos(rightWindow, s2, d2Dis);
		Object[] leftPos = getWindowPos(leftWindow, s2, d2Dis);
		Object[] rightPos = getWindowPos(rightWindow, s2, d2Dis);
		int disLeft = INT_MAX;
		int disRight = INT_MAX;
		int ovlDis = INT_MAX;
 		int lenOverlap = 0;
 		int lLenOverlap = 0;
 		int rLenOverlap = 0;
 		int oldWindowSize = windowSize;
 		int oldThreshold = THRESHOLD;
 		
		// if all leftPos[i] are -1, disLeft will be kept to be INT_MAX.
 		for (int i=0; i<leftPos.length; i++) {
 			int lPos = ((Integer)(leftPos[i])).intValue();
 			int tLenOverlap = s2.length() - lPos; 
			windowSize = tLenOverlap;	//use ovlWindowSize to calculate overlap distance
			THRESHOLD = THRESHOLD_OVL;
	 		int tmpDis = INT_MAX;
			if (tLenOverlap > s1.length()) {	//if s1 is included in s2
				windowSize = s1.length();
				tmpDis = this.getD2Distance(s1.substring(0, s1.length()), s2.substring(lPos));
			} else {
				tmpDis = this.getD2Distance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if (tmpDis < disLeft) {
				disLeft = tmpDis;
				lLenOverlap = tLenOverlap;
			}
			if (tmpDis == -1) {
				System.out.println("LenOverlap is: " + tLenOverlap);
				System.out.println("s1 is: " + s1 + "\ns2 is: " + s2);
			}
			windowSize = oldWindowSize;	//recover the windowsize
			THRESHOLD = oldThreshold;
 		}
		
		// if all rightPos[i] are -1, disRight will be kept to be INT_MAX.
 		for (int i=0; i<rightPos.length; i++) {
 			int rPos = ((Integer)(rightPos[i])).intValue();
			int tLenOverlap = rPos + windowSize; 
			int lenInS1 = s1.length()-tLenOverlap;

			windowSize = tLenOverlap;	//use ovlWindowSize to calculate overlap distance
			THRESHOLD = THRESHOLD_OVL;
	 		int tmpDis = INT_MAX;
			if (lenInS1 < 0) {	//if s1 is included in s2
				windowSize = s1.length();
				tmpDis = this.getD2Distance(s1.substring(0), s2.substring(0, tLenOverlap));
			} else {
				tmpDis = this.getD2Distance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
			}
			if (tmpDis < disRight) {
				disRight = tmpDis;
				rLenOverlap = tLenOverlap;
			}
			windowSize = oldWindowSize;	//recover the windowsize
			THRESHOLD = oldThreshold;
 		}

		/*
 		// if both leftPos[0] and leftPos[1] are -1, disLeft will be kept to be INT_MAX.
   		if (leftPos[1] != -1) {	// first consider the right-most position in s2
			lenOverlap = s2.length() - leftPos[1]; 
	 		
			windowSize = lenOverlap;	//use ovlWindowSize to calculate overlap distance
			THRESHOLD = THRESHOLD_OVL;
			if (lenOverlap > s1.length()) {	//if s1 is included in s2
				windowSize = s1.length();
				disLeft = this.getD2Distance(s1.substring(0, s1.length()), s2.substring(leftPos[1]));
			} else {
				disLeft = this.getD2Distance(s1.substring(0, lenOverlap), s2.substring(leftPos[1]));
			}
			windowSize = oldWindowSize;	//recover the windowsize
			THRESHOLD = oldThreshold;
		} else if (leftPos[0] != -1) {
			lenOverlap = s2.length() - leftPos[0]; 

			windowSize = lenOverlap;	//use ovlWindowSize to calculate overlap distance
			THRESHOLD = THRESHOLD_OVL;
			if (lenOverlap > s1.length()) {	//if s1 is included in s2
				windowSize = s1.length();
				disLeft = this.getD2Distance(s1.substring(0, s1.length()), s2.substring(leftPos[0]));
			} else {
				disLeft = this.getD2Distance(s1.substring(0, lenOverlap), s2.substring(leftPos[0]));
			}
			windowSize = oldWindowSize;	//recover the windowsize
			THRESHOLD = oldThreshold;
		}
		
		// if rightPos[0] is -1, disRight will be kept to be INT_MAX.
		if (rightPos[0] != -1) {
			lenOverlap = rightPos[0] + windowSize; 
			int lenInS1 = s1.length()-lenOverlap;

			windowSize = lenOverlap;	//use ovlWindowSize to calculate overlap distance
			THRESHOLD = THRESHOLD_OVL;
			if (lenInS1 < 0) {	//if s1 is included in s2
				windowSize = s1.length();
				disRight = this.getD2Distance(s1.substring(0), s2.substring(0, lenOverlap));
			} else {
				disRight = this.getD2Distance(s1.substring(lenInS1), s2.substring(0, lenOverlap));
			}
			windowSize = oldWindowSize;	//recover the windowsize
			THRESHOLD = oldThreshold;
		} */
		
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
	 * Get the window position in string 's2' such that the d2 distance between the two input string is 'd2Dis'. 
	 * If multiple positions are found, the function will return all of them.
	 * If no position is found, the function will return an empty array.
	 * 
	 * 'w' has the length of windowSize, it only has one window which is itself.
	 * 's2' is a string. The function finds starting position of the window on
	 * 's2' which makes the distance of 'w' and 's2' is 'd2Dis'.
	 * 
	 * @param w String the first string, s2 String the second string, d2Dis int d2 distance.
	 * @return an Integer array with two elements.
	 */
	private Object[] getWindowPos(String w, String s2, int d2Dis) {
		ArrayList<Integer> aPos = new ArrayList<Integer> ();
		//int[] positions = new int[] {-1, -1};
		//int index = 0;	//index of positions
		
		getVectorsForS2(s2);	// put frequency of all the words on all the windows in string s2 into 'v2Array'
		int[] v1 = new int[words.length];	// store frequency of words on one window of s1
		v1 = getFreqInWindow(w);

		// calculate sed for w and all the windows in s2 and find those with the value of d2Dis.
		for (int j=0; j<v2Array.length; j++) {
			int sed = 0;
			// calculate sed for two windows in w and s2
			for (int k=0; k<words.length; k++) {
				sed = sed + (v1[k] - v2Array[j][k])* (v1[k] - v2Array[j][k]);
			}
			
			if (sed == d2Dis) {
				aPos.add(Integer.valueOf(j));
				/*if (index > 1) {	//we only need put two values in the variable 'positions'
					index = 1;
				}
				positions[index++] = j;*/
			}
		}
		Object[] positions = aPos.toArray();
		return positions;
	}

	/**
	 * Assign values to the class variable words.
	 * 'words' store all the possible words on the alphabet with length of 'boundOfWord'.
	 */
	private void initWords() {
		int wLen = (int)java.lang.Math.pow(alphabet.length, boundOfWord);
		words = new String[wLen];
		
		int index = 0;
		int len = alphabet.length;
		for (int i=0; i<len; i++) {
			for (int j=0; j<len; j++) {
				for (int k=0; k<len; k++) {	// if boundOfWord changes, we have to change here.
					for (int k1=0; k1<len; k1++) {
						for (int k2=0; k2<len; k2++) {
							StringBuffer tmpWord = new StringBuffer();
							tmpWord.append(alphabet[i]).append(alphabet[j]).append(alphabet[k]).append(alphabet[k1]).append(alphabet[k2]);
							words[index] = tmpWord.toString();
							index++;
							
						}
					}
				}
			}
		}
	}
	
	/**
	 * Assign values to the class variable 'v2Array'.
	 * 'v2Array' stores all the frequency of words for all the windows in the input string.
	 * 
	 *  @param s String
	 */
	private void getVectorsForS2(String s) {
		initV2(s.length());
		int l = s.length() - windowSize + 1;
		for (int i=0; i<l; i++) {
			v2Array[i] = getFreqInWindow(s.substring(i, i+windowSize));
		}
	}
	
	/**
	 * Initialize the class variable 'v2Array'
	 */
	private void initV2(int len) {
		v2Array = new int[len-windowSize+1][words.length];
	}
	
	/**
	 * Calculate all the frequency of words for the input string. Return the results as an int array.
	 * 
	 * @param s String
	 * @return an array to store the frequency of words
	 */
	private int[] getFreqInWindow(String s) {
		int[] freqArray = new int[words.length];
		for (int i=0; i<words.length; i++) {
			int sPosition = 0;
			int num = 0;
			while (sPosition < s.length()) {
				int tmpIndex = s.indexOf(words[i], sPosition);
				if (tmpIndex != -1 ) {
					sPosition = tmpIndex + 1;
					num++;
				} else {
					break;
				}
			}
			freqArray[i] = num;
		}
		
		return freqArray;
	}
	
	public static void main(String args[]) {
		D2 d= new D2();
		String s1 = "AGATGGTTCCTTTCTTTACCTGCTCTGTAAAGA";
		String s2 = "AAGGAGGTGAACATCGATGACGACGAACACTGT";
		int dis = d.getD2Distance(s1, s2);
		int[] dis1 = d.getOVLDistance(s1, s2, dis);
		System.out.println("d2 distance is: " + dis);
		System.out.println("length of overlap is: " + dis1[0]);
		System.out.println("overlap distance is: " + dis1[1]);
	}
}
