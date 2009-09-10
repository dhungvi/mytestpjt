import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Stack;

import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedEdgeIterator;

import neobio.alignment.BasicScoringScheme;
import neobio.alignment.IncompatibleScoringSchemeException;
import neobio.alignment.NeedlemanWunsch;
import neobio.alignment.SmithWaterman;

/**
 * This class implements the D2 algorithm. Feb 15, 2009.
 * 
 * Note: we need change the method 'initWords' to generate more words if boundOfWord changes.
 */

public class D2 {
	protected final int INT_MAX = 2147483647;
	protected int windowSize;	// the size of window
	protected final int boundOfWord; 	// the upper bound and the lower bound have the same value
	protected int THRESHOLD;	// THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
						// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
	protected int THRESHOLD_OVL;	// THRESHOLD for overlap distance
	protected int InclusionThreshold;	// use this value to define overlap distance of two inclusion subsequence.
						// this value is used in getOVLDistance for judging inclusion(s1 includes s2, or versa).
						// When there is no error in est, we can set it to be zero;
						// When error occur, if the average overlap length is len, we can set it to be (1-(len-4)/len)*100; here 4 means 
						//	allowing 2 different bases in two ests. That is, if the different bases<=2, we assume them to be inclusion.
	protected int alignmentThreshold; //It is used in NewD2 class. It's the threshold for alignment. That is, all the alignment with
											// the distance which is bigger than the value will be seen as infinity. 
	private char[] alphabet = new char[]{'A', 'T', 'C', 'G'};	// alphabet of the two compared strings

	//private String[] words;	// store all the possible words.
	class ValObj {
		int v1, v2;
	}
	private Map<String, ValObj> words = new HashMap<String, ValObj>();
	
	public D2(Properties props) {
		windowSize = Integer.parseInt(props.getProperty("windowSize"));
		boundOfWord = Integer.parseInt(props.getProperty("boundOfWord"));
		THRESHOLD = Integer.parseInt(props.getProperty("THRESHOLD"));
		THRESHOLD_OVL = Integer.parseInt(props.getProperty("boundOfWord"));
		InclusionThreshold = Integer.parseInt(props.getProperty("InclusionThreshold"));
		alignmentThreshold = Integer.parseInt(props.getProperty("alignmentThreshold"));
		initWords();	// initialize the variable 'words'
	}
	
	public int getWindowSize() {
		return windowSize;
	}
	
	/**
	 * Get the window position in string 's2' such that the d2 distance between the two input string is minimal. 
	 * If multiple positions are found, the function will return all of them.
	 * If no position is found, the function will return an empty array.	 
	 * 
	 * @param s1 String the first string with the length of windonwSize, s2 String the second string.
	 * @return an Object array.
	 */
	public Object[] getD2Sed(String s1, String s2) {
		ArrayList<Integer> aPos = new ArrayList<Integer> ();
		
		int minSed = INT_MAX;	//minSed is initialized to the maximum value of int
		/*int[] v1 = getFreqInWindow(s1);
		int[] v2 = getFreqInWindow(s2.substring(0, windowSize));
		*/
		int initSed = getInitSed(s1, s2.substring(0, windowSize));
		
		int l = s2.length() - windowSize + 1;
		int[] disArray = new int[l];	//store value of d2Dis for all the windows in s2.

		/*int sed = 0;
		for (int j=0; j<v2.length; j++) {
			for (int k=0; k<words.length; k++) {
				sed = sed + (v1[k] - v2[k]) * (v1[k] - v2[k]);
			}
		}*/
		
		if (initSed > THRESHOLD) {
			disArray[0] = INT_MAX;
		} else {
			disArray[0] = initSed;
		}
		if (initSed < minSed) {
			minSed = initSed;
		} 
		
		for (int i=1; i<l; i++) {
			String firstWord = s2.substring(i-1, i-1+boundOfWord);
			String lastWord = s2.substring(i+windowSize-boundOfWord, i+windowSize);
			/*for (int j=0; j<words.length; j++) {
				if (words[j].equalsIgnoreCase(firstWord)) {
					v2[j]--;
				}
				
				if (words[j].equalsIgnoreCase(lastWord)) {
					v2[j]++;
				}
			}
			for (int j=0; j<v2.length; j++) {
				for (int k=0; k<words.length; k++) {
					sed = sed + (v1[k] - v2[k])* (v1[k] - v2[k]);
				}
			}*/
			
			int sed = initSed;
			
			ValObj vo = words.get(firstWord);
			if (vo == null) {
				System.out.println("firstWord = " + firstWord + " includes unexpected character!");
			} else {
				int orgFirstValue = (int)Math.pow((vo.v1 - vo.v2), 2);
				vo.v2 -= 1;
				sed = sed - orgFirstValue + (int)Math.pow((vo.v1 - vo.v2), 2);
			}
			
			
			vo = words.get(lastWord);
			if (vo == null) {
				System.out.println("lastWord = " + lastWord + " includes unexpected character!");
			} else {
				int orgLastValue = (int)Math.pow((vo.v1 - vo.v2), 2);
				vo.v2 += 1;
				sed = sed - orgLastValue + (int)Math.pow((vo.v1 - vo.v2), 2);
			}
			
			if (sed > THRESHOLD) {
				disArray[i] = INT_MAX;
			}else {
				disArray[i] = sed;
			}
			if (sed < minSed) {
				minSed = sed;
			} 
			initSed = sed;
		}
		
		// get positions with the value of minSed
		if (minSed < INT_MAX) {
			for (int i=0; i<disArray.length; i++) {
				if (minSed == disArray[i]) {
					aPos.add(Integer.valueOf(i));
				}
			}
		}

		Object[] positions = aPos.toArray();
		return positions;
	}

	/**
	 * Get the overlap distance of the two strings. 
	 * This function 
	 * 	1) tries to find the position with the minimal d2 distance;
	 * 	2) uses alignment to get the similarity value of two substrings. And it sets 
	 * 		the similarity value to be the overlap distance of s1 and s2.
	 * 
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
	protected int[] getOVLDistance(String tS1, String tS2) {
		String s1 = "";
		String s2 = "";
		int flag = 1;	//1 - no switch for tS1 and tS2; -1 - switch.
		/*
		 * put the shorter string to s1 and the longer one to s2 in 
		 * order to identify inclusion. Now we just need to identify the
		 * situation when s1 is included in s2. 
		 */
		if (tS1.length() > tS2.length()) {
			s1 = tS2;
			s2 = tS1;
			flag = -1; //tS1 and tS2 are switched
		} else {
			s1 = tS1;
			s2 = tS2;
		}
		
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
				tmpDis = this.getDistance(s1.substring(0, s1.length()), s2.substring(lPos, lPos+s1.length()));
				tLenOverlap = s1.length();
			} else {
				tmpDis = this.getDistance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if (tmpDis < disLeft){ // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
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
				tmpDis = this.getDistance(s1.substring(0), s2.substring(tLenOverlap-s1.length(), tLenOverlap));
				tLenOverlap = s1.length();
			} else {
				tmpDis = this.getDistance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
			}
			if (tmpDis < disRight) {// && (tLenOverlap > rLenOverlap). do we need to use two conditions or just one?
				disRight = tmpDis;
				rLenOverlap = tLenOverlap;
			}
 		}

		
		// compare disLeft and disRight, select the one with smaller value.
		if (disLeft < disRight) {
			ovlDis = -1 * disLeft * flag;	//minus represents that s2 is to the left of s1
			lenOverlap = -1 * lLenOverlap * flag;
		} else {
			ovlDis = disRight * flag;	//s2 is to the right of s1
			lenOverlap = rLenOverlap * flag;
		} 
			
		/*if s1 is included in s2, we will ignore s2 by assigning INT_MAX to the overlap distance.
		 * We do not need to consider that s1 includes s2 because we have switched them at the beginning 
		 * of this function if s1 is longer than s2.
		 * 
		 * Note: if ovlWindowSize > length of s1, disRight and disLeft would not be equal to zero.
		 */
		if ((Math.abs(lenOverlap) == s1.length()) && 
				(s1.length() <= s2.length()) && 
				((disRight <= InclusionThreshold) || (Math.abs(disLeft) <= InclusionThreshold))) {
			lenOverlap = 0;
			ovlDis = INT_MAX;
		}
		
		returnValues[0] = lenOverlap;
		returnValues[1] = ovlDis;
		return returnValues;
	}

	/*
	 * judge if s1 is included in s2
	 * @return true or false
	 */
	protected boolean checkInclusion(String s1, String s2) {
		String leftWindow = s1.substring(0, windowSize);
		String rightWindow = s1.substring(s1.length()-windowSize);
		Object[] leftPos = getWindowPos(leftWindow, s2);
		Object[] rightPos = getWindowPos(rightWindow, s2);
		int disLeft = INT_MAX;
		int disRight = INT_MAX;
 		int lenOverlap = 0;
 		int lLenOverlap = 0;
 		int rLenOverlap = 0;
  		
		// if all leftPos[i] are -1, disLeft will be kept to be INT_MAX.
 		for (int i=0; i<leftPos.length; i++) {
 			int lPos = ((Integer)(leftPos[i])).intValue();
 			int tLenOverlap = s2.length() - lPos; 
	 		int tmpDis = INT_MAX;
			if (tLenOverlap > s1.length()) {	//if s1 is included in s2
				tmpDis = this.getDistance(s1.substring(0, s1.length()), s2.substring(lPos, lPos+s1.length()));
				tLenOverlap = s1.length();
			} else {
				tmpDis = this.getDistance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if (tmpDis < disLeft){ // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
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
				tmpDis = this.getDistance(s1.substring(0), s2.substring(tLenOverlap-s1.length(), tLenOverlap));
				tLenOverlap = s1.length();
			} else {
				tmpDis = this.getDistance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
			}
			if (tmpDis < disRight) {// && (tLenOverlap > rLenOverlap). do we need to use two conditions or just one?
				disRight = tmpDis;
				rLenOverlap = tLenOverlap;
			}
 		}

		
		// compare disLeft and disRight, select the one with smaller value.
		if (disLeft < disRight) {
			lenOverlap = -1 * lLenOverlap;
		} else {
			lenOverlap = rLenOverlap;
		} 
			
		/*if s1 is included in s2, return true, else return false
		 */
		if ((Math.abs(lenOverlap) == s1.length()) && 
				(s1.length() <= s2.length()) && 
				((disRight <= InclusionThreshold) || (Math.abs(disLeft) <= InclusionThreshold))) {
			return true;
		} else {
			return false;
		}
	}

	/**
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
		return getD2Sed(w, s2);
	}

	/*
	 * Calculate the distance of two strings.
	 * dis = (1 - similarityScore/lengthOfLongerString)*a, actually, in our case, s1 has the same length as s2. 
	 * Now, we set a=100. So the return value would be [0, 100]
	 * @param s1, s2
	 * @return int distance.
	 */
	public int getDistance(String s1, String s2) {
		int score = getSimlarityScore(s1, s2);
		//int score = getLocalSimlarityScore(s1, s2);
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
		int gap = -2;
		NeedlemanWunsch algorithm = new NeedlemanWunsch();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		
		int score = INT_MAX;
		try {
			score = algorithm.getScore();
			//System.out.println(algorithm.getPairwiseAlignment().getGappedSequence1());
			//System.out.println(algorithm.getPairwiseAlignment().getGappedSequence2());
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

	/*
	 * Use Needleman-Wunsch algorithm to get global alignment.
	 * @param s1, s2
	 * @return string[], [0] and [1] are the two aligned sequences, [2] is the pairwise alignment.
	 */
	public String[] getGlobalAlignment(String s1, String s2) {
		int match = 1;
		int mismatch = -1;
		int gap = -2;
		NeedlemanWunsch algorithm = new NeedlemanWunsch();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		String[] strs = new String[3];
		
		try {
			strs[0] = algorithm.getPairwiseAlignment().getGappedSequence1();
			strs[1] = algorithm.getPairwiseAlignment().getGappedSequence2();
			strs[2] = algorithm.getPairwiseAlignment().toString();
			System.out.println(algorithm.getPairwiseAlignment());
		} catch (IncompatibleScoringSchemeException e) {
			// TODO Auto-generated catch block
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return strs;
	}

	/*
	 * Use Smith-Waterman algorithm to calculate similarity score of two string.
	 * @param s1, s2
	 * @return int similarity score(>=0), if the value is less than 0, it's set to be zero.
	 */
	public int getLocalSimlarityScore(String s1, String s2) {
		int match = 1;
		int mismatch = -1;
		int gap = -2;
		SmithWaterman algorithm = new SmithWaterman();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		
		int score = INT_MAX;
		try {
			score = algorithm.getScore();
			//System.out.println(algorithm.getPairwiseAlignment().getGappedSequence1());
			//System.out.println(algorithm.getPairwiseAlignment().getGappedSequence2());
			//System.out.println(algorithm.getPairwiseAlignment());
			//String tmp = algorithm.getPairwiseAlignment().toString();
			//System.out.println(tmp);
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

	/*
	 * Use Smith-Waterman algorithm to calculate A score of two string.
	 * @param s1, s2; s1-provided sequence, s2-assembled sequence
	 * @return int A-score.
	 * A-score = (2*sequence length) - (15*no.of insertions) -
     *     (15*no.of deletions) - (5*no. of substitutions); 
     *     here we set sequence length=length of s2
	 */
/*	public int getAScore(String s1, String s2) {
		int match = 1;
		int mismatch = -1;
		int gap = -1;
		SmithWaterman algorithm = new SmithWaterman();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		
		int score = INT_MAX;
		try {
			//System.out.println(algorithm.getPairwiseAlignment());

			//score1=len of same - no. of gap - no. of mismatch
			int score1 = algorithm.getScore(); 
			
			match = 3;
			mismatch = -2;
			gap = -2;
			scoring = new BasicScoringScheme(match, mismatch, gap);
			algorithm.setScoringScheme(scoring);
			//System.out.println(algorithm.getPairwiseAlignment());

			//score2=(3*len of same) - (2*no. of gap) - (2*no. of mismatch), and this two scoring should 
			// come from the same alignment sequence because in every case the penalty of gap and mismatch
			// are same.
			String seq1 = algorithm.getPairwiseAlignment().getGappedSequence1();
			String seq2 = algorithm.getPairwiseAlignment().getGappedSequence2();
			int score2 = algorithm.getScore();
			
			int lenOfSame = score2 - 2*score1;
			int nGap = 0;
			if (seq1.indexOf("-") != -1) {
				nGap = seq1.split("-").length - 1;
			}
			if (seq2.indexOf("-") != -1) {
				nGap = (seq2.split("-").length - 1) + nGap;
			}
			
			int nMismatch = lenOfSame - nGap - score1;
			
			//if seq2 is substring of s2, we will consider all the other characters in s2 as mismatch.
			String tSeq2 = seq2.replaceAll("-", "");
			nMismatch = nMismatch + (s2.length() - tSeq2.length());
			
			score = 2*s2.length() - 15*nGap - 5*nMismatch;
		} catch (IncompatibleScoringSchemeException e) {
			// TODO Auto-generated catch block
			System.err.println(e.getMessage());
			System.exit(1);
		}

		return score;
	}
*/
	/*
	 * Use Smith-Waterman algorithm to get local alignment.
	 * @param s1, s2
	 * @return string[], [0] and [1] are the two aligned sequences, [2] is the pairwise alignment.
	 */
	public String[] getLocalAlignment(String s1, String s2) {
		int match = 1;
		int mismatch = -1;
		int gap = -2;
		SmithWaterman algorithm = new SmithWaterman();
		BasicScoringScheme scoring = new BasicScoringScheme(match, mismatch, gap);
		algorithm.setScoringScheme(scoring);
		algorithm.loadSequences(s1, s2);
		String[] strs = new String[3];
		
		try {
			strs[0] = algorithm.getPairwiseAlignment().getGappedSequence1();
			strs[1] = algorithm.getPairwiseAlignment().getGappedSequence2();
			strs[2] = algorithm.getPairwiseAlignment().toString();
			//System.out.println(algorithm.getPairwiseAlignment());
		} catch (IncompatibleScoringSchemeException e) {
			// TODO Auto-generated catch block
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return strs;
	}

	/**
	 * Assign values to the class variable words.
	 * 'words' store all the possible words on the alphabet with length of 'boundOfWord'.
	 */
	private void initWords() {
		int wLen = (int)java.lang.Math.pow(alphabet.length, boundOfWord);
		//words = new String[wLen];
		
		int index = 0;
		int len = alphabet.length;
		for (int i=0; i<len; i++) {
			for (int j=0; j<len; j++) {
				for (int k=0; k<len; k++) {	// if boundOfWord changes, we have to change here.
					for (int k1=0; k1<len; k1++) {
						for (int k2=0; k2<len; k2++) {
							for (int k3=0; k3<len; k3++) {
							StringBuffer tmpWord = new StringBuffer();
							tmpWord.append(alphabet[i]).append(alphabet[j]).append(alphabet[k]).append(alphabet[k1]).append(alphabet[k2]).append(alphabet[k3]);
							//words[index] = tmpWord.toString();
							String key = tmpWord.toString();
							ValObj value = new ValObj();
							words.put(key, value);
							index++;
							}
						}
					}
				}
			}
		}
	}
	


	private int getInitSed(String s1, String s2sub) {
		int sed = 0;
		for (Iterator<String> i = words.keySet().iterator(); i.hasNext(); ) {
			String key = i.next();
			
			int sPosition = 0;
			int num1 = 0, num2 = 0;
			while (sPosition < s1.length()) {
				int tmpIndex = s1.indexOf(key, sPosition);
				if (tmpIndex != -1 ) {
					sPosition = tmpIndex + 1;
					num1++;
				} else {
					break;
				}
			}
			sPosition = 0;
			while (sPosition < s2sub.length()) {
				int tmpIndex = s2sub.indexOf(key, sPosition);
				if (tmpIndex != -1 ) {
					sPosition = tmpIndex + 1;
					num2++;
				} else {
					break;
				}
			}
			//freqArray[i] = num;
			ValObj vo = words.get(key);
			vo.v1 = num1;
			vo.v2 = num2;
			int tempInt = num1 - num2;
			sed += Math.pow(tempInt, 2);
		}
		return sed;
	}


	public static void main(String args[]) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}
		
		D2 d2 = new D2(props);
		String s1 = "TCTCCAGGGGATGGGTAAGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCGTTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGAAGAGCCCCGGACTCGCCTTGTAACAAGGCGCGGGA";
		String s2 = "CGGGCCTTCGTTTTACGAAAACAGGTGCGCGAAGCCTGCAAATTTGACGGGGATGGGTATGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCGTTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGAAGAGCCCCGGACTCGCCTTGTAACAAGGCGCGGGA";
		//System.out.println(d2.getSimlarityScore(s1,s2));
		System.out.println(d2.getLocalAlignment(s1,s2));
		//String s="CGGGCCTTCGTTTTACGAAAACAGGTGC";
		//String s3 = "TTTACGA-AAACA";
		//s = s.replace(s3.replace("-", ""), s3);
		//System.out.println(d2.getGlobalAlignment(s1, s2));
		//System.out.println(d2.getAScore(s1, s2));
		//System.out.println(d2.checkInclusion(s1, s1));
		//int offset = s1.indexOf(s2);
		//System.out.println("offset=" + offset);
		//double dis = d2.getLocalSimlarityScore(s1,s2);
		//if ((dis/s2.length()) < 0.95) { //if tStr is not included in retStr, attach it to retStr.
		//	System.out.println("dis=" + dis + "; len=" + s2.length());
		//}
		//d2.getOVLDistance(s1, s2);
		//System.out.println(d2.checkInclusion(s1, s2));
	}
	
	//only used for test by main
	protected static Properties getProperties(String fName) throws IOException {
		Properties props = new Properties();
		File f = new File(fName);
        
        if (!f.exists()) {
        	return props;
        }
        
        props.load(new FileInputStream(f)); 
        return props;
    }
}
