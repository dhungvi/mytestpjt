import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;


/**
 * This class includes all the methods related to overlap distance. 
 * 
 */

public class OvlDistance {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected final int INT_MIN = Integer.MIN_VALUE;
	protected int windowSize;	// the size of window
	protected int InclusionThreshold;	// use this value to define overlap distance of two inclusion subsequence.
	// this value is used in getOVLDistance for judging inclusion(s1 includes s2, or versa).
	// When there is no error in est, we can set it to be zero;
	// When error occur, if the average overlap length is len, we can set it to be (1-(len-4)/len)*100; here 4 means 
	//	allowing 2 different bases in two ests. That is, if the different bases<=2, we assume them to be inclusion.
	// the distance which is bigger than the value will be seen as infinity. 
	protected D2 d2;
	protected Alignment alignment;


	public OvlDistance(Properties props) {
		windowSize = Integer.parseInt(props.getProperty("windowSize"));
		InclusionThreshold = Integer.parseInt(props.getProperty("InclusionThreshold"));
		d2 = new D2(props);
		alignment = new Alignment(props);
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
	 * 		If s1 and s2 has inclusion, the overlap distance is INT_MIN.
	 * 
	 * @param s1 String the first string, s2 String the second string, d2Dis int d2 distance.
	 * @return the first element is the overlap length, the second is the overlap distance.
	 * If s2 is to the left of s1, the length are negative, the distance is zero or negative.
	 * If s2 is to the right of s1, the length are positive, the distance is zero or positive.
	 * If no overlap is found, the distance is INT_MAX.
	 * If s2 is included in s1, the distance is INT_MIN.
	 * If s1 is included in s2, the distance is INT_MIN.
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
		
		BestWindowMatches best = d2.matchEndWindows(s1, s2);
		int[] leftPos = best.bestLeftStart;
		int[] rightPos = best.bestRightStart;
		
		int disLeft = INT_MAX;
		int disRight = INT_MAX;
		int ovlDis = INT_MAX;
		int lenOverlap = 0;
		int lLenOverlap = 0;
		int rLenOverlap = 0;

		for (int i=0; i<leftPos.length; i++) {
			int lPos = leftPos[i];
			int tLenOverlap = s2.length() - lPos; 
			int tmpDis = INT_MAX;
			if (tLenOverlap > s1.length()) {	//if s1 is included in s2
				tmpDis = alignment.getDistance(s1.substring(0, s1.length()), s2.substring(lPos, lPos+s1.length()));
				tLenOverlap = s1.length();
			} else {
				tmpDis = alignment.getDistance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if (tmpDis < disLeft){ // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
				disLeft = tmpDis;
				lLenOverlap = tLenOverlap;
			}
		}

		for (int i=0; i<rightPos.length; i++) {
			int rPos = rightPos[i];
			int tLenOverlap = rPos + windowSize; 
			int lenInS1 = s1.length()-tLenOverlap;

			int tmpDis = INT_MAX;
			if (lenInS1 < 0) {	//if s1 is included in s2
				tmpDis = alignment.getDistance(s1.substring(0), s2.substring(tLenOverlap-s1.length(), tLenOverlap));
				tLenOverlap = s1.length();
			} else {
				tmpDis = alignment.getDistance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
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
			lenOverlap = Math.min(s1.length(), s2.length());
			ovlDis = INT_MIN;
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
		BestWindowMatches best = d2.matchEndWindows(s1, s2);
		int[] leftPos = best.bestLeftStart;
		int[] rightPos = best.bestRightStart;
		int disLeft = INT_MAX;
		int disRight = INT_MAX;
		int lenOverlap = 0;
		int lLenOverlap = 0;
		int rLenOverlap = 0;

		// if all leftPos[i] are -1, disLeft will be kept to be INT_MAX.
		for (int i=0; i<leftPos.length; i++) {
			int lPos = leftPos[i];
			int tLenOverlap = s2.length() - lPos; 
			int tmpDis = INT_MAX;
			if (tLenOverlap > s1.length()) {	//if s1 is included in s2
				tmpDis = alignment.getDistance(s1.substring(0, s1.length()), s2.substring(lPos, lPos+s1.length()));
				tLenOverlap = s1.length();
			} else {
				tmpDis = alignment.getDistance(s1.substring(0, tLenOverlap), s2.substring(lPos));
			}
			if (tmpDis < disLeft){ // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
				disLeft = tmpDis;
				lLenOverlap = tLenOverlap;
			}
		}

		// if all rightPos[i] are -1, disRight will be kept to be INT_MAX.
		for (int i=0; i<rightPos.length; i++) {
			int rPos = rightPos[i];
			int tLenOverlap = rPos + windowSize; 
			int lenInS1 = s1.length()-tLenOverlap;

			int tmpDis = INT_MAX;
			if (lenInS1 < 0) {	//if s1 is included in s2
				tmpDis = alignment.getDistance(s1.substring(0), s2.substring(tLenOverlap-s1.length(), tLenOverlap));
				tLenOverlap = s1.length();
			} else {
				tmpDis = alignment.getDistance(s1.substring(lenInS1), s2.substring(0, tLenOverlap));
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


	public static void main(String args[]) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
			return;
		}

		OvlDistance ovl = new OvlDistance(props);
		String s1 = "GCAacGCTgTGCAGGAATCTAgTACGCTTTGAGGTGCAAAACAGCCCTTCAtCGCCGGAGATGCCGCTATgGGGAGCGGAGACGCCCAGCANCCTGGGGCCGTCATATCTTGTGgNCCCATGGCCCGCTATGGAGCGGACTAaCTTNCAGATTACcAGCAGCAaGCCTGGTGCtGCgGCGGACGGCAGCaTCGGCGGCTgGGGcgATGGtGGCGaGCAAGGTCGCTgACAGTCGCTGCCAGAACCGCCAGTCATTAAGCAcCAGAAAcTATCCGCTCGTGGGaAACGcAGTGCCCCAGCAGCAGTTATcNGTGcAGCTCtCGcCACGCAGGGAGCGCCGGAtaCCAGCCCCaGCAggTACAGCTATAAATTTTCAGCAACAGAtTCTCCcTTTtGGGCTTTCTtTNaTTTGTTGCTTACCGtCTAGTCGCCGcCTGGAAGcCTGCGCTAATCGGCCTGACACCtTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGAGGcN";
		String s2 = "CAcCCCTaCAgCGCCGGcGtNGCCGgTATGGGGAGCGGAGcCGCCCAGCAGCCTGGGGCCGacATATCTTGTGCCCCCATGgtCCCGCTATGGAGCgGGACACCTTACAGATTACAGCAGCATGtTCTGGTGCTGCCGCGtACGGCAGCCTCGGCGGCTGGGGTCNTGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAtATTATCCGCTCGTGcGAACGCAGTGCCCCAcCaAGCaGTTATCCGTGGAGCTCGCGCACGNAGGGtAGCGCCGGAGGCCAGCCCCAGCaCCTACAGCTATAAATTTTCAGCAAaAaATTCTgCATTTTGGGCTTCTTTGGTTTGtTGCTTACCgcCTAGtCGCCgTCTGGAAGGCTGCGCTAATCGGtCTGACACCCTTGGgGCGTTGCCTCTCAtTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCcCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCaGAAGTCTCaCTTCTGCCaTg";
		int[] a = ovl.getOVLDistance(s1, s2);
		for (int i=0; i<a.length; i++) {
			System.out.println(a[i]);
		}
		System.out.println(ovl.checkInclusion(s1, s2));
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


