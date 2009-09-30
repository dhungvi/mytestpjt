import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;


/**
 * This class implements the D2 algorithm. 
 * 
 */

public class D2 {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected int windowSize;	// the size of window
	protected final int wordSize; 	// the upper bound and the lower bound have the same value
	protected int THRESHOLD;	// THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
	// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
	protected int numWords;     // number of different words;
	protected int wordFilter;


	public D2(Properties props) {
		windowSize = Integer.parseInt(props.getProperty("windowSize"));
		wordSize = Integer.parseInt(props.getProperty("boundOfWord"));
		THRESHOLD = Integer.parseInt(props.getProperty("THRESHOLD"));
		wordFilter = (1 << 2*wordSize) - 1;
		numWords = 1 << 2*wordSize;
	}

	public int getWindowSize() {
		return windowSize;
	}

	public int getWordSize() {
		return wordSize;
	}

	private int encodeBase(char c) {
		switch (c) {
		case 'A' :
		case 'a' : return 0;
		case 'C' :
		case 'c' : return 1;
		case 'G' :
		case 'g' : return 2;
		case 'T' :
		case 't' : return 3;
		case 'n' :
		case 'N' : return -1;
		}
		return -2;
	}

	private int[] createWindowHash(String s, int leftCoord) {
		int[] H = new int[numWords];

		int currentWordCode = 0;
		int currentWordSize = 0;
		for (int i=0; i < windowSize; i++) {
			int c = encodeBase(s.charAt(i + leftCoord));
			if (c < 0) {
				currentWordCode = 0;
				currentWordSize = 0;
			}
			else {
				currentWordCode = ((currentWordCode << 2) | c) & wordFilter;
				currentWordSize = Math.min(currentWordSize+1, wordSize);
			}
			if (currentWordSize == wordSize) {
				H[currentWordCode]++;
			}
		}	       	   
		return H;
	}


	// Returns the word starting at base leftCoord o
	private int encodeWord(String s, int leftCoord) {
		int code = 0;
		for (int i=0; i < wordSize; i++) {
			int c = encodeBase(s.charAt(i + leftCoord));
			if (c < 0)
				return -1*i - 1;
			else
				code = ((code << 2) | c) & wordFilter;
		}
		return code;
	}

	public BestWindowMatches matchEndWindows(String s1, String s2) {
		int[] H1_left = createWindowHash(s1, 0);
		int[] H1_right = createWindowHash(s1, s1.length() - windowSize);
		int[] H2 = createWindowHash(s2, 0);

		int d2_left = 0;
		int d2_right = 0;
		for (int i=0; i < numWords; i++) {
			if (H1_left[i] != H2[i]) {
				d2_left += Math.pow(H1_left[i] - H2[i], 2);
			}
			if (H1_right[i] != H2[i]) {
				d2_right += Math.pow(H1_right[i] - H2[i], 2);
			}
		}

		int[] bestLeftWindow = new int[s2.length() - windowSize + 1];
		bestLeftWindow[0] = 0;
		int numBestLeft = 1;
		int bestLeftScore = d2_left;

		int[] bestRightWindow = new int[s2.length() - windowSize + 1];
		bestRightWindow[0] = 0;
		int numBestRight = 1;
		int bestRightScore = d2_right;

		for (int i=0; i < s2.length() - windowSize; i++) {
			int firstWord = encodeWord(s2, i);
			int lastWord = encodeWord(s2, i + windowSize - wordSize + 1);

			if (firstWord != lastWord) {
				if (firstWord >= 0 && lastWord >= 0) {
					// This is what the adjustment to d2 should be:
					//d2 = d2 - (int)Math.pow(H1_left[firstWord] - H2[firstWord], 2) + (int)Math.pow(H1_left[firstWord] - (H2[firstWord]-1), 2) - 
					//(int)Math.pow(H1_left[lastWord] - H2[lastWord], 2) + (int)Math.pow(H1_left[lastWord] - (H2[lastWord]+1), 2);
					// Hazelhurst proves the following is equivilent, IF I am understanding correctly -- must be checked.
					d2_left += (((H2[lastWord] - H1_left[lastWord]) - (H2[firstWord] - H1_left[firstWord]) + 1) << 1);
					d2_right += (((H2[lastWord] - H1_right[lastWord]) - (H2[firstWord] - H1_right[firstWord]) + 1) << 1);
					H2[firstWord]--;
					H2[lastWord]++;		   
				}
				else {
					if (firstWord >= 0) {
						d2_left = (int) (d2_left - Math.pow(H1_left[firstWord] - H2[firstWord], 2) + Math.pow(H1_left[firstWord] - (H2[firstWord]-1), 2));
						d2_right = (int) (d2_right - Math.pow(H1_right[firstWord] - H2[firstWord], 2) + Math.pow(H1_right[firstWord] - (H2[firstWord]-1), 2));
						H2[firstWord]--;
					}
					if (lastWord >= 0) {
						d2_left = (int) (d2_left - Math.pow(H1_left[lastWord] - H2[lastWord], 2) + Math.pow(H1_left[lastWord] - (H2[lastWord]+1), 2));
						d2_right = (int) (d2_right - Math.pow(H1_right[lastWord] - H2[lastWord], 2) + Math.pow(H1_right[lastWord] - (H2[lastWord]+1), 2));
						H2[lastWord]++;
					}
				}

			}
			if (d2_left == bestLeftScore) {
				bestLeftWindow[numBestLeft++] = i+1;
			}
			else if (d2_left < bestLeftScore) {
				bestLeftScore = d2_left;
				bestLeftWindow[0] = i+1;
				numBestLeft = 1;
			}

			if (d2_right == bestRightScore) {
				bestRightWindow[numBestRight++] = i+1;
			}
			else if (d2_right < bestRightScore) {
				bestRightScore = d2_right;
				bestRightWindow[0] = i+1;
				numBestRight = 1;
			}
		}
		
		if (bestLeftScore > THRESHOLD) {
			bestLeftWindow = null;
			numBestLeft = 0;
			bestLeftScore = INT_MAX;
		}
		if (d2_right > THRESHOLD) {
			bestRightWindow = null;
			numBestRight = 0;
			bestRightScore = INT_MAX;
		}
		return new BestWindowMatches(bestLeftWindow, numBestLeft, bestLeftScore, bestRightWindow, numBestRight, bestRightScore);
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
		//String s1 = "TTACGaCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCtGAaGtGCATCCaCTGtCGTGACAACCACGGTGCCGTGGCCGAGaAGgCCCTGCgGcCGCGCCAAGTAGCtGGACATTTGGAcTTGGTTGTgGCGCGCTGCAGCTGAACGGGCCGaCCGTTCGAGCGGTGGCGTTCCtCTACGCAGTAGCgGCGCGcACGGGCACCATgGgAAGTCGCATGGTTTTCATGTT";
		//String s2 = "AGCtGGACATTTGGAcTTGGTTGTgGCGCGCTGCAGCTGAACGGGCCGaCCGTTCGAGCGGTGGCGTTCCtCTACGCAGTAGCgGCGCGcACGGGCACCATgGgAAGTCGCATGGTTTTCATGTTTTACGaCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCtGAaGtGCATCCaCTGtCGTGACAACCACGGTGCCGTGGCCGAGaAGgCCCTGCgGcCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTtCTGGCTAGAgCTGACACTAAaAAATGGGCAAgTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGCCGGTGGATCAGCaGgCAGCAACTACCATACaTCAGTgGCGAGCGaTGCAGCAGGAGCTCGGCGCgATTTgGTCTCATCCTTATATcAGAACTAgCAGaAGCAGtGCGTTCcAGtcGTATGGAaGCGTAGCATGAGCAAATGgCCcCc";
		String s1 = "GCAacGCTgTGCAGGAATCTAgTACGCTTTGAGGTGCAAAACAGCCCTTCAtCGCCGGAGATGCCGCTATgGGGAGCGGAGACGCCCAGCANCCTGGGGCCGTCATATCTTGTGgNCCCATGGCCCGCTATGGAGCGGACTAaCTTNCAGATTACcAGCAGCAaGCCTGGTGCtGCgGCGGACGGCAGCaTCGGCGGCTgGGGcgATGGtGGCGaGCAAGGTCGCTgACAGTCGCTGCCAGAACCGCCAGTCATTAAGCAcCAGAAAcTATCCGCTCGTGGGaAACGcAGTGCCCCAGCAGCAGTTATcNGTGcAGCTCtCGcCACGCAGGGAGCGCCGGAtaCCAGCCCCaGCAggTACAGCTATAAATTTTCAGCAACAGAtTCTCCcTTTtGGGCTTTCTtTNaTTTGTTGCTTACCGtCTAGTCGCCGcCTGGAAGcCTGCGCTAATCGGCCTGACACCtTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGAGGcN";
		String s2 = "aGGaAtCANTCAgTAAGCATGGCTATcGCAGGAATCTATACGCTTTNANGTgCAACAACAGCCgTTCAGtGCcaCGtAGATGCCGCTgTGGGGAGCGGAGACGCgCAGCAGCCNGGGGCCGTCATATtTgGcGCcCCgCATGGCCCCTATGGAGCTGaATACCTTACAGATTACAGCcGCATGCTCTGGTGCTGCCGCGGACtGGCAGCCaCGGCGGCTtGGGTCATGGCGGCGAcGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCatGCTCGtTGaGAAtGCAGTGCCaCAttAGCNGTTATCCGTGGAGCTCGCGgACGCANtGtGCGCCGGAGGCCAGCCCCAGCACCtACAGtTATAAATTTTtAtNtACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCG";
		BestWindowMatches best = d2.matchEndWindows(s1, s2);

		System.out.print("bestLeftStart: ");
		for (int i = 0; i < best.numBestLeftWindows; i++)
		    System.out.print(" " + best.bestLeftStart[i]);
		System.out.println("");
		System.out.println("bestLeftD2 = " + best.bestLeftD2);

		System.out.print("bestRightStart: ");
		for (int i = 0; i < best.numBestRightWindows; i++)
		    System.out.print(" " + best.bestRightStart[i]);
		System.out.println("");
		System.out.println("bestRightD2 = " + best.bestRightD2);
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

class BestWindowMatches {    
	public int[] bestLeftStart;
	public int numBestLeftWindows;
	public int bestLeftD2;

	public int[] bestRightStart;
	public int numBestRightWindows;
	public int bestRightD2;

	public BestWindowMatches(int[] leftStart, int numBestLeft, int leftD2, int[] rightStart, int numBestRight, int rightD2) {
		bestLeftStart = new int[numBestLeft];
		for (int i=0; i < numBestLeft; i++)
			bestLeftStart[i] = leftStart[i];
		numBestLeftWindows = numBestLeft;
		bestLeftD2 = leftD2;

		bestRightStart = new int[numBestRight];
		for (int i=0; i < numBestRight; i++)
			bestRightStart[i] = rightStart[i];
		numBestRightWindows = numBestRight;
		bestRightD2 = rightD2;
	}
}



