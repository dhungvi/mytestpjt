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
	protected int THRESHOLD;	// THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
	// if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).


	// d2 parameters
	protected int d2WordSize; 	
	protected int d2NumWords;     // number of different words;
	protected int d2WordFilter;

	// heuristic parameters
	protected int heuristicWordSize;
	protected int heuristicNumWords;
	protected int heuristicWordFilter;
	protected int u;
	protected int uv_skip;
	protected int t;
	protected int tv_max;


	public D2(Properties props) {
		windowSize = Integer.parseInt(props.getProperty("windowSize"));
		THRESHOLD = Integer.parseInt(props.getProperty("THRESHOLD"));

		// d2 parameters
		d2WordSize = Integer.parseInt(props.getProperty("boundOfWord"));
		d2WordFilter = (1 << (d2WordSize << 1)) - 1; // 2^(2*d2WordSize) - 1
		d2NumWords = 1 << (d2WordSize << 1); // 2^(2*d2WordSize)

		// heuristic parameters
		heuristicWordSize = Integer.parseInt(props.getProperty("HeuristicWordSize"));
		heuristicWordFilter = (1 << (heuristicWordSize << 1)) - 1;
		heuristicNumWords = 1 << (heuristicWordSize << 1);

		u = Integer.parseInt(props.getProperty("u"));
		uv_skip = Integer.parseInt(props.getProperty("uv_skip"));
		t = Integer.parseInt(props.getProperty("t"));
		tv_max = Integer.parseInt(props.getProperty("tv_max"));
	}

	public int getWindowSize() {
		return windowSize;
	}

	public int getd2WordSize() {
		return d2WordSize;
	}

	public int getHeuristicWordSize() {
		return heuristicWordSize;
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

	private int[] createWindowHash(String s, int leftCoord, int windowSize, int wordSize, int wordFilter, int numWords) {
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
	private int encodeWord(String s, int leftCoord, int wordSize, int wordFilter) {
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
		if (!uv_tv_Heuristic(s1, s2))
			return new BestWindowMatches(new int[1], 0, 0, new int [1], 0, 0);

		//System.out.println("HERE");

		int[] H1_left = createWindowHash(s1, 0, windowSize, d2WordSize, d2WordFilter, d2NumWords);
		int[] H1_right = createWindowHash(s1, s1.length() - windowSize, windowSize, d2WordSize, d2WordFilter, d2NumWords);
		int[] H2 = createWindowHash(s2, 0, windowSize, d2WordSize, d2WordFilter, d2NumWords);

		int d2_left = 0;
		int d2_right = 0;
		for (int i=0; i < d2NumWords; i++) {
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
			int firstWord = encodeWord(s2, i, d2WordSize, d2WordFilter);
			int lastWord = encodeWord(s2, i + windowSize - d2WordSize + 1, d2WordSize, d2WordFilter);

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

	// Hazelhusrt uv heuristic: v = object word size
	private boolean uv_tv_Heuristic(String s1, String s2) {

		// The u/v heuristic
		// Look at every (uv_skip) word on s1 and count the number of instances on s2.
		// Return false if the value is less than u.
		int[] H = createWindowHash(s2, 0, s2.length(), heuristicWordSize, heuristicWordFilter, heuristicNumWords);
		int total = 0;
		for (int i=0; total < u && i <= s1.length() - heuristicWordSize; i += uv_skip) {
			int code = encodeWord(s1, i, heuristicWordSize, heuristicWordFilter);
			if (code >= 0) {
				total += H[code];
			}
		}
		if (total < u)
			return false;

		// the t/v heursitc
		// Must find at least t words on s2 that occur within 100 bases of eachother on s1.
		int[] arr = new int[tv_max];   // Assuming this gets initilized automatically
		int current_position = 0;
		int current_code = encodeWord(s1, 0, heuristicWordSize, heuristicWordFilter);
		while (current_code < 0 && current_position <=  s1.length() - heuristicWordSize) {
			current_position += -current_code;
			current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);  // Shift over past the N
		}
		if (current_code >= 0) {
			total = H[current_code];
			arr[current_position % tv_max] = total;
		}
		else
			total = 0;

		while (current_position < s1.length() - heuristicWordSize ) {
			if (total >= t)
				return true;

			int next_char = encodeBase(s1.charAt(current_position+heuristicWordSize));
			if (next_char >= 0)  {
				current_code = ((current_code << 2) | next_char) & heuristicWordFilter;
				current_position++;
			}
			else {
				current_position += 2;
				current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);
				while (current_code < 0) {
					current_position += -current_code;
					if (current_position > s1.length() - heuristicWordSize)
						break;
					current_code = encodeWord(s1, current_position, heuristicWordSize, heuristicWordFilter);
				}
			}

			int current_index = current_position % tv_max;
			total = total - arr[current_index] + H[current_code];
			arr[current_index] = H[current_code];	    
		}
		return false;
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


