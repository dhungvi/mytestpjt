import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

import neobio.alignment.BasicScoringScheme;
import neobio.alignment.IncompatibleScoringSchemeException;
import neobio.alignment.NeedlemanWunsch;
import neobio.alignment.SmithWaterman;

/**
 * This class implements the D2 algorithm. Feb 15, 2009.
 * 
 * Note: we need change the method 'initWords' to generate more words if boundOfWord changes.
 */

public class Alignment {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected int alignmentThreshold; //It is used in NewD2 class. It's the threshold for alignment. That is, all the alignment with
											// the distance which is bigger than the value will be seen as infinity. 

	
	public Alignment(Properties props) {
		alignmentThreshold = Integer.parseInt(props.getProperty("alignmentThreshold"));
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


	public static void main(String args[]) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}
		
		Alignment al = new Alignment(props);
		String s1 = "GCAacGCTgTGCAGGAATCTAgTTTTTACGCTTTGAGGTGCAAACCCCCCCACAGCCCTTCAtCGCCGGAGTTTTTTTTTTTATGCCGCTATgGGGAGCGGAGACGCCCAGCANCCTGGGGCCGTCATATCTTGTGgNCCCATGGCCCGCTATGGAGCGGACTAaCTTNCAGATTACcAGCAGCAaGCCTGGTGCtGCgGCGGACGGCAGCaTCGGCGGCTgGGGcgATGGtGGCGaGCAAGGTCGCTgACAGTCGCTGCCAGAACCGCCAGTCATTAAGCAcCAGAAAcTATCCGCTCGTGGGaAACGcAGTGCCCCAGCAGCAGTTATcNGTGcAGCTCtCGcCACGCAGGGAGCGCCGGAtaCCAGCCCCaGCAggTACAGCTATAAATTTTCAGCAACAGAtTCTCCcTTTtGGGCTTTCTtTNaTTTGTTGCTTACCGtCTAGTCGCCGcCTGGAAGcCTGCGCTAATCGGCCTGACACCtTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGAGGcN";
		String s2 = "CAcCCCTaCAgCGCCGGcGtNGCCGgTATGGGNNNNNNNGAGCGGAGcCGCCGGGGGGGGCAGCAGCCTTTGGGGCCGacNNATATCTTGTGCCCCCATGgtCCCGCTATGGAGCgGGACACCTTACAGATTACAGCAGCATGtTCTGGTGCTGCCGCGtACGGCAGCCTCGGCGGCTGGGGTCNTGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAtATTATCCGCTCGTGcGAACGCAGTGCCCCAcCaAGCaGTTATCCGTGGAGCTCGCGCACGNAGGGtAGCGCCGGAGGCCAGCCCCAGCaCCTACAGCTATAAATTTTCAGCAAaAaATTCTgCATTTTGGGCTTCTTTGGTTTGtTGCTTACCgcCTAGtCGCCgTCTGGAAGGCTGCGCTAATCGGtCTGACACCCTTGGgGCGTTGCCTCTCAtTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCcCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCaGAAGTCTCaCTTCTGCCaTg";
		s1 = s1.toUpperCase();
		s2 = s2.toUpperCase();
		//System.out.println(d2.getSimlarityScore(s1,s2));
		System.out.println(al.getLocalAlignment(s1,s2)[0] + "\n" + al.getLocalAlignment(s1,s2)[1]);
		System.out.println(al.getLocalSimlarityScore(s1, s2));
		//String s="CGGGCCTTCGTTTTACGAAAACAGGTGC";
		//String s3 = "TTTACGA-AAACA";
		//s = s.replace(s3.replace("-", ""), s3);
		//System.out.println(al.getGlobalAlignment(s1, s2));
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
