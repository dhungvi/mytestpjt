import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;


/*
 * This class analyze the results of assembly.
 * There are several input files to this class:
 * 	cap3.out: the contigs from cap3, one contig in one line, optional;
 *  cap3Singleton.out: the singletons from cap3, one contig in one line, optional;
 *  esta.out: the contigs from my tool, one contig in one line;
 *  estaSingleton.out: the singletons from my tool, one contig in one line.
 *  
 * The class evaluate the result using the following parameters:
 *  number of contigs, ideally it is 1;
 *  number of singletons, ideally it is 0;
 *  length of the longest contig;
 *  number of difference between the longest contig and the corresponding sub-sequence in the provided sequence;
 *  A-score, ideally it is 2*(length of the provided sequence);
 */
public class ResultAnalysis {
	D2 d2;
	String oriSeq = "";
	//String genSeq = "";
	ArrayList<String> cap3Contigs;
	ArrayList<String> eastContigs;
	ArrayList<String> cap3Singletons;
	ArrayList<String> eastSingletons;
	String cap3Time;
	String estaTime;
	String cap3NumUsedEsts;
	String eastNumUsedEsts;
	//Alignment alignment;
	
	public ResultAnalysis(Properties props) {
		//alignment = new Alignment(props);
		//read the original gene file
		ArrayList<String> tList1 = readFile(props.getProperty("GeneFileName"));
		if (tList1 != null) {
			oriSeq = tList1.get(0);
		}

		//contigs from cap3
		cap3Contigs = readFile(props.getProperty("Cap3Contigs"));

		//sigletons from cap3
		cap3Singletons = readFile(props.getProperty("Cap3Singletons"));

		//time used by cap3
		ArrayList<String> strs = readFile(props.getProperty("Cap3Time"));
		if (strs != null)	cap3Time = strs.get(0);
		
		//number of ESTs that are used by cap3
		strs = readFile(props.getProperty("Cap3NumUsedEsts"));
		if (strs != null)	cap3NumUsedEsts = strs.get(0);
		
		//contigs from esta
		eastContigs = readFile(props.getProperty("EastContigs"));

		//sigletons from esta
		eastSingletons = readFile(props.getProperty("EastSingletons"));

		//time used by esta
		strs = readFile(props.getProperty("EastTime"));
		if (strs != null)	estaTime = strs.get(0);

		//number of ESTs that are used by east
		strs = readFile(props.getProperty("EastNumUsedEsts"));
		if (strs != null)	eastNumUsedEsts = strs.get(0);
	}

	/*
	 * read strings from a file, one string in one line in the input file. put 
	 * all the strings into an ArrayList, one string is one element.
	 */
	private ArrayList<String> readFile(String fileName) {
		ArrayList<String> ret = new ArrayList<String> ();
		try{ 
			/*
			 * get the provided sequence
			 */
			File oriFile = (new File(fileName));
			if (!oriFile.exists()) {
				System.out.println("The file " + fileName + " does not exist!");
				return null;
			}
			BufferedReader in = new BufferedReader(new FileReader(oriFile));
			String str = in.readLine();
			while (str != null) {
				ret.add(str.trim());
				str = in.readLine();
			}
			in.close();	
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
		
		return ret;
	}
	
	/*
	 * compare the original sequence with the assembly results, put all the analysis into 
	 * a file named "analysis.out". The output file only has two lines, the first line is
	 * a comment, the second line includes ten columns: 
	 * cap3-number of contigs; esta-number of contigs;
	 * cap3-number of singletons; esta-number of singletons;
	 * cap3-length of the longest contig;
	 * esta-length of the longest contig;
	 * cap3- A-score;
	 * esta- A-score;
	 * 
	 */
	public void analysis() {
		int c[] = new int[13];
		if (cap3Contigs != null) {
			c[1] = cap3Contigs.size();

			//get the longest contig in cap3
			String lStr = "";
			int len = 0;
			for (int i=0; i<c[1]; i++) {
				String tStr = cap3Contigs.get(i);
				int tLen = tStr.length();
				if (tLen > len) {
					len = tLen;
					lStr = tStr;
				}
			}
			c[5] = len;
			//c[9] = getDiff(lStr);
			c[7] = getAScore(lStr, oriSeq);
		}
		
		if (eastContigs != null) {
			c[2] = eastContigs.size();

			//get the longest contig in esta
			String lStr = "";
			int len = 0;
			for (int i=0; i<c[2]; i++) {
				String tStr = eastContigs.get(i);
				int tLen = tStr.length();
				if (tLen > len) {
					len = tLen;
					lStr = tStr;
				}
			}
			c[6] = len;
			//c[10] = getDiff(lStr);
			c[8] = getAScore(lStr, oriSeq);
		}
		
		if (cap3Singletons != null) {
			c[3] = cap3Singletons.size();
		}
		if (eastSingletons != null) {
			c[4] = eastSingletons.size();
		}
		
		if (cap3Time != null) {
			c[9] = Integer.valueOf(cap3Time).intValue();
		}
		
		if (estaTime != null) {
			c[10] = Integer.valueOf(estaTime).intValue();
		}
		
		if (cap3NumUsedEsts != null) {
			c[11] = Integer.valueOf(cap3NumUsedEsts).intValue();
		}
		
		if (eastNumUsedEsts != null) {
			c[12] = Integer.valueOf(eastNumUsedEsts).intValue();
		}
		
		/*
		 * write c1-c10 into the output file "analysis.out" in fasta format.
		 */
		try{ 
			File outFile = new File("OutputFile");
			boolean bExists = outFile.exists();
			if (bExists) {
				outFile.delete();
			}
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile, true));
			out.write(">CAP3 vs. EAST\n");
			out.write(">Cap3NumContigs  EastNumContigs  Cap3NumSing  EastNumSin  Cap3Len  EastLen  Cap3AScore  EastAScore  Cap3Time  EastTime  Cap3NumUsedEsts  EastNumUsedEsts\n");
			for (int i=1; i<c.length; i++) {
				out.write(String.valueOf(c[i]));
				out.write("\t");
			}
			out.flush();
			out.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}
	
	/*
	 * Use modified Needleman-Wunsch algorithm to get alignment score.
	 * @param s1, s2
	 * @return A-score
	 */
	public int getAScore(String consensus, String gene) {
/*		int match = 2;
		int gap = -13;	 //gap penalty
		int mis = -3;	 //mismatch penalty
		int n = consensus.length();  //consensus sequence
		int m = gene.length();  //gene
		int[][] A = new int[n+1][m+1];
		for (int j=1; j<m; j++) {
			A[0][j] = 0;
		}

		for (int i=1; i<=n; i++) {
			if (i == n) {
				gap = 0;
			}
			
			A[i][0] = i * gap;
			
			for (int j=1; j<=m; j++) {
				int t1 = A[i-1][j-1] + (consensus.charAt(i-1)==gene.charAt(j-1) ? match : mis);
				int t2 = A[i-1][j] + gap;
				int t3 = A[i][j-1] + gap;
				A[i][j] = Math.max(t1, Math.max(t2, t3));
				//System.out.print(A[i][j]+"\t");
			}
			//System.out.println();
		}

		
		return A[n][m];
*/	
	//we have to use linear space alignment to comput A-score to avoid Java out-of-memory error	
	    Substitution sub = new ForAScore();
	    AScore ascore = (new AScore (sub, 13, consensus, gene));
	    return ascore.getScore();
	
	}

	protected static Properties getProperties(String fName) throws IOException {
		Properties props = new Properties();
		File f = new File(fName);
        
        if (!f.exists()) {
        	return props;
        }
        
        props.load(new FileInputStream(f)); 
        return props;
    }

	/*
	 * analyze the difference between the oriSeq and the input string.
	 * In the input file, the first line is the original sequence and the second line is the consensus.
	 * 
	 * @return number of difference
	 */
/*	public int getDiff(String genSeq) {
		//we will compare the consensus from startPoss[0] of the original sequence to the generated sequence 
		// starting from startPos[1].
		int[] startPoss = getStartCompPos(genSeq);
		
		int num = startPoss[1]; //number of difference
		int idxOfOri = startPoss[0];
		for (int i=startPoss[1]; i<genSeq.length(); i++) {
			if (oriSeq.length()-1 < idxOfOri) {
				num++;
				continue;
			}
			if (oriSeq.charAt(idxOfOri++) != genSeq.charAt(i)) {
				//System.out.println("Index " + i + " is different!");
				num++;
			} 
		}
		//System.out.println("Start to compare from the position: " + startPoss[0] + "," + startPoss[1]);
		//System.out.println("The number of difference is: " + num);
		//System.out.println("The total number of compared base is: " + (genSeq.length()));
		
		return num;
	}
*/	
	/*
	 * Get starting position in the original gene where we start comparison between 
	 * the original and the calculated sequence.
	 */
/*	private int[] getStartCompPos(String genSeq) {
		int ret[] = new int[2];
		String parts[] = alignment.getLocalAlignment(oriSeq, genSeq);
		ret[0] = oriSeq.indexOf(parts[0].replace("-", ""));
		ret[1] = genSeq.indexOf(parts[1].replace("-", ""));

		return ret;
		
	}
*/	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Properties props = null;
		try {
			if (args.length > 0) {
				props = getProperties(args[0]);
			} else {
				props = getProperties("analysis.properties");
			}
		} catch (IOException e) {
			System.err.println("Get analysis.properties failed, " + e);
	    	return;
		}
		
		ResultAnalysis ana = new ResultAnalysis(props);
		//String s1 = "GAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCACCCGGCCGCCTACACCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCCGCCATCGCCGACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCTTATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCACAGTCGAAGGTGTTGACGAGTTCTGTAGGTTCATACATGGCAGTCGTTGAAGTTGCGTTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCATGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGATGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTG";
		//String s2 = "TAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCACCCGGCCGCCTACACCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCCGCCATCGCCGACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCTTATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCACAGTCGAAGGTGTTGACGAGTTCTGTAGGTTCATACATGGCAGTCGTTGAAGTTGCGTTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCATGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGATGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCT";
		//System.out.println(ana.getDiff("CGGGGATGGGTATGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAG-TTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCGTTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGA"));
		ana.analysis();
		//String consensus = "CCCGATAGTCTTTGACCGTTCTATTTTTAGGGACCCGTNCCTAGTGGACGTTAAACTCTNGAGGCAACCTTTGACCGTTCTATTGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCTCCCGGCCGCCTACGCCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCTGCCATCGCCGACACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACCTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCAAATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCGCAGTCGAAGGTGTTGACGAGCTCTGTAGGTTCATACGTGGCAATCGTTGAAGTTGCATTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCGTGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGGTGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCTCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGAAAAACAAAT";
		//String gene = "GAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTGATCGTCTTGCGCCACCCTCTCGCTCGCTACGCATGCGCACCCGGCCGCCTACACCGTGGCCAACACGGCCAGTTGCAGAGGCCGACATAACAGATCTGCAAACATCGGAATGAGAGACCCGGGAGACGTTCCGCCATCGCCGACCACCATGATCAGTCGCTCTTTGCTCTGGCCGATCTGCCTAGACCATGGAGCTATGTTCGGACCTATCTGACTGTGATCGATATTATCCTGAGCAGATCGCGACTACTAGCGAGCGGAGTAGACGGCCGCTGGCGCTCTTATCATTCTAACGGAAGAATGCCCACAGAGCGTTGGGGAGTTTAAGCCAGCTTGAGGACGGAGATCGCCATGATGCGGACGTTGCCGAGTGCCGGAGCGAAACGCCGCGGTGGAGAGCATCTGCAAAATGCGAGAACTCGCACAGTCGAAGGTGTTGACGAGTTCTGTAGGTTCATACATGGCAGTCGTTGAAGTTGCGTTTGTCGTGCTCATGTCTTTCAGTGTTTGAGACGCACATGTTGCATCGATGTTGGAGGCGCTGTATGCCCCAGCTCGCCCAGGAGCTCATAGGAGTGGTGCCCTGAATGCCCTGAGGCTCATTTGGTTGTGTTATGTCTGAGGTTCGCCATTAGCCGCCGTGAACGGAGCGACAGGCTCAACGCGCATGGCTATCACGCCCGCGAGCCCTCTCCATTGCGGTTCACACCAAGGCAGTGTAGGCATGGTTTGTGAGACTTTCAGACGCATGACTGCACTGAGACAGACGTACACGCACGCAGACTGGCAATATGCATACGGCACACCCTGCACAGGCCCAGACGTTGCAGACGAGCGCCGCATTGACATTGCGCCGCACGCTGGGGGGGGAGTCGCCGTCAGGCACACGAAACGTTCATTGTCCTATCATAGGGCGCGACAGCTGTCATGACCTGGCGTGTGACCGTGCATATCGTGTGACTGAGTGATTTGGATGTGGTTTTGGCGTGTGGGAAGAGGGGAACGTGCCCGAGGTTGGTTGGCTCGACACGGTGTATTGCATGGGAGGCCGGAGTACGGGTGATCCGAGAGAGGATATAGGTCGCGTCTGAGAGCGGCTAGAGCCGTCTCTCATTCCGATCCTAGGCGAGATAGTCTTTGACCGTTCTATATGAGGGACCCGTCCTAGTGGACGTTAAACTCTGAGGCAACATGTAAGCCGGAACGTG";
		//System.out.println(ana.getAScore(consensus, gene));
//		ana.getDiff(consensus);
	}

}
