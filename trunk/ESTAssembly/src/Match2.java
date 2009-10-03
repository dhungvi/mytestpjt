/*
 * Note: this java file is found from the Internet because we need a linear space alignment to avoid Java
 * out-of-memory error.
 * It is used by :
 * 		getAScore in ResultAnalysis.java to compute A-Score.
 * 		processMoreConsensusWithInclusion in Reconstruction.java to make local alignment between two consensus 
 * 			sequences.
 * It is not very convenient to set scoring system in this file. We have to use a two-dimensional array in 
 * subclasses of "Substitution" class to store the scoring system.
*/


// Implementation of some algorithms for pairwise alignment from
// Durbin et al: Biological Sequence Analysis, CUP 1998, chapter 2.
// Peter Sestoft, sestoft@itu.dk 1999-09-25, 2003-04-20 version 1.4
// Reference:  http://www.itu.dk/people/sestoft/bsa.html

// License: Anybody can use this code for any purpose, including
// teaching, research, and commercial purposes, provided proper
// reference is made to its origin.  Neither the author nor the Royal
// Veterinary and Agricultural University, Copenhagen, Denmark, can
// take any responsibility for the consequences of using this code.

// Compile with:
//      javac Match2.java
// Run with:
//      java Match2 HEAGAWGHEE PAWHEAE

//  Class hierarchies
//  -----------------
//  Align                   general pairwise alignment
//     AlignSimple          alignment with simple gap costs
//        NW                global alignment with simple gap costs
//        SW                local alignment with simple gap costs
//        RM                repeated matches with simple gap costs
//        OM                overlap matches with simple gap costs 
//     AlignAffine          alignment with affine gap costs (FSA model)
//        NWAffine          global alignment with affine gap costs
//     AlignSmart           alignment using smart linear-space algorithm
//        NWSmart           global alignment using linear space
//        SWSmart           local alignment using linear space
//     AlignSmartAffine     alignment w affine gap costs in linear space
//        SWSmartAffine     local alignment w affine gap costs in linear space
//  Traceback               traceback pointers
//     Traceback2           traceback for simple gap costs
//     Traceback3           traceback for affine gap costs
//  Substitution            substitution matrices with fast lookup
//     Blosum50             the BLOSUM50 substitution matrix
//  Output                  general text output
//     SystemOut            output to the console (in the application)
//     TextAreaOut          output to a TextArea (in the applet)

// Notational conventions: 
//   i in {0..n} indexes columns and sequence seq1
//   j in {0..m} indexes rows    and sequence seq2
//   k in {0..2} indexes states (in affine alignment)

// The class of substitution (scoring) matrices

// Modified by Jenna on 09/10/03 to make it useful for computing A-Score.
abstract class Substitution {
	public int[][] score;

	void buildscore(String residues, int[][] residuescores) {
		// Allow lowercase and uppercase residues (ASCII code <= 127):
		score = new int[127][127];
		for (int i=0; i<residues.length(); i++) {
			char res1 = residues.charAt(i);
			for (int j=0; j<=i; j++) {
				char res2 = residues.charAt(j);
				score[res1][res2] = score[res2][res1] 
				                                = score[res1][res2+32] = score[res2+32][res1] 
				                                                                        = score[res1+32][res2] = score[res2][res1+32] 
				                                                                                                             = score[res1+32][res2+32] = score[res2+32][res1+32] 
				                                                                                                                                                        = residuescores[i][j];
			}
		}
	}

	abstract public String getResidues();
}


// The substitution matrix for gene (Jenna, 09.10.03)

class ForGene extends Substitution {

	private String residues = "AGCTN";

	public String getResidues() 
	{ return residues; }

	private int[][] residuescores = 
	{ 		/* A  G  C  T  N*/
			/* A */ {  1              },
			/* G */ { -1, 1           },
			/* C */ { -1,-1, 1        },
			/* T */ { -1,-1, -1, 1    },
			/* N */ { -1,-1, -1, -1, 2},
	};

	public ForGene() 
	{ buildscore(residues, residuescores); }
}


//The substitution matrix for computing A-Score (Jenna, 09.10.03)

class ForAScore extends Substitution {

	private String residues = "AGCTN";

	public String getResidues() 
	{ return residues; }

	private int[][] residuescores = 
	{ 		/* A  G  C  T  N*/
			/* A */ {  2              },
			/* G */ { -3, 2           },
			/* C */ { -3,-3, 2        },
			/* T */ { -3,-3, -3, 2    },
			/* N */ { -3,-3, -3, -3, 2},
	};

	public ForAScore() 
	{ buildscore(residues, residuescores); }
}


// Pairwise sequence alignment 

abstract class Match2 {
	Substitution sub;             // substitution matrix
	int d;                        // gap cost
	String seq1, seq2;            // the sequences
	int n, m;                     // their lengths
	Traceback B0;                 // the starting point of the traceback

	final static int NegInf = Integer.MIN_VALUE/2; // negative infinity

	public Match2(Substitution sub, int d, String seq1, String seq2) {
		this.sub = sub;
		this.seq1 = strip(seq1); this.seq2 = strip(seq2);
		this.d = d;
		this.n = this.seq1.length(); this.m = this.seq2.length();
	}

	public String strip(String s) {
		boolean[] valid = new boolean[127];
		String residues = sub.getResidues();
		for (int i=0; i<residues.length(); i++) {
			char c = residues.charAt(i);
			if (c < 96) 
				valid[c] = valid[c+32] = true;
			else
				valid[c-32] = valid[c] = true;
		}
		StringBuffer res = new StringBuffer(s.length());
		for (int i=0; i<s.length(); i++)
			if (valid[s.charAt(i)])
				res.append(s.charAt(i));
		return res.toString();
	}

	// Return two-element array containing an alignment with maximal score

	public String[] getMatch() {
		StringBuffer res1 = new StringBuffer();
		StringBuffer res2 = new StringBuffer();
		Traceback tb = B0;
		int i = tb.i, j = tb.j;
		while ((tb = next(tb)) != null) {
			if (i == tb.i) 
				res1.append('-');
			else
				res1.append(seq1.charAt(i-1)); 
			if (j == tb.j) 
				res2.append('-');
			else
				res2.append(seq2.charAt(j-1)); 
			i = tb.i; j = tb.j;
		}        
		String[] res = { res1.reverse().toString(), res2.reverse().toString() };
		return res;
	}

	public String fmtscore(int val) {
		if (val < NegInf/2) 
			return "-Inf";
		else
			return Integer.toString(val);
	}

	// Get the next state in the traceback
	public Traceback next(Traceback tb)
	{ return tb; }                // dummy implementation for the `smart' algs.

	// Return the score of the best alignment
	public abstract int getScore();

	// Auxiliary functions
	static int max(int x1, int x2) 
	{ return (x1 > x2 ? x1 : x2); }

	static int max(int x1, int x2, int x3) 
	{ return max(x1, max(x2, x3)); }

	static int max(int x1, int x2, int x3, int x4) 
	{ return max(max(x1, x2), max(x3, x4)); }

	static String padLeft(String s, int width) {
		int filler = width - s.length();
		if (filler > 0) {           // and therefore width > 0
			StringBuffer res = new StringBuffer(width);
			for (int i=0; i<filler; i++)
				res.append(' ');
			return res.append(s).toString();
		} else
			return s;
	}
}


// Alignment with simple gap costs

abstract class AlignSimple extends Match2 {
	int[][] F;                    // the matrix used to compute the alignment
	Traceback2[][] B;             // the traceback matrix

	public AlignSimple(Substitution sub, int d, String seq1, String seq2) {
		super(sub, d, seq1, seq2);
		F = new int[n+1][m+1];
		B = new Traceback2[n+1][m+1];
	}

	public Traceback next(Traceback tb) {
		Traceback2 tb2 = (Traceback2)tb;
		return B[tb2.i][tb2.j];
	}

	public int getScore() 
	{ return F[B0.i][B0.j]; }

}


// Traceback objects

abstract class Traceback {
	int i, j;                     // absolute coordinates
}


// Traceback2 objects for simple gap costs

class Traceback2 extends Traceback {
	public Traceback2(int i, int j)
	{ this.i = i; this.j = j; }
}


// Global alignment with the Needleman-Wunsch algorithm (simple gap costs)

class NW extends AlignSimple {

	public NW(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		int n = this.n, m = this.m;
		int[][] score = sub.score;
		for (int i=1; i<=n; i++) {
			F[i][0] = -d * i;
			B[i][0] = new Traceback2(i-1, 0);
		}
		for (int j=1; j<=m; j++) {
			F[0][j] = -d * j;
			B[0][j] = new Traceback2(0, j-1);
		}
		for (int i=1; i<=n; i++)
			for (int j=1; j<=m; j++) {
				int s = score[seq1.charAt(i-1)][seq2.charAt(j-1)];
				int val = max(F[i-1][j-1]+s, F[i-1][j]-d, F[i][j-1]-d);
				F[i][j] = val;
				if (val == F[i-1][j-1]+s)
					B[i][j] = new Traceback2(i-1, j-1);
				else if (val == F[i-1][j]-d)
					B[i][j] = new Traceback2(i-1, j);
				else if (val == F[i][j-1]-d)
					B[i][j] = new Traceback2(i, j-1);
				else
					throw new Error("NW 1");
			}
		B0 = new Traceback2(n, m);
	}
}


// Local alignment with the Smith-Waterman algorithm (simple gap costs)

class SW extends AlignSimple {
	public SW(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		int n = this.n, m = this.m;
		int[][] score = sub.score;
		int maxi = n, maxj = m;
		int maxval = NegInf;
		for (int i=1; i<=n; i++)
			for (int j=1; j<=m; j++) {
				int s = score[seq1.charAt(i-1)][seq2.charAt(j-1)];
				int val = max(0, F[i-1][j-1]+s, F[i-1][j]-d, F[i][j-1]-d);
				F[i][j] = val;
				if (val == 0)
					B[i][j] = null;
				else if (val == F[i-1][j-1]+s)
					B[i][j] = new Traceback2(i-1, j-1);
				else if (val == F[i-1][j]-d)
					B[i][j] = new Traceback2(i-1, j);
				else if (val == F[i][j-1]-d)
					B[i][j] = new Traceback2(i, j-1);
				else
					throw new Error("SW 1");
				if (val > maxval) {
					maxval = val;
					maxi = i; maxj = j;
				}
			}
		B0 = new Traceback2(maxi, maxj);
	}
}



// Alignment with simple gap costs; smart linear-space algorithm

abstract class AlignSmart extends Match2 {
	int[][] F;

	public AlignSmart(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		F = new int[2][m+1];
	}

	static void swap01(Object[] A) 
	{ Object tmp = A[1]; A[1] = A[0]; A[0] = tmp; }
}


// Global alignment (simple gap costs, smart linear-space algorithm)

class NWSmart extends AlignSmart {
	int u;     // Halfway through seq1
	int[][] c; // Best alignment from (0,0) to (i,j) passes through (u, c[i][j]) 

	public NWSmart(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		int n = this.n, m = this.m;
		u = n/2;
		c = new int[2][m+1];
		int[][] score = sub.score;
		for (int j=0; j<=m; j++)
			F[1][j] = -d * j;
		for (int i=1; i<=n; i++) {
			swap01(F); swap01(c);
			// F[1] represents (new) column i and F[0] represents (old) column i-1
			F[1][0] = -d * i;
			for (int j=1; j<=m; j++) {
				int s = score[seq1.charAt(i-1)][seq2.charAt(j-1)];
				int val = max(F[0][j-1]+s, F[0][j]-d, F[1][j-1]-d);
				F[1][j] = val;
				if (i == u)
					c[1][j] = j;
				else 
					if (val == F[0][j-1]+s)
						c[1][j] = c[0][j-1];
					else if (val == F[0][j]-d)
						c[1][j] = c[0][j];
					else if (val == F[1][j-1]-d)
						c[1][j] = c[1][j-1];
					else
						throw new Error("NWSmart 1");
			}
		}
	}

	public int getV() 
	{ return c[1][m]; }

	public String[] getMatch() {
		int v = getV();
		if (n > 1 && m > 1) {
			NWSmart al1, al2;
			al1 = new NWSmart(sub, d, seq1.substring(0, u), seq2.substring(0, v));
			al2 = new NWSmart(sub, d, seq1.substring(u),    seq2.substring(v));
			String[] match1 = al1.getMatch();
			String[] match2 = al2.getMatch();
			String[] match = { match1[0] + match2[0], match1[1] + match2[1] };
			return match;
		} else {
			NW al = new NW(sub, d, seq1, seq2);
			return al.getMatch();
		}
	}

	public int getScore() 
	{ return F[1][m]; }
}

//Global alignment (simple gap costs, smart linear-space algorithm)

class AScore extends AlignSmart {
	int u;     // Halfway through seq1
	int[][] c; // Best alignment from (0,0) to (i,j) passes through (u, c[i][j]) 

	public AScore(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		int n = this.n, m = this.m;
		u = n/2;
		c = new int[2][m+1];
		int[][] score = sub.score;
		for (int j=0; j<=m; j++)
			F[1][j] = 0;
		for (int i=1; i<=n; i++) {
			swap01(F); swap01(c);
			if (i == n) {
				d = 0;
			}
			// F[1] represents (new) column i and F[0] represents (old) column i-1
			F[1][0] = -d * i;
			for (int j=1; j<=m; j++) {
				int s = score[seq1.charAt(i-1)][seq2.charAt(j-1)];
				int val = max(F[0][j-1]+s, F[0][j]-d, F[1][j-1]-d);
				F[1][j] = val;
				if (i == u)
					c[1][j] = j;
				else 
					if (val == F[0][j-1]+s)
						c[1][j] = c[0][j-1];
					else if (val == F[0][j]-d)
						c[1][j] = c[0][j];
					else if (val == F[1][j-1]-d)
						c[1][j] = c[1][j-1];
					else
						throw new Error("AScore 1");
			}
		}
	}

	public int getV() 
	{ return c[1][m]; }

	public String[] getMatch() {
		int v = getV();
		if (n > 1 && m > 1) {
			NWSmart al1, al2;
			al1 = new NWSmart(sub, d, seq1.substring(0, u), seq2.substring(0, v));
			al2 = new NWSmart(sub, d, seq1.substring(u),    seq2.substring(v));
			String[] match1 = al1.getMatch();
			String[] match2 = al2.getMatch();
			String[] match = { match1[0] + match2[0], match1[1] + match2[1] };
			return match;
		} else {
			NW al = new NW(sub, d, seq1, seq2);
			return al.getMatch();
		}
	}

	public int getScore() 
	{ return F[1][m]; }
}


// Local alignment with the Smith-Waterman algorithm (simple gap
// costs, smart linear space algorithm)

class SWSmart extends AlignSmart {
	Traceback2[][] start; // Best alignment ending at (i,j) begins at start[i][j]
	int maxval;           // Score of best alignment
	int start1, start2;   // Best alignment begins at (start1, start2)
	int end1, end2;       // Best alignment ends at (end1, end2)

	public SWSmart(Substitution sub, int d, String sq1, String sq2) {
		super(sub, d, sq1, sq2);
		int n = this.n, m = this.m;
		int[][] score = sub.score;
		start = new Traceback2[2][m+1];
		maxval = NegInf;
		// Initialize first column (i=0):
		for (int j=0; j<=m; j++)
			start[1][j] = new Traceback2(0, j);
		for (int i=1; i<=n; i++) {
			swap01(F); swap01(start);
			// F[1] represents (new) column i and F[0] represents (old) column i-1
			// Initialize first row (j=0):
			start[1][0] = new Traceback2(i, 0);
			for (int j=1; j<=m; j++) {
				int s = score[seq1.charAt(i-1)][seq2.charAt(j-1)];
				int val = max(0, F[0][j-1]+s, F[0][j]-d, F[1][j-1]-d);
				F[1][j] = val;
				if (val == 0)           // Best alignment starts (and ends) here
					start[1][j] = new Traceback2(i, j);
				else if (val == F[0][j-1]+s)
					start[1][j] = start[0][j-1];
				else if (val == F[0][j]-d)
					start[1][j] = start[0][j];
				else if (val == F[1][j-1]-d)
					start[1][j] = start[1][j-1];
				else
					throw new Error("SWSmart 1");
				if (val > maxval) {
					maxval = val;
					Traceback2 sij = start[1][j];
					start1 = sij.i; start2 = sij.j;
					end1 = i; end2 = j;
				}
			}
		}
	}

	public int getScore() 
	{ return maxval; }

	public String[] getMatch() {
		String subseq1 = seq1.substring(start1, end1);
		String subseq2 = seq2.substring(start2, end2);
		// The optimal local alignment between seq1 and seq2 is the
		// optimal global alignment between subseq1 and subseq2:
		return (new NWSmart(sub, d, subseq1, subseq2)).getMatch();
	}
}





