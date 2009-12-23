import com.mhhe.clrs2e.MergeSort;

public class Debugger {

	/*
	 * print elements in 'dGraph'
	 */
	static void printDgraph(int[][] d){
		for (int i=0; i<d.length; i++) {
			System.out.println(d[i][0] + "\t" + d[i][1] + "\t" + d[i][2] + "\t" + d[i][3]);
		}
		System.out.println();
	}	
	
	/*
	 * print information of all the assumed left-end nodes.
	 *  	starting position of the node;
	 * 		whether or not they are real left ends;
	 * 		If they are false left ends, print the overlap length they have with other nodes.
	 */
	static String printLeftEndInfo(int leftEnd, Graph g) {
		String ret = "";
		int sp = Integer.parseInt(g.getNameOfNode(leftEnd));	//actual starting position of the node
		ret = "Node " + leftEnd + " starts from " + sp + "\n";
		int ln = g.getLenOfNode(leftEnd);
		int flag = 0;
		for (int k=0; k<g.getSizeofGraph(); k++) {
			int tmpSp = Integer.parseInt(g.getNameOfNode(k));
			int tmpLn = g.getLenOfNode(k);
			if ((sp > tmpSp) && 
					(sp < (tmpSp+tmpLn-1)) &&
					((sp+ln-1) > (tmpSp+tmpLn-1))) {
				//if overlap length is less than windowsize, we consider they're not overlapping.
				//if not, we see them as overlapping,so this is not a real left end.
				if ((tmpSp+tmpLn-sp)>=(g.ovl.d2.getWindowSize())) {
					ret = ret + "Node " + leftEnd + " is not a real left-most node. ";
					ret = ret + "Overlap length with node " + k + " is " + (tmpSp+tmpLn-sp) + "\n"; //(sp+ln-tmpSp));
					flag = 1;
				}
			}
		}
		
		if (flag == 0) {
			ret = ret + "Node " + leftEnd + " is a real left-most node.\n";
		} 
		return ret;
	}
	
	/* 
	 * Print the assembled starting position and the actual position for all the ests
	 */
	static void printSPos(Graph g, int[] sPosDebug) {
		System.out.println("Calculated s_i	Actual s_i	sPos");
		for (int i=0; i<sPosDebug.length; i++) {
			System.out.println(sPosDebug[i] + "	" + g.getNameOfNode(i));
		}
	}



	/* 
	 * Calculate inversions for all the calculated positions of ESTs
	 */
	static void calcInversion(Graph g, int[] sPosDebug) {
			//Firstly, sort the array sPos
			StartPos2[] resultArray = new StartPos2[sPosDebug.length]; //store the starting positions of ests
			for (int i=0; i<sPosDebug.length; i++) {
				resultArray[i] = new StartPos2(sPosDebug[i], g.getNameOfNode(i));
			}
			MergeSort merge = new MergeSort();
			merge.sort(resultArray);
			
			int[] inversionArray = new int[resultArray.length];
			for (int j=0; j<resultArray.length; j++) {
				inversionArray[j] = resultArray[j].realStartPos;
			}
			System.out.print( "The assembly has " );
            System.out.println(Inversions.countInversions(inversionArray) + " inversions.");
	}

	static class StartPos2 implements Comparable<StartPos2> {
		int pos;
		int realStartPos;
		public StartPos2(int p, String s) {
			pos = p;
			realStartPos = (Integer.valueOf(s)).intValue();
		}
		
		public int compareTo(StartPos2 other) {
			//Returns 0 if the argument is equal to this; 			
			//a value less than 0 if the argument is greater than this; 
			//and a value greater than 0 if the argument is less than this. 
			if (this.pos == other.pos) {
				return 0;
			} else if (this.pos > other.pos) {
				return 1;
			} else {
				return -1;
			}
		}
	}	


}

