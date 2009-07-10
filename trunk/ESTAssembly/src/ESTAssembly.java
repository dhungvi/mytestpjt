import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Properties;

import com.mhhe.clrs2e.*;

public class ESTAssembly {
	int INT_MAX = 2147483647;
	int INT_MIN = -2147483647;
	String oriFileName;	//store the original gene data. It's used in "printEsts" method.
	String mstFile;	//store the MST. It is generated from PEACE.
	String inFileName;
	String resultFileName;
	String consensusFileName;
	ArrayList<String> ests;	//store all the ests. It is generated in 'readEstFile' function.
	NewGraph g;	//graph to store all the ests. It is generated in 'readEstFile' function.
	WeightedAdjacencyListGraph mstForG;	//minimum spanning tree generated from 'g'
										//It is generated in 'createAlignArray' function.
	/*
	 * store the position of aligned nodes
	 * 1st-dimension: index of nodes in graph;
	 * 2rd-dimension: the first is the index of node on the left, 
	 * 					the second is the overlap length with + or -,
	 * 					the third is the index of node on the right, 
	 * 					the forth is the overlap length with + or -.
	 * 					For the first and third one, if no node is found, the value is -1;
	 * 					For the second and fourth one, if no node is found, the value is infinite.	 
	 */
	int alignArray[][];
	int[] sPos;	//starting positions of all the nodes, initialized in "processAlignArray" function.
				//the index in the array is the index of the node, the value is its starting position. 
	
	public ESTAssembly(Properties props) {
		oriFileName = props.getProperty("SourceFile");
		mstFile = props.getProperty("MSTFile");
		inFileName = props.getProperty("OutFile");
		resultFileName = props.getProperty("ResultFile");
		consensusFileName = props.getProperty("ConsensusFile");
		ests = new ArrayList<String> ();
		alignArray = null;
		sPos = null;
		//g= new Graph(props);
	}
	
	/*
	 * Read ests from the input file;
	 * Generate a Graph object, all the ests are considered to be one node in the graph;
	 * No edge in the graph. Edges will be added in "createAlignArray" function.
	 * 
	 * FASTA format:
	 * A sequence in FASTA format begins with a single-line description, followed by lines of 
	 * sequence data. The description line is distinguished from the sequence data by a greater-than
	 * (">") symbol in the first column. The word following the ">" symbol is the identifier of the 
	 * sequence, and the rest of the line is the description (both are optional). There should be no 
	 * space between the ">" and the first letter of the identifier. It is recommended that all lines 
	 * of text be shorter than 80 characters. The sequence ends if another line starting with a ">" 
	 * appears; this indicates the start of another sequence.
	 * 
	 * EST file comment format:
	 * >g001_001169_001679: first one is number of gene, second is index of tarting position, third is 
	 * 						index of ending position. Index starts from 0.
	 */
	protected void readEstFile() {
		boolean bExists = false;
		
		try{ 
			File f = (new File(inFileName));
			bExists = f.exists();
			if (!bExists) {
				System.out.println("File does not exist!");
				return;
			}

			BufferedReader in = new BufferedReader(new FileReader(f));
			String str = in.readLine();
			while (str != null) {
				str = str.trim();
				// first line is comment line which begins from '>'
				if (str.charAt(0) == '>') {	//comment line begins from '>'
					String[] paras = str.split("_");
					ests.add(paras[1]);
					
					//get est in the next lines
					str = in.readLine();
					StringBuffer estStr = new StringBuffer();
					while (str != null) {
						str = str.trim();
						if (str.compareTo("") != 0)	{
							if (str.charAt(0) != '>') {
								estStr.append(str.trim());
							} else  {
								ests.add(estStr.toString());
								break;
							}
						}
						str = in.readLine();
					}
					if (str == null) {
						ests.add(estStr.toString());
					}
				} 
			}
			in.close();			
		}catch(IOException e){ 
			System.out.println(e.toString());
			return;
		}
		
		//generate a graph from the input file
		int i=0;
		while (i<ests.size()) {
			//the sequence is upper-case
			g.addNode(new Node(ests.get(i), ests.get(i+1).toUpperCase()));
			i = i+2;
		}
	}
	
	public void createAlignArray() {
		//alignArray = g.alignNodes();	
		//mstForG = g.genMST();
		mstForG = g.readMST(mstFile);
		System.out.println("End to generate MST.");
		System.out.println("Start to generate 6-tuples.");
		alignArray = g.get2CloseNodesFromMST(mstForG);
	}

	/*
	 * Find the left-most node(alignNodes[][0]==-1);
	 * Use alignNodes to construct a tree whose values are maximal, that is, 
	 * 		the tree has maximal overlap length;
	 * Calculate the starting position of each node.
	 */
	public void processAlignArray() {
		ArrayList <Integer> leftMostNodes = new ArrayList<Integer> ();
		//get all the nodes which has no left nodes to them
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] == -1) {
				leftMostNodes.add(Integer.valueOf(i));	//store index of the node
			}
		}
		//find the two most closest nodes for those leftmostnodes(because some of them may be falsely left-most)
		int numOfLeftMostNodes = leftMostNodes.size();
		if (numOfLeftMostNodes > 1) {
			for (int i=0; i<leftMostNodes.size(); i++) {
				int cNode = leftMostNodes.get(i).intValue();
				//int[] lNode = g.get2CloseNodes(cNode);
				int[] lNode = g.get2CloseNodesFromGrand(mstForG, cNode, alignArray[cNode]);
				alignArray[cNode][0] = lNode[0];
				alignArray[cNode][1] = lNode[1];
				alignArray[cNode][2] = lNode[2];
				alignArray[cNode][3] = lNode[3];
			}
		} 
		
		//find the two most closest nodes for those rightmostnodes(because we need these information to make a directed maximum spanning tree)
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][2] == -1) {
				//int[] rNode = g.get2CloseNodes(i);
				int[] rNode = g.get2CloseNodesFromGrand(mstForG, i, alignArray[i]);
				alignArray[i][0] = rNode[0];
				alignArray[i][1] = rNode[1];
				alignArray[i][2] = rNode[2];
				alignArray[i][3] = rNode[3];			
			}
		} 

		/* construct a directed graph from alignArray, 
		 *  	the second dimension has two elements:
		 *  			index of starting node,
		 *  			index of ending node, 
		 *  			weight between them (positive value).
		 *  	if there is no edge, weight=INT_MAX.
		 */
		int len = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] != -1) {
				len++;
			}
			if (alignArray[i][2] != -1) {
				len++;
			}
		}
		
		int[][] dGraph = new int[len][3];
		int indexOfDGraph = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] != -1) {
				dGraph[indexOfDGraph][0] = alignArray[i][0];
				dGraph[indexOfDGraph][1] = i;
				dGraph[indexOfDGraph][2] = Math.abs(alignArray[i][1]);
				indexOfDGraph++;
			}
			
			if (alignArray[i][2] != -1) {
				dGraph[indexOfDGraph][0] = i;
				dGraph[indexOfDGraph][1] = alignArray[i][2];
				dGraph[indexOfDGraph][2] = alignArray[i][3];
				indexOfDGraph++;
			}

		}
		
		/*
		 * Remove those false left ends. For example:
		 * Node 2 has the set in alignArray: [-1, 0, 4, 8]
		 * but Node 8 has this set: [7, 5, 2, 6], this means ovlDis(8,2)=6, 
		 * 		so node 2 is not left end because node 8 is to its left.
		 */
		//Get all the nodes which has the value of -1 in alignArray[x][0]
		ArrayList <Integer> tmpLeftNodes = new ArrayList<Integer> ();
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] == -1) {
				tmpLeftNodes.add(Integer.valueOf(i));	//store index of the node
			}
		}
		//remove false left ends
		leftMostNodes.clear();
		for (int i=0; i<tmpLeftNodes.size(); i++) {
			int tEnd = tmpLeftNodes.get(i).intValue();
			int f = 0;
			//if the left end appears in second element of dGraph, that means some 
			//node is on its left, so it is not a real left end.
			for (int j=0; j<dGraph.length; j++) {
				if (dGraph[j][1] == tEnd) {	// false left end
					f = 1;
					break;
				}
			}
			if (f == 0) {
				leftMostNodes.add(Integer.valueOf(tEnd));
			}
		}
		
		//print mst
		System.out.println("Minimum Spanning Tree:");
		System.out.println(mstForG);
		//print dGraph
		System.out.println("dGraph:");
		printDgraph(dGraph);
		
		System.out.println("There are " + leftMostNodes.size() + " left-most nodes.");
		
		/*
		 *  1. print information of all the assumed left-end nodes.
		 *  	starting position of the node;
		 * 		whether or not they are real left ends;
		 * 		If they are false left ends, print the overlap length they have with other nodes.
		 *  2. For each left-end node, starting from it to calculate positions for each node.
		 *  	Because the Prim algorithm starts from index 0 to generate MST, we have to
		 *  		put left-end node to index 0 in order to get the MST we want. If Prim does 
		 *  		not start from the left-end node, the directed tree will be unconnected.
		 */
		sPos = new int[alignArray.length];	//store starting positions of all the nodes
		WeightedAdjacencyListGraph primMST = null;
		for (int i=0; i<leftMostNodes.size(); i++) {
			int leftEnd = leftMostNodes.get(i).intValue();
			printLeftEndInfo(leftEnd);

			// Calculate starting positions using maximum spanning tree starting from this left-end node.
			for (int t=0; t<dGraph.length; t++) {
				//exchange index of node 0 and the left-end node so as Prim starts from the left end.
				if (dGraph[t][0] == 0) {
					dGraph[t][0] = leftEnd;
				} else if (dGraph[t][0] == leftEnd) {
					dGraph[t][0] = 0;
				}
				if (dGraph[t][1] == 0) {
					dGraph[t][1] = leftEnd;
				} else if (dGraph[t][1] == leftEnd) {
					dGraph[t][1] = 0;
				}
			}
			primMST = constructMaxTree(alignArray.length, dGraph);
			
			//put leftEnd node to index 0 in array sPos to be consistent with dGraph and primMST
			sPos[leftEnd] = sPos[0];
			sPos[0] = Integer.parseInt(g.getNameOfNode(leftEnd)); //starting position of the node
			//get starting positions for the nodes in primMST
			getStartPos(0, leftEnd, primMST);
			//exchange sPos[0] and sPos[leftEnd] to recover index 0 in sPos
			int tmp = sPos[0];
			sPos[0] = sPos[leftEnd];
			sPos[leftEnd] = tmp;

			//re-exchange index of node 0 and the left-end node to recover dGraph to its original values. 
			for (int t=0; t<dGraph.length; t++) {
				if (dGraph[t][0] == 0) {
					dGraph[t][0] = leftEnd;
				} else if (dGraph[t][0] == leftEnd) {
					dGraph[t][0] = 0;
				}
				if (dGraph[t][1] == 0) {
					dGraph[t][1] = leftEnd;
				} else if (dGraph[t][1] == leftEnd) {
					dGraph[t][1] = 0;
				}
			}
		}
		//print mst
		//System.out.println("Minimum Spanning Tree for starting positions:");
		//System.out.println(primMST);
		
	}
	
	/*
	 * print elements in 'dGraph'
	 */
	protected void printDgraph(int[][] d){
		for (int i=0; i<d.length; i++) {
			System.out.println(d[i][0] + "\t" + d[i][1] + "\t" + d[i][2]);
		}
		System.out.println();
	}
	/*
	 * print information of all the assumed left-end nodes.
	 *  	starting position of the node;
	 * 		whether or not they are real left ends;
	 * 		If they are false left ends, print the overlap length they have with other nodes.
	 */
	protected void printLeftEndInfo(int leftEnd) {
		int sp = Integer.parseInt(g.getNameOfNode(leftEnd));	//actual starting position of the node
		System.out.println("Node " + leftEnd + " starts from " + sp);
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
				if ((tmpSp+tmpLn-sp)>=(g.d2.getWindowSize())) {
					System.out.print("Node " + leftEnd + " is not a real left-most node. ");
					System.out.println("Overlap length with node " + k + " is " + (tmpSp+tmpLn-sp)); //(sp+ln-tmpSp));
					flag = 1;
				}
			}
		}
		
		if (flag == 0) {
			System.out.println("Node " + leftEnd + " is a real left-most node.");
		} else {
			System.out.println();
		}

	}
	
	/* 
	 * Calculate starting positions for each node
	 */
	protected void getStartPos(int parentNode, int leftEnd, WeightedAdjacencyListGraph tree) {
		WeightedEdgeIterator ite = (WeightedEdgeIterator) tree.edgeIterator(parentNode);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			int index = v.getIndex();
			int overlapLen = (int) Math.abs(ite.getWeight());
			if (parentNode == 0) { // it's left end node actually
				sPos[index] = sPos[parentNode] + g.getLenOfNode(leftEnd) - overlapLen;
			} else if (parentNode == leftEnd) { // it's node 0 actually
				sPos[index] = sPos[parentNode] + g.getLenOfNode(0) - overlapLen;
			} else {
				sPos[index] = sPos[parentNode] + g.getLenOfNode(parentNode) - overlapLen;
			}
			getStartPos(index, leftEnd, tree);
		}
	}
	
	/* 
	 * Print the assembled starting position and the actual position for all the ests
	 */
	public void printSPos() {
		System.out.println("Calculated s_i	Actual s_i");
		for (int i=0; i<sPos.length; i++) {
			System.out.println(sPos[i] + "	" + g.getNameOfNode(i));
		}
	}

	/* 
	 * Print all the ests to the file named as the value of "resultFile"
	 */
	public void printEsts() {
		try{ 
			File outFile = new File(resultFileName);
			boolean bExists = outFile.exists();
			if (bExists) {
				outFile.delete();
			}
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile, true));

			/*
			 * print the original sequence
			 */
			File oriFile = (new File(oriFileName));
			if (!oriFile.exists()) {
				System.out.println("SourceFile does not exist!");
				return;
			}
			BufferedReader in = new BufferedReader(new FileReader(oriFile));
			out.write(in.readLine());
			out.write("\n");
			in.close();	

			/*
			 * print all the ests
			 */
			//Firstly, sort the array sPos
			StartPos[] resultArray = new StartPos[sPos.length]; //store the starting positions of ests
			for (int i=0; i<sPos.length; i++) {
				resultArray[i] = new StartPos(sPos[i], g.getSeqOfNode(i));
			}
			MergeSort merge = new MergeSort();
			merge.sort(resultArray);
			//Secondly, print the ests
			for (int i=0; i<resultArray.length; i++) {
				for (int j=0; j<resultArray[i].pos; j++) {
					out.write(" ");
				}
				out.write(resultArray[i].seq);
				out.write("\n");
			}

			out.flush();
			out.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}

	static class StartPos implements Comparable<StartPos> {
		int pos;
		String seq;
		public StartPos(int p, String s) {
			pos = p;
			seq = s;
		}
		
		public int compareTo(StartPos other) {
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

	/* 
	 * Calculate inversions for all the calculated positions of ESTs
	 */
	public void calcInversion() {
			//Firstly, sort the array sPos
			StartPos2[] resultArray = new StartPos2[sPos.length]; //store the starting positions of ests
			for (int i=0; i<sPos.length; i++) {
				resultArray[i] = new StartPos2(sPos[i], g.getNameOfNode(i));
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
	
	/*
	 * Get consensus from the ESTs.
	 * According to the starting positions of ESTs, calculate the consensus base for each position.
	 */
	public String getConsensus() {
		//sort the array sPos
		StartPos[] resultArray = new StartPos[sPos.length]; //store the starting positions of ests
		for (int i=0; i<sPos.length; i++) {
			resultArray[i] = new StartPos(sPos[i], g.getSeqOfNode(i));
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);

		// Create an array which stores all the bases with the same position for all the ESTs.
		int lenOfArray = resultArray[resultArray.length-1].pos+resultArray[resultArray.length-1].seq.length(); 
		ArrayList<Character> [] tmpArraylists = new ArrayList [lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			tmpArraylists[i] = new ArrayList<Character>();
		}
		
		// Put all the bases into the array
		for (int i=0; i<resultArray.length; i++) {
			int tmpSPos = resultArray[i].pos;
			String tmpStr = resultArray[i].seq;
			for (int j=0; j<tmpStr.length(); j++) {
				int p = tmpSPos+j;
				if (p < lenOfArray) {
					tmpArraylists[p].add(tmpStr.charAt(j));
				}
			}
		}
		
		// Calculate consensus base for each position, and put them into an char array
		char[] consensus = new char[lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			consensus[i] = getConsensusBase(tmpArraylists[i]);
		}
		
		String retStr = String.valueOf(consensus);
		//System.out.println(String.valueOf(consensus));
		return retStr;
	}
	
	/*
	 * Calculate the consensus base from several characters
	 */
	public char getConsensusBase(ArrayList<Character> lst) {
		int numA = 0;
		int numG = 0;
		int numC = 0;
		int numT = 0;
		int numDash = 0;
		
		for (int i=0; i<lst.size(); i++) {
			switch (lst.get(i).charValue()) {
				case 'A': 
					numA++;
					break;
				case 'G':
					numG++;
					break;
				case 'C':
					numC++;
					break;
				case 'T':
					numT++;
					break;
				case '-':
					numDash++;
					break;
			}
		}
		
		int max = Math.max(Math.max(Math.max(Math.max(numA, numG), numC), numT), numDash);
		
		if (max == 0) {
			return ' ';
		}
		if (numA == max) {
			return 'A';
		} else if (numG == max) {
			return 'G';
		} else if (numC == max) {
			return 'C';
		} else if (numT == max){
			return 'T';
		} else {
			return '-';
		}
	}
	
	/*
	 * print the original sequence and consensus into a file which is specified in the property file.
	 */
	public void printConsensus() {
		try{ 
			File outFile = new File(consensusFileName);
			boolean bExists = outFile.exists();
			if (bExists) {
				outFile.delete();
			}
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile, true));

			/*
			 * print the original sequence
			 */
			File oriFile = (new File(oriFileName));
			if (!oriFile.exists()) {
				System.out.println("SourceFile does not exist!");
				return;
			}
			BufferedReader in = new BufferedReader(new FileReader(oriFile));
			out.write(in.readLine());
			out.write("\n");
			in.close();	

			/*
			 * print the consensus
			 */
			out.write(getConsensus());
			out.write("\n");
			
			out.flush();
			out.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}
	
	/*
	 * Construct a directed Maximum spanning tree.
	 * 	Assign the corresponding negative value to the weight of each edge (e.g., assign -5 to 5);
	 * 	Construct Minimum Spanning Tree.
	 *  return the tree.
	 * 
	 *  @param nOfNodes number of nodes
	 *  @param g a directed graph, the second dimension has three elements:
	 *  	index of starting node, index of ending node, weight between them.
	 *  	if there is no edge, weight=0.
	 */
	protected WeightedAdjacencyListGraph constructMaxTree(int nOfNodes, int[][] g) {
		// Make a directed graph.
		WeightedAdjacencyListGraph dGraph =
		    new WeightedAdjacencyListGraph(nOfNodes, true);
		for (int i=0; i<nOfNodes; i++) {
			dGraph.addVertex(i, Integer.toString(i));
		}
		for (int j=0; j<g.length; j++) {
			if (g[j][2] != 0) {	//there is an edge between the nodes
				dGraph.addEdge(g[j][0], g[j][1], -g[j][2]);
			}
		}
		WeightedAdjacencyListGraph mst = (new Prim()).computeMST(dGraph);
		return mst;
	}
	
	public int[][] getResults() {
		createAlignArray();
		return alignArray;
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

	
	public static void main(String[] args) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}

		ESTAssembly assemble = new ESTAssembly(props);

		assemble.readEstFile();
		assemble.createAlignArray();
		assemble.processAlignArray();
		assemble.printSPos();
	}
}
