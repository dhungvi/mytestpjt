// reconstruction from all the ESTs

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Properties;
import java.util.TreeSet;

import com.mhhe.clrs2e.MergeSort;
import com.mhhe.clrs2e.Prim;
import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

public class Reconstruction {
	Graph g;
	ArrayList<SixTuple> alignArray;
	int[] sPos;	//starting positions of all the nodes
				//the index in the array is the index of the node, the value is its starting position. 
	int[] sPosDebug;	//starting positions of all the nodes
						//it's used for debugging. All the left ends will be assigned to their actual value in order to calculate inversions later. 
	Alignment alignment;
	InclusionNodes incNodes;
	ArrayList<SixTuple> leftMostNodes;
	String consensusFileName;
	String singletonFileName;
	String printStr;
	ArrayList<String> allConsensus;	//store all the generated sequences
	ArrayList<String> numOfNodes;	//store number of used nodes, corresponds to each element in allConsensus.
	ArrayList<String> allSingletons;	//store all the singletons.
	ArrayList<String> firstEsts;	//includes all the left end sequence corresponding to each element in allConsensus.
	
	public Reconstruction(Properties props, Graph graph, ArrayList<SixTuple> align, ArrayList<SixTuple> leftEnds, InclusionNodes inc) {
		g = graph;
		alignArray = align;
		sPos = null;
		sPosDebug = null;
		incNodes = inc;
		leftMostNodes = leftEnds;
		consensusFileName = props.getProperty("ConsensusFile");
		singletonFileName = props.getProperty("SingletonFile");
		alignment = new Alignment(props);
		
		printStr = "";
		allConsensus= new ArrayList<String> ();	
		numOfNodes = new ArrayList<String> ();	
		allSingletons= new ArrayList<String> ();	
		firstEsts = new ArrayList<String> ();

	}
	
	public void getConsensus() {
		// Get the components of the time
	    long time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to reconstruct.");
		printConsensus();
		System.out.println("End to reconstruct.");
		System.out.println("The time used to reconstruct is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		//Debugger.printSPos(g, sPosDebug);
		//Debugger.calcInversion(g, sPosDebug);
	}
	
	/*
	 * Print the original sequence and multiple consensus into a file which is specified in the property file.
	 */
	private void printConsensus() {
		try{ 
			/*
			 * print consensus sequences
			 */
			File outFile1 = new File(consensusFileName);
			boolean bExists = outFile1.exists();
			if (bExists) {
				outFile1.delete();
			}
			BufferedWriter out1 = new BufferedWriter(new FileWriter(outFile1, true));
			ArrayList<String>[] results = this.reconstruct();
			ArrayList<String> consensus = results[0];
			int index = 1;
			for (int i=0; i<consensus.size(); i++) {
				String str = consensus.get(i);
				if (str.indexOf("\n") != -1) { //there is "\n" in the sequence
					String[] tStrs = str.split("\n");
					for (int j=0; j<tStrs.length; j++) {
						out1.write(">contig " + (index++) + "\n");
						out1.write(tStrs[j]);
						out1.write("\n");
					}
				} else {
					out1.write(">contig " + (index++) + "\n");
					out1.write(consensus.get(i));
					out1.write("\n");
				}
			}
			out1.flush();
			out1.close();
			
			/*
			 * print the singletons
			 */
			File outFile2 = new File(singletonFileName);
			bExists = outFile2.exists();
			if (bExists) {
				outFile2.delete();
			}
			BufferedWriter out2 = new BufferedWriter(new FileWriter(outFile2, true));
			ArrayList<String> singletons = results[1];
			for (int i=0; i<singletons.size(); i++) {
				out2.write(singletons.get(i));
				out2.write("\n");
			}
			out2.flush();
			out2.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}
	
	/*
	 * reconstruct the sequence
	 * 
	 * 1. Generate dGraph. 
	 * 2. For each left-end node, starting from it to calculate positions for each node. 
	 * In order to get starting positions, it constructs a MST. The weight of the Minimum Spanning tree is 
	 * the overlap distance instead of overlap length. 
	 * Then reconstruct the sequence from the set of ESTs.
	 * 
	 * @return two arraylist: The assembled sequences, the singletons.
	 */
	private ArrayList<String>[] reconstruct() {
		int[][]dGraph = genDGraph();
		
		sPosDebug = new int[g.graphNodes.size()];	
		for (int i=0; i<leftMostNodes.size(); i++) { //start for
			genConsensusFromOneLeftEnd(leftMostNodes.get(i).curNode, dGraph);
		} //end for
		
		//if there are more than one consensus, process them.
		int tmpSize = allConsensus.size();
		if (tmpSize > 1) { 
			ArrayList<String> s = processMoreConsensus();
			allConsensus = s;
			
			//for debug
			printStr = printStr + "The consensus from above " + tmpSize + " sequences:\n";
			for (int p=0; p<s.size(); p++) {
				printStr = printStr + s.get(p) + "\n";
			}
		}
		
		//print debug information about the generated consensus.
		System.out.println("*********************consensus:********************");
		System.out.println(printStr);
		
		ArrayList<String>[] rets = new ArrayList[2];
		rets[0] = allConsensus;
		rets[1] = allSingletons;
		return rets;
	}
	
	private int[][] genDGraph() {
		/*
		 * Calculate the length of dGraph.
		 */

		int len = 0;
		for (int i=0; i<alignArray.size(); i++) {
			if (alignArray.get(i).leftNode != -1) {
				len++;
			}
			if (alignArray.get(i).rightNode != -1) {
				len++;
			}
		}
		/*
		 * generate dGraph.
		 */
		int[][] dGraph = new int[len][4];
		int indexOfDGraph = 0;
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode != -1) {
				dGraph[indexOfDGraph][0] = curTuple.leftNode;
				dGraph[indexOfDGraph][1] = curTuple.curNode;
				dGraph[indexOfDGraph][2] = Math.abs(curTuple.lDis);	//distance
				dGraph[indexOfDGraph][3] = Math.abs(curTuple.lOvlLen);	//overlap length
				indexOfDGraph++;
			}
			
			if (curTuple.rightNode != -1) {
				dGraph[indexOfDGraph][0] = curTuple.curNode;
				dGraph[indexOfDGraph][1] = curTuple.rightNode;
				dGraph[indexOfDGraph][2] = curTuple.rDis;	//distance
				dGraph[indexOfDGraph][3] = curTuple.rOvlLen;	//overlap length
				indexOfDGraph++;
			}
		}

		//print dGraph
		//System.out.println("dGraph:");
		//Debugger.printDgraph(dGraph);
		
		return dGraph;
	}


	/*
	 *  1. Print information of the left-end node. It include:
	 *  	starting position of the node;
	 * 		whether or not they are real left ends;
	 * 		If they are false left ends, print the overlap length they have with other nodes.
	 *  2. For the left-end node, starting from it to calculate positions for each node.
	 *  	Because the Prim algorithm starts from index 0 to generate MST, we have to
	 *  		put left-end node to index 0 in order to get the MST we want. If Prim does 
	 *  		not start from the left-end node, the directed tree will be unconnected.

	 */
	private void genConsensusFromOneLeftEnd(int leftEnd, int[][] dGraph) {
		sPos = new int[g.graphNodes.size()];	//store starting positions of all the nodes

		//Debugger.printLeftEndInfo(leftEnd, g);
		
		// Calculate starting positions using minimum spanning tree starting from this left-end node.
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

		WeightedAdjacencyListGraph primMST = constructMinTree(g.graphNodes.size(), dGraph); //the first param is the total number of ESTs.

		//put leftEnd node to index 0 in array sPos to be consistent with dGraph and primMST
		sPos[leftEnd] = sPos[0];
			sPosDebug[leftEnd] = sPosDebug[0];
		sPos[0] = 0; //starting position of the left end is assigned to be 0.
			sPosDebug[0] = Integer.parseInt(g.getNameOfNode(leftEnd));
		//get starting positions for the nodes in primMST
		getStartPos(0, leftEnd, primMST, dGraph);
			getStartPosDebug(0, leftEnd, primMST, dGraph);
		//exchange sPos[0] and sPos[leftEnd] to recover index 0 in sPos
		int tmp = sPos[0];
			int tmp1 = sPosDebug[0];
		sPos[0] = sPos[leftEnd];
			sPosDebug[0] = sPosDebug[leftEnd];
		sPos[leftEnd] = tmp;
			sPosDebug[leftEnd] = tmp1;

		
		//reconstruct this sequence 
		//sort the array sPos (in ascending order)
		ArrayList<StartPos> tmpArray = new ArrayList<StartPos> ();
		for (int j=0; j<sPos.length; j++) {
			if ((j == leftEnd) || (sPos[j] != 0)) {
				tmpArray.add(new StartPos(sPos[j], j));
			}
		}

		if (tmpArray.size() == 1) { //singleton
			allSingletons.add(g.getCommentOfNode(leftEnd) + "\n" + g.getSeqOfNode(leftEnd));
			return;
		} 
		
		firstEsts.add(g.getSeqOfNode(leftEnd));

		String[] tStr = reconstructSeq(tmpArray);
		allConsensus.add(tStr[0]);
		numOfNodes.add(tStr[1]);
		printStr = printStr + tStr[0] + "\n";

		
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
	
	
	 /* 
	 * This method is used when there are more than one consensus in the assembly.
	 * Then number of consensus is equal to the number of left ends. And the different consensus does not mean that they are 
	 * not overlapped(or they correspond to different part of the gene). Some left ends may include each other, and consensus 
	 * from them overlap with each other.
	 * This method is designed to remove all the dependent consensus and extract all the independent ones, which means, we intend 
	 * to find all the consensus which represent different part of the gene.
	 * 
	 * 			allConsensus an arraylist which includes all the generated consensus from the calling method;
	 * 		 	firstEsts an arraylist which includes all the left end sequence corresponding to all the consensus;
	 * 			numOfNodes an arraylist which stores the number of used nodes for each consensus.
	 * @return the combined consensus.
	 */
	 private ArrayList<String> processMoreConsensus() {
		//String retStr = "";
		ArrayList<String> allOutputContigs= new ArrayList<String> ();	//store all the generated sequences
		
		int sizeOfs = allConsensus.size();
		Consensus[] resultArray = new Consensus[sizeOfs]; //store the starting positions of ests
		for (int i=0; i<sizeOfs; i++) {
			resultArray[i] = new Consensus(Integer.parseInt(numOfNodes.get(i)), allConsensus.get(i), firstEsts.get(i));
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		ArrayList<Consensus> allConsensus = new ArrayList<Consensus> ();
		for (int i=sizeOfs-1; i>=0; i--) {
			allConsensus.add(resultArray[i]);
		}
		
		while (true) {
			String s1 = allConsensus.get(0).firstEst;
			ArrayList<Consensus> includeStrs = new ArrayList<Consensus>();
			includeStrs.add(allConsensus.get(0));
			ArrayList<Consensus> excludeStrs = new ArrayList<Consensus>();
			for (int i=1; i<allConsensus.size(); i++) {
				boolean b = g.ovl.checkInclusion(allConsensus.get(i).firstEst, s1); //if resultArray[i].firstEst is included in s1
				if (b) {
					includeStrs.add(allConsensus.get(i));
				} else {
					excludeStrs.add(allConsensus.get(i));
				}
			}

			allOutputContigs.add(processMoreConsensusWithInclusion(includeStrs));
			if (excludeStrs.size() == 0) {
				break;
			} else {
				allConsensus = excludeStrs;
			}
		}


		return allOutputContigs;
	 }
	 


	/*
	 * This method is called by "processMoreConsensus".
	 * This method is used to process those consensus which starts from the left ends that include each other.
	 * These input consensus originate from more than one left ends which include each other. 
	 * 
	 * Here we combine all the consensus into one.
	 * We combine the one with the longest first EST and the one using the most number of ESTs together to 
	 * form a new one and return it.
	 */
	private String processMoreConsensusWithInclusion(ArrayList<Consensus> includeStrs) {
		int maxLen = 0; //the maximal length of the first EST.
		int idxMaxLen = 0;
		int maxNumNodes = 0;
		int idxMaxNumNodes = 0;
		for (int i=0; i<includeStrs.size(); i++) {
			int tLen = includeStrs.get(i).lenOfFirstEst;
			if (tLen > maxLen) {
				maxLen = tLen;
				idxMaxLen = i;
			}
			
			int num = includeStrs.get(i).numOfUsedNodes;
			if (num > maxNumNodes) {
				maxNumNodes = num;
				idxMaxNumNodes = i;
			}
		}
		
		if (includeStrs.size() == 0) {
			return "";
		} else if (includeStrs.size() == 1) {
			return includeStrs.get(0).seq;
		} else {
			String s1 = includeStrs.get(idxMaxLen).seq;
			String s2 = includeStrs.get(idxMaxNumNodes).seq;
			String[] strs = alignment.getLocalAlignment(s1, s2);
			int offset = s1.indexOf(strs[0].replace("-", ""));
			return (s1.substring(0, offset) + s2);
		}
	}
	
	/*
	 * reconstruct a sequence which starts from a left end.
	 * return: ret[0]-the consensus, ret[1]-the number of nodes used for reconstruction.
	 */
	public String[] reconstructSeq(ArrayList<StartPos> a) {
		String[] ret = new String[2];
		int sizeOfa = a.size();
		if (sizeOfa == 0) {
			return null;
		} else if (sizeOfa == 1) {
			 ret[0] = g.getSeqOfNode(a.get(0).index);
			 ret[1] = Integer.toString(1);
			 return ret;
		}
		
		StartPos[] tmpResultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			tmpResultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(tmpResultArray);
		
		TreeSet<UsedNode> addedNodes = addInclusionNodes(tmpResultArray);  //add all those related inclusion nodes into it for reconstruction.
		System.out.println(addedNodes.size() + " nodes are used to reconstruct the sequence.\n");
		ret[1] = Integer.toString(addedNodes.size());

		StartPos[] resultArray = new StartPos[addedNodes.size()]; 
		Object[] r1 = addedNodes.toArray();
		for (int i=0; i<resultArray.length; i++) {
			UsedNode tmpNode = (UsedNode) r1[i];
			resultArray[i] = new StartPos(tmpNode.pos, tmpNode.index);
		}
		merge.sort(resultArray);

/*		System.out.println("These nodes are:");
		for (int r=0; r<resultArray.size(); r++) {
			System.out.print(resultArray.get(r).index+"  ");
		}
		System.out.println();
*/		
		ArrayList<SingleBase> bases = new ArrayList<SingleBase> ();
		String tConsensus = g.getSeqOfNode(resultArray[0].index);
		String curSeq = "";
		int len = resultArray.length - 1;
		for (int i=1; i<=len; i++) {
			curSeq = g.getSeqOfNode(resultArray[i].index);
			String[] strs = alignment.getLocalAlignment(tConsensus, curSeq);
			tConsensus = tConsensus.replace(strs[0].replace("-", ""), strs[0]);
			int offset = tConsensus.indexOf(strs[0]);
			
			String tSeq = curSeq.replace(strs[1].replace("-", ""), strs[1]);
			curSeq = tSeq.substring(tSeq.indexOf(strs[1]));
			
			if (i == 1) {
				int len1 = tConsensus.length();
				int len2 = curSeq.length();
				int end = Math.max(offset+len2, len1);
				for (int j=0; j<end; j++) {
					if ((j < len1) && (j-offset >= 0) && (j-offset < len2)) { //overlap part
						bases.add(new SingleBase(tConsensus.charAt(j), curSeq.charAt(j-offset)));
					} else if ((j-offset < 0) || (j-offset >= len2)) {
						bases.add(new SingleBase(tConsensus.charAt(j)));
					} else if (j >= len1) {
						bases.add(new SingleBase(curSeq.charAt(j-offset)));
					} 
				}
			
			} else {
				int len1 = tConsensus.length();
				int len2 = curSeq.length();
				int end = Math.max(offset+len2, len1);
				for (int j=offset; j<end; j++) {
					if ((j < len1) && (j-offset < len2)) { //overlap part
						char c1 = tConsensus.charAt(j);
						char c2 = curSeq.charAt(j-offset);
						if (c1 != '-') {
							bases.get(j).addOneBase(c2);
						} else {
							bases.add(j, new SingleBase(c1, c2));
						}
					} else if (j >= len1) {
						bases.add(new SingleBase(curSeq.charAt(j-offset)));
					} 
				}
			}
			
			tConsensus = getCurConsensus(bases);
		}
		
		ret[0]= tConsensus.replace("P", "");
		return ret;
	}

	private String getCurConsensus(ArrayList<SingleBase> bases) {
		int len = bases.size();
		StringBuffer tStr = new StringBuffer();
		for (int i=0; i<len; i++) {
			tStr.append(bases.get(i).getCurBase());
		}
		return tStr.toString();
	}
	
	/*
	 * Add all the inclusion nodes into the input arraylist.
	 * For each element in the arraylist, put its corresponding node just after it. 
	 */
	private TreeSet<UsedNode> addInclusionNodes(StartPos[] input) {
		TreeSet<UsedNode> retList = new TreeSet<UsedNode>();
		int size = input.length;
		
		for (int i=0; i<size; i++) {
			int curIdx = input[i].index;
			int pos = input[i].pos;
			retList.add(new UsedNode(curIdx, pos));
			
			int[] chdIdx = incNodes.containPNode(curIdx, g.graphNodes.size()); //inclusion children index of the curIdx if exist.
			if (chdIdx != null) {
				for (int j=0; j<chdIdx.length; j++) {
					retList.add(new UsedNode(chdIdx[j], pos+1));
				}
			}
		}
		return retList;
	}
	
	/*
	 * Construct a directed Miminum spanning tree.
	 * 
	 *  @param nOfNodes number of nodes
	 *  @param g a directed graph, the second dimension has three elements:
	 *  	index of starting node, index of ending node, weight between them.
	 */
	private WeightedAdjacencyListGraph constructMinTree(int nOfNodes, int[][] g) {
		// Make a directed graph.
		WeightedAdjacencyListGraph dGraph =
		    new WeightedAdjacencyListGraph(nOfNodes, true);
		for (int i=0; i<nOfNodes; i++) {
			dGraph.addVertex(i, Integer.toString(i));
		}
		for (int j=0; j<g.length; j++) {
			if (g[j][3] != 0) {	//there is an edge between the nodes
				dGraph.addEdge(g[j][0], g[j][1], g[j][2]);
			}
		}
		WeightedAdjacencyListGraph mst = (new Prim()).computeMST(dGraph);
		return mst;
	}

	
	/* 
	 * Calculate starting positions for each node. 
	 */
	
	private void getStartPos(int parentNode, int leftEnd, WeightedAdjacencyListGraph tree, int[][] d) {
		WeightedEdgeIterator ite = (WeightedEdgeIterator) tree.edgeIterator(parentNode);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			
			int index = v.getIndex();
			
			int overlapLen = 0;
			for (int i=0; i<d.length; i++) {
					if ((d[i][0] == parentNode) && (d[i][1] == index)) {
						overlapLen = d[i][3];
						break;
					}
			}
			
			if (parentNode == 0) { // it's left end node actually
				sPos[index] = sPos[parentNode] + g.getLenOfNode(leftEnd) - overlapLen;
			} else if (parentNode == leftEnd) { // it's node 0 actually
				sPos[index] = sPos[parentNode] + g.getLenOfNode(0) - overlapLen;
			} else {
				sPos[index] = sPos[parentNode] + g.getLenOfNode(parentNode) - overlapLen;
			}
			getStartPos(index, leftEnd, tree, d);
		}
	}

	/* 
	 * Used for debugging. 
	 */
	
	private void getStartPosDebug(int parentNode, int leftEnd, WeightedAdjacencyListGraph tree, int[][] d) {
		WeightedEdgeIterator ite = (WeightedEdgeIterator) tree.edgeIterator(parentNode);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			
			int index = v.getIndex();
			
			int overlapLen = 0;
			for (int i=0; i<d.length; i++) {
					if ((d[i][0] == parentNode) && (d[i][1] == index)) {
						overlapLen = d[i][3];
						break;
					}
			}
			
			if (parentNode == 0) { // it's left end node actually
				sPosDebug[index] = sPosDebug[parentNode] + g.getLenOfNode(leftEnd) - overlapLen;
			} else if (parentNode == leftEnd) { // it's node 0 actually
				sPosDebug[index] = sPosDebug[parentNode] + g.getLenOfNode(0) - overlapLen;
			} else {
				sPosDebug[index] = sPosDebug[parentNode] + g.getLenOfNode(parentNode) - overlapLen;
			}
			getStartPosDebug(index, leftEnd, tree, d);
		}
	}


	class UsedNode implements Comparable<UsedNode> {
		int index; //index of the node
		int pos;;
		public UsedNode(int idx, int p) {
			index = idx;
			pos = p;
		}
		
		public int compareTo(UsedNode other) {
			//Returns 0 if the argument is equal to this; 			
			//a value less than 0 if the argument is greater than this; 
			//and a value greater than 0 if the argument is less than this. 
			if (this.index == other.index) {
				return 0;
			} else if (this.index > other.index) {
				return 1;
			} else {
				return -1;
			}
		}
	}
	

	class StartPos implements Comparable<StartPos> {
		int pos;
		int index; //index of the node
		public StartPos(int p, int idx) {
			pos = p;
			index = idx;
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

	class Consensus implements Comparable<Consensus> {
		 int lenOfFirstEst;
		 int numOfUsedNodes;
		 String seq;
		 String firstEst;
		 public Consensus(int n, String s1, String s2) {
			 lenOfFirstEst = s2.length();
			 numOfUsedNodes = n;
			 seq = s1;
			 firstEst = s2;
		 }

		 public int compareTo(Consensus other) {
			 //Returns 0 if the argument is equal to this; 			
			 //a value less than 0 if the argument is greater than this; 
			 //and a value greater than 0 if the argument is less than this. 
			 if (this.lenOfFirstEst == other.lenOfFirstEst) {
				 return 0;
			 } else if (this.lenOfFirstEst > other.lenOfFirstEst) {
				 return 1;
			 } else {
				 return -1;
			 }
		 }
	}	
}



