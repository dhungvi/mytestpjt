import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Properties;

import com.mhhe.clrs2e.*;

public class ESTAssembly {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected final int INT_MIN = Integer.MIN_VALUE;
	String oriFileName;	//store the original gene data. 
	String mstFile;	//store the MST. It is generated from PEACE.
	String inFileName;
	String resultFileName;
	String consensusFileName;
	String singletonFileName;
	ArrayList<String> ests;	//store all the ests. It is generated in 'readEstFile' function.
	Graph g;	//graph to store all the ests. It is generated in 'readEstFile' function.
	WeightedAdjacencyListGraph mstForG;	//minimum spanning tree generated from 'g'
										//It is generated in 'createAlignArray' function.
	ArrayList<SixTuple> leftMostNodes;
	InclusionNodes incNodes;
	Alignment alignment;

	/*
	 * store the position of aligned nodes
	 * index of nodes in graph;
	 * the first is the index of node on the left, 
	 * the second is the overlap length with + or -.
	 * the third is the distance with + or -.
	 * the fourth is the index of node on the right, 
	 * the fifth is the overlap length with + or -.
	 * the sixth is the distance with + or -.
	 */
	ArrayList<SixTuple> alignArray;
	int[] sPos;	//starting positions of all the nodes, initialized in "processAlignArray" function.
				//the index in the array is the index of the node, the value is its starting position. 
	int[] sPosDebug;	//starting positions of all the nodes, initialized in "processAlignArray" function.
	//it's used for debugging. All the left ends will be assigned to their actual value in order to calculate inversions later. 
	
	public ESTAssembly(Properties props) {
		oriFileName = props.getProperty("SourceFile");
		mstFile = props.getProperty("MSTFile");
		inFileName = props.getProperty("EstFile");
		//resultFileName = props.getProperty("ResultFile");
		consensusFileName = props.getProperty("ConsensusFile");
		singletonFileName = props.getProperty("SingletonFile");
		ests = new ArrayList<String> ();
		alignArray = null;
		sPos = null;
		sPosDebug = null;
		g = new Graph(props);
		leftMostNodes = new ArrayList<SixTuple> ();
		incNodes = new InclusionNodes();
		alignment = new Alignment(props);
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
					//ests.add("0"); //starting position
					ests.add(str); //comment
					
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
			g.addNode(new Node(ests.get(i), ests.get(i+1), ests.get(i+2).toUpperCase()));
			i = i+3;
		}
	}
	
	public void createAlignArray() {
		mstForG = g.readMST(mstFile);
		System.out.println("End to generate MST.");
		System.out.println("Start to generate 6-tuples.");
		alignArray = g.get2CloseNodesFromMST(mstForG, incNodes);
		System.out.println("There are " + incNodes.getSize() + " nodes in the inclusion list.");
		//incNodes.printAllNodes();
		System.out.println("There are " + alignArray.size() + " nodes in alignArray.");
	}

	/* 
	 * Get the assumed left ends(alignNodes[][0]==-1), check them to find all the real left ends.
	 */
	public void processAlignArray() {
		ArrayList<SixTuple> rightMostNodes = new ArrayList<SixTuple> ();
		//get all the nodes which has no left nodes to them
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode == -1) {
				leftMostNodes.add(curTuple);	//store sixtuple of the node
			}
			if (curTuple.rightNode == -1) {
				rightMostNodes.add(curTuple);	//store sixtuple of the node
			}
		}

		//Re-calculating six-tuples for those assumed left ends and put the new six-tuple into alignArray.
		int numOfLeftMostNodes = leftMostNodes.size();
		System.out.println("The number of original left ends is " + numOfLeftMostNodes);
		System.out.println("The number of original right ends is " + rightMostNodes.size());
		if (numOfLeftMostNodes > 1) {
			for (int i=0; i<leftMostNodes.size(); i++) {
				SixTuple curTuple = leftMostNodes.get(i);
				int cNode = curTuple.curNode;
				SixTuple lNode = g.get2CloseNodesFromGrand(mstForG, cNode, curTuple, incNodes);
				curTuple.leftNode = lNode.leftNode;
				curTuple.lOvlLen = lNode.lOvlLen;	//overlap length
				curTuple.lDis = lNode.lDis;	//distance
				if (curTuple.rDis > lNode.rDis) { //get a smaller distance
					curTuple.rightNode = lNode.rightNode;
					curTuple.rOvlLen = lNode.rOvlLen;
					curTuple.rDis = lNode.rDis;
				}
			}
		} 
		for (int i=0; i<rightMostNodes.size(); i++) {
			SixTuple curTuple = rightMostNodes.get(i);
			int cNode = curTuple.curNode;
			SixTuple lNode = g.get2CloseNodesFromGrand(mstForG, cNode, curTuple, incNodes);
			if (lNode.rightNode != -1) {
				curTuple.rightNode = lNode.rightNode;
				curTuple.rOvlLen = lNode.rOvlLen;
				curTuple.rDis = lNode.rDis;
			}
		}
		
		leftMostNodes.clear();
		rightMostNodes.clear();
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode == -1) {
				leftMostNodes.add(curTuple);	//store sixtuple of the node
			}
			if (curTuple.rightNode == -1) {
				rightMostNodes.add(curTuple);	//store sixtuple of the node
			}
		}
		System.out.println("There are " + leftMostNodes.size() + " left-most nodes after running 3 levels.");
		System.out.println("There are " + rightMostNodes.size() + " right-most nodes after running 3 levels.");
		
		/* Recalculate 6-tuples for all the current left nodes in order to remove all the false left ends.
		 * Specifically, for those assumed left ends,start to calculate from fourth level until meeting one node which 
		 * makes six-tuple[0] != -1 or until the level we specified in the property file, then return the six-tuple.
		 * If we fail to find any node, we consider it a real left end.
		 */
		for (int i=0; i<leftMostNodes.size(); i++) {
			SixTuple curTuple = leftMostNodes.get(i);
			int tEnd = curTuple.curNode; //index of the node
			SixTuple tmpTuple = g.checkLeftEndFromMST(mstForG, tEnd, curTuple, incNodes);
			if (tmpTuple != null) {
				curTuple.leftNode = tmpTuple.leftNode;
				curTuple.lOvlLen = tmpTuple.lOvlLen;	//overlap length
				curTuple.lDis = tmpTuple.lDis;	//distance
				curTuple.rightNode = tmpTuple.rightNode;
				curTuple.rOvlLen = tmpTuple.rOvlLen;
				curTuple.rDis = tmpTuple.rDis;
			}
		}
		for (int i=0; i<rightMostNodes.size(); i++) {
			SixTuple curTuple = rightMostNodes.get(i);
			int tEnd = curTuple.curNode; //index of the node
			SixTuple tmpTuple = g.checkRightEndFromMST(mstForG, tEnd, curTuple, incNodes);
			if (tmpTuple != null) {
				curTuple.rightNode = tmpTuple.rightNode;
				curTuple.rOvlLen = tmpTuple.rOvlLen;
				curTuple.rDis = tmpTuple.rDis;
			}
		}
		
		leftMostNodes.clear();
		rightMostNodes.clear();
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode == -1) {
				leftMostNodes.add(curTuple);	//store sixtuple of the node
			}
			if (curTuple.rightNode == -1) {
				rightMostNodes.add(curTuple);	//store sixtuple of the node
			}
		}
		System.out.println("\nThere are " + leftMostNodes.size() + " left-most nodes after checking left ends.");
		System.out.println("\nThere are " + rightMostNodes.size() + " right-most nodes after checking right ends.");
		
		
		
		
		/* construct a temporary directed graph from alignArray, 
		 *  	the second dimension has two elements:
		 *  			index of starting node,
		 *  			index of ending node, 
		 *  			weight between them (positive value, weight is abs(their distance)).
		 *  	if there is no edge, weight=INT_MAX.
		 */
		int tLen = 0;
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode != -1) {
				tLen++;
			}
			if (curTuple.rightNode != -1) {
				tLen++;
			}
		}
		
		int[][] tmpDGraph = new int[tLen][4];
		int tmpIndex = 0;
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			int curIdx = curTuple.curNode;
			if (curTuple.leftNode != -1) {
				tmpDGraph[tmpIndex][0] = curTuple.leftNode;
				tmpDGraph[tmpIndex][1] = curIdx;
				tmpDGraph[tmpIndex][2] = Math.abs(curTuple.lDis);	//distance
				tmpDGraph[tmpIndex][3] = Math.abs(curTuple.lOvlLen);	//overlap length
				tmpIndex++;
			}
			
			if (curTuple.rightNode != -1) {
				tmpDGraph[tmpIndex][0] = curIdx;
				tmpDGraph[tmpIndex][1] = curTuple.rightNode;
				tmpDGraph[tmpIndex][2] = curTuple.rDis;	//distance
				tmpDGraph[tmpIndex][3] = curTuple.rOvlLen;	//overlap length
				tmpIndex++;
			}
		}
		
		/*
		 * Remove those false left ends. For example:
		 * Node 2 has the set in alignArray: [-1, 0, 5, 4, 8, 9]
		 * but Node 8 has this set: [7, 5, 3, 2, 7, 6], this means ovlDis(8,2)=6, 
		 * 		so node 2 is not left end because node 8 is to its left.
		 */
		//Get all the nodes which has the value of -1 in alignArray[x][0]
		
		ArrayList<SixTuple> tmpLeftNodes = new ArrayList<SixTuple> ();
		for (int i=0; i<alignArray.size(); i++) {
			SixTuple curTuple = alignArray.get(i);
			if (curTuple.leftNode == -1) {
				tmpLeftNodes.add(curTuple);	//store index of the node
			}
		}
		//remove false left ends
		leftMostNodes.clear();
		for (int i=0; i<tmpLeftNodes.size(); i++) {
			int tEnd = tmpLeftNodes.get(i).curNode;
			int f = 0;
			//if the left end appears in second element of dGraph, that means some 
			//node is on its left, so it is not a real left end.
			for (int j=0; j<tmpDGraph.length; j++) {
				if (tmpDGraph[j][1] == tEnd) {	// false left end
					f = 1;
					break;
				}
			}
			if (f == 0) {
				leftMostNodes.add(tmpLeftNodes.get(i));
			}
		}
		System.out.println("\nThere are " + leftMostNodes.size() + " left-most nodes after processing.");
	}

	

	/*
	 * print elements in 'dGraph'
	 */
	private void printDgraph(int[][] d){
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
	private void printLeftEndInfo(int leftEnd) {
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
				if ((tmpSp+tmpLn-sp)>=(g.ovl.d2.getWindowSize())) {
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
	 * Print the assembled starting position and the actual position for all the ests
	 */
	public void printSPos() {
		System.out.println("Calculated s_i	Actual s_i	sPos");
		for (int i=0; i<sPosDebug.length; i++) {
			System.out.println(sPosDebug[i] + "	" + g.getNameOfNode(i) + "	" + sPos[i]);
		}
	}


	static class StartPos implements Comparable<StartPos> {
		int pos;
		int index; //index of the node
		String seq;
		public StartPos(int p, String s, int idx) {
			pos = p;
			seq = s;
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

	/* 
	 * Calculate inversions for all the calculated positions of ESTs
	 */
	public void calcInversion() {
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
	public ArrayList<String>[] reconstruct() {
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
		System.out.println("dGraph:");
		printDgraph(dGraph);
		
		/*
		 *  1. Print information of all the left-end nodes. It include:
		 *  	starting position of the node;
		 * 		whether or not they are real left ends;
		 * 		If they are false left ends, print the overlap length they have with other nodes.
		 *  2. For each left-end node, starting from it to calculate positions for each node.
		 *  	Because the Prim algorithm starts from index 0 to generate MST, we have to
		 *  		put left-end node to index 0 in order to get the MST we want. If Prim does 
		 *  		not start from the left-end node, the directed tree will be unconnected.
		 *     Then reconstruct the sequence from the set of ESTs.
		 *     Return all the generated sequences which are separated by newline character.
		 */
		WeightedAdjacencyListGraph primMST = null;
		String printStr = "";
		ArrayList<String> allConsensus= new ArrayList<String> ();	//store all the generated sequences
		ArrayList<String> numOfNodes = new ArrayList<String> ();	//store number of used nodes, corresponds to the element in allConsensus.
		ArrayList<String> allSingletons= new ArrayList<String> ();	//store all the singletons
		ArrayList<String> firstEsts = new ArrayList<String> ();
		
		sPosDebug = new int[g.graphNodes.size()];	
		for (int i=0; i<leftMostNodes.size(); i++) { //start for
			sPos = new int[g.graphNodes.size()];	//store starting positions of all the nodes

			int leftEnd = leftMostNodes.get(i).curNode;
			printLeftEndInfo(leftEnd);
			
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
			primMST = constructMinTree(g.graphNodes.size(), dGraph); //the first param is the total number of ESTs.

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
					tmpArray.add(new StartPos(sPos[j], g.getSeqOfNode(j), j));
				}
			}

			if (tmpArray.size() == 1) { //singleton
				allSingletons.add(g.getCommentOfNode(leftEnd) + "\n" + g.getSeqOfNode(leftEnd));
				continue;
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
		} //end for
		
		//if there are more than one consensus, process them.
		int tmpSize = allConsensus.size();
		if (tmpSize > 1) { 
			printStr = printStr + "The consensus from above " + tmpSize + " sequences:\n";
			ArrayList<String> s = processMoreConsensus(allConsensus, firstEsts, numOfNodes);
			for (int p=0; p<s.size(); p++) {
				printStr = printStr + s.get(p) + "\n";
			}
			allConsensus = s;
		}
		
		//print debug information about the generated consensus.
		System.out.println("*********************consensus:********************");
		System.out.println(printStr);
		
		ArrayList<String>[] rets = new ArrayList[2];
		rets[0] = allConsensus;
		rets[1] = allSingletons;
		return rets;
	}
	
	 /* 
	 * This method is used when there are more than one consensus in the assembly.
	 * Then number of consensus is equal to the number of left ends. And the different consensus does not mean that they are 
	 * not overlapped(or they correspond to different part of the gene). Some left ends may include each other, and consensus 
	 * from them overlap with each other.
	 * This method is designed to remove all the dependent consensus and extract all the independent ones, which means, we intend 
	 * to find all the consensus which represent different part of the gene.
	 * 
	 * @param s an arraylist which includes all the generated consensus from the calling method;
	 * 		 	firstEsts an arraylist which includes all the left end sequence corresponding to all the consensus;
	 * 			numOfNodes an arraylist which stores the number of used nodes for each consensus.
	 * @return the combined consensus.
	 */
	 private ArrayList<String> processMoreConsensus(ArrayList<String> s, ArrayList<String> firstEsts, ArrayList<String> numOfNodes) {
		//String retStr = "";
		ArrayList<String> allOutputContigs= new ArrayList<String> ();	//store all the generated sequences
		
		int sizeOfs = s.size();
		Consensus[] resultArray = new Consensus[sizeOfs]; //store the starting positions of ests
		for (int i=0; i<sizeOfs; i++) {
			resultArray[i] = new Consensus(firstEsts.get(i).length(), Integer.parseInt(numOfNodes.get(i)), s.get(i), firstEsts.get(i));
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
			//retStr = retStr  + "\n" + processMoreConsensusWithInclusion(includeStrs);
			allOutputContigs.add(processMoreConsensusWithInclusion(includeStrs));
			if (excludeStrs.size() == 0) {
				break;
			} else {
				allConsensus = excludeStrs;
			}
		}


		return allOutputContigs;
	 }
	 
	 static class Consensus implements Comparable<Consensus> {
		 int lenOfFirstEst;
		 int numOfUsedNodes;
		 String seq;
		 String firstEst;
		 public Consensus(int p, int n, String s1, String s2) {
			 lenOfFirstEst = p;
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


	/*
	 * This method is called by "processMoreConsensus".
	 * This method is used to process those consensus which starts from the left ends that include each other.
	 * These input consensus originate from more than one left ends which include each other. 
	 * 
	 * Here we combine all the consensus into one.
	 * If there are two in includeStrs, we return the longer one. If there are more than or equal to three, 
	 * we combine them by calling "reconstructSeq".
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
		} /*else if (includeStrs.size() == 2) {
			int i1 = includeStrs.get(0).numOfUsedNodes;
			int i2 = includeStrs.get(1).numOfUsedNodes;
			if (i1 >= i2) {
				return includeStrs.get(0).seq;
			} else {
				return includeStrs.get(1).seq;
			}
		}*/ else {
/*			ArrayList<StartPos> input = new ArrayList<StartPos>();
			for (int i=0; i<includeStrs.size(); i++) {
				input.add(new StartPos(maxLen-includeStrs.get(i).lenOfFirstEst, includeStrs.get(i).seq, -1));
			}
			return reconstructSeq(input);
*/		
			String s1 = includeStrs.get(idxMaxLen).seq;
			String s2 = includeStrs.get(idxMaxNumNodes).seq;
			String[] strs = alignment.getLocalAlignment(s1, s2);
			int offset = s1.indexOf(strs[0].replace("-", ""));
			return (s1.substring(0, offset) + s2);
		}
		
		
/*		//return the longest consensus
		int idxMaxLen = 0; //the index of the longest consensus.
		int maxLen = 0;
		for (int i=0; i<includeStrs.size(); i++) {
			int tLen = includeStrs.get(i).seq.length();
			if (tLen > maxLen) {
				maxLen = tLen;
				idxMaxLen = i;
			}
		}
		return includeStrs.get(idxMaxLen).seq;
*/		
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
	 * reconstruct a sequence which starts from a left end.
	 * return: 0-the consensus, 1-the number of nodes used for reconstruction.
	 */
	public String[] reconstructSeq(ArrayList<StartPos> a) {
		String[] ret = new String[2];
		int sizeOfa = a.size();
		if (sizeOfa == 0) {
			return null;
		} else if (sizeOfa == 1) {
			 ret[0] = a.get(0).seq;
			 ret[1] = Integer.toString(1);
			 return ret;
		}
		StartPos[] tmpResultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			tmpResultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(tmpResultArray);
		
		ArrayList<StartPos> resultArray = new ArrayList<StartPos> ();
		for (int i=0; i<tmpResultArray.length; i++) {
			resultArray.add(tmpResultArray[i]);
		}
		resultArray = addInclusionNodes(resultArray);  //add all those related inclusion nodes into it for reconstruction.
		
		System.out.println(resultArray.size() + " nodes are used to reconstruct the sequence.\n");
		ret[1] = Integer.toString(resultArray.size());
/*		System.out.println("These nodes are:");
		for (int r=0; r<resultArray.size(); r++) {
			System.out.print(resultArray.get(r).index+"  ");
		}
		System.out.println();
*/		
		ArrayList<SingleBase> bases = new ArrayList<SingleBase> ();
		String tConsensus = resultArray.get(0).seq;
		String curSeq = "";
		int len = resultArray.size() - 1;
		for (int i=1; i<=len; i++) {
			curSeq = resultArray.get(i).seq;
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

	/*
	 * Add all the inclusion nodes into the input arraylist.
	 * For each element in the arraylist, put its corresponding node just after it. 
	 */
	private ArrayList<StartPos> addInclusionNodes(ArrayList<StartPos> input) {
		ArrayList<StartPos> retList = input;
		int size = retList.size();
		ArrayList<StartPos> tmpList = new ArrayList<StartPos>();
		for (int i=0; i<size; i++) {
			tmpList.add(retList.get(i));
		}
		
		int addedNum = 1;
		for (int i=0; i<size; i++) {
			int curIdx = tmpList.get(i).index;
			int[] chdIdx = incNodes.containPNode(curIdx, g.graphNodes.size()); //inclusion children index of the curIdx if exist.
			if (chdIdx != null) {
				for (int j=0; j<chdIdx.length; j++) {
					retList.add(i+addedNum, new StartPos(0, g.getSeqOfNode(chdIdx[j]),0));
					addedNum++;
				}
			}
		}
		return retList;
	}
	
	
	/*
	 * Calculate the consensus base from several characters
	 */
	public char getConsensusBase(ArrayList<Character> lst) {
		int numA = 0;
		int numG = 0;
		int numC = 0;
		int numT = 0;
		int numN = 0;
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
				case 'N':
					numN++;
					break;
				case '-':
					numDash++;
					break;
			}
		}
		
		int max = Math.max(Math.max(Math.max(Math.max(Math.max(numA, numG), numC), numT), numN), numDash);
		
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
		} else if (numN == max){
			return 'N';
		} else {
			return '-';
		}
	}
	
	/*
	 * Print the original sequence and multiple consensus into a file which is specified in the property file.
	 */
	public void printConsensus() {
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
	 * Construct a directed Miminum spanning tree.
	 * 
	 *  @param nOfNodes number of nodes
	 *  @param g a directed graph, the second dimension has three elements:
	 *  	index of starting node, index of ending node, weight between them.
	 */
	protected WeightedAdjacencyListGraph constructMinTree(int nOfNodes, int[][] g) {
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
	
	protected void getStartPos(int parentNode, int leftEnd, WeightedAdjacencyListGraph tree, int[][] d) {
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
	
	protected void getStartPosDebug(int parentNode, int leftEnd, WeightedAdjacencyListGraph tree, int[][] d) {
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

		// Get the components of the time
	    long time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to generate MST.");
		assemble.createAlignArray();
		System.out.println("End to generate 6-tuples.");
		System.out.println("The time used to generate 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to process 6-tuples.");
		assemble.processAlignArray();
		System.out.println("End to process 6-tuples.");
		System.out.println("The time used to process 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to reconstruct.");
		assemble.printConsensus();
		System.out.println("End to reconstruct.");
		System.out.println("The time used to reconstruct is " + (new GregorianCalendar().getTimeInMillis()-time1));

		assemble.printSPos();
		assemble.calcInversion();
		
		assemble.incNodes.printAllNodes();
//		ArrayList<StartPos> input = new ArrayList<StartPos>();
//		input.add(new StartPos(0, "CGGGCCTTCGTTTTACGAAAACAGGTGCGCGAAGCCTGCAAATTTGACGGGGATGGGTATGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCGTTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGAAGAGCCCCGGACTCGCCTTGTAACAAGGCGCGGGA", -1));
//		input.add(new StartPos(1, "GAAAACAGGTGCGCGAAGCCTGCAAATTTGACGGGGATGGGTATGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCATTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGAAGAGCCCCGGACTCGCCTTGTAGCAAAGCGCAGGC", -1));
//		input.add(new StartPos(2, "AGGTGCGCGAAGCCTGCAAATTTGACGGGGATGGGTATGCTGTTCAGTGCAGCAGCAGAGGCGGAGCATTGGCAAGGAGCGGCTGCTGGCGCGAAGCAGCAACAGGAGCAGTCAATAAGCATGGCTATGCAGGAATCTAATACGCTTTGAGGTGCAACAACAGCCCTTCAGCGCCGGAGATGCCGCTATGGGGAGCGGAGACGCCCAGCAGCCTGGGGCCGTCATATCTTGTGCCCCCATGGCCCGCTATGGAGCTGGACTACCTTACAGATTACAGCAGCATGCTCTGGTGCTGCCGCGGACGGCAGCCTCGGCGGCTGGGGTCATGGCGGCGAGCAAGGTCGCTGACAGTCGCTGCCAGAACCGCCAGTCATTAAGCATCAGAAATTATCCGCTCGTGGGAACGCAGTGCCCCAGCAGCAGTTATCCGTGGAGCTCGCGCACGCAGGGAGCGCCGGAGGCCAGCCCCAGCACCTACAGCTATAAATTTTCAGCAACAGATTCTCCATTTTGGGCTTTCTTTGGTTTGTTGCTTACCGGCTAGTCGCCGTCTGGAAGGCTGCGCTAATCGGCCTGACACCCTTGGCGCGTTGCCTCTCACTGTGCAGGCAAAGATGGGTCGTCGCCCCGCCCGCTGCTACCGGCAGTCTAAGGGTAAGCCCTACCCGAAGTCTCGCTTCTGCCGTGGTGTGCCCGACCCCAAGATCCGCATTTACGATGCGGGTATGAAGAGGGCCGATGTGGACACCTTCCCCTGCTGCGTCCACCTGGCCAGGTGCGTTAGTGACAGGACAATGGGGTGAAGTGGTGGCAAGCTGGGGGGCGGGTTTCTGCTTAGCACGGGGACTTGCGTGGGGAGCGGACTGGCAGGATAGGGGCAACAGCCGGGGTTGGGGCGGTCTTGGAGGCTGTGATCTCATTGCCATTGGACCGGCGTCGGGCCTCTCAACATCATTCGGAGGTGCCGGAATGAGCGCACGGAGCAGCGCGGAGCGCCACCACTGCGGTGGCGCCATGGCTTGGAGCTGGAGCGGGGCCAACGTGTGCCTGTGGCGGGCAGCGGCAGCCTTGCTGGCAGCAGCGCGGGCAGTAGCCATTGCGCCCGCAAGCATGGACCGGGCGGAGGGCCTGGGGCGCTTGATATGATTGGAGCTGGGCGGCCGCGCTCGCTGCTCGGAACTACCGCGCCTCGTGCCTGGACAGCAATGAACAGCGGAGCGGGCTACTGAACGCAAGGGGGACGGGACGCATTCGCAGTGCCTACCGCCCGGCGTTTGGGACTGGGGTTGAGACGGGAACGGGCTGCATTCAGCTGTGCTGACGTGCTGTTCCTGTTGTTGCTGTGCTGTACAAACAGTTGGGAGAAGGAGAACGTGACCAGTGAGGCGCTGGAGGCTGCCCGTGTGGCGGCTAACAAGTACATGGTGAAGAACGCCGGAAAGGAGGCGTTCCACCTGCGCGTGCGCGTGCACCCCTTCCACGTGCTGCGCATCAACAAGGCAAGCAGGGTCTGGGCAAGCAGGGTGCTTGTCGGTACCTTGAAAAGATAGGCCCTTCGAGCTGGGATGCATTTGCGCACTCTTTCCGCTGTTGCGGACCATCAGCATGGGACAATGTGGTCATGCTCATTACGGTGGTTGTTTTGGCGTGCATCCCAGTTAGATGGTCATGCTGTGCTGACCCGGCTTCATCGCGCCTTTGTGCTCGCAGATGCTTTCGTGCGCAGGCGCTGATCGCCTGCAGACCGGTATGCGTGGTGCTTTCGGCAAGCCCAACGGCGTCTGCGCTCGTGTGCAGATCGGCCAGGTGAGGCACGCGTGCCGGTGGCAGCTGGACACGAGGAATAGGGGATGGGGCGCTGCCACCGCAAAGGGTGCTAAGTACATCCCAAGTTCGGTCTCAATCTGGGCAGGAGGCGTGGGGTGGATCGACATGGACTGTCGCGCATTGGGGGCACTGAGGCCGCACAACACGGAGCAGGAGCATTTGGATTGTGGTTGTGGCGCGCTGCAGCTGAACGGGCCGGCCGTTCGAGCGGTGGCGTTCCGCTACGCAGTAGCGGCGCGGACGGGCACCATGGAAAGTCGCATGGTTTTCATGTTTTACGGCTGACCAGTCGGTTCATGCTCCTTGAATTGCCCAGGTGCTGATGAGCATCCGCTGCCGTGACAACCACGGTGCCGTGGCCGAGGAGGCCCTGCGCCGCGCCAAGTTCAAGTTCCCCGGCCGGCAGAAGATCATCCGGTCCAACAACTGGTGAGTTGACAGCGCTGCTGGCTAGATCTGACACTAAGAAATGGGCAAGTGACTGGCGTCGGGCAGGCGGGTTCCTGAGGGCCGGTGGATCAGCGGGCAGCAACTACCATACTCAGTAGCGAGCGGTGCAGCAGGAGCTCGGCGCGATTTTGTCTCATCCTTATATGAGAACTAGCAGTAGCAGGGCGTTCTAGCTGTATGGACGCGTAGCATGAGCAAATGCCCGCCTTTTCCCCTCGTGGCAACCATGAGCGAACCTCTTGCGCCTGATCGCGTATTGCTGTTTTCCTGTCCTCGCCCTCTTCCACCTGTACAGGGGCTTCACCAACCTGTCCCGCAAGGACTTCAAGATCTTCCGCGAGGAGGGCCGCCTGATCAACGACGGCTCGCACGTCAAGGTCATCACCAACAAGGGCCCCCTGGCCGAGCGCGAGCCCGACCACATCTTCGACTACCCCGCCTTCAAGCACCACACACCCCTGCACAAGGACGAGTAAACAGCTGGCGCTAGCAGCGTCGGTTGCGGAGGGCGGCAGCAGCATCGGGCGGCGCGCCGCAGCGTGCGCGGCGTGAGCTGTGCAGGCGCGTGGCGTGGCGTCGTTGCGGGGCCGCTCCGCCTGCACCCGTTCTGTGCTCGTGTGCGTAGTGTGCGGCTAGGCAAGGAGGCTGTTGGTTGGCGTGAAGTGACAGGCCCGGCTGCGTGTCTTTTGCGGGGGCGGTCTGCTGCCATTGCAGAGGCGCTTAGTGTCGGCCGTGGGCTGGCCCCAGGCCACGACTGAACTGTCGGGCGTGGTAGGGTCGCTGCTGTGTTCCGCGCATGTTCCCCTGGTGGTGCATGTTGGGGCTCATGGCGGCATACCCCACCGCGAGTCCAACGTTGGGTGGGGAAGAGCCCCGGACTCGCCTTGTAGCAAAGCGCAGGC", -1));
//		System.out.println(assemble.reconstructSeq(input));
	}
}
