import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.Properties;

import com.mhhe.clrs2e.*;

public class ESTAssembly {
	int INT_MAX = 2147483647;
	int INT_MIN = -2147483647;
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
	ArrayList <Integer> leftMostNodes;

	/*
	 * store the position of aligned nodes
	 * 1st-dimension: index of nodes in graph;
	 * 2rd-dimension: 	the first is the index of node on the left, 
	 * 					the second is the overlap length with + or -.
	 * 					the third is the distance with + or -.
	 * 					the fourth is the index of node on the right, 
	 * 					the fifth is the overlap length with + or -.
	 * 					the sixth is the distance with + or -.
	 */
	SixTuple[] alignArray;
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
		leftMostNodes = new ArrayList<Integer> ();
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
		alignArray = g.get2CloseNodesFromMST(mstForG);
	}

	/* 
	 * Get the assumed left ends(alignNodes[][0]==-1), check them to find all the real left ends.
	 */
	public void processAlignArray() {
		//get all the nodes which has no left nodes to them
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode == -1) {
				leftMostNodes.add(Integer.valueOf(i));	//store index of the node
			}
		}

		//Re-calculating six-tuples for those assumed left ends and put the new six-tuple into alignArray.
		int numOfLeftMostNodes = leftMostNodes.size();
		System.out.println("The number of original left ends is " + numOfLeftMostNodes);
		if (numOfLeftMostNodes > 1) {
			for (int i=0; i<leftMostNodes.size(); i++) {
				int cNode = leftMostNodes.get(i).intValue();
				SixTuple lNode = g.get2CloseNodesFromGrand(mstForG, cNode, alignArray[cNode]);
				alignArray[cNode].leftNode = lNode.leftNode;
				alignArray[cNode].lOvlLen = lNode.lOvlLen;	//overlap length
				alignArray[cNode].lDis = lNode.lDis;	//distance
				
				if (alignArray[cNode].rDis > lNode.rDis) { //get a smaller distance
					alignArray[cNode].rightNode = lNode.rightNode;
					alignArray[cNode].rOvlLen = lNode.rOvlLen;
					alignArray[cNode].rDis = lNode.rDis;
				}
			}
		} 
		
		
		leftMostNodes.clear();
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode == -1) {
				leftMostNodes.add(Integer.valueOf(i));	//store index of the node
			}
		}
		System.out.println("There are " + leftMostNodes.size() + " left-most nodes after running 3 levels.");
		
		/* Recalculate 6-tuples for all the current left nodes in order to remove all the false left ends.
		 * Specifically, for those assumed left ends,start to calculate from fourth level until meeting one node which 
		 * makes six-tuple[0] != -1 or until the level we specified in the property file, then return the six-tuple.
		 * If we fail to find any node, we consider it a real left end.
		 */
		for (int i=0; i<leftMostNodes.size(); i++) {
			int tEnd = leftMostNodes.get(i).intValue(); //index of the node
			SixTuple tmpTuple = g.checkLeftEndFromMST(mstForG, tEnd, alignArray[tEnd]);
			if (tmpTuple != null) {
				alignArray[tEnd].leftNode = tmpTuple.leftNode;
				alignArray[tEnd].lOvlLen = tmpTuple.lOvlLen;	//overlap length
				alignArray[tEnd].lDis = tmpTuple.lDis;	//distance
				alignArray[tEnd].rightNode = tmpTuple.rightNode;
				alignArray[tEnd].rOvlLen = tmpTuple.rOvlLen;
				alignArray[tEnd].rDis = tmpTuple.rDis;
			}
		}
		
		leftMostNodes.clear();
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode == -1) {
				leftMostNodes.add(Integer.valueOf(i));	//store index of the node
			}
		}
		System.out.println("\nThere are " + leftMostNodes.size() + " left-most nodes after checking left ends.");
		
		
		
		
		/* construct a temporary directed graph from alignArray, 
		 *  	the second dimension has two elements:
		 *  			index of starting node,
		 *  			index of ending node, 
		 *  			weight between them (positive value, weight is abs(their distance)).
		 *  	if there is no edge, weight=INT_MAX.
		 */
		int tLen = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode != -1) {
				tLen++;
			}
			if (alignArray[i].rightNode != -1) {
				tLen++;
			}
		}
		
		int[][] tmpDGraph = new int[tLen][4];
		int tmpIndex = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode != -1) {
				tmpDGraph[tmpIndex][0] = alignArray[i].leftNode;
				tmpDGraph[tmpIndex][1] = i;
				tmpDGraph[tmpIndex][2] = Math.abs(alignArray[i].lDis);	//distance
				tmpDGraph[tmpIndex][3] = Math.abs(alignArray[i].lOvlLen);	//overlap length
				tmpIndex++;
			}
			
			if (alignArray[i].rightNode != -1) {
				tmpDGraph[tmpIndex][0] = i;
				tmpDGraph[tmpIndex][1] = alignArray[i].rightNode;
				tmpDGraph[tmpIndex][2] = alignArray[i].rDis;	//distance
				tmpDGraph[tmpIndex][3] = alignArray[i].rOvlLen;	//overlap length
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
		
		ArrayList <Integer> tmpLeftNodes = new ArrayList<Integer> ();
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode == -1) {
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
			for (int j=0; j<tmpDGraph.length; j++) {
				if (tmpDGraph[j][1] == tEnd) {	// false left end
					f = 1;
					break;
				}
			}
			if (f == 0) {
				leftMostNodes.add(Integer.valueOf(tEnd));
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
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode != -1) {
				len++;
			}
			if (alignArray[i].rightNode != -1) {
				len++;
			}
		}
		/*
		 * generate dGraph.
		 */
		int[][] dGraph = new int[len][4];
		int indexOfDGraph = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i].leftNode != -1) {
				dGraph[indexOfDGraph][0] = alignArray[i].leftNode;
				dGraph[indexOfDGraph][1] = i;
				dGraph[indexOfDGraph][2] = Math.abs(alignArray[i].lDis);	//distance
				dGraph[indexOfDGraph][3] = Math.abs(alignArray[i].lOvlLen);	//overlap length
				indexOfDGraph++;
			}
			
			if (alignArray[i].rightNode != -1) {
				dGraph[indexOfDGraph][0] = i;
				dGraph[indexOfDGraph][1] = alignArray[i].rightNode;
				dGraph[indexOfDGraph][2] = alignArray[i].rDis;	//distance
				dGraph[indexOfDGraph][3] = alignArray[i].rOvlLen;	//overlap length
				indexOfDGraph++;
			}
		}

		//print dGraph
		//System.out.println("dGraph:");
		//printDgraph(dGraph);
		
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
		ArrayList<String> allSingletons= new ArrayList<String> ();	//store all the singletons
		ArrayList<String> lastEsts = new ArrayList<String> ();
		ArrayList<String> firstEsts = new ArrayList<String> ();
		
		sPosDebug = new int[alignArray.length];	
		for (int i=0; i<leftMostNodes.size(); i++) { //start for
			sPos = new int[alignArray.length];	//store starting positions of all the nodes

			int leftEnd = leftMostNodes.get(i).intValue();
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
			primMST = constructMinTree(alignArray.length, dGraph);

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
					tmpArray.add(new StartPos(sPos[j], g.getSeqOfNode(j)));
				}
			}
			System.out.println(tmpArray.size() + " nodes are used to reconstruct the sequence.\n");
			if (tmpArray.size() == 1) { //singleton
				allSingletons.add(g.getCommentOfNode(leftEnd) + "\n" + g.getSeqOfNode(leftEnd));
				continue;
			} 
			
			lastEsts.add(tmpArray.get(tmpArray.size()-1).seq);
			firstEsts.add(g.getSeqOfNode(leftEnd));

			String tStr = reconstructSeq(tmpArray, 0);
			allConsensus.add(tStr);
			printStr = printStr + tStr + "\n";

			
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
			ArrayList<String> s = processMoreConsensus(allConsensus, firstEsts, lastEsts);
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
	 * not overlapped(or they correspond to different part of the gene). Recall that we treat inclusion as infinite distance, 
	 * so some left ends may include each other, and consensus from them overlap with each other.
	 * This method is designed to remove all the dependent consensus and extract all the independent ones, which means, we intend 
	 * to find all the consensus which represent different part of the gene.. 
	 * 
	 * @param s an arraylist which includes all the generated consensus from the calling method;
	 * 		 	firstEsts an arraylist which includes all the left end sequence corresponding to all the consensus;
	 * 			lastEsts an arraylist which includes all last est sequence corresponding to all the consensus.
	 * @return the combined consensus.
	 */
	 private ArrayList<String> processMoreConsensus(ArrayList<String> s, ArrayList<String> firstEsts, ArrayList<String> lastEsts) {
		//String retStr = "";
		ArrayList<String> allOutputContigs= new ArrayList<String> ();	//store all the generated sequences
		
		int sizeOfs = s.size();
		Consensus[] resultArray = new Consensus[sizeOfs]; //store the starting positions of ests
		for (int i=0; i<sizeOfs; i++) {
			resultArray[i] = new Consensus(firstEsts.get(i).length(), s.get(i), firstEsts.get(i), lastEsts.get(i));
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
				boolean b = g.d2.checkInclusion(allConsensus.get(i).firstEst, s1); //if resultArray[i].firstEst is included in s1
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
		 String seq;
		 String firstEst;
		 String lastEst;
		 public Consensus(int p, String s1, String s2, String s3) {
			 lenOfFirstEst = p;
			 seq = s1;
			 firstEst = s2;
			 lastEst = s3;
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
	 * These input consensus originate from more than one left ends which include each other. Although the 
	 * consensus are different at their beginning, they often end with the same characters.
	 * For example, they will look like:
	 * AGGCTCTCCCCAAGTCCACTAGTTCAGACGGGACAATATAACGGACTGCATGGCAGCGCATGTCGAGCTCCACGCGCATCTACACTCACCTCGCATGGACTGCACAAT
	 *                       TTCAGACGGGACAATATAACGGACTGCATGGCAGCGCATGTCGAGCTCCACGCGCATCTACACTCACCTCGCATGGACTGCACAAT
	 *               TCCACTAGTTCAGACGGGACAATATAACGGACTGCATGGCAGCGCATGTCGAGCTCCACGCGCATCTACACTCACCTCGCATGGACTGCACAAT
	 *           CAAGTCCACTAGTTCAGACGGGACAATATAACGGACTGCATGGCAGCGCATGTCGAGCTCCACGCGCATCTACACTCACCTCGCATGGACTGCACAAT
	 * 
	 * So in this method, we reverse all the consensus, and get one consensus from them. 
	 * 
	 * In some situation, the consensus may not ends with the same character even if their left ends have inclusion.
	 * In this method, we judge if two consensus ends with the same Est. If it does, we process them. If not, 
	 * we will judge if the consensus are part of the longer one. If it does, we will return the longer one. If not,
	 * we will add it to the return string and return it.
	 */
	private String processMoreConsensusWithInclusion(ArrayList<Consensus> includeStrs) {
		//ArrayList<String> s, ArrayList<String> lastEsts
		int maxLen = 0; //the maximal length of all the strings.
		int indexOfMax = -1;
		for (int i=0; i<includeStrs.size(); i++) {
			int tLen = includeStrs.get(i).seq.length();
			if (tLen > maxLen) {
				maxLen = tLen;
				indexOfMax = i;
			}
		}
		String maxLastStr = includeStrs.get(indexOfMax).lastEst;
		ArrayList<String> tmpStr1 = new ArrayList<String>(); //will be processed
		ArrayList<String> tmpStr2 = new ArrayList<String>(); //won't be processed
		tmpStr1.add(includeStrs.get(indexOfMax).seq); //tmpStr1 has at least one element
		for (int i=0; i<includeStrs.size(); i++) {
			if (i != indexOfMax) {
				if (includeStrs.get(i).lastEst.compareToIgnoreCase(maxLastStr) == 0) {
					tmpStr1.add(includeStrs.get(i).seq);
				} else {
					tmpStr2.add(includeStrs.get(i).seq);
				}
			}
		}

		String retStr = "";
		int size = tmpStr1.size();
		if (size > 2) {
			String[] strs = new String[size];
			for (int i=0; i<size; i++) {
				strs[i] = tmpStr1.get(i);
			}

			// Calculate consensus base for each position, and put them into an char array
			char[] consensus = new char[maxLen];
			for (int i=0; i<maxLen; i++) {
				ArrayList<Character> tmpArraylist = new ArrayList<Character>();
				for (int j=0; j<size; j++) {
					int index = strs[j].length() - (i+1);
					if (index >= 0) {
						tmpArraylist.add(strs[j].charAt(index));
					}
				}
				consensus[maxLen-i-1] = getConsensusBase(tmpArraylist);
			}
			retStr = retStr + String.valueOf(consensus);
		} else if ((size == 1)|| (size == 2)) {
			retStr = retStr + tmpStr1.get(0);
		}
		
		String tmpRetStr = retStr;
		for (int i=0; i<tmpStr2.size(); i++) {
			//retStr = retStr + "\n" + tmpStr2.get(i);
			String tStr = tmpStr2.get(i);
			double dis = g.d2.getLocalSimlarityScore(tmpRetStr, tStr);
			if ((dis/tStr.length()) < 0.95) { //if tStr is not included in retStr, attach it to retStr.
				//System.out.println(tmpRetStr);
				//System.out.println(tStr);
				retStr = retStr + "\n" + tStr;
			}
		}
		
		return retStr;
	}
	
	/*
	 * reconstruct a sequence which starts from a left end.
	 */
	private String reconstructSeq(ArrayList<StartPos> a, int breakPoint) {
		int sizeOfa = a.size();
		if (sizeOfa == 0) {
			return null;
		} else if (sizeOfa == 1) {
			return a.get(0).seq;
		}
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);


		//local align every two adjacent nodes, e.g, node 1 and 2, node 2 and 3..., node i-1 and node i.
		//put them into two arraylists, the local alignment sequences of node i and i+1 are put into arraylist 1 and 2 respectively.
		ArrayList<String> arrList1 = new ArrayList<String> ();
		ArrayList<String> arrList2 = new ArrayList<String> ();
		for (int i=0; i<resultArray.length-1; i++) {
			String[] strs = g.d2.getLocalAlignment(resultArray[i].seq, resultArray[i+1].seq);
			arrList1.add(strs[0]);
			arrList2.add(strs[1]);
		}
		
		int size = arrList1.size();
		//store the starting positions of all the elements in arrList1
		int[] sposArr1 = new int[size];
		//store the starting positions of all the elements in arrList2
		int[] sposArr2 = new int[size];
		//calculate the starting positions and put them into sposArr1 and sposArr2.
		sposArr1[0] = 0;
		sposArr2[0] = 0;
		int bPoint = size-1; //record the point when the the program breaks the loop. default is when it goes through all the loops.
		int flag = 0; //identify if the we need execute reconstruction again for the program breaks the loop.
		ArrayList<StartPos> b = new ArrayList<StartPos> (); //the input parameter for this function if flag=1.
		for (int i=breakPoint; i<size-1; i++) {
			String s1 = arrList2.get(i);
			String s2 = arrList1.get(i+1);
			String s1Tmp = s1.replace("-", "");
			String s2Tmp = s2.replace("-", "");
			
			// if s1 and s2 do not overlap, record i to bPoint, and break this loop.
			// get the assembled sequence from all the i<bPoint, add "\n" to the end.
			// and then restart the reconstruction from i+1.
			String tmpStr = resultArray[i+1].seq;
			int tmpPos1 = tmpStr.indexOf(s1Tmp);
			int tmpPos2 = tmpStr.indexOf(s2Tmp);
			if (tmpPos1+s1.length() <= tmpPos2+30) {
				bPoint = i;
				flag = 1;
				//put the left data to variable b as the input parameter
				for (int j=i+1; j<sizeOfa; j++) {
					b.add(resultArray[j]);
				}
				break;
			}
			// in general case, s1 should be on the left of s2.
			// but rarely, s1 may be on the right of s2. If it is the case, we would exchange s1 and s2, 
			//     calculate the offset, and then give the offset a negative value.
			int posDirection = 1;
			tmpPos1 = tmpStr.indexOf(s1Tmp); //tmpStr may have been changed, so calculate it again.
			tmpPos2 = tmpStr.indexOf(s2Tmp);
			if (tmpPos1 > tmpPos2) {
				posDirection = -1; //set direction to be -1.
				String tmp = s1;
				s1 = s2;
				s2 = tmp;
			}

			if (posDirection != -1) {
				// get a correct s1 if there is "-" in s1 in order to get a correct position for s2.
				boolean b1 = (s1.indexOf("-") != -1);
				if (b1) {
					/* 
					 * align the six acquired sequences in sposArr1 and sposArr2 before s1(including s1), and put  
					 * the correct sequence of s1 into arrList2. Then recalculate the offset of subS1 (we would get a new s1).
					 */
					if (i == 1) { //get the previous four sequences
						String[] alignedStrs = new String[4];
						alignedStrs[0] = arrList1.get(0);
						alignedStrs[1] = arrList2.get(0);
						alignedStrs[2] = arrList1.get(1);
						alignedStrs[3] = arrList2.get(1);
						int[] startPos = new int[4];
						startPos[0] = sposArr1[0];
						startPos[1] = sposArr2[0];
						startPos[2] = sposArr1[1];
						startPos[3] = sposArr2[1];
						s1 = getCorrectS1(alignedStrs, startPos);
					} else if (i >= 2) { //get the previous six sequences
						String[] alignedStrs = new String[6];
						alignedStrs[0] = arrList1.get(i-2);
						alignedStrs[1] = arrList2.get(i-2);
						alignedStrs[2] = arrList1.get(i-1);
						alignedStrs[3] = arrList2.get(i-1);
						alignedStrs[4] = arrList1.get(i);
						alignedStrs[5] = arrList2.get(i);
						int[] startPos = new int[6];
						startPos[0] = sposArr1[i-2];
						startPos[1] = sposArr2[i-2];
						startPos[2] = sposArr1[i-1];
						startPos[3] = sposArr2[i-1];
						startPos[4] = sposArr1[i];
						startPos[5] = sposArr2[i];
						s1 = getCorrectS1(alignedStrs, startPos);
					}
					arrList2.set(i, s1);
				}
			}
			
			/*calculate starting position for s2
			 * 1. Find the position of subS1 in s1, we name it as offset.
			 * 2. Get starting position of s2.
			 * 	Note that subS1 and subS2 may not be the overlap of s1 and s2, they can be shorter than the overlap.
			 * 	For example, 
			 * 		GAGACAAGACAAGGCTCTCCCCAAGTCCACTAGTTCAGACGGGACA
			 *                            CTGTCCACTAGTTCAGACGGGACAATATAACGGACTGCATGGCAGC
			 *  subS1 will be: GTCCACTAGTTCAGACGGGACA, not AAGTCCACTAGTTCAGACGGGACA.
			 *  
			 *  So we will find the position of subS2 in s2, and replace s2 with the new one, that is: 
			 *  GTCCACTAGTTCAGACGGGACAATATAACGGACTGCATGGCAGC. Correspondingly, we will replace s3 with the substring
			 *  of the old s3 which corresponds to the new s2.
			 *  
			 *  3. If there is "-" in the new s2, we will correct it. Because "-" in s2 means there may be some error
			 *  in s3, which will affect the starting position of the following sequence.
			 */
			//get the position of subS1 in s1.
			s1Tmp = s1.replace("-", "");
			s2Tmp = s2.replace("-", "");
			String[] tStrs = g.d2.getLocalAlignment(s1Tmp, s2Tmp);
			String subS1 = tStrs[0].replace("-", "");
			int offset = s1.indexOf(subS1);
			if (offset == -1) { 
				for (int idx=0; idx<s1.length(); idx++) {
					if (s1.charAt(idx) != '-') {
						int tmpIdx = s1.substring(idx).replace("-", "").indexOf(subS1);
						if (tmpIdx == 0) {
							offset = idx;
							break;
						}
					}
				}
			} 
			//get the starting position for s2
			String subS2 = tStrs[1].replace("-", "");
			String s3 = arrList2.get(i+1);
			int tOffset = s2.indexOf(subS2);
			if (tOffset == -1) { 
				for (int idx=0; idx<s2.length(); idx++) {
					if (s2.charAt(idx) != '-') {
						int tmpIdx = s2.substring(idx).replace("-", "").indexOf(subS2);
						if (tmpIdx == 0) {
							if (posDirection == -1) {
								offset = offset - idx;
							} else {
								tOffset = idx;
							}
							break;
						}
					}
				}
			} else {
				if (posDirection == -1) {
					offset = offset - tOffset;
				} 
			}
			if (posDirection != -1) {
				//set s2 and s3 to be correct ones
				s2 = s2.substring(tOffset);
				s3 = s3.substring(tOffset);
				arrList1.set(i+1, s2);
				arrList2.set(i+1, s3);
			}
			
			//set starting positions
			offset = offset * posDirection;
			sposArr1[i+1] = sposArr2[i]+offset;
			sposArr2[i+1] = sposArr1[i+1];

			// get a correct s2 if there is "-" in s2.
			s2 = arrList1.get(i+1); // re-get s2 to avoid the effect of exchange and changes on it. 
			boolean b2 = (s2.indexOf("-") != -1);
			if (b2) {
				/* 
				 * align the seven acquired sequences in sposArr1 and sposArr2 before s2(including s1), and put  
				 * the correct sequence of s2 into arrList1.
				 */
				if (i == 1) { //get the previous five sequences
					String[] alignedStrs = new String[5];
					alignedStrs[0] = arrList1.get(0);
					alignedStrs[1] = arrList2.get(0);
					alignedStrs[2] = arrList1.get(1);
					alignedStrs[3] = arrList2.get(1);
					alignedStrs[4] = arrList1.get(2);
					int[] startPos = new int[5];
					startPos[0] = sposArr1[0];
					startPos[1] = sposArr2[0];
					startPos[2] = sposArr1[1];
					startPos[3] = sposArr2[1];
					startPos[4] = sposArr1[2];
					s2 = getCorrectS1(alignedStrs, startPos);
				} else if (i >= 2) { //get the previous seven sequences
					String[] alignedStrs = new String[7];
					alignedStrs[0] = arrList1.get(i-2);
					alignedStrs[1] = arrList2.get(i-2);
					alignedStrs[2] = arrList1.get(i-1);
					alignedStrs[3] = arrList2.get(i-1);
					alignedStrs[4] = arrList1.get(i);
					alignedStrs[5] = arrList2.get(i);
					alignedStrs[6] = arrList1.get(i+1);
					int[] startPos = new int[7];
					startPos[0] = sposArr1[i-2];
					startPos[1] = sposArr2[i-2];
					startPos[2] = sposArr1[i-1];
					startPos[3] = sposArr2[i-1];
					startPos[4] = sposArr1[i];
					startPos[5] = sposArr2[i];
					startPos[6] = sposArr1[i+1];
					s2 = getCorrectS1(alignedStrs, startPos);
				}
				arrList1.set(i+1, s2);
				arrList2.set(i+1, s2);
			}
		} //end for
		
		
		// Create an array which stores all the bases with the same position for all the ESTs.
		int lenOfArray = resultArray[bPoint+1].pos
						- resultArray[0].pos
						+ resultArray[bPoint+1].seq.length(); 
		ArrayList<Character> [] tmpArraylists = new ArrayList [lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			tmpArraylists[i] = new ArrayList<Character>();
		}
		// Put all the bases into the array
		for (int i=breakPoint; i<=bPoint; i++) {
			int tmpSPos1 = sposArr1[i];
			String tmpStr1 = arrList1.get(i);
			int tmpSPos2 = sposArr2[i];
			String tmpStr2 = arrList2.get(i);
			
			//int f=0;
			for (int j=0; j<tmpStr1.length(); j++) {
				int p = tmpSPos1+j;
				if ((p < lenOfArray) && (p >= 0)) {
					tmpArraylists[p].add(tmpStr1.charAt(j));
				}
				//if (p==14575) f=1;
			}
/*			if (f==1) {
				System.out.println("tmpSPos1="+tmpSPos1);
				System.out.println(tmpStr1.charAt(14575-tmpSPos1));
				System.out.println(tmpStr1);
			}
			f=0;
*/			for (int j=0; j<tmpStr2.length(); j++) {
				int p = tmpSPos2+j;
				if ((p < lenOfArray) && (p >= 0)){
					tmpArraylists[p].add(tmpStr2.charAt(j));
				}
				//if (p==14575) f=1;
			}
/*			if (f==1){
				System.out.println("tmpSPos2="+tmpSPos2);
				System.out.println(tmpStr2.charAt(14575-tmpSPos2));
				System.out.println(tmpStr2+"\n\n");
			}
*/		}
		
		/* Calculate consensus base for each position, and put them into an char array.
		 * Do not include those bases with less than 3/5 prepared-bases to remove errors near both ends.
		 */
		char[] consensus = new char[lenOfArray];
		int start = 0;
		int end = lenOfArray-1;
		for (int i=0; i<=end; i++) {
			if (tmpArraylists[i].size() >= 1) {
				start = i;
				break;
			}
		}
		for (int i=end; i>=start; i--) {
			if (tmpArraylists[i].size() >= 1) {
				end = i;
				break;
			}
		}
		
		//System.out.println("start=" + start + "; end=" + end);
		for (int i=start; i<=end; i++) {
			consensus[i-start] = getConsensusBase(tmpArraylists[i]);
		}
		
		String retStr = String.valueOf(consensus).substring(0, end-start+1);
		retStr = retStr.replace("-", "");
		retStr = retStr.replace(" ", "");
		
		if (flag == 1) {//execute the function again
			retStr = retStr + "\n" + reconstructSeq(b, bPoint+1);
		} 
		return retStr;
	}
	
	/*
	 * Get the correct sequence for the last element in the input string array.
	 * @param strs an String array which includes all the strings.
	 * @param pos an int array which records the starting positions of all the elements in strs.
	 * @return the consensus of the last element in strs.
	 * 
	 * Get the consensus according to the two input parameters, then extract the subsequence which corresponds
	 * to the last element in strs from the consensus.
	 */
	private String getCorrectS1(String[] strs, int[] pos) {
		int minPos = INT_MAX;
		int maxPos = 0;
		int len = strs.length;
		for (int i=0; i<len; i++) {
			if (pos[i] < minPos) {
				minPos = pos[i];
			}
			if (pos[i] > maxPos) {
				maxPos = pos[i];
			}
		}
		
		// Create an array which stores all the bases with the same position.
		int lenOfArray = pos[len-1] - minPos + strs[len-1].length();
		ArrayList<Character> [] tmpArraylists = new ArrayList [lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			tmpArraylists[i] = new ArrayList<Character>();
		}

		// Put all the bases into the array
		for (int i=0; i<len; i++) {
			int tmpSPos = pos[i] - minPos;
			String tmpStr = strs[i];
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
		retStr = retStr.substring(pos[len-1]-minPos);

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
		
	}
}
