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
	protected final int INT_MAX = Integer.MAX_VALUE;

	Graph g;
	Properties props;
	ArrayList<SixTuple> alignArray;
	int[] sPos;	//starting positions of all the nodes
				//the index in the array is the index of the node, the value is its starting position. 
	int[] sPosDebug;	//starting positions of all the nodes
						//it's used for debugging. All the left ends will be assigned to their actual value in order to calculate inversions later. 
	Alignment alignment;
	InclusionNodes incNodes;
	//ArrayList<SixTuple> leftMostNodes;
	String consensusFileName;
	String singletonFileName;
	String numOfUsedESTsFileName;
	String printStr;
	ArrayList<String> numOfNodes;	//store number of used nodes, corresponds to each element in allConsensus.
	ArrayList<String> firstEsts;	//includes all the left end sequence corresponding to each element in allConsensus.
	int[] usedNodes;	//index of the array is the index of the node, 1-used, 0-not used. It is used to identify singletons.
	int checkType; //0-not check AS, 1-typeI, 2-typeII.
	
	public Reconstruction(Properties p, Graph graph, ArrayList<SixTuple> align, InclusionNodes inc, int type) {
		g = graph;
		alignArray = align;
		sPos = null;
		sPosDebug = null;
		incNodes = inc;
		props = p;
		checkType = type;
		consensusFileName = props.getProperty("ConsensusFile");
		singletonFileName = props.getProperty("SingletonFile");
		numOfUsedESTsFileName = props.getProperty("NumOfUsedESTs");
		alignment = new Alignment(props);
		
		printStr = "";
		numOfNodes = new ArrayList<String> ();	
		firstEsts = new ArrayList<String> ();
		
		usedNodes = new int[g.graphNodes.size()];
		
		int[] chdNodes = incNodes.getAllChdNodes(); //all the children nodes in inclusion list should not be considered as singletons.
		for (int i=0; i<chdNodes.length; i++) {
			usedNodes[chdNodes[i]] = 1;
		}
	}
	
	public void getConsensus(ArrayList<SixTuple> leftEnds) {
		// Get the components of the time
	    long time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to reconstruct.");
		printConsensus(leftEnds);
		System.out.println("End to reconstruct.");
		System.out.println("The time used to reconstruct is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		//Debugger.printSPos(g, sPosDebug);
		//Debugger.calcInversion(g, sPosDebug);
	}
	
	/*
	 * Print the original sequence and multiple consensus into a file which is specified in the property file.
	 */
	private void printConsensus(ArrayList<SixTuple> leftEnds) {
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
			ArrayList<String> consensus = this.reconstruct(leftEnds);
			int index = 1;
			for (int i=0; i<consensus.size(); i++) {
				out1.write(">contig " + (index++) + "\n");
				out1.write(consensus.get(i));
				out1.write("\n");
			}
			out1.flush();
			out1.close();
			
			/*
			 * print singletons and numOfUsedESTs
			 */
			File outFile2 = new File(singletonFileName);
			bExists = outFile2.exists();
			if (bExists) {
				outFile2.delete();
			}
			BufferedWriter out2 = new BufferedWriter(new FileWriter(outFile2, true));
			int num = 0;
			for (int i=0; i<usedNodes.length; i++) {
				if (usedNodes[i] == 0) {	//singleton
					out2.write(g.getCommentOfNode(i) + "\n");
					out2.write(g.getSeqOfNode(i) + "\n");
				} else {
					num++;
				}
			}
			out2.flush();
			out2.close();
			
			File outFile3 = new File(numOfUsedESTsFileName);
			bExists = outFile3.exists();
			if (bExists) {
				outFile3.delete();
			}
			BufferedWriter out3 = new BufferedWriter(new FileWriter(outFile3, true));
			out3.write(">number of used ESTs by EAST\n");
			out3.write(Integer.toString(num));
			out3.flush();
			out3.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}
	
	/*
	 * reconstruct the sequence
	 * 
	 * @return an arraylist: The assembled sequences.
	 */
	private ArrayList<String> reconstruct(ArrayList<SixTuple> leftEnds) {
		sPosDebug = new int[g.graphNodes.size()];	
		
		ArrayList<String> ret = processLeftEnds(leftEnds);

		//for debug
		printStr = printStr + ret.size() + " consensus from above " + leftEnds.size() + " left ends:\n";
		for (int p=0; p<ret.size(); p++) {
			printStr = printStr + ret.get(p) + "\n";
		}
		
		//print debug information about the generated consensus.
		System.out.println("*********************consensus:********************");
		System.out.println(printStr);
		
		return ret;
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
	 *  If flag = 0, return number of used nodes for this left end.
	 *  If flag = 1, return consensus from this left end. 
	 *  	1+str: there is typeI AS;
	 *  	2+str: there is typeII AS;
	 *  	0+str: no AS.
	 *  
	 * 1. Generate dGraph. 
	 * 2. For each left-end node, starting from it to calculate positions for each node. 
	 * In order to get starting positions, it constructs a MST. The weight of the Minimum Spanning tree is 
	 * the overlap distance instead of overlap length. 
	 * Then reconstruct the sequence from the set of ESTs.
	 * 
	 *  	Because the Prim algorithm starts from index 0 to generate MST, we have to
	 *  		put left-end node to index 0 in order to get the MST we want. If Prim does 
	 *  		not start from the left-end node, the directed tree will be unconnected.
	 *  
	 */
	private String getInfoOfLeftEnd(int leftEnd, int[][] dGraph, int flag) {
		String ret = null;
		sPos = new int[g.graphNodes.size()];	//store starting positions of all the nodes

		
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
			//sPosDebug[leftEnd] = sPosDebug[0];
		sPos[0] = 0; //starting position of the left end is assigned to be 0.
			//sPosDebug[0] = Integer.parseInt(g.getNameOfNode(leftEnd));
		//get starting positions for the nodes in primMST
		getStartPos(0, leftEnd, primMST, dGraph);
			//getStartPosDebug(0, leftEnd, primMST, dGraph);
		//exchange sPos[0] and sPos[leftEnd] to recover index 0 in sPos
		int tmp = sPos[0];
			int tmp1 = sPosDebug[0];
		sPos[0] = sPos[leftEnd];
			//sPosDebug[0] = sPosDebug[leftEnd];
		sPos[leftEnd] = tmp;
			//sPosDebug[leftEnd] = tmp1;

		
		//sort the array sPos (in ascending order)
		ArrayList<StartPos> tmpArray = new ArrayList<StartPos> ();
		for (int j=0; j<sPos.length; j++) {
			if ((j == leftEnd) || (sPos[j] != 0)) {
				tmpArray.add(new StartPos(sPos[j], j));
			}
		}

		
		if (flag == 0) { //get number of used nodes
			int num = getNumUsedNodes(tmpArray);
			ret = String.valueOf(num);
		} else if (flag == 1) { //get consensus
			printStr = printStr + Debugger.printLeftEndInfo(leftEnd, g); //get information of this left end which is used to start reconstruction.
			
			if (checkType == 1) {
				int[] asPos = checkAS1(tmpArray);
				if (asPos[0] != -1) {
					printStr = printStr + "alternative splicing type I exists!\n";
					String asVariants = reconstructAS(tmpArray, asPos);
					printStr = printStr + asVariants + "\n\n";
					ret = "1" + asVariants;
				} else {
					String[] tStr = reconstructSeq(tmpArray);
					if (tStr != null) {
						//this is an estimated value, not an exact one because some nodes may repeat (both from nodes2 in inclusion list and also exist in the ordinary list)
						//the value in numOfUsedEsts file is an exact value.
						printStr = printStr + String.valueOf(tStr[1]) + " nodes are used to reconstruct the sequence.(estimated. real number<= the number)\n";
						printStr = printStr + tStr[0] + "\n\n";
						ret = "0" + tStr[0]; 
					}
				}
			} else if (checkType == 2) {
				int[] asPos2 = checkAS2(tmpArray);
				if (asPos2[0] != -1){
					printStr = printStr + "alternative splicing type II exists!\n";
					String asVariants = reconstructAS2(tmpArray, asPos2); //this will be wrote into "consensus.out"
					printStr = printStr + asVariants + "\n\n";
					ret = "2" + asVariants;
				} else {
					String[] tStr = reconstructSeq(tmpArray);
					if (tStr != null) {
						//this is an estimated value, not an exact one because some nodes may repeat (both from nodes2 in inclusion list and also exist in the ordinary list)
						//the value in numOfUsedEsts file is an exact value.
						printStr = printStr + String.valueOf(tStr[1]) + " nodes are used to reconstruct the sequence.(estimated. real number<= the number)\n";
						printStr = printStr + tStr[0] + "\n\n";
						ret = "0" + tStr[0]; 
					}
				}
			} else {
				String[] tStr = reconstructSeq(tmpArray);
				if (tStr != null) {
					//this is an estimated value, not an exact one because some nodes may repeat (both from nodes2 in inclusion list and also exist in the ordinary list)
					//the value in numOfUsedEsts file is an exact value.
					printStr = printStr + String.valueOf(tStr[1]) + " nodes are used to reconstruct the sequence.(estimated. real number<= the number)\n";
					printStr = printStr + tStr[0] + "\n\n";
					ret = "0" + tStr[0]; 
				}
			}
		}
		
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
		return ret;
	}
	
	// write a string to a file, no comment line.
	private void writeToFile(String str, String fileName) {
		try{ 
			File outFile = new File(fileName);
			boolean bExists = outFile.exists();
			if (bExists) {
				outFile.delete();
			}
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
			out.write(str);
			out.flush();
			out.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 

	}
	
	 /* 
	 * This method is designed to group left ends into independent sets. All the dependent left ends (include one by one) are put 
	 * into one group. 
	 * In each group, we find the left end with the longest length and the one which will use the most number of ESTs to reconstruct,
	 * and combine them together to form a consensus by calling the function "processLeftEndsWithInclusion". 

	 * @return all the consensus sequences.
	 */
	 private ArrayList<String> processLeftEnds(ArrayList<SixTuple> leftEnds) {
		ArrayList<String> allOutputContigs= new ArrayList<String> ();	//store all the generated sequences
		
		int[][]dGraph = genDGraph();
		int sizeOfs = leftEnds.size();
		LeftEnd[] resultArray = new LeftEnd[sizeOfs]; //store the starting positions of ests
		for (int i=0; i<sizeOfs; i++) { //start for
			int idx = leftEnds.get(i).curNode;
		//System.out.println(Debugger.printLeftEndInfo(idx, g));  //print the information of all the assumed left ends.
			String num = getInfoOfLeftEnd(idx, dGraph, 0);
			resultArray[i] = new LeftEnd(idx, Integer.parseInt(num), g.getSeqOfNode(idx));
		} //end for

		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		ArrayList<LeftEnd> allLeftEnds = new ArrayList<LeftEnd> ();
		for (int i=sizeOfs-1; i>=0; i--) {
			allLeftEnds.add(resultArray[i]);
		}
		
		while (true) {
			String s1 = allLeftEnds.get(0).seq;
			ArrayList<LeftEnd> includedEnds = new ArrayList<LeftEnd>();
			includedEnds.add(allLeftEnds.get(0));
			ArrayList<LeftEnd> excludedEnds = new ArrayList<LeftEnd>();
			for (int i=1; i<allLeftEnds.size(); i++) {
				boolean b = g.ovl.checkInclusion(allLeftEnds.get(i).seq, s1); //if resultArray[i].firstEst is included in s1
				if (b) {
					includedEnds.add(allLeftEnds.get(i));
				} else {
					excludedEnds.add(allLeftEnds.get(i));
				}
			}

			if ((checkType == 2) && (includedEnds.size() != sizeOfs)) {
				checkType = 0;
			}
			String contig = processLeftEndsWithInclusion(includedEnds);
			if ((contig != null) && (contig.trim().compareTo("") != 0)) {
				allOutputContigs.add(contig);
			}
			if (excludedEnds.size() == 0) {
				break;
			} else {
				allLeftEnds = excludedEnds;
			}
		}
		return allOutputContigs;
	 }
	 


	/*
	 * This method is called by "processLeftEnds".
	 * This method is used to process those left ends that include each other.
	 * 
	 * We combine the one with the longest first EST and the one using the most number of ESTs together to 
	 * form a new one and return it.
	 */
	 private String processLeftEndsWithInclusion(ArrayList<LeftEnd> includeStrs) {
		 if (includeStrs.size() == 0) {
			 return "";
		 } else if (includeStrs.size() == 1) {
			 int[][]dGraph = genDGraph();
			 String tStr = getInfoOfLeftEnd(includeStrs.get(0).index, dGraph, 1);
			 int type = Integer.parseInt(tStr.substring(0,1));
			 tStr = tStr.substring(1); 
			 if (type == 2) { //typeII AS
				 writeToFile(tStr, "ast2Consensus.out"); //write to a file, one consensus in one line, no comment line.
			 } else if (type == 1) { //typeI AS
				 writeToFile(tStr, "asConsensus.out"); //write to a file, one consensus in one line, no comment line.
			 }
			 return tStr;
		 } else {
			 int maxLen = 0; //the maximal length of the first EST.
			 int idxMaxLen = 0;
			 int maxNumNodes = 0;
			 int idxMaxNumNodes = 0;
			 
			 for (int i=0; i<includeStrs.size(); i++) {
				 int tLen = includeStrs.get(i).lenOfSeq;
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

			 String s1 = includeStrs.get(idxMaxLen).seq;
			 usedNodes[includeStrs.get(idxMaxLen).index] = 1;	//mark that the node is used.
			 int[][]dGraph = genDGraph();
			 String tStr = getInfoOfLeftEnd(includeStrs.get(idxMaxNumNodes).index, dGraph, 1);
			 //if there is alternative splicing, the consensus sequences will be separated by "\n".
			 //if not, s2 only has one sequence.
			 int type = Integer.parseInt(tStr.substring(0,1));
			 tStr = tStr.substring(1);
			 String retStr = "";
			 if (type == 2) { //typeII AS
				 retStr = tStr;
				 writeToFile(tStr, "ast2Consensus.out"); //write to a file, one consensus in one line, no comment line.
			 } else { //typeI AS, or not AS
				 //We have to use linear space alignment to avoid Java out-of-memory error.
				 Substitution sub = new ForGene();
				 String s2 = tStr;
				 SWSmart smart = (new SWSmart (sub, 2, s1, s2));
				 String[] strs = smart.getMatch();
				 int offset = s1.indexOf(strs[0].replace("-", ""));
				 retStr += s1.substring(0, offset) + s2 + "\n";
				 
				 if (type == 1) {
					 writeToFile(retStr, "asConsensus.out"); //write to a file, one consensus in one line, no comment line.
				 }
			 }

			 return retStr;
		 }
	 }

	/*
	 * get the number of used nodes when making a consensus. It has the same value as ret[1] from "reconstructSeq".
	 */
	private int getNumUsedNodes(ArrayList<StartPos> a) {
		int ret = 0;
		int sizeOfa = a.size();
		if (sizeOfa == 0) {
			return 0;
		} else if (sizeOfa == 1) {
			return 1;
		}

		StartPos[] tmpResultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			tmpResultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(tmpResultArray);

		TreeSet<UsedNode> addedNodes = addInclusionNodes(tmpResultArray);  //add all those related inclusion nodes into it for reconstruction.
		ret = addedNodes.size();
		
		return ret;
	}

	/*
	 * reconstruct a sequence which starts from a left end.
	 * return: ret[0]-the consensus, ret[1]-the total number of nodes used for reconstruction.
	 */
	public String[] reconstructSeq(ArrayList<StartPos> a) {
		String[] ret = new String[2];
		int sizeOfa = a.size();
		if (sizeOfa == 0) {
			return null;
		} else if (sizeOfa == 1) {
			 //ret[0] = g.getSeqOfNode(a.get(0).index);
			 //ret[1] = Integer.toString(1);
			 //usedNodes[a.get(0).index] = 1;	//mark that the node is used.
			 return null; //this node should be a singleton if it is not used in other place. it shouldn't appear as a consensus sequence.
		}
		
		StartPos[] tmpResultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			tmpResultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(tmpResultArray);
		
		TreeSet<UsedNode> addedNodes = addInclusionNodes(tmpResultArray);  //add all those related inclusion nodes into it for reconstruction.
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
		usedNodes[resultArray[0].index] = 1;	//mark that the node is used.
		String curSeq = "";
		int len = resultArray.length - 1;
		for (int i=1; i<=len; i++) {
			curSeq = g.getSeqOfNode(resultArray[i].index);
			usedNodes[resultArray[i].index] = 1;	//mark that the node is used.
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

	private String reconstructAS(ArrayList<StartPos> a, int[] asPos) {
		int sizeOfa = a.size();
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		//separate the ESTs which is related and not related to alternative splicing
		ArrayList<StartPos> noAsSet = new ArrayList<StartPos> ();
		ArrayList<StartPos> asSet = new ArrayList<StartPos> ();
		for (int i=0; i<sizeOfa; i++) {
			int tmp = resultArray[i].pos;
			if (tmp <= (asPos[0])) {
				noAsSet.add(resultArray[i]);
			} 
//			if (tmp > asPos[1]) {
			if (tmp > (asPos[0])) {
				asSet.add(resultArray[i]);
			}
		}
		String frontPart = reconstructSeq(noAsSet)[0];
		
		//put all the ESTs which is related to alternative splicing into a file for reconstruction
		try{ 
			/*
			 * print consensus sequences
			 */
			File outFile1 = new File("asEstFile.fa");
			boolean bExists = outFile1.exists();
			if (bExists) {
				outFile1.delete();
			}
			BufferedWriter out1 = new BufferedWriter(new FileWriter(outFile1, true));
			for (int i=0; i<asSet.size(); i++) {
				int idx = asSet.get(i).index;
				
				out1.write(">" + g.getCommentOfNode(idx) + "\n");
				out1.write(g.getSeqOfNode(idx) + "\n");
			}
			out1.flush();
			out1.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 

		
		return frontPart;
	}
	
	private String reconstructAS2(ArrayList<StartPos> a, int[] asPos) {
		int sizeOfa = a.size();
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		ArrayList<StartPos> asSet = new ArrayList<StartPos> ();
		for (int i=0; i<sizeOfa; i++) {
			int tmp = resultArray[i].pos;
//			if (tmp > asPos[1]) {
			if (tmp > asPos[0]) {
				asSet.add(resultArray[i]);
			}
		}
		String firstStr = reconstructSeq(asSet)[0];
		String secondStr = reconstructSeq(a)[0];
		
		
		return firstStr+"\n"+secondStr;
	}

	private String concatConsensus(String s1, String s2) {
		Substitution sub = new ForGene();
		SWSmart smart = (new SWSmart (sub, 2, s1, s2));
		String[] strs = smart.getMatch();
		int offset = s1.indexOf(strs[0].replace("-", ""));
		int offset2 = s2.indexOf(strs[1].replace("-", "")) + strs[1].length();
		
		return s1.substring(0, offset) + strs[0].replace("-", "") + s2.substring(offset2) + "\n";
	}

	/*
	 * check if there is typeI alternative splicing.
	 * Type I: AB, AC, AD
	 * Type II: AB, ACB
	 * If there is, return the range of starting positions where alternative splicing starts to occur.
	 */
/*	public int[] checkAS1(ArrayList<StartPos> a) {
		System.out.println("\ncheck type I alternative splicing");

		int pos[] = new int[2];
		pos[0] = -1;
		pos[1] = -1;
		int sizeOfa = a.size();
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		int asWindow = 150;
		
		int startPoint = 0;
		int endPoint = startPoint + asWindow;
		boolean flag1 = false;
		boolean flag2 = false;
		ArrayList<String> seqInWindow = new ArrayList<String> ();
		for (int i=0; i<sizeOfa; i++) {
			if (resultArray[i].pos < endPoint) {
				System.out.println(g.getCommentOfNode(resultArray[i].index));
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			} else {
				//process seqInWindow. Get the number of overlapping pairs by calculating each pair of sequences.
				int numOfNonOvl = getNonOvlNum(seqInWindow);
				int threshold = seqInWindow.size()/2;
				System.out.println("start from " + startPoint + ", end at " + endPoint);
				System.out.println(numOfNonOvl + " non overlap, " + seqInWindow.size() + " Ests.");
				if (flag1) {
					if (numOfNonOvl > threshold) { //alternative splicing
						flag2 = true; //if the last loop is assumed AS, flag2 will not be true.
						pos[0] = startPoint;
						pos[1] = endPoint;
						break;
					} else {
//						pos[0] = -1;
//						pos[1] = -1;
						flag1 = false;
					}
				} else {
					if (numOfNonOvl > threshold) { //assumed alternative splicing
						flag1 = true;
					}
				}
				
				//reset startPoint and endPoint
				startPoint = endPoint;
				endPoint = startPoint + asWindow;
				
				//clear seqInWindow
				seqInWindow.clear();
				//put last sequence into seqInWindow
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			}
		}

		if (!flag2) { //if no alternative splicing occurs, continue to process the last seqInWindow to check if AS exists. 
			//process last seqInWindow. Get the number of overlapping pairs by calculating each pair of sequences.
			int numOfNonOvl = getNonOvlNum(seqInWindow);
			int threshold = seqInWindow.size() * 2;
			System.out.println("start from " + startPoint + ", end at " + endPoint);
			System.out.println(numOfNonOvl + " non overlap, " + seqInWindow.size() + " Ests.");
			if (numOfNonOvl > threshold) { //alternative splicing
				pos[0] = startPoint;
				pos[1] = endPoint;
			}
		}

		if (!flag2) {
			pos[0] = -1;
			pos[1] = -1;
		} 
		return pos;
	}
*/

	
	public int[] checkAS1(ArrayList<StartPos> a) {
		System.out.println("\ncheck type I alternative splicing");

		int pos[] = new int[2];
		pos[0] = -1;
		pos[1] = -1;
		int sizeOfa = a.size();
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		int asWindow = 150;
		
		int startPoint = 0;
		int endPoint = startPoint + asWindow;
		ArrayList<String> seqInWindow = new ArrayList<String> ();
		for (int i=0; i<sizeOfa; i++) {
			if (resultArray[i].pos < endPoint) {
				System.out.println(g.getCommentOfNode(resultArray[i].index));
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			} else {
				//process seqInWindow. Get the number of overlapping pairs by calculating each pair of sequences.
				int numOfNonOvl = getNonOvlNum(seqInWindow);
				int threshold = seqInWindow.size() * seqInWindow.size() * 3/16;
				System.out.println("start from " + startPoint + ", end at " + endPoint);
				System.out.println(numOfNonOvl + " non overlap, " + seqInWindow.size() + " Ests.");
				if (numOfNonOvl > threshold) { //alternative splicing
					pos[0] = startPoint;
					pos[1] = endPoint;
					break;
				} 
				
				//reset startPoint and endPoint
				startPoint = endPoint;
				endPoint = startPoint + asWindow;
				
				//clear seqInWindow
				seqInWindow.clear();
				//put last sequence into seqInWindow
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			}
		}

		return pos;
	}

	
	/*
	 * check if there is typeII alternative splicing.
	 * Type I: AB, AC, AD
	 * Type II: AB, ACB
	 * If there is, return the range of starting positions where alternative splicing starts to occur.
	 */
	public int[] checkAS2(ArrayList<StartPos> a) {
		System.out.println("\ncheck type II alternative splicing");

		int pos[] = new int[2];
		pos[0] = -1;
		pos[1] = -1;
		int sizeOfa = a.size();
		StartPos[] resultArray = new StartPos[sizeOfa]; //store the starting positions of ests
		for (int i=0; i<sizeOfa; i++) {
			resultArray[i] = a.get(i);
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);
		
		int asWindow = 150;
		
		int startPoint = 0;
		int endPoint = startPoint + asWindow;
		boolean flag1 = false;
		boolean flag2 = false;
		ArrayList<String> seqInWindow = new ArrayList<String> ();
		int totalLen = resultArray[sizeOfa-1].pos + g.getSeqOfNode(resultArray[sizeOfa-1].index).length();
		int numOfWindow = (int)(totalLen/asWindow);
		int threshold = (int) ((sizeOfa/numOfWindow) * 1.5);
		System.out.println("Threshold="+threshold+"; totalLen=" + totalLen + "; numOfWindow=" + numOfWindow);
		
		for (int i=0; i<sizeOfa; i++) {
			if (resultArray[i].pos < endPoint) {
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			} else {
				//process seqInWindow. Get the number of overlapping pairs by calculating each pair of sequences.
				int num = seqInWindow.size();
				System.out.println("start from " + startPoint + ", end at " + endPoint);
				System.out.println("	" + num + " Ests in this window.");
				if (flag1) {
					if (num > threshold) { //alternative splicing
						flag2 = true; //if the last loop is assumed AS, flag2 will not be true.
						break;
					} else {
						pos[0] = -1;
						pos[1] = -1;
						flag1 = false;
					}
				} else {
					if (num > threshold) { //assumed alternative splicing
						pos[0] = startPoint;
						pos[1] = endPoint;
						flag1 = true;
					}
				}
				
				//reset startPoint and endPoint
				startPoint = endPoint;
				endPoint = startPoint + asWindow;
				
				//clear seqInWindow
				seqInWindow.clear();
				//put last sequence into seqInWindow
				seqInWindow.add(g.getSeqOfNode(resultArray[i].index));
			}
		}
		
		if (!flag2) {
			pos[0] = -1;
			pos[1] = -1;
		} 
		return pos;
	}

	private int getNonOvlNum(ArrayList<String> seqInWindow) {
		int size = seqInWindow.size();
		int numOfNonOvl = 0;
		for (int i=0; i<size; i++) {
			String s1 = seqInWindow.get(i);
			for (int j=i+1; j<size; j++) {
				String s2 = seqInWindow.get(j);
				int distance = Math.abs((g.ovl.getOVLDistance(s1, s2))[1]);
				if (distance == INT_MAX) { //do not overlap or include.
					numOfNonOvl++;
				}
			}
		}
		return numOfNonOvl;
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

	class LeftEnd implements Comparable<LeftEnd> {
		 int index;
		 int lenOfSeq;
		 int numOfUsedNodes;
		 String seq;
		 public LeftEnd(int idx, int n, String s) {
			 index = idx;
			 lenOfSeq = s.length();
			 numOfUsedNodes = n;
			 seq = s;
		 }

		 public int compareTo(LeftEnd other) {
			 //Returns 0 if the argument is equal to this; 			
			 //a value less than 0 if the argument is greater than this; 
			 //and a value greater than 0 if the argument is less than this. 
			 if (this.lenOfSeq == other.lenOfSeq) {
				 return 0;
			 } else if (this.lenOfSeq > other.lenOfSeq) {
				 return 1;
			 } else {
				 return -1;
			 }
		 }
	}	
}



