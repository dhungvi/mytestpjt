import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.Properties;

import com.mhhe.clrs2e.MergeSort;
import com.mhhe.clrs2e.Prim;
import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

/*
 * This class is a subclass of ESTAssembly. It assigns NewGraph to member variable g.
 */
public class NewESTAssembly extends ESTAssembly{
	public NewESTAssembly(Properties props) {
		super(props);
		g = new NewGraph(props);
	}

	/* Overload the method of parent class
	 * Find the left-most node(alignNodes[][0]==-1);
	 * Use alignNodes to construct a tree whose values are maximal, that is, 
	 * 		the tree has maximal overlap length;
	 * Calculate the starting position of each node.
	 * 
	 * This function is different from the same-name method in its parent class in that 
	 * in order to get starting positions it constructs a MST instead of Maximum Spanning
	 * tree like what the parent method does. The weight of the Minimum Spanning tree is 
	 * the overlap distance instead of overlap length.
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
		System.out.println("The number of original left ends is " + numOfLeftMostNodes);
		if (numOfLeftMostNodes > 1) {
			for (int i=0; i<leftMostNodes.size(); i++) {
				int cNode = leftMostNodes.get(i).intValue();
				//int[] lNode = g.get2CloseNodes(cNode);
				int[] lNode = g.get2CloseNodesFromGrand(mstForG, cNode, alignArray[cNode]);
				alignArray[cNode][0] = lNode[0];
				alignArray[cNode][1] = lNode[1];	//overlap length
				alignArray[cNode][2] = lNode[2];	//distance
				
				if (alignArray[cNode][5] > lNode[5]) { //get a smaller distance
					alignArray[cNode][3] = lNode[3];
					alignArray[cNode][4] = lNode[4];
					alignArray[cNode][5] = lNode[5];
				}
			}
		} 
		
		//find the two most closest nodes for those rightmostnodes(because we need these information to make a directed maximum spanning tree)
		/*for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][3] == -1) {
				//int[] rNode = g.get2CloseNodes(i);
				int[] rNode = g.get2CloseNodesFromGrand(mstForG, i, alignArray[i]);
				if (Math.abs(alignArray[i][2]) > Math.abs(rNode[2])) { //get a smaller distance
					alignArray[i][0] = rNode[0];
					alignArray[i][1] = rNode[1];
					alignArray[i][2] = rNode[2];
				}
				alignArray[i][3] = rNode[3];			
				alignArray[i][4] = rNode[4];
				alignArray[i][5] = rNode[5];			
			}
		} */

		/* construct a directed graph from alignArray, 
		 *  	the second dimension has two elements:
		 *  			index of starting node,
		 *  			index of ending node, 
		 *  			weight between them (positive value, weight is abs(their distance)).
		 *  	if there is no edge, weight=INT_MAX.
		 */
		int len = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] != -1) {
				len++;
			}
			if (alignArray[i][3] != -1) {
				len++;
			}
		}
		
		int[][] dGraph = new int[len][4];
		int indexOfDGraph = 0;
		for (int i=0; i<alignArray.length; i++) {
			if (alignArray[i][0] != -1) {
				dGraph[indexOfDGraph][0] = alignArray[i][0];
				dGraph[indexOfDGraph][1] = i;
				dGraph[indexOfDGraph][2] = Math.abs(alignArray[i][2]);	//distance
				dGraph[indexOfDGraph][3] = Math.abs(alignArray[i][1]);	//overlap length
				indexOfDGraph++;
			}
			
			if (alignArray[i][3] != -1) {
				dGraph[indexOfDGraph][0] = i;
				dGraph[indexOfDGraph][1] = alignArray[i][3];
				dGraph[indexOfDGraph][2] = alignArray[i][5];	//distance
				dGraph[indexOfDGraph][3] = alignArray[i][4];	//overlap length
				indexOfDGraph++;
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
		System.out.println("Original minimum Spanning Tree:");
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
		//int[][] tmpGraph = dGraph;	//keep the original values of dGraph to pass to constructMinTree as parameter.
									//because dGraph will be changed in the following for loop.
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
			primMST = constructMinTree(alignArray.length, dGraph);
			//primMST = constructMaxTree(alignArray.length, dGraph);
			
			//put leftEnd node to index 0 in array sPos to be consistent with dGraph and primMST
			sPos[leftEnd] = sPos[0];
			sPos[0] = Integer.parseInt(g.getNameOfNode(leftEnd)); //starting position of the node
			//get starting positions for the nodes in primMST
			getStartPos(0, leftEnd, primMST, dGraph);
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
			//print mst
			//System.out.println("Minimum Spanning Tree for starting positions:");
			//System.out.println(primMST);
		}
		//this.printConsensus(leftMostNodes);
	}

	protected WeightedAdjacencyListGraph constructMaxTree(int nOfNodes, int[][] g) {
		// Make a directed graph.
		WeightedAdjacencyListGraph dGraph =
		    new WeightedAdjacencyListGraph(nOfNodes, true);
		for (int i=0; i<nOfNodes; i++) {
			dGraph.addVertex(i, Integer.toString(i));
		}
		for (int j=0; j<g.length; j++) {
			dGraph.addEdge(g[j][0], g[j][1], -g[j][2]);
		}
		WeightedAdjacencyListGraph mst = (new Prim()).computeMST(dGraph);
		return mst;
	}

	/*
	 * print elements in 'dGraph'
	 */
	protected void printDgraph(int[][] d){
		for (int i=0; i<d.length; i++) {
			System.out.println(d[i][0] + "\t" + d[i][1] + "\t" + d[i][2] + "\t" + d[i][3]);
		}
		System.out.println();
	}
	
	/*
	 * reconstruct the sequence
	 */
	public String reconstruct() {
		//sort the array sPos
		StartPos[] resultArray = new StartPos[sPos.length]; //store the starting positions of ests
		for (int i=0; i<sPos.length; i++) {
			resultArray[i] = new StartPos(sPos[i], g.getSeqOfNode(i));
		}
		MergeSort merge = new MergeSort();
		merge.sort(resultArray);

		//local align every two adjacent nodes, e.g, node 1 and 2, node 2 and 3..., node i-1 and node i.
		//put them into two arraylists, the local alingment sequences of node i and i+1 are put into arraylist 1 and 2 respectively.
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
		for (int i=0; i<size-1; i++) {
			String s1 = arrList2.get(i);
			String s2 = arrList1.get(i+1);
			String s1Tmp = s1.replace("-", "");
			String s2Tmp = s2.replace("-", "");
			String subS1 = (g.d2.getLocalAlignment(s1Tmp, s2Tmp))[0];
			int pos = s1.indexOf("-");
			int offset = s1Tmp.indexOf(subS1);
			while (pos != -1) {
				String tmp = s1.substring(pos+1).replace("-", "");
				s1Tmp = s1.substring(0, pos+1) + tmp;
				if (s1Tmp.indexOf(subS1) != -1) {
					offset++;
				}
				pos = s1.indexOf("-", pos+1);
			}
			sposArr1[i+1] = sposArr2[i]+offset;
			sposArr2[i+1] = sposArr1[i+1];
		}
		
				
		// Create an array which stores all the bases with the same position for all the ESTs.
		int lenOfArray = resultArray[resultArray.length-1].pos+resultArray[resultArray.length-1].seq.length(); 
		ArrayList<Character> [] tmpArraylists = new ArrayList [lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			tmpArraylists[i] = new ArrayList<Character>();
		}
		// Put all the bases into the array
		for (int i=0; i<size; i++) {
			int tmpSPos1 = sposArr1[i];
			String tmpStr1 = arrList1.get(i);
			int tmpSPos2 = sposArr2[i];
			String tmpStr2 = arrList2.get(i);
			for (int j=0; j<tmpStr1.length(); j++) {
				int p = tmpSPos1+j;
				if (p < lenOfArray) {
					tmpArraylists[p].add(tmpStr1.charAt(j));
				}
			}
			for (int j=0; j<tmpStr2.length(); j++) {
				int p = tmpSPos2+j;
				if (p < lenOfArray) {
					tmpArraylists[p].add(tmpStr2.charAt(j));
				}
			}
		}
		
		// Calculate consensus base for each position, and put them into an char array
		char[] consensus = new char[lenOfArray];
		for (int i=0; i<lenOfArray; i++) {
			consensus[i] = getConsensusBase(tmpArraylists[i]);
		}
		
		String retStr = String.valueOf(consensus);
		retStr = retStr.replace("-", "");
		return retStr;
	}

	/*
	 * overload the same-name method of parent class.
	 * print the original sequence and multiple consensus into a file which is specified in the property file.
	 * The printed consensuses start from all the assumed starting position.
	 */
	/*
	public void printConsensus(ArrayList <Integer> leftEnds) {
		try{ 
			File outFile = new File(consensusFileName);
			boolean bExists = outFile.exists();
			if (bExists) {
				outFile.delete();
			}
			BufferedWriter out = new BufferedWriter(new FileWriter(outFile, true));

			
			// print the original sequence
			File oriFile = (new File(oriFileName));
			if (!oriFile.exists()) {
				System.out.println("SourceFile does not exist!");
				return;
			}
			BufferedReader in = new BufferedReader(new FileReader(oriFile));
			out.write(in.readLine());
			out.write("\n");
			in.close();	

			
			// print all the consensuses
			String tmpStr = getConsensus();
			for (int i=0; i<leftEnds.size()-1; i++) {
				int pos1 = leftEnds.get(i).intValue();
				int pos2 = leftEnds.get(i+1).intValue();
				String s = tmpStr.substring(pos1, pos2+1);
				for (int j=0; j<pos1; j++) {
					out.write(" ");
				}
				out.write(s);
				out.write("\n");
			}
			int pos2 = leftEnds.get(leftEnds.size()-1);
			String s = tmpStr.substring(pos2);
			for (int j=0; j<pos2; j++) {
				out.write(" ");
			}
			out.write(s);
			out.write("\n");
			
			out.flush();
			out.close();
		}catch(IOException e){ 
			System.out.println(e.toString());
		} 
	}*/
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
			String oriStr = in.readLine();
			in.close();	

			/*
			 * print the consensus
			 */
			out.write(oriStr);
			out.write("\n");
			String consensus = this.reconstruct();
			out.write(consensus);
			out.write("\n");
			
			//out.write((g.d2.getLocalAlignment(oriStr, consensus))[2]);
			out.flush();
			out.close();
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
	 * The same-name method as the parent class. But it adds one more parameter d.
	 * Calculate starting positions for each node
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
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Properties props = null;
		try {
			props = getProperties("config.properties");
		} catch (IOException e) {
			System.err.println("Get config.properties failed, " + e);
	    	return;
		}

		NewESTAssembly assemble = new NewESTAssembly(props);

		assemble.readEstFile();

		// Get the components of the time
	    long time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to generate MST.");
		assemble.createAlignArray();
		System.out.println("End to generate 6-tuples.");
		System.out.println("The time used to process 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to process 6-tuples.");
		assemble.processAlignArray();
		System.out.println("End to process 6-tuples.");
		System.out.println("The time used to process 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1));
		
		assemble.printSPos();
		//assemble.printConsensus();
		//assemble.printEsts();
		assemble.calcInversion();
		
		time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to reconstruct.");
		assemble.printConsensus();
		System.out.println("End to reconstruct.");
		System.out.println("The time used to reconstruct is " + (new GregorianCalendar().getTimeInMillis()-time1));

	}

}
