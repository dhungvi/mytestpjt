// get six-tuple for each node.

import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.Properties;


public class SixTuplesGeneration {
	Graph g;
	InclusionNodes incNodes;
	ArrayList<SixTuple> leftMostNodes;
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
	
	public SixTuplesGeneration(Properties props, Graph graph, InclusionNodes inc) {
		g = graph;
		incNodes = inc;
		alignArray = null;
		leftMostNodes = new ArrayList<SixTuple> ();
		init();
	}
	
	public ArrayList<SixTuple> getAlignArray() {
		return alignArray;
	}
	
	public ArrayList<SixTuple> getLeftEnds() {
		return leftMostNodes;
	}
	
	private void init() {
		// Get the components of the time
	    long time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to generate 6-tuples.");
		createAlignArray();
		System.out.println("End to generate 6-tuples.");
		System.out.println("The time used to generate 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1) + "\n");
		System.out.println("There are " + incNodes.getSize() + " nodes in the inclusion list.");
		System.out.println("There are " + alignArray.size() + " nodes in alignArray.\n");
		
		time1 = new GregorianCalendar().getTimeInMillis();
		System.out.println("Start to process 6-tuples.");
		processAlignArray();
		System.out.println("End to process 6-tuples.");
		System.out.println("The time used to process 6-tuples is " + (new GregorianCalendar().getTimeInMillis()-time1) + "\n");

	}
	
	private void createAlignArray() {
		alignArray = g.get2CloseNodesFromMST();
	}

	/* 
	 * Get the assumed left ends(alignNodes[][0]==-1), check them to find all the real left ends.
	 */
	private void processAlignArray() {
		ArrayList<SixTuple> rightMostNodes = new ArrayList<SixTuple> ();
		//get all the nodes which has no left nodes or no right nodes to them
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
				SixTuple lNode = g.get2CloseNodesFromGrand(cNode, curTuple);
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
			SixTuple lNode = g.get2CloseNodesFromGrand(cNode, curTuple);
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
		System.out.println("\nThere are " + leftMostNodes.size() + " left-most nodes after running 3 levels.");
		System.out.println("There are " + rightMostNodes.size() + " right-most nodes after running 3 levels.");
		
		/* Recalculate 6-tuples for all the current left nodes in order to remove all the false left ends.
		 * Specifically, for those assumed left ends,start to calculate from fourth level until meeting one node which 
		 * makes six-tuple[0] != -1 or until the level we specified in the property file, then return the six-tuple.
		 * If we fail to find any node, we consider it a real left end.
		 */
		for (int i=0; i<leftMostNodes.size(); i++) {
			SixTuple curTuple = leftMostNodes.get(i);
			int tEnd = curTuple.curNode; //index of the node
			SixTuple tmpTuple = g.checkLeftEndFromMST(tEnd, curTuple);
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
			SixTuple tmpTuple = g.checkRightEndFromMST(tEnd, curTuple);
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
		System.out.println("There are " + rightMostNodes.size() + " right-most nodes after checking right ends.");
		
		/* 
		 * Remove those false left ends.
		 * 
		 * construct a temporary directed graph from alignArray, 
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
		System.out.println("\nThere are " + leftMostNodes.size() + " left-most nodes after processing false left ends.");
		System.out.println("There are " + rightMostNodes.size() + " right-most nodes after processing false left ends.");
	}


}