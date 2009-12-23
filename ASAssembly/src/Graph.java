import java.util.ArrayList;
import java.util.Properties;
import java.util.Stack;

import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

public class Graph {
	protected final int INT_MAX = Integer.MAX_VALUE;
	protected final int INT_MIN = Integer.MIN_VALUE;
	private int numOfLevels;

	ArrayList<Node> graphNodes;
	OvlDistance ovl;
	CalculatedOvlDistance calDist;
	WeightedAdjacencyListGraph mst;	//minimum spanning tree generated from peace
	InclusionNodes inc;
	
	public Graph(Properties p, InclusionNodes in) {
		numOfLevels = Integer.parseInt(p.getProperty("NumOfLevels"));
		graphNodes = new ArrayList<Node> ();
		ovl = new OvlDistance(p);
		calDist = new CalculatedOvlDistance();
		inc = in;
	}
	
	public void setMst(WeightedAdjacencyListGraph m) {
		mst = m;
	}
	
	public void addNode(Node s) {
		graphNodes.add(s);
	}
	
	public void removeNode (int index) {
		graphNodes.remove(index);
	}
	
	public int getSizeofGraph() {
		return graphNodes.size();
	}
	
	/*
	 * get length of the node with index i
	 */
	public int getLenOfNode(int i) {
		return graphNodes.get(i).getLen();
	}

	/*
	 * get ID of the node with index i
	 */
	public String getNameOfNode(int i) {
		return graphNodes.get(i).getName();
	}
	
	/*
	 * get sequence of the node with index i
	 */
	public String getSeqOfNode(int i) {
		return graphNodes.get(i).getSeq();
	}
	
	/*
	 * get comment of the node with index i
	 */
	public String getCommentOfNode(int i) {
		return graphNodes.get(i).getComment();
	}
	


	private ArrayList<SixTuple> handleInclusion() {
		ArrayList<SixTuple> nodes = new ArrayList<SixTuple>();
		int nOfNodes = mst.getCardV();
		
		for (int i=0; i<nOfNodes; i++) {
			WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(i);
			boolean flag = true;
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index = v.getIndex();
				
				String curSeq = graphNodes.get(i).getNodeStr();
				String comSeq = graphNodes.get(index).getNodeStr();

				int[] ovlDis = calDist.searchDistance(i, index);
				if (ovlDis == null) {
					ovlDis = (ovl).getOVLDistance(curSeq, comSeq);
					//add to CalculatedOvlDistance
					calDist.addDistance(i, index, ovlDis[1], ovlDis[0]);
				}
				
				if (curSeq.length() <= comSeq.length()) {
					if (ovlDis[1] == INT_MIN) { //has inclusion
						inc.addNode(i, index); //get rid of i and put it into inclusion list.
						flag = false;
						break;
					} 
				} 
			}
			if (flag) {
				nodes.add(new SixTuple(i));
			}
		}
		return nodes;
		
	}
	/**
	 * Get two closest nodes which is on the left and on the right to every node 
	 * from the input minimum spanning tree, and store the data into an array list.
	 * This function removes inclusion nodes from the tree, and then finds parent and 
	 * children for all the left nodes, and calculates overlap distance between the node and 
	 * others, select the two nodes with minimum left and right overlap distance.
	 * 
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @return an array list of SixTuple which stores two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the tree. 
	 */
	public ArrayList<SixTuple> get2CloseNodesFromMST() {
		ArrayList<SixTuple> alignedNodes = handleInclusion();
		int nOfNodes = alignedNodes.size();
		
		for (int i=0; i<nOfNodes; i++) {
			int curIdx = alignedNodes.get(i).curNode;
			if (inc.containInclusionNode(curIdx)) continue;
			int leftNode = -1;
			int rightNode = -1;
			int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
			int minRight = INT_MAX;	//minimum right distance
			int overlapLeft = 0;
			int overlapRight = 0;
			
			WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(curIdx);
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index = v.getIndex();
				if (inc.containInclusionNode(index)) continue; 

				int[] ovlDis = calDist.searchDistance(curIdx, index);
				if (ovlDis == null) {

					ovlDis = (ovl).getOVLDistance(graphNodes.get(curIdx).getNodeStr(), 
							graphNodes.get(index).getNodeStr());
					//add to CalculatedOvlDistance
					calDist.addDistance(curIdx, index, ovlDis[1], ovlDis[0]);
				}
				
				if (ovlDis[1] != INT_MAX) {	// there is overlap between them
					if (ovlDis[0] < 0) {
						if (ovlDis[1] > maxLeft){
							maxLeft = ovlDis[1];
							overlapLeft = ovlDis[0];
							leftNode = index;
						} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
							if (Math.abs(ovlDis[0]) > Math.abs(overlapLeft)) {
								overlapLeft = ovlDis[0];
								leftNode = index;
							}
						}
					}
					if (ovlDis[0] > 0) {
						if (ovlDis[1] < minRight) {
							minRight = ovlDis[1];
							overlapRight = ovlDis[0];
							rightNode = index;
						} else if (ovlDis[1] == minRight) {	//if they are equal, find that one with maximal overlap
							if (Math.abs(ovlDis[0]) > Math.abs(overlapRight)) {
								overlapRight = ovlDis[0];
								rightNode = index;
							}
						}
					}
				}
			}
			//leftNode, index of node on the left
			//overlapLeft, overlap length
			//maxLeft, overlap distance
			//rightNode, index of node on the right
			//overlapRight, overlap length
			//minRight, overlap distance
			alignedNodes.get(i).setSixTuple(leftNode, overlapLeft, maxLeft, rightNode, overlapRight, minRight);
		}
		return alignedNodes;
	}

	
	/**
	 * Get two closest nodes which is on the left and on the right to the 'index' node
	 * from the input minimum spanning tree, and store the data into an array.
	 * To get the six-tuple for the node, the function calculates distance from 
	 * this node to all the other nodes which are less than or equal to three levels from it. 
	 * 
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param index The index of current node in the tree and graph.
	 * @param sixTuple The SixTuple for this node with the index
	 * @return SixTuple for the current node which has the input "index". 
	 *						the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the distance with + or -.
	 * 						the fourth is the index of node on the right, the fifth is the overlap length with + or -.
	 * 						the sixth is the distance with + or -.
	 * 						For the first and fourth one, if no node is found, the value is -1;
	 * 						For the second and fifth one, if no node is found, the value is 0;
	 * 						For the third and sixth one, if no node is found, the value is INT_MIN or INT_MAX.
	 */
	public SixTuple get2CloseNodesFromGrand(int index, SixTuple sixTuple) {
		SixTuple closeNode = new SixTuple();
		
		int leftNode = -1;
		int rightNode = -1;
		//int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
		//int minRight = INT_MAX;	//minimum right distance
		//int overlapLeft = 0;
		//int overlapRight = 0;
		int maxLeft = sixTuple.lDis;
		int minRight = sixTuple.rDis;
		int overlapLeft = sixTuple.lOvlLen;
		int overlapRight = sixTuple.rOvlLen;
		
		//put all the nodes within three levels from the current node into the stack "allNodes". Do not 
		// include parents and children because they have been processed.
		Stack<Integer> allNodes = new Stack<Integer> ();
		Stack<Integer> partNodes = new Stack<Integer> ();
		partNodes.push(Integer.valueOf(index));
		
		for (int i=0; i<3; i++) {
			Stack<Integer> tmpStack = new Stack<Integer> ();
			
			while (!partNodes.empty()) {
				int tmpIndex = partNodes.pop();
				WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(tmpIndex);
				while (ite.hasNext()) {
					Vertex v = (Vertex) ite.next();
					int tIndex = v.getIndex();
					tmpStack.push(Integer.valueOf(tIndex));
					
					if ((i != 0) && (tIndex != index)){ //Do not include parents and children because they have been processed.
						allNodes.push(Integer.valueOf(tIndex));
					}
				}
			}
			partNodes = tmpStack;
		}

		//find two closest nodes
		String s1 = graphNodes.get(index).getNodeStr();
		while (!allNodes.isEmpty()) {
			int tmpIndex = allNodes.pop();
			if (inc.containInclusionNode(tmpIndex)) continue; 
			String s2 = graphNodes.get(tmpIndex).getNodeStr();

			int[] ovlDis = calDist.searchDistance(index, tmpIndex);
			if (ovlDis == null) {
				ovlDis = (ovl).getOVLDistance(s1, s2);

				//add to CalculatedOvlDistance
				calDist.addDistance(index, tmpIndex, ovlDis[1], ovlDis[0]);
			}
			
			if (ovlDis[1] == INT_MIN) {	// there is inclusion between them
				if (s1.length() >= s2.length()) {
					inc.addNode2(tmpIndex, index);
				} else {
					inc.addNode2(index, tmpIndex);
				}
			} else if (ovlDis[1] != INT_MAX) {	// there is overlap between them
				if (ovlDis[0] < 0) {
					if (ovlDis[1] > maxLeft){
						maxLeft = ovlDis[1];
						overlapLeft = ovlDis[0];
						leftNode = tmpIndex;
					} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapLeft)) {
							overlapLeft = ovlDis[0];
							leftNode = tmpIndex;
						}
					}
				}
				if (ovlDis[0] > 0) {
					if (ovlDis[1] < minRight) {
						minRight = ovlDis[1];
						overlapRight = ovlDis[0];
						rightNode = tmpIndex;
					} else if (ovlDis[1] == minRight) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapRight)) {
							overlapRight = ovlDis[0];
							rightNode = tmpIndex;
						}
					}
				}
			}
		}
		closeNode.leftNode = leftNode;	//index of node on the left
		closeNode.lOvlLen = overlapLeft;	//overlap length
		closeNode.lDis = maxLeft;	//overlap distance
		
		closeNode.rightNode = rightNode;	//index of node on the right
		closeNode.rOvlLen = overlapRight;	//overlap length
		closeNode.rDis = minRight;	//overlap distance
		
		//System.out.println("Get 6-tuple for node "+index);
		return closeNode;
	}

	/**
	 * Recalculate 6-tuples for an assumed left end in order to make sure it is a real left end.
	 * Specifically, for the assumed left end,start to calculate from fourth level until meeting one node which 
	 * makes six-tuple.leftEnd != -1, or until the level we specified in the property file, then return the six-tuple.
	 * If we fail to find any node, we consider it a real left end and six-tuple.leftEnd will be set to be -1.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param index The index of current node in the tree and graph.
	 * @param sixTuple The SixTuple for this node with the index
	 * @return SixTuple which stores two closest nodes which is to the left and to the right if it is not left end;
	 * if it does, six-tuple.leftEnd = -1.
	 */
	public SixTuple checkLeftEndFromMST(int index, SixTuple sixTuple) {
		if (numOfLevels == 0) { //keep checking until leaves
			numOfLevels = INT_MAX;
		}
		
		Stack<Integer>[] nodes = new Stack[2];
		nodes[0] = new Stack<Integer>();
		nodes[1] = new Stack<Integer>();
		nodes[0].push(Integer.valueOf(index));
		nodes[1].push(Integer.valueOf(-1));
		
		for (int i=0; i<3; i++) { //skip over all the nodes within 3 levels
			nodes = getNodesFromMST(nodes);
		}
		
		for (int foundLevel=4; foundLevel<=numOfLevels; foundLevel++) { //start from level 4
			nodes = getNodesFromMST(nodes);
			
//			System.out.println("GetNodeFromMST for foundLevel=" + foundLevel
//					+ "; index=" + index);
//			System.out.println("\tnumber of nodes = " + nodes[0].size());
		
			if (nodes[0].size() == 0) {
				break;
			} else {
				SixTuple closeNode = findAdjacentNode(nodes[0], index, sixTuple);
				if (closeNode.leftNode != -1) {
//					System.out.println("findAdjacentNode for index=" + index + "; adjNode=" + closeNode.leftNode);
					return closeNode;
				}
			}
		}
		
//		System.out.println("Fail to find AdjacentNode for index=" + index);
		return null;
	}
	
	public SixTuple checkRightEndFromMST(int index, SixTuple sixTuple) {
		if (numOfLevels == 0) { //keep checking until leaves
			numOfLevels = INT_MAX;
		}
		
		Stack<Integer>[] nodes = new Stack[2];
		nodes[0] = new Stack<Integer>();
		nodes[1] = new Stack<Integer>();
		nodes[0].push(Integer.valueOf(index));
		nodes[1].push(Integer.valueOf(-1));
		
		for (int i=0; i<3; i++) { //skip over all the nodes within 3 levels
			nodes = getNodesFromMST(nodes);
		}
		
		for (int foundLevel=4; foundLevel<=numOfLevels; foundLevel++) { //start from level 4
			nodes = getNodesFromMST(nodes);
			
			if (nodes[0].size() == 0) {
				break;
			} else {
				SixTuple closeNode = findAdjacentNode(nodes[0], index, sixTuple);
				if (closeNode.rightNode != -1) {
					return closeNode;
				}
			}
		}
		
		return null;
	}

	/*
	 * find the most adjacent node to the current node from allNodes.
	 * @param allNodes Store indices of all the nodes which will be compared to the current node.
	 * @param index The index of current node.
	 * @sixTuple The sixTuple for the current node.
	 */
	private SixTuple findAdjacentNode(Stack<Integer> nodes, int index, SixTuple sixTuple) {
		Stack<Integer> allNodes = new Stack<Integer> (); //put all the values of nodes into another stack so that we won't change nodes(it's a pointer) because it may be used by the calling method.
		for (int i=0; i<nodes.size(); i++) { //get(0) will get the bottom element of the stack
			allNodes.push(nodes.get(i));
		}
		
		SixTuple closeNode = new SixTuple();
		
		int leftNode = -1;
		int rightNode = -1;
		int maxLeft = sixTuple.lDis;
		int minRight = sixTuple.rDis;
		int overlapLeft = sixTuple.lOvlLen;
		int overlapRight = sixTuple.rOvlLen;

		//find two closest nodes
		String s1 = graphNodes.get(index).getNodeStr();
		while (!allNodes.empty()) {
			int tmpIndex = (Integer) allNodes.pop();
			if (inc.containInclusionNode(tmpIndex)) continue; 
			if (tmpIndex == index) continue;
			String s2 = graphNodes.get(tmpIndex).getNodeStr();

			int[] ovlDis = calDist.searchDistance(index, tmpIndex);
			if (ovlDis == null) {
				ovlDis = (ovl).getOVLDistance(s1, s2);

				//add to CalculatedOvlDistance
				calDist.addDistance(index, tmpIndex, ovlDis[1], ovlDis[0]);
			}

			
			if (ovlDis[1] == INT_MIN) {	// there is inclusion between them
				if (s1.length() >= s2.length()) {
					inc.addNode2(tmpIndex, index);
				} else {
					inc.addNode2(index, tmpIndex);
				}
			} else if (ovlDis[1] != INT_MAX) {	// there is overlap between them
				if (ovlDis[0] < 0) {
					if (ovlDis[1] > maxLeft){
						maxLeft = ovlDis[1];
						overlapLeft = ovlDis[0];
						leftNode = tmpIndex;
					} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapLeft)) {
							overlapLeft = ovlDis[0];
							leftNode = tmpIndex;
						}
					}
				}
				if (ovlDis[0] > 0) {
					if (ovlDis[1] < minRight) {
						minRight = ovlDis[1];
						overlapRight = ovlDis[0];
						rightNode = tmpIndex;
					} else if (ovlDis[1] == minRight) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapRight)) {
							overlapRight = ovlDis[0];
							rightNode = tmpIndex;
						}
					}
				}
			}
		}
		closeNode.leftNode = leftNode;	//index of node on the left
		closeNode.lOvlLen = overlapLeft;	//overlap length
		closeNode.lDis = maxLeft;	//overlap distance
		
		closeNode.rightNode = rightNode;	//index of node on the right
		closeNode.rOvlLen = overlapRight;	//overlap length
		closeNode.rDis = minRight;	//overlap distance
		
		return closeNode;
	}

	/**
	 * Get all the children for the input nodes and put them into a stack.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param nodes The stack array, nodes[0], the indexes of all the nodes for which we will find their children;
	 * 		nodes[1], the indexes of the parent node of every element in nodes[0].
	 * @return a stack array which stores all the found nodes' indexes. ret[0], the indexes of all the found nodes
	 * 		(that is, the children of input nodes); ret[1], the indexes of the parent node of every element in nodes[0].
	 */
	private Stack<Integer>[] getNodesFromMST(Stack<Integer>[] nodes) {
		Stack<Integer>[] ret = new Stack[2];
		ret[0] = new Stack<Integer>();
		ret[1] = new Stack<Integer>();
		while (!(nodes[0].empty())) {
			int curIndex = (Integer) nodes[0].pop();
			int parentIndex = (Integer) nodes[1].pop();
			
			WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(curIndex);
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index2 = v.getIndex();
				
				if (index2 != parentIndex) {
					ret[0].push(Integer.valueOf(index2));
					ret[1].push(Integer.valueOf(curIndex));
				}
			}
		}
		return ret;
	}

	
	
	public static void main(String args[]) {
		/*
		Graph g= new Graph();
		//g.getNSmallValues(data, 5);
		g.addNode(new Node("0","ATCGTGCAAATTT"));
		g.addNode(new Node("1","GTGCAAATTTGGG"));
		g.addNode(new Node("2","CAAATTTGGGCAT"));
		g.addNode(new Node("3","ATTTGGGCATCGGA"));
		g.addNode(new Node("4","CGGATTCAACCTG"));
		g.addNode(new Node("5","AACCTGAGT"));
		g.addNode(new Node("6","CCTGAGTTCGTCA"));
		g.addNode(new Node("7","TCGTCAAGTCAGT"));
		g.addNode(new Node("8","AAGTCAGTTCCG"));
		
		//g.addNode(new Node("TCCACTAGT"));
		//g.addNode(new Node("TCAGACGG"));
		//g.addNode(new Node("GACAATA"));
		
		//int a[][] = g.alignNodes();
		WeightedAdjacencyListGraph mst = g.genMST();
		g.get2CloseNodesFromGrand(mst, 2);
	*/
	}

}
