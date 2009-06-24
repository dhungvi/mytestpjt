import java.util.ArrayList;
import java.util.Properties;

import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

/*
 * This class is subclass of Graph. It assigns NewD2 to member variable d2.
 */
public class NewGraph extends Graph{
	
	public NewGraph(Properties p) {
		super(p);
		//graphNodes = new ArrayList<Node> ();
		d2 = new NewD2(p);
	}

	/**
	 * Overload the method of parent class
	 * Get two closest nodes which is on the left and on the right to every node 
	 * from the input minimum spanning tree, and store the data into an array.
	 * This function finds parent and children of current node, and calculates
	 * overlap distance between the node and others, select the two nodes with minimum 
	 * left and right overlap distance.
	 * 
	 * This function is different from the same-name method in its parent class in that 
	 * the second dimension of returned array has six fields instead of four in the parent
	 * method. It adds distance value into the array instead of overlap length.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @return an array which stores two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the tree. 
	 * 		1st-dimension: index of nodes in the graph 'graphNodes';
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the distance with + or -.
	 * 						the fourth is the index of node on the right, the fifth is the overlap length with + or -.
	 * 						the sixth is the distance with + or -.
	 * 						For the first and fourth one, if no node is found, the value is -1;
	 * 						For the second and fifth one, if no node is found, the value is 0;
	 * 						For the third and sixth one, if no node is found, the value is INT_MIN or INT_MAX.
	 */
	public int[][] get2CloseNodesFromMST(WeightedAdjacencyListGraph mst) {
		int nOfNodes = mst.getCardV();
		/*
		 * store the position of aligned nodes
		 * 1st-dimension: index of nodes in graph;
		 * 2rd-dimension: the first is the index of node on the left, 
		 * 					the second is the overlap length with + or -,
		 * 					the third is the index of node on the right, 
		 * 					the forth is the overlap length with + or -.
		 */
		int[][] alignedNodes = new int[nOfNodes][6];
		
		for (int i=0; i<nOfNodes; i++) {
			int leftNode = -1;
			int rightNode = -1;
			int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
			int minRight = INT_MAX;	//minimum right distance
			int overlapLeft = 0;
			int overlapRight = 0;
			
			WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(i);
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index = v.getIndex();
				int[] ovlDis = ((NewD2)d2).getOVLDistance(graphNodes.get(i).getNodeStr(), 
												graphNodes.get(index).getNodeStr());
				//int[] ovlDis = d2.getOVLDistance(graphNodes.get(i).getNodeStr(), 
				//		graphNodes.get(index).getNodeStr(),
				//		(int)ite.getWeight());
				//if (ovlDis[0] != 0) {	// there is overlap between them
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
			alignedNodes[i][0] = leftNode;	//index of node on the left
			alignedNodes[i][1] = overlapLeft;	//overlap length
			alignedNodes[i][2] = maxLeft;	//overlap distance
			
			alignedNodes[i][3] = rightNode;	//index of node on the right
			alignedNodes[i][4] = overlapRight;	//overlap length
			alignedNodes[i][5] = minRight;	//overlap distance
			
			//System.out.println("Get 6-tuple for node "+i);
		}
		return alignedNodes;
	}

	
	/**
	 * Overload the method of parent class
	 * Get two closest nodes which is on the left and on the right to the 'index' node
	 * from the input minimum spanning tree, and store the data into an array.
	 * To get the six-tuple for the node, the function calculates distance from 
	 * this node to: its parent, its children, its grandparent, its grand-children, its
	 * grand-grand-parent and its grand-grand-children. 
	 * 
	 * This function is different from the same-name method in its parent class in that 
	 * the second dimension of returned array has six fields instead of four in the parent
	 * method. It adds distance value into the array instead of overlap length.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param index The index of current node in the tree and graph.
	 * @param sixTuple The six tuple for this node with the index
	 * @return an array which stores two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the tree. 
	 * 		1st-dimension: index of nodes in the graph 'graphNodes';
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the distance with + or -.
	 * 						the fourth is the index of node on the right, the fifth is the overlap length with + or -.
	 * 						the sixth is the distance with + or -.
	 * 						For the first and fourth one, if no node is found, the value is -1;
	 * 						For others, if no node is found, the value is infinite.
	 */
	public int[] get2CloseNodesFromGrand(WeightedAdjacencyListGraph mst, int index, int[] sixTuple) {
		/*
		 * store the position of aligned nodes
		 * 1st-dimension: index of nodes in graph;
		 * 2rd-dimension: the first is the index of node on the left, 
		 * 					the second is the overlap length with + or -,
		 * 					the third is the index of node on the right, 
		 * 					the forth is the overlap length with + or -.
		 */
		int[] closeNode = new int[6];
		
		int leftNode = -1;
		int rightNode = -1;
		//int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
		//int minRight = INT_MAX;	//minimum right distance
		//int overlapLeft = 0;
		//int overlapRight = 0;
		int maxLeft = sixTuple[2];
		int minRight = sixTuple[5];
		int overlapLeft = sixTuple[1];
		int overlapRight = sixTuple[4];
		
		//put all the parents,grandparents,children,grand-children nodes into an arraylist. Do not 
		// include parents and children because they have been processed.
		ArrayList<Vertex> allNodes = new ArrayList<Vertex> ();
		WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(index);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			//allNodes.add(v);
			int index2 = v.getIndex();
			WeightedEdgeIterator ite2 = (WeightedEdgeIterator) mst.edgeIterator(index2);
			while (ite2.hasNext()) {
				Vertex v2 = (Vertex) ite2.next();
				if (v2.getIndex() != index) {
					allNodes.add(v2);

					int index3 = v2.getIndex();
					WeightedEdgeIterator ite3 = (WeightedEdgeIterator) mst.edgeIterator(index3);
					while (ite3.hasNext()) {
						Vertex v3 = (Vertex) ite3.next();
						if (v3.getIndex() != index2) {
							allNodes.add(v3);
						}
					}
				}
			}
		}
		
		//find two closest nodes
		String s1 = graphNodes.get(index).getNodeStr();
		while (!allNodes.isEmpty()) {
			Vertex tmpV = allNodes.remove(0);
			int tmpIndex = tmpV.getIndex();
			String s2 = graphNodes.get(tmpIndex).getNodeStr();
			int[] ovlDis = ((NewD2)d2).getOVLDistance(s1, s2);
			
			//if (ovlDis[0] != 0) {	// there is overlap between them
			if (ovlDis[1] != INT_MAX) {	// there is overlap between them
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
		closeNode[0] = leftNode;	//index of node on the left
		closeNode[1] = overlapLeft;	//overlap length
		closeNode[2] = maxLeft;	//overlap distance
		
		closeNode[3] = rightNode;	//index of node on the right
		closeNode[4] = overlapRight;	//overlap length
		closeNode[5] = minRight;	//overlap distance
		
		//System.out.println("Get 6-tuple for node "+index);
		return closeNode;
	}

	/**
	 * Recalculate 6-tuples for an assumed left end in order to remove all the false left ends.
	 * Specifically, for the assumed left end,start to calculate from fourth level until meeting one node which 
	 * makes six-tuple[][0] != -1, then return the six-tuple.
	 * If we fail to find any node after running out of all the nodes in the MST, we consider it a real left end
	 * and six-tuple[][0] == -1.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param index The index of current node in the tree and graph.
	 * @param sixTuple The six tuple for this node with the index
	 * @return an array which stores two closest nodes which is to the left and to the right if it is not left end;
	 * if it does, six-tuple[][0] = -1.
	 * 		1st-dimension: index of nodes in the graph 'graphNodes';
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the index of node on the right, the fourth is the overlap length with + or -.
	 * 						For the first and third one, if no node is found, the value is -1;
	 * 						For the second and fourth one, if no node is found, the value is infinite.
	 */
	public int[] checkLeftEndFromMST(WeightedAdjacencyListGraph mst, int index, int[] sixTuple) {
		for (int foundLevel=4; ; foundLevel++) {
			ArrayList<Vertex> allNodes = new ArrayList<Vertex> ();
			allNodes = getNodesFromMST(mst, index, foundLevel, 0, -1, allNodes);

			System.out.println("GetNodeFromMST for foundLevel=" + foundLevel
					+ "; index=" + index);
			System.out.println("\tnumber of nodes = " + allNodes.size());

			if (allNodes.size() == 0) {
				break;
			} else {
				int[] closeNode = findAdjacentNode(allNodes, index, sixTuple);
				if (closeNode[0] != -1) {
					System.out.println("findAdjacentNode for index=" + index);
					return closeNode;
				}
			}
		}
		
		System.out.println("Fail to findAdjacentNode for index=" + index);
		return null;
	}
	
	/*
	 * find the most adajcent node to the current node from allNodes.
	 * @param allNodes Store all the nodes which will be compared to the current node.
	 * @param index The index of current node.
	 * @sixTuple The sixTuple for the current node.
	 */
	private int[] findAdjacentNode(ArrayList<Vertex> allNodes, int index, int[] sixTuple) {
		/*
		 * store the position of aligned nodes
		 * 1st-dimension: index of nodes in graph;
		 * 2rd-dimension: the first is the index of node on the left, 
		 * 					the second is the overlap length with + or -,
		 * 					the third is the index of node on the right, 
		 * 					the forth is the overlap length with + or -.
		 */
		int[] closeNode = new int[6];
		
		int leftNode = -1;
		int rightNode = -1;
		int maxLeft = sixTuple[2];
		int minRight = sixTuple[5];
		int overlapLeft = sixTuple[1];
		int overlapRight = sixTuple[4];

		//find two closest nodes
		String s1 = graphNodes.get(index).getNodeStr();
		while (!allNodes.isEmpty()) {
			Vertex tmpV = allNodes.remove(0);
			int tmpIndex = tmpV.getIndex();
			String s2 = graphNodes.get(tmpIndex).getNodeStr();
			int[] ovlDis = ((NewD2)d2).getOVLDistance(s1, s2);
			
			if (ovlDis[1] != INT_MAX) {	// there is overlap between them
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
		closeNode[0] = leftNode;	//index of node on the left
		closeNode[1] = overlapLeft;	//overlap length
		closeNode[2] = maxLeft;	//overlap distance
		
		closeNode[3] = rightNode;	//index of node on the right
		closeNode[4] = overlapRight;	//overlap length
		closeNode[5] = minRight;	//overlap distance
		
		return closeNode;
	}

	/**
	 * Get all the nodes which is at the "foundLevel" from the node "curIndex".
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param curIndex The index of current node from which we want to find the nodes.
	 * @param foundLevel We will find the nodes which is at this level.
	 * @param curLevel Current level which is used by the recursion to record at which level we are.
	 * @param parentIndex The parent node of the current node. For the first calling, it is -1; it has value during recursion.
	 * @param nodes The arraylist which stores all the found nodes.
	 * @return an array which stores all the found nodes.
	 */
	private ArrayList<Vertex> getNodesFromMST(WeightedAdjacencyListGraph mst, int curIndex, int foundLevel, int curLevel, int parentIndex, ArrayList<Vertex> nodes) {
		
		ArrayList<Vertex> allNodes = nodes;
		WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(curIndex);
		
		if (foundLevel < curLevel) {
			System.out.println("getNodesFromMST has an error: foundLevel is less than curLevel!");
			return null;
		}
		
		if (foundLevel == curLevel+1){
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index2 = v.getIndex();
				if (parentIndex == -1) {
					allNodes.add(v);
				} else if (parentIndex != index2) {
					allNodes.add(v);
				}
			}
		} else { //foundLevel - curLevel >= 1
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index2 = v.getIndex();
				if ((parentIndex == -1) || (parentIndex != index2)){
					getNodesFromMST(mst, index2, foundLevel, curLevel+1, curIndex, nodes);
				} 
			}
		}
		
		return allNodes;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
