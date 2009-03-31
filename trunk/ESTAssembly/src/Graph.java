import java.util.ArrayList;
import java.util.Properties;

import com.mhhe.clrs2e.Prim;
import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

public class Graph {
	protected final int INT_MAX = 2147483647;
	protected final int INT_MIN = -2147483647;
	private final int num = 5;	//specify the number of the closest nodes the program finds for one node;
								//function 'getNClosestNodes' uses 'num' as parameter.

	ArrayList<Node> graphNodes;
	D2 d2;
	
	public Graph(Properties p) {
		graphNodes = new ArrayList<Node> ();
		d2 = new D2(p);
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
	 * generate a minimum spanning tree from the graph 'graphNodes'
	 */
	public WeightedAdjacencyListGraph genMST() {
		/*
		 * 	1st-dimension: index of nodes in graph;
		 * 	2nd-dimension: index of n closest nodes;
		 * 	3rd-dimension: the first is the value of distance, the second is the index of node. It's sorted 
		 * 	in ascending order according to the first value.
		 */
		int[][][] nodes = getNClosestNodes();
		int nOfNodes = nodes.length;
		// Make a undirected graph.
		WeightedAdjacencyListGraph udGraph =
		    new WeightedAdjacencyListGraph(nOfNodes, false);
		for (int i=0; i<nOfNodes; i++) {	//i is the index of the node in graph
			udGraph.addVertex(i, Integer.toString(i));
		}
		for (int i=0; i<nOfNodes; i++) {
			for (int j=0; j<nodes[i].length; j++) {
				if (nodes[i][j][0] != INT_MAX) {	//there is an edge between the nodes
					udGraph.addEdge(i, nodes[i][j][1], nodes[i][j][0]);
				}
				
			}
		}
		WeightedAdjacencyListGraph mst = (new Prim()).computeMST(udGraph);
		return mst;
	}
	
	/**
	 * Get two closest nodes which is on the left and on the right to every node 
	 * from the input minimum spanning tree, and store the data into an array.
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @return an array which stores two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the tree. 
	 * 		1st-dimension: index of nodes in the graph 'graphNodes';
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the index of node on the right, the fourth is the overlap length with + or -.
	 * 						For the first and third one, if no node is found, the value is -1;
	 * 						For the second and fourth one, if no node is found, the value is infinite.
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
		int[][] alignedNodes = new int[nOfNodes][4];
		
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
				int[] ovlDis = d2.getOVLDistance(graphNodes.get(i).getNodeStr(), 
												graphNodes.get(index).getNodeStr(), 
												(int)ite.getWeight());
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
			alignedNodes[i][2] = rightNode;	//index of node on the right
			alignedNodes[i][3] = overlapRight;	//overlap length
		}
		return alignedNodes;
	}
	
	/**
	 * Get two closest nodes which is on the left and on the right to the 'index' node
	 * from the input minimum spanning tree, and store the data into an array.
	 * To get the four-element set for the node, the function calculates distance from 
	 * this node to: its parent, its children, its grandparent, its grand-children, its
	 * grand-grand-parent and its grand-grand-children. 
	 * 
	 * @param mst a Minimum Spanning Tree.
	 * @param index The index of current node in the tree and graph.
	 * @return an array which stores two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the tree. 
	 * 		1st-dimension: index of nodes in the graph 'graphNodes';
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the index of node on the right, the fourth is the overlap length with + or -.
	 * 						For the first and third one, if no node is found, the value is -1;
	 * 						For the second and fourth one, if no node is found, the value is infinite.
	 */
	public int[] get2CloseNodesFromGrand(WeightedAdjacencyListGraph mst, int index) {
		/*
		 * store the position of aligned nodes
		 * 1st-dimension: index of nodes in graph;
		 * 2rd-dimension: the first is the index of node on the left, 
		 * 					the second is the overlap length with + or -,
		 * 					the third is the index of node on the right, 
		 * 					the forth is the overlap length with + or -.
		 */
		int[] closeNodes = new int[4];
		
		int leftNode = -1;
		int rightNode = -1;
		int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
		int minRight = INT_MAX;	//minimum right distance
		int overlapLeft = 0;
		int overlapRight = 0;
		
		//put all the parents,grandparents,children,grand-children nodes into an arraylist.
		ArrayList<Vertex> allNodes = new ArrayList<Vertex> ();
		ArrayList<Vertex> debugNodesL1 = new ArrayList<Vertex> ();
		ArrayList<Vertex> debugNodesL2 = new ArrayList<Vertex> ();
		ArrayList<Vertex> debugNodesL3 = new ArrayList<Vertex> ();
		ArrayList<Vertex> debugNodesL4 = new ArrayList<Vertex> ();
		WeightedEdgeIterator ite = (WeightedEdgeIterator) mst.edgeIterator(index);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			allNodes.add(v);
			//debugNodesL1.add(v);
			int index2 = v.getIndex();
			WeightedEdgeIterator ite2 = (WeightedEdgeIterator) mst.edgeIterator(index2);
			while (ite2.hasNext()) {
				Vertex v2 = (Vertex) ite2.next();
				if (v2.getIndex() != index) {
					allNodes.add(v2);
					//debugNodesL2.add(v2);

					int index3 = v2.getIndex();
					WeightedEdgeIterator ite3 = (WeightedEdgeIterator) mst.edgeIterator(index3);
					while (ite3.hasNext()) {
						Vertex v3 = (Vertex) ite3.next();
						if (v3.getIndex() != index2) {
							allNodes.add(v3);
							//debugNodesL3.add(v3);

							int index4 = v3.getIndex();
							WeightedEdgeIterator ite4 = (WeightedEdgeIterator) mst.edgeIterator(index4);
							while (ite4.hasNext()) {
								Vertex v4 = (Vertex) ite4.next();
								if (v4.getIndex() != index3) {
									//allNodes.add(v4);
									//debugNodesL4.add(v4);
								}
							}
						}
					}
				}
			}
		}
		// debug: print all the nodes in different levels
		/*System.out.println("Node " + index + "\n1st Level: ");
		while (!debugNodesL1.isEmpty()) {
			System.out.print(debugNodesL1.remove(0).getIndex()+"\t");
		}
		System.out.println("\n2nd Level: ");
		while (!debugNodesL2.isEmpty()) {
			System.out.print(debugNodesL2.remove(0).getIndex()+"\t");
		}
		System.out.println("\n3rd Level: ");
		while (!debugNodesL3.isEmpty()) {
			System.out.print(debugNodesL3.remove(0).getIndex()+"\t");
		}
		System.out.println("\n4th Level: ");
		while (!debugNodesL4.isEmpty()) {
			System.out.print(debugNodesL4.remove(0).getIndex()+"\t");
		}
		System.out.println("\n");*/
		
		//calculate two closest nodes
		while (!allNodes.isEmpty()) {
			Vertex tmpV = allNodes.remove(0);
			int tmpIndex = tmpV.getIndex();
			String s1 = graphNodes.get(index).getNodeStr();
			String s2 = graphNodes.get(tmpIndex).getNodeStr();
			int d2Dis = d2.getD2Distance(s1, s2);
			int[] ovlDis = d2.getOVLDistance(s1, s2, d2Dis);
			
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
		closeNodes[0] = leftNode;	//index of node on the left
		closeNodes[1] = overlapLeft;	//overlap length
		closeNodes[2] = rightNode;	//index of node on the right
		closeNodes[3] = overlapRight;	//overlap length
		
		return closeNodes;
	}
	
	/**
	 * Get two closest node on the left and closet node on the right for the input node.
	 * This function computes the distance from the 'index' node to all the other nodes in the graph.
	 * 
	 * @param index the index of the node in the graph 'graphNodes'.
	 * @return an array which store two closest nodes which is to the left and to the right
	 * respectively for the input node. 
	 * 		the first is the index of node on the left, the second is the overlap length with + or -.
	 * 		the third is the index of node on the right, the fourth is the overlap length with + or -.
	 * 		For the first and third one, if no node is found, the value is -1;
	 * 		For the second and fourth one, if no node is found, the value is infinite.
	 */

	public int[] get2CloseNodes(int index) {
		int lenOfGraph = graphNodes.size();
		/*
		 * store the two closest nodes to the node with index 'index'
		 * 		the first is the index of node on the left, 
		 * 		the second is the overlap length with + or -,
		 * 		the third is the index of node on the right, 
		 * 		the forth is the overlap length with + or -.
		 */
		int[] closeNodes = new int[4];
		
		int disToAllNodes[] = new int[lenOfGraph];
		Node currentNode = graphNodes.get(index);
		for (int i=0; i<lenOfGraph; i++) {
			if (i != index) {
				Node tmpNode = graphNodes.get(i);
				disToAllNodes[i] = d2.getD2Distance(currentNode.getNodeStr(), tmpNode.getNodeStr());
			} else {
				disToAllNodes[i] = INT_MAX;
			}
		}
		
		int[][] nodes = getNSmallValues(disToAllNodes);

		int leftNode = -1;
		int rightNode = -1;
		int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
		int minRight = INT_MAX;	//minimum right distance
		int overlapLeft = 0;
		int overlapRight = 0;

		for (int j=0; j<nodes.length; j++) {
			int[] ovlDis = d2.getOVLDistance(currentNode.getNodeStr(), 
					graphNodes.get(nodes[j][1]).getNodeStr(), 
					nodes[j][0]);
			//if (ovlDis[0] != 0) {	// there is overlap between them
			if (ovlDis[1] != INT_MAX) {	// there is overlap between them
				if (ovlDis[0] < 0) {
					if (ovlDis[1] > maxLeft){
						maxLeft = ovlDis[1];
						overlapLeft = ovlDis[0];
						leftNode = nodes[j][1];
					} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapLeft)) {
							overlapLeft = ovlDis[0];
							leftNode = nodes[j][1];
						}
					}
				}
				if (ovlDis[0] > 0) {
					if (ovlDis[1] < minRight) {
						minRight = ovlDis[1];
						overlapRight = ovlDis[0];
						rightNode = nodes[j][1];
					} else if (ovlDis[1] == minRight) {	//if they are equal, find that one with maximal overlap
						if (Math.abs(ovlDis[0]) > Math.abs(overlapRight)) {
							overlapRight = ovlDis[0];
							rightNode = nodes[j][1];
						}
					}
				}
			}
			closeNodes[0] = leftNode;	//index of node on the left
			closeNodes[1] = overlapLeft;	//overlap length
			closeNodes[2] = rightNode;	//index of node on the right
			closeNodes[3] = overlapRight;	//overlap length
		}
		return closeNodes;
	}

	/**
	 * Get closest node on the left and closet node on the right for every node in the graph,
	 * and store the data into an array.
	 * 
	 * @param 
	 * @return an array which store two closest nodes which is to the left and to the right 
	 * respectively for all the nodes in the graph. 
	 * 		1st-dimension: index of nodes in graph;
	 * 		2rd-dimension: the first is the index of node on the left, the second is the overlap length with + or -.
	 * 						the third is the index of node on the right, the fourth is the overlap length with + or -.
	 * 						For the first and third one, if no node is found, the value is -1;
	 * 						For the second and fourth one, if no node is found, the value is infinite.
	 */
	public int[][] alignNodes() {
		int lenOfGraph = graphNodes.size();
		/*
		 * store the position of aligned nodes
		 * 1st-dimension: index of nodes in graph;
		 * 2rd-dimension: the first is the index of node on the left, 
		 * 					the second is the overlap length with + or -,
		 * 					the third is the index of node on the right, 
		 * 					the forth is the overlap length with + or -.
		 */
		int[][] alignedNodes = new int[lenOfGraph][4];
		int[][][] nodes = getNClosestNodes();
		
		for (int i=0; i<nodes.length; i++) {
			int leftNode = -1;
			int rightNode = -1;
			int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
			int minRight = INT_MAX;	//minimum right distance
			int overlapLeft = 0;
			int overlapRight = 0;
			
			for (int j=0; j<num; j++) {
				int[] ovlDis = d2.getOVLDistance(graphNodes.get(i).getNodeStr(), 
												graphNodes.get(nodes[i][j][1]).getNodeStr(), 
												nodes[i][j][0]);
				//if (ovlDis[0] != 0) {	// there is overlap between them
				if (ovlDis[1] != INT_MAX) {	// there is overlap between them
					if (ovlDis[0] < 0) {
						if (ovlDis[1] > maxLeft){
							maxLeft = ovlDis[1];
							overlapLeft = ovlDis[0];
							leftNode = nodes[i][j][1];
						} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
							if (Math.abs(ovlDis[0]) > Math.abs(overlapLeft)) {
								overlapLeft = ovlDis[0];
								leftNode = nodes[i][j][1];
							}
						}
					}
					if (ovlDis[0] > 0) {
						if (ovlDis[1] < minRight) {
							minRight = ovlDis[1];
							overlapRight = ovlDis[0];
							rightNode = nodes[i][j][1];
						} else if (ovlDis[1] == minRight) {	//if they are equal, find that one with maximal overlap
							if (Math.abs(ovlDis[0]) > Math.abs(overlapRight)) {
								overlapRight = ovlDis[0];
								rightNode = nodes[i][j][1];
							}
						}
					}
				}
			}
			alignedNodes[i][0] = leftNode;	//index of node on the left
			alignedNodes[i][1] = overlapLeft;	//overlap length
			alignedNodes[i][2] = rightNode;	//index of node on the right
			alignedNodes[i][3] = overlapRight;	//overlap length
		}
		return alignedNodes;
	}
	
	/**
	 * Get 'num' closest nodes for all the nodes in the graph.
	 * Note: some nodes may include same elements in their set of NClosestNodes.
	 * 
	 * @param num int.
	 * @return an array which store 'num' closest nodes for all the nodes in the graph. 
	 * 	1st-dimension: index of nodes in graph;
	 * 	2nd-dimension: index of n closest nodes;
	 * 	3rd-dimension: the first is the value of distance, the second is the index of node. It's sorted 
	 * 	in ascending order according to the first value. 
	 */
	private int[][][] getNClosestNodes() {
		int lenOfGraph = graphNodes.size();

		/*
		 * 1st-dimension: index of nodes in graph;
		 * 2nd-dimension: index of n closest nodes;
		 * 3rd-dimension: the first is the value of distance, the second is the index of node. It's sorted 
		 * in ascending order according to the first value.
		 */
		int nClosetNodes[][][] = new int[lenOfGraph][num][2];	
		int disToAllNodes[] = new int[lenOfGraph];
		
		for (int i=0; i<lenOfGraph; i++) {
			Node currentNode = graphNodes.get(i);
			
			for (int j=0; j<lenOfGraph; j++) {
				if (j != i) {
					Node tmpNode = graphNodes.get(j);
					disToAllNodes[j] = d2.getD2Distance(currentNode.getNodeStr(), tmpNode.getNodeStr());
				} else {
					disToAllNodes[j] = INT_MAX;
				}
			}
			
			nClosetNodes[i] = getNSmallValues(disToAllNodes);
		}
		return nClosetNodes;
	}
	
	/**
	 * Get the 'num' smallest values in the array 'data'. 
	 * 
	 * @param data int[], num int.
	 * @return an array with the 'num' smallest values in 'data'. 
	 * the array is sorted in ascending order according to values[][0].
	 * 		values[][0] is the value in 'data';
	 * 		values[][1] is the index of this value in 'data'.
	 */
	private int[][] getNSmallValues(int[] data) {
		int values[][] = new int[num][2];
		int index = 0;
		for (int i=0; i<num; i++) {
			int min = INT_MAX;
			for (int j=0; j<data.length; j++) {
				if (data[j] < min) {
					min = data[j];
					index = j;
				}
			}
			values[i][0] = min;
			values[i][1] = index;
			data[index] = INT_MAX;
		}
		return values;
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
