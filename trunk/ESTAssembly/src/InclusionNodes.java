import java.util.ArrayList;
import java.util.Iterator;
import java.util.Stack;
import java.util.TreeSet;

import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

/*
 * record all the nodes which are included in other nodes.
 * We get rid of them from the input MST and put them into the class. Then 
 * we use them again during reconstruction the consensus.
 */
public class InclusionNodes {
	TreeSet<InclusionNode> nodes; //Used in handleInclusion in Graph.java. The compared member variable has to be unique.
	ArrayList<PairNode> nodes2; //Used in get2CloseNodesFromGrand and findAdjacentNode in Graph.java.
	WeightedAdjacencyListGraph directedGraph;
	TreeSet<PNode> pNodes; //store all the parent nodes in nodes and nodes2.

	
	public InclusionNodes() {
		nodes = new TreeSet<InclusionNode>();
		nodes2 = new ArrayList<PairNode>();
		pNodes = new TreeSet<PNode>();
		directedGraph = null;
	}
	
	// chd is included in parent.
	public void addNode(int chd, int parent) {
		nodes.add(new InclusionNode(chd, parent));
	}
	
	public void addNode2(int chd, int parent) {
		nodes2.add(new PairNode(chd, parent));
	}

	public int getSize() {
		return nodes.size();
	}
	
	/*
	 * check if the idx is in the "nodes".
	 */
	public boolean containInclusionNode(int idx) {
		return nodes.contains(new InclusionNode(idx,0));
	}

	/*
	 * If pIdx is in nodes or nodes2, return the corresponding idxChd;
	 * else return -1.
	 * Called by addInclusionNodes in ESTAssembly.java.
	 */
	public int[] containPNode(int pIdx, int numOfTotalNodes) {
		if (directedGraph == null) {
			makeDirectedGraph(numOfTotalNodes);
		}
		
		if (!pNodes.contains(new PNode(pIdx))) {
			return null;
		}
		
		Stack<Integer> allNodes = new Stack<Integer> ();
//		Stack<Integer> partNodes = new Stack<Integer> ();
//		partNodes.push(Integer.valueOf(pIdx));
//		
//		while (!partNodes.empty()) {
//			Stack<Integer> tmpStack = new Stack<Integer> ();
//			while (!partNodes.empty()) {
//				int tmpIndex = partNodes.pop();
//				WeightedEdgeIterator ite = (WeightedEdgeIterator) directedGraph.edgeIterator(tmpIndex);
//				while (ite.hasNext()) {
//					Vertex v = (Vertex) ite.next();
//					int tIndex = v.getIndex();
//					tmpStack.push(Integer.valueOf(tIndex));
//					allNodes.push(Integer.valueOf(tIndex));
//				}
//			}
//			partNodes = tmpStack;
//		}
		
		WeightedEdgeIterator ite = (WeightedEdgeIterator) directedGraph.edgeIterator(pIdx);
		while (ite.hasNext()) {
			Vertex v = (Vertex) ite.next();
			int tIndex = v.getIndex();
			allNodes.push(Integer.valueOf(tIndex));
		}

		int[] values = new int[allNodes.size()];
		for (int i=0; i<allNodes.size(); i++) {
			values[i] = allNodes.get(i);
		}

		return values;		
	}
	
	private void makeDirectedGraph(int numOfTotalNodes) {
		// Make a directed graph.
		int nOfNodes = numOfTotalNodes;
		directedGraph = new WeightedAdjacencyListGraph(nOfNodes, true);
		
		for (int i=0; i<nOfNodes; i++) {
			directedGraph.addVertex(i, Integer.toString(i));
		}

		Iterator<InclusionNode> ite = nodes.iterator();
		while (ite.hasNext()) {
			InclusionNode n1 = ite.next();
			directedGraph.addEdge(n1.idxP, n1.idxChd, 1);
			pNodes.add(new PNode(n1.idxP));
		}
		
		Iterator<PairNode> ite2 = nodes2.iterator();
		while (ite2.hasNext()) {
			PairNode n1 = ite2.next();
			directedGraph.addEdge(n1.idxP, n1.idxChd, 1);
			pNodes.add(new PNode(n1.idxP));
		}
	}
	
	public void printAllNodes() {
		System.out.println("childIndex\tparentIndex");
		Iterator<InclusionNode> ite = nodes.iterator();
		while (ite.hasNext()) {
			InclusionNode n1 = ite.next();
			System.out.println(n1.idxChd+"\t"+n1.idxP);
		}
		Iterator<PairNode> ite2 = nodes2.iterator();
		while (ite2.hasNext()) {
			PairNode n1 = ite2.next();
			System.out.println(n1.idxChd+"\t"+n1.idxP);
		}
	}
	
	public static void main(String args[]) {
		InclusionNodes in = new InclusionNodes();
		in.addNode(4,3);
		in.addNode(5,3);
		in.addNode(6,4);
		in.addNode(7,4);
		in.addNode2(5,1);
		in.addNode2(7,1);
		
		//System.out.println(in.containInclusionNode(7));
		int[] a = in.containPNode(3, 8);
		System.out.println(a.length);
		//System.out.println(in.containInclusionNode(3));
	}
}

class InclusionNode implements Comparable<InclusionNode>{
	int idxChd; //index of the node.
	int idxP;  //index of the node which include curNode.
	
	InclusionNode(int c, int p) {
		idxChd = c;
		idxP = p;
	}
	
	public int compareTo(InclusionNode other) {
		//Returns 0 if the argument is equal to this; 			
		//a value less than 0 if the argument is greater than this; 
		//and a value greater than 0 if the argument is less than this. 
		if (this.idxChd == other.idxChd) {
			return 0;
		} else if (this.idxChd > other.idxChd) {
			return 1;
		} else {
			return -1;
		}
	}
}

class PairNode { 
	int idxP; //index of the node.
	int idxChd;  //index of the node which is included by pNode.
	
	PairNode(int c, int p) {
		idxP = p;
		idxChd = c;
	}
}

class PNode implements Comparable<PNode>{
	int idxP;  //index of the node which include curNode.
	
	PNode(int p) {
		idxP = p;
	}
	
	public int compareTo(PNode other) {
		//Returns 0 if the argument is equal to this; 			
		//a value less than 0 if the argument is greater than this; 
		//and a value greater than 0 if the argument is less than this. 
		if (this.idxP == other.idxP) {
			return 0;
		} else if (this.idxP > other.idxP) {
			return 1;
		} else {
			return -1;
		}
	}
}