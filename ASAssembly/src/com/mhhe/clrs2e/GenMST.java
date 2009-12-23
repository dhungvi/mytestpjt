package com.mhhe.clrs2e;

public class GenMST {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		WeightedAdjacencyListGraph wm = new WeightedAdjacencyListGraph(5, false);
		wm.addVertex(0, "node 1");
		wm.addVertex(1, "node 2");
		wm.addVertex(2, "node 3");
		wm.addVertex(3, "node 4");
		wm.addVertex(4, "node 5");
		wm.addEdge(0, 1, 2);
		wm.addEdge(0, 2, 4);
		wm.addEdge(0, 3, 5);
		wm.addEdge(0, 4, 6);
		wm.addEdge(1, 2, 7);
		wm.addEdge(4, 2, 3);
		
		Prim p = new Prim();
		WeightedAdjacencyListGraph p2 = p.computeMST(wm);
		
		
	}

}
