import com.mhhe.clrs2e.Prim;
import com.mhhe.clrs2e.Vertex;
import com.mhhe.clrs2e.WeightedAdjacencyListGraph;
import com.mhhe.clrs2e.WeightedEdgeIterator;

public class LinearTree {
	static final int INT_MAX = Integer.MAX_VALUE;
	
	public static WeightedAdjacencyListGraph constructLinearTree(int nOfNodes, int[][] g) {
		// Make a directed graph.
		WeightedAdjacencyListGraph dGraph = new WeightedAdjacencyListGraph(nOfNodes, true);
		WeightedAdjacencyListGraph lst = new WeightedAdjacencyListGraph(nOfNodes, true);
		for (int i=0; i<nOfNodes; i++) {
			dGraph.addVertex(i, Integer.toString(i));
			lst.addVertex(i, Integer.toString(i));
		}
		for (int j=0; j<g.length; j++) {
			if (g[j][3] != 0) {	//there is an edge between the nodes
				double weight = g[j][2] + (1.0/g[j][3]);
				dGraph.addEdge(g[j][0], g[j][1], weight);
			}
		}
		
		int curIdx = 0;
		while (true) {
			int nextIdx = -1;
			double minWeight = INT_MAX;
			WeightedEdgeIterator ite = (WeightedEdgeIterator) dGraph.edgeIterator(curIdx);
			boolean flag = false;
			
			while (ite.hasNext()) {
				Vertex v = (Vertex) ite.next();
				int index2 = v.getIndex();
				if (ite.getWeight() < minWeight) {
					minWeight = ite.getWeight();
					nextIdx = index2;
					flag = true;
				}
			}
			
			if (flag) {
				lst.addEdge(curIdx, nextIdx, minWeight);
				curIdx = nextIdx;
			} else {
				break;
			}
		}
		
		return lst;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
