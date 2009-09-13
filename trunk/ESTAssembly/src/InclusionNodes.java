import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

/*
 * record all the nodes which are included in other nodes.
 * We get rid of them from the input MST and put them into the class. Then 
 * we use them again during reconstruction the consensus.
 */
public class InclusionNodes {
	TreeSet<InclusionNode> nodes; 
	//TreeSet<PNode> prtNodes; 
	
	
	public InclusionNodes() {
		nodes = new TreeSet<InclusionNode>();
		//prtNodes = new TreeSet<PNode>();
	}
	
	// chd is included in parent.
	public void addNode(int chd, int parent) {
		nodes.add(new InclusionNode(chd, parent));
		//prtNodes.add(new PNode(parent, chd));
	}
	
	public int getSize() {
		return nodes.size();
	}
	
	/*
	 * check if the idx is in the nodes.
	 */
	public boolean containInclusionNode(int idx) {
		return nodes.contains(new InclusionNode(idx,0));
	}

	/*
	 * If pIdx is in prtNodes, return the corresponding idxChd;
	 * else return -1.
	 */
	public int[] containPNode(int pIdx) {
		ArrayList<Integer> ret = new ArrayList<Integer> ();
		Iterator<InclusionNode> ite = nodes.iterator();
		while (ite.hasNext()) {
			InclusionNode n1 = ite.next();
			if (n1.idxP == pIdx) {
				ret.add(Integer.valueOf(n1.idxChd));
			} 
		}
		
		int[] values = new int[ret.size()];
		for (int i=0; i<ret.size(); i++) {
			values[i] = ret.get(i);
		}
		return values;		
	}
	
	public void printAllNodes() {
		System.out.println("childIndex\tparentIndex");
		Iterator<InclusionNode> ite = nodes.iterator();
		while (ite.hasNext()) {
			InclusionNode n1 = ite.next();
			System.out.println(n1.idxChd+"\t"+n1.idxP);
		}
	}
	
	public static void main(String args[]) {
		InclusionNodes in = new InclusionNodes();
		in.addNode(3,5);
		in.addNode(5,7);
		in.addNode(7,7);
		in.addNode(6,7);
		
		System.out.println(in.containInclusionNode(7));
		int[] a = in.containPNode(7);
		if (a.length == 0) {
			System.out.println("false");
		}
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

/*
 * similar to InclusionNode except that the comparable variable is pNodes instead of curNode.
 */
/*class PNode implements Comparable<PNode>{ 
	int idxP; //index of the node.
	int idxChd;  //index of the node which is included by pNode.
	
	PNode(int p, int c) {
		idxP = p;
		idxChd = c;
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
}*/
