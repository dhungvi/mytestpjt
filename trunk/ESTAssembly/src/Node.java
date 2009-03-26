
public class Node {
	String sequence;	//the bases of the est
	String name;	//ID of the est,currently it's the starting position of the node.
	
	public Node(String n, String s) {
		name = n;
		sequence = s;
	}
	
	public String getNodeStr() {
		return sequence;
	}
	
	/*
	 * get length of the node
	 */
	public int getLen() {
		return sequence.length();
	}
	
	/*
	 * get ID of the node
	 */
	public String getName() {
		return name;
	}
}
