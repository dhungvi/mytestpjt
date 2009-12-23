
public class Node {
	String sequence;	//the bases of the est
	String name;	//ID of the est,currently it's the starting position of the node.
	String comment; //comment to the est in the input est file
	
	public Node(String n, String c, String s) {
		name = n;
		comment = c;
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

	/*
	 * get sequence of the node
	 */
	public String getSeq() {
		return sequence;
	}

	/*
	 * get sequence of the node
	 */
	public String getComment() {
		return comment;
	}
}
