#ifndef ZY_Node_HH
#define ZY_Node_HH

class Node {
public:
	std::string sequence; //the bases of the est
	std::string name; //ID of the est,currently it's the starting position of the node.
	std::string comment; //comment to the est in the input est file

	inline Node(std::string n, std::string c, std::string s) {
		name = n;
		comment = c;
		sequence = s;
	}

	inline std::string getNodeStr() {
		return sequence;
	}

	/*
	 * get length of the node
	 */
	inline int getLen() {
		return sequence.length();
	}

	/*
	 * get ID of the node
	 */
	inline std::string getName() {
		return name;
	}

	/*
	 * get sequence of the node
	 */
	inline std::string getSeq() {
		return sequence;
	}

	/*
	 * get sequence of the node
	 */
	inline std::string getComment() {
		return comment;
	}
};

#endif
