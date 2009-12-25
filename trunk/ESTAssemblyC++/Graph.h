#ifndef ZY_Graph_HH
#define ZY_Graph_HH

#include <vector>
#include <string>
#include <limits.h>
#include <stack>
#include <map>
#include "Param.h"
#include "Node.h"
#include "OvlDistance.h"
#include "CalculatedOvlDistance.h"
#include "InclusionNodes.h"
#include "SixTuple.h"
#include "Prim.h"

class Graph {
public:
	int numOfLevels;
	std::vector<Node> graphNodes;
	OvlDistance ovl;
	CalculatedOvlDistance calDist;
	DefGraph mst; //minimum spanning tree generated from peace
	InclusionNodes* inc;

	Graph() {numOfLevels=0; inc=NULL;};
	Graph(InclusionNodes* in);
	std::vector<SixTuple*> get2CloseNodesFromMST();
	SixTuple get2CloseNodesFromGrand(int index, const SixTuple sixTuple);
	SixTuple checkLeftEndFromMST(int index, const SixTuple sixTuple);
	SixTuple checkRightEndFromMST(int index, const SixTuple sixTuple);
	inline void setMst(DefGraph& m) {
		mst = m;
	}
private:
	std::vector<SixTuple*> handleInclusion();
	SixTuple findAdjacentNode(std::stack<int> nodes, int index, const SixTuple sixTuple);
	std::vector<std::stack<int> > getNodesFromMST(std::vector<std::stack<int> > nodes);

public:
	inline void addNode(Node s) {
		graphNodes.push_back(s);
	}

	inline void removeNode(int index) {
		graphNodes.erase(graphNodes.begin() + index);
	}

	inline int getSizeofGraph() {
		return graphNodes.size();
	}

	/*
	 * get length of the node with index i
	 */
	inline int getLenOfNode(int i) {
		return graphNodes.at(i).getLen();
	}

	/*
	 * get ID of the node with index i
	 */
	inline std::string getNameOfNode(int i) {
		return graphNodes.at(i).getName();
	}

	/*
	 * get sequence of the node with index i
	 */
	inline std::string getSeqOfNode(int i) {
		return graphNodes.at(i).getSeq();
	}

	/*
	 * get comment of the node with index i
	 */
	inline std::string getCommentOfNode(int i) {
		return graphNodes.at(i).getComment();
	}


};

#endif
