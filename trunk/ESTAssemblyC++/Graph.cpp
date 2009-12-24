#include <math.h>
#include <iostream>
#include "Graph.h"

using namespace std;

Graph::Graph(InclusionNodes* in) {
	numOfLevels = NUMOFLEVELS;
	inc = in;
}

vector<SixTuple*> Graph::handleInclusion() {
	vector<SixTuple*> nodes;
	int nOfNodes = mst.size();

	for (int i=0; i<nOfNodes; i++) {
		bool flag = true;
		for (EdgeIterator it=mst[i].begin(); it!=mst[i].end(); it++) {
			int index = (*it).node;

			string curSeq = graphNodes[i].getNodeStr();
			string comSeq = graphNodes[index].getNodeStr();
			vector<int> ovlDis = calDist.searchDistance(i, index);

			if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) {
				ovlDis = ovl.getOVLDistance(curSeq, comSeq);
				//add to CalculatedOvlDistance
				calDist.addDistance(i, index, ovlDis[1], ovlDis[0]);
			}

			if (curSeq.length() <= comSeq.length()) {
				if (ovlDis[1] == INT_MIN) { //has inclusion
					inc->addNode1(i, index); //get rid of i and put it into inclusion list.
					flag = false;
					break;
				}
			}

		}
		if (flag) {
			nodes.push_back(new SixTuple(i));
		}
	}
	return nodes;
}

/**
 * Get two closest nodes which is on the left and on the right to every node
 * from the input minimum spanning tree, and store the data into an array list.
 * This function removes inclusion nodes from the tree, and then finds parent and
 * children for all the left nodes, and calculates overlap distance between the node and
 * others, select the two nodes with minimum left and right overlap distance.
 *
 *
 * @param mst a Minimum Spanning Tree.
 * @return an array list of SixTuple which stores two closest nodes which is to the left and to the right
 * respectively for all the nodes in the tree.
 */
vector<SixTuple*> Graph::get2CloseNodesFromMST() {
	vector<SixTuple*> alignedNodes = handleInclusion();
	int nOfNodes = alignedNodes.size();

	for (int i=0; i<nOfNodes; i++) {
		int curIdx = alignedNodes[i]->curNode;
		if (inc->containInclusionNode(curIdx)) continue;
		int leftNode = -1;
		int rightNode = -1;
		int maxLeft = INT_MIN;	//maximum left distance because left distance is negative
		int minRight = INT_MAX;	//minimum right distance
		int overlapLeft = 0;
		int overlapRight = 0;

		for (EdgeIterator it=mst[curIdx].begin(); it!=mst[curIdx].end(); it++) {
			int index = (*it).node;
			if (inc->containInclusionNode(index)) continue;

			vector<int> ovlDis = calDist.searchDistance(curIdx, index);
			if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) {
				ovlDis = ovl.getOVLDistance(graphNodes[curIdx].getNodeStr(),
						graphNodes[index].getNodeStr());
				//add to CalculatedOvlDistance
				calDist.addDistance(curIdx, index, ovlDis[1], ovlDis[0]);
			}

			if (ovlDis[1] != INT_MAX) {	// there is overlap between them
				if (ovlDis[0] < 0) {
					if (ovlDis[1] > maxLeft){
						maxLeft = ovlDis[1];
						overlapLeft = ovlDis[0];
						leftNode = index;
					} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
						if (abs(ovlDis[0]) > abs(overlapLeft)) {
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
						if (abs(ovlDis[0]) > abs(overlapRight)) {
							overlapRight = ovlDis[0];
							rightNode = index;
						}
					}
				}
			}
		}
		//leftNode, index of node on the left
		//overlapLeft, overlap length
		//maxLeft, overlap distance
		//rightNode, index of node on the right
		//overlapRight, overlap length
		//minRight, overlap distance
		alignedNodes[i]->setSixTuple(leftNode, overlapLeft, maxLeft, rightNode, overlapRight, minRight, false);
	}
	return alignedNodes;
}

/**
 * Get two closest nodes which is on the left and on the right to the 'index' node
 * from the input minimum spanning tree, and store the data into an array.
 * To get the six-tuple for the node, the function calculates distance from
 * this node to all the other nodes which are less than or equal to three levels from it.
 *
 *
 * @param mst a Minimum Spanning Tree.
 * @param index The index of current node in the tree and graph.
 * @param sixTuple The SixTuple for this node with the index
 * @return SixTuple for the current node which has the input "index".
 *						the first is the index of node on the left, the second is the overlap length with + or -.
 * 						the third is the distance with + or -.
 * 						the fourth is the index of node on the right, the fifth is the overlap length with + or -.
 * 						the sixth is the distance with + or -.
 * 						For the first and fourth one, if no node is found, the value is -1;
 * 						For the second and fifth one, if no node is found, the value is 0;
 * 						For the third and sixth one, if no node is found, the value is INT_MIN or INT_MAX.
 */
SixTuple Graph::get2CloseNodesFromGrand(int index, const SixTuple sixTuple) {
	SixTuple closeNode;

	int leftNode = -1;
	int rightNode = -1;
	int maxLeft = sixTuple.lDis;
	int minRight = sixTuple.rDis;
	int overlapLeft = sixTuple.lOvlLen;
	int overlapRight = sixTuple.rOvlLen;

	//put all the nodes within three levels from the current node into the stack "allNodes". Do not
	// include parents and children because they have been processed.
	stack<int> allNodes;
	stack<int> partNodes;
	partNodes.push(index);

	for (int i=0; i<3; i++) {
		stack<int> tmpStack;

		while (!partNodes.empty()) {
			int tmpIndex = partNodes.top();
			partNodes.pop();
			for (EdgeIterator it=mst[tmpIndex].begin(); it!=mst[tmpIndex].end(); it++) {
				int tIndex = (*it).node;
				tmpStack.push(tIndex);

				if ((i != 0) && (tIndex != index)){ //Do not include parents and children because they have been processed.
					allNodes.push(tIndex);
				}
			}
		}
		partNodes = tmpStack;

	}

	//find two closest nodes
	string s1 = graphNodes[index].getNodeStr();
	while (!allNodes.empty()) {
		int tmpIndex = allNodes.top();
		allNodes.pop();
		if (inc->containInclusionNode(tmpIndex)) continue;
		string s2 = graphNodes[tmpIndex].getNodeStr();

		vector<int> ovlDis = calDist.searchDistance(index, tmpIndex);
		if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) {
			ovlDis = ovl.getOVLDistance(s1, s2);

			//add to CalculatedOvlDistance
			calDist.addDistance(index, tmpIndex, ovlDis[1], ovlDis[0]);
		}

		if (ovlDis[1] == INT_MIN) {	// there is inclusion between them
			if (s1.length() >= s2.length()) {
				inc->addNode2(tmpIndex, index);
			} else {
				inc->addNode2(index, tmpIndex);
			}
		} else if (ovlDis[1] != INT_MAX) {	// there is overlap between them
			if (ovlDis[0] < 0) {
				if (ovlDis[1] > maxLeft){
					maxLeft = ovlDis[1];
					overlapLeft = ovlDis[0];
					leftNode = tmpIndex;
				} else if (ovlDis[1] == maxLeft) {	//if they are equal, find that one with maximal overlap
					if (abs(ovlDis[0]) > abs(overlapLeft)) {
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
					if (abs(ovlDis[0]) > abs(overlapRight)) {
						overlapRight = ovlDis[0];
						rightNode = tmpIndex;
					}
				}
			}
		}
	}
	closeNode.leftNode = leftNode;	//index of node on the left
	closeNode.lOvlLen = overlapLeft;	//overlap length
	closeNode.lDis = maxLeft;	//overlap distance

	closeNode.rightNode = rightNode;	//index of node on the right
	closeNode.rOvlLen = overlapRight;	//overlap length
	closeNode.rDis = minRight;	//overlap distance
	closeNode.isNull = false;

	return closeNode;
}

/**
 * Recalculate 6-tuples for an assumed left end in order to make sure it is a real left end.
 * Specifically, for the assumed left end,start to calculate from fourth level until meeting one node which
 * makes six-tuple.leftEnd != -1, or until the level we specified in the property file, then return the six-tuple.
 * If we fail to find any node, we consider it a real left end and six-tuple.leftEnd will be set to be -1.
 *
 * @param mst a Minimum Spanning Tree.
 * @param index The index of current node in the tree and graph.
 * @param sixTuple The SixTuple for this node with the index
 * @return SixTuple which stores two closest nodes which is to the left and to the right if it is not left end;
 * if it does, six-tuple.leftEnd = -1.
 */
SixTuple Graph::checkLeftEndFromMST(int index, const SixTuple sixTuple) {
	if (numOfLevels == 0) { //keep checking until leaves
		numOfLevels = INT_MAX;
	}

	vector<stack<int> > nodes (2);
	nodes[0].push(index);
	nodes[1].push(-1);

	for (int i=0; i<3; i++) { //skip over all the nodes within 3 levels
		nodes = getNodesFromMST(nodes);
	}

	SixTuple closeNode;

	for (int foundLevel=4; foundLevel<=numOfLevels; foundLevel++) { //start from level 4
		nodes = getNodesFromMST(nodes);

		if (nodes[0].size() == 0) {
			break;
		} else {
			closeNode = findAdjacentNode(nodes[0], index, sixTuple);
			if (closeNode.leftNode != -1) {
				closeNode.isNull = false;
				return closeNode;
			}
		}
	}

	closeNode.isNull = true;
	return closeNode;
}

SixTuple Graph::checkRightEndFromMST(int index, const SixTuple sixTuple) {
	if (numOfLevels == 0) { //keep checking until leaves
		numOfLevels = INT_MAX;
	}

	vector<stack<int> > nodes (2);
	nodes[0].push(index);
	nodes[1].push(-1);

	for (int i=0; i<3; i++) { //skip over all the nodes within 3 levels
		nodes = getNodesFromMST(nodes);
	}

	SixTuple closeNode;

	for (int foundLevel=4; foundLevel<=numOfLevels; foundLevel++) { //start from level 4
		nodes = getNodesFromMST(nodes);

		if (nodes[0].size() == 0) {
			break;
		} else {
			closeNode = findAdjacentNode(nodes[0], index, sixTuple);
			if (closeNode.rightNode != -1) {
				closeNode.isNull = false;
				return closeNode;
			}
		}
	}

	closeNode.isNull = true;
	return closeNode;
}

/*
 * find the most adjacent node to the current node from allNodes.
 * @param allNodes Store indices of all the nodes which will be compared to the current node.
 * @param index The index of current node.
 * @sixTuple The sixTuple for the current node.
 */
SixTuple Graph::findAdjacentNode(stack<int> nodes, int index, const SixTuple sixTuple) {
	stack<int> allNodes; //put all the values of nodes into another stack so that we won't change nodes(it's a pointer) because it may be used by the calling method.
	int stackSize = nodes.size();
	vector<int> tmpStack(stackSize);
	for (int i=stackSize-1; i>=0; i--) {
		tmpStack[i] = nodes.top();
		nodes.pop();
	}
	for (int i=0; i<stackSize; i++) { //get(0) will get the bottom element of the stack
		allNodes.push(tmpStack[i]);
	}

	SixTuple closeNode;

	int leftNode = -1;
	int rightNode = -1;
	int maxLeft = sixTuple.lDis;
	int minRight = sixTuple.rDis;
	int overlapLeft = sixTuple.lOvlLen;
	int overlapRight = sixTuple.rOvlLen;

	//find two closest nodes
	string s1 = graphNodes[index].getNodeStr();
	while (!allNodes.empty()) {
		int tmpIndex = allNodes.top();
		allNodes.pop();
		if (inc->containInclusionNode(tmpIndex)) continue;
		if (tmpIndex == index) continue;
		string s2 = graphNodes[tmpIndex].getNodeStr();

		vector<int> ovlDis = calDist.searchDistance(index, tmpIndex);
		if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) {
			ovlDis = ovl.getOVLDistance(s1, s2);

			//add to CalculatedOvlDistance
			calDist.addDistance(index, tmpIndex, ovlDis[1], ovlDis[0]);
		}

		if (ovlDis[1] == INT_MIN) { // there is inclusion between them
			if (s1.length() >= s2.length()) {
				inc->addNode2(tmpIndex, index);
			} else {
				inc->addNode2(index, tmpIndex);
			}
		} else if (ovlDis[1] != INT_MAX) { // there is overlap between them
			if (ovlDis[0] < 0) {
				if (ovlDis[1] > maxLeft) {
					maxLeft = ovlDis[1];
					overlapLeft = ovlDis[0];
					leftNode = tmpIndex;
				} else if (ovlDis[1] == maxLeft) { //if they are equal, find that one with maximal overlap
					if (abs(ovlDis[0]) > abs(overlapLeft)) {
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
				} else if (ovlDis[1] == minRight) { //if they are equal, find that one with maximal overlap
					if (abs(ovlDis[0]) > abs(overlapRight)) {
						overlapRight = ovlDis[0];
						rightNode = tmpIndex;
					}
				}
			}
		}
	}
	closeNode.leftNode = leftNode; //index of node on the left
	closeNode.lOvlLen = overlapLeft; //overlap length
	closeNode.lDis = maxLeft; //overlap distance

	closeNode.rightNode = rightNode; //index of node on the right
	closeNode.rOvlLen = overlapRight; //overlap length
	closeNode.rDis = minRight; //overlap distance
	closeNode.isNull = false;

	return closeNode;
}

/**
 * Get all the children for the input nodes and put them into a stack.
 *
 * @param mst a Minimum Spanning Tree.
 * @param nodes The stack array, nodes[0], the indexes of all the nodes for which we will find their children;
 * 		nodes[1], the indexes of the parent node of every element in nodes[0].
 * @return a stack array which stores all the found nodes' indexes. ret[0], the indexes of all the found nodes
 * 		(that is, the children of input nodes); ret[1], the indexes of the parent node of every element in nodes[0].
 */
vector<stack<int> > Graph::getNodesFromMST(vector<stack<int> > nodes) {
	vector<stack<int> > ret (2);
	while (!(nodes[0].empty())) {
		int curIndex = nodes[0].top();
		nodes[0].pop();
		int parentIndex = nodes[1].top();
		nodes[1].pop();

		for (EdgeIterator it=mst[curIndex].begin(); it!=mst[curIndex].end(); it++) {
			int index2 = (*it).node;

			if (index2 != parentIndex) {
				ret[0].push(index2);
				ret[1].push(curIndex);
			}
		}
	}
	return ret;
}

