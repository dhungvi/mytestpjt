#include "SixTuplesGeneration.h"

using namespace std;

SixTuplesGeneration::SixTuplesGeneration(Graph* graph, InclusionNodes* inc) {
	g = graph;
	incNodes = inc;
	init();
}

void SixTuplesGeneration::init() {
	// Get the components of the time
	cout << "Start to generate 6-tuples." << endl;
	createAlignArray();
	cout << "End to generate 6-tuples." << endl;
	cout << "There are " << incNodes->getSize() << " nodes in the inclusion list." << endl;
	cout << "There are " << alignArray.size() << " nodes in alignArray.\n" << endl;

	cout << "Start to process 6-tuples." << endl;
	processAlignArray();
	cout << "End to process 6-tuples.\n" << endl;
}

void SixTuplesGeneration::createAlignArray() {
	alignArray = g->get2CloseNodesFromMST();
}

void SixTuplesGeneration::processAlignArray() {
	vector<SixTuple*> rightMostNodes;
	map<int, int> ends; //used to avoid duplicate calculation for right most nodes
	//get all the nodes which has no left nodes or no right nodes to them
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode == -1) {
			leftMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
		if (curTuple->rightNode == -1) {
			rightMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
	}

	//Re-calculating six-tuples for those assumed left ends and put the new six-tuple into alignArray.
	int numOfLeftMostNodes = leftMostNodes.size();
	cout << "The number of original left ends is " << numOfLeftMostNodes << endl;
	cout << "The number of original right ends is " << rightMostNodes.size() << endl;
	if (numOfLeftMostNodes > 1) {
		for (int i=0; i<numOfLeftMostNodes; i++) {
			SixTuple* curTuple = leftMostNodes[i];
			int cNode = curTuple->curNode;
			ends.insert(pair<int, int> (cNode, 0));
			SixTuple lNode = g->get2CloseNodesFromGrand(cNode, *curTuple);
			curTuple->leftNode = lNode.leftNode;
			curTuple->lOvlLen = lNode.lOvlLen;	//overlap length
			curTuple->lDis = lNode.lDis;	//distance
			if (curTuple->rDis > lNode.rDis) { //get a smaller distance
				curTuple->rightNode = lNode.rightNode;
				curTuple->rOvlLen = lNode.rOvlLen;
				curTuple->rDis = lNode.rDis;
			}
		}
	}

	for (int i=0; i<rightMostNodes.size(); i++) {
		SixTuple* curTuple = rightMostNodes[i];
		int cNode = curTuple->curNode;
		pair<map<int,int>::iterator,bool> insertRet = ends.insert(pair<int, int> (cNode, 0));
		if (insertRet.second == true) { //if cNode has not been calculated as a left end
			SixTuple lNode = g->get2CloseNodesFromGrand(cNode, *curTuple);
			if (lNode.rightNode != -1) {
				curTuple->rightNode = lNode.rightNode;
				curTuple->rOvlLen = lNode.rOvlLen;
				curTuple->rDis = lNode.rDis;
			}
		}
	}

	leftMostNodes.clear();
	rightMostNodes.clear();
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode == -1) {
			leftMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
		if (curTuple->rightNode == -1) {
			rightMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
	}
	cout << "\nThere are " << leftMostNodes.size() << " left-most nodes after running 3 levels." << endl;
	cout << "There are " << rightMostNodes.size() << " right-most nodes after running 3 levels." << endl;
	ends.clear();
	/* Recalculate 6-tuples for all the current left nodes in order to remove all the false left ends.
	 * Specifically, for those assumed left ends,start to calculate from fourth level until meeting one node which
	 * makes six-tuple[0] != -1 or until the level we specified in the property file, then return the six-tuple.
	 * If we fail to find any node, we consider it a real left end.
	 */
	for (int i=0; i<leftMostNodes.size(); i++) {
		SixTuple* curTuple = leftMostNodes[i];
		int tEnd = curTuple->curNode; //index of the node
		ends.insert(pair<int, int> (tEnd, 0));
		SixTuple tmpTuple = g->checkLeftEndFromMST(tEnd, *curTuple);
		if (!tmpTuple.isNull) {
			curTuple->leftNode = tmpTuple.leftNode;
			curTuple->lOvlLen = tmpTuple.lOvlLen;	//overlap length
			curTuple->lDis = tmpTuple.lDis;	//distance
			curTuple->rightNode = tmpTuple.rightNode;
			curTuple->rOvlLen = tmpTuple.rOvlLen;
			curTuple->rDis = tmpTuple.rDis;
			curTuple->isNull = tmpTuple.isNull;
		}
	}
	for (int i=0; i<rightMostNodes.size(); i++) {
		SixTuple* curTuple = rightMostNodes[i];
		int tEnd = curTuple->curNode; //index of the node
		pair<map<int,int>::iterator,bool> insertRet = ends.insert(pair<int, int> (tEnd, 0));
		if (insertRet.second == true) { //if tEnd has not been calculated as a left end
			SixTuple tmpTuple = g->checkRightEndFromMST(tEnd, *curTuple);
			if (!tmpTuple.isNull) {
				curTuple->rightNode = tmpTuple.rightNode;
				curTuple->rOvlLen = tmpTuple.rOvlLen;
				curTuple->rDis = tmpTuple.rDis;
				curTuple->isNull = tmpTuple.isNull;
			}
		}
	}

	leftMostNodes.clear();
	rightMostNodes.clear();
	ends.clear();
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode == -1) {
			leftMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
		if (curTuple->rightNode == -1) {
			rightMostNodes.push_back(curTuple);	//store sixtuple of the node
		}
	}
	cout << "\nThere are " << leftMostNodes.size() << " left-most nodes after checking left ends." << endl;
	cout << "There are " << rightMostNodes.size() << " right-most nodes after checking right ends." << endl;

	/*
	 * Remove those false left ends.
	 *
	 * construct a temporary directed graph from alignArray,
	 *  	the second dimension has two elements:
	 *  			index of starting node,
	 *  			index of ending node,
	 *  			weight between them (positive value, weight is abs(their distance)).
	 *  	if there is no edge, weight=INT_MAX.
	 */
	int tLen = 0;
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode != -1) {
			tLen++;
		}
		if (curTuple->rightNode != -1) {
			tLen++;
		}
	}

	vector<vector<int>  > tmpDGraph = vector<vector<int> > (tLen, vector<int> (4));
	int tmpIndex = 0;
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		int curIdx = curTuple->curNode;
		if (curTuple->leftNode != -1) {
			tmpDGraph[tmpIndex][0] = curTuple->leftNode;
			tmpDGraph[tmpIndex][1] = curIdx;
			tmpDGraph[tmpIndex][2] = abs(curTuple->lDis);	//distance
			tmpDGraph[tmpIndex][3] = abs(curTuple->lOvlLen);	//overlap length
			tmpIndex++;
		}

		if (curTuple->rightNode != -1) {
			tmpDGraph[tmpIndex][0] = curIdx;
			tmpDGraph[tmpIndex][1] = curTuple->rightNode;
			tmpDGraph[tmpIndex][2] = curTuple->rDis;	//distance
			tmpDGraph[tmpIndex][3] = curTuple->rOvlLen;	//overlap length
			tmpIndex++;
		}
	}

	/*
	 * Remove those false left ends. For example:
	 * Node 2 has the set in alignArray: [-1, 0, 5, 4, 8, 9]
	 * but Node 8 has this set: [7, 5, 3, 2, 7, 6], this means ovlDis(8,2)=6,
	 * 		so node 2 is not left end because node 8 is to its left.
	 */
	//Get all the nodes which has the value of -1 in alignArray[x][0]

	vector<SixTuple*> tmpLeftNodes;
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode == -1) {
			tmpLeftNodes.push_back(curTuple);	//store index of the node
		}
	}
	//remove false left ends
	leftMostNodes.clear();
	for (int i=0; i<tmpLeftNodes.size(); i++) {
		int tEnd = tmpLeftNodes[i]->curNode;
		int f = 0;
		//if the left end appears in second element of dGraph, that means some
		//node is on its left, so it is not a real left end.
		for (int j=0; j<tmpDGraph.size(); j++) {
			if (tmpDGraph[j][1] == tEnd) {	// false left end
				f = 1;
				break;
			}
		}
		if (f == 0) {
			leftMostNodes.push_back(tmpLeftNodes[i]);
		}
	}
	cout << "\nThere are " << leftMostNodes.size() << " left-most nodes after processing false left ends." << endl;
	cout << "There are " << rightMostNodes.size() << " right-most nodes after processing false left ends." << endl;
}
