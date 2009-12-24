#ifndef ZY_Reconstruction_HH
#define ZY_Reconstruction_HH
// reconstruction from all the ESTs

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include "SixTuple.h"
#include "Prim.h"
#include "Alignment.h"
#include "InclusionNodes.h"
#include "Graph.h"
#include "Param.h"
#include "SingleBase.h"

class UsedNode {
public:
	int index; //index of the node
	int pos;
	UsedNode(int idx, int p);
	friend bool operator< (const UsedNode& k1, const UsedNode& k2);
};

class StartPos {
public:
	int pos;
	int index; //index of the node
	StartPos() {};
	StartPos(int p, int idx);
	friend bool operator< (const StartPos& k1, const StartPos& k2);
};

class LeftEnd {
public:
	int index;
	int lenOfSeq;
	int numOfUsedNodes;
	std::string seq;
	LeftEnd() {};
	LeftEnd(int idx, int n, std::string s);
	friend bool operator< (const LeftEnd& k1, const LeftEnd& k2);
};

class Reconstruction {
public:
	Graph* g;
	std::vector<SixTuple*> alignArray;
	std::vector<int> sPos;	//starting positions of all the nodes
				//the index in the array is the index of the node, the value is its starting position.
	std::vector<int> sPosDebug;	//starting positions of all the nodes
						//it's used for debugging. All the left ends will be assigned to their actual value in order to calculate inversions later.
	Alignment alignment;
	InclusionNodes* incNodes;
	std::vector<SixTuple*> leftMostNodes;
	std::string consensusFileName;
	std::string singletonFileName;
	std::string numOfUsedESTsFileName;
	std::string printStr;
	std::vector<std::string> numOfNodes;	//store number of used nodes, corresponds to each element in allConsensus.
	std::vector<std::string> firstEsts;	//includes all the left end sequence corresponding to each element in allConsensus.
	std::vector<int> usedNodes;	//index of the array is the index of the node, 1-used, 0-not used. It is used to identify singletons.

	Reconstruction();
	Reconstruction(Graph* graph, std::vector<SixTuple*> align, std::vector<SixTuple*> leftEnds, InclusionNodes* inc, const std::string& con, const std::string& sing, const std::string& numF);
	void getConsensus();
	std::vector<std::string> reconstructSeq(std::vector<StartPos>& a);

private:
	void printConsensus();
	std::vector<std::string> reconstruct();
	std::vector<std::vector<int> > genDGraph();
	std::string getInfoOfLeftEnd(int leftEnd, std::vector<std::vector<int> >& dGraph, int flag);
	std::vector<std::string> processLeftEnds();
	std::string processLeftEndsWithInclusion(std::vector<LeftEnd>& includeStrs);
	int getNumUsedNodes(std::vector<StartPos>& a);
	std::string getCurConsensus(std::vector<SingleBase*> bases);
	std::vector<UsedNode> addInclusionNodes(std::vector<StartPos>& input);
	DefGraph constructMinTree(int nOfNodes, const std::vector<std::vector<int> >& input, int source);
	void getStartPos(int parentNode, DefGraph& tree, std::vector<std::vector<int> >& d);
	std::string replace(std::string str, const std::string& old, const std::string& newstr);
	std::string printLeftEndInfo(int leftEnd);
};



#endif
