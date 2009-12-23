#include "Reconstruction.h"
#include <sstream>
#include <cstdlib>

using namespace std;

UsedNode::UsedNode(int idx, int p) {
	index = idx;
	pos = p;
}

bool operator< (const UsedNode& k1, const UsedNode& k2) {
	if (k1.index < k2.index)
		return true;
	else if ((k1.index == k2.index) && (k1.pos < k2.pos))
		return true;
	else
		return false;
}

StartPos::StartPos(int p, int idx) {
	pos = p;
	index = idx;
}

bool operator< (const StartPos& k1, const StartPos& k2) {
	return (k1.pos < k2.pos);
}

LeftEnd::LeftEnd(int idx, int n, string s) {
	 index = idx;
	 lenOfSeq = s.length();
	 numOfUsedNodes = n;
	 seq = s;
}

bool operator< (const LeftEnd& k1, const LeftEnd& k2) {
	return (k1.lenOfSeq < k2.lenOfSeq);
}

Reconstruction::Reconstruction() {
	g = NULL;
	incNodes = NULL;
	consensusFileName = CONSENSUSFILE;
	singletonFileName = SINGLETONFILE;
	numOfUsedESTsFileName = NUMOFUSEDEST;
}

Reconstruction::Reconstruction(Graph* graph, vector<SixTuple*> align, vector<SixTuple*> leftEnds, InclusionNodes* inc) {
	g = graph;
	alignArray = align;
	incNodes = inc;
	leftMostNodes = leftEnds;
	consensusFileName = CONSENSUSFILE;
	singletonFileName = SINGLETONFILE;
	numOfUsedESTsFileName = NUMOFUSEDEST;

	usedNodes = vector<int> (g->graphNodes.size(), 0);

	vector<int> chdNodes = incNodes->getAllChdNodes(); //all the children nodes in inclusion list should not be considered as singletons.
	for (int i=0; i<chdNodes.size(); i++) {
		usedNodes[chdNodes[i]] = 1;
	}
}

void Reconstruction::getConsensus() {
	cout << "Start to reconstruct." << endl;
	printConsensus();
	cout << "End to reconstruct." << endl;
}

/*
 * Print the original sequence and multiple consensus into a file which is specified in the property file.
 */
void Reconstruction::printConsensus() {
	/*
	 * print consensus sequences
	 */
	 ofstream outFile1;
	 outFile1.open(consensusFileName.c_str(),ios::trunc);

	 vector<string> consensus = reconstruct();
	 int index = 1;
	 for (int i=0; i<consensus.size(); i++) {
		 string str = consensus[i];
		 size_t found = str.find('\n');

		 if (found != string::npos) { //there is "\n" in the sequence
			 outFile1 << ">contig " << (index++) << endl;

			 int start = 0;
			 while (found != string::npos) {
				 outFile1 << str.substr(start, (int)found-start) << endl;
				 start = found + 1;
				 found = str.find('\n', found+1);
			 }
			 outFile1 << str.substr(start) << endl;
		 } else {
			 outFile1 << ">contig " << (index++) << endl;
			 outFile1 << consensus[i] << endl;
		 }
	 }
	 outFile1.close();


	 /*
	  * print singletons and numOfUsedESTs
	  */
	 ofstream outFile2;
	 outFile2.open(singletonFileName.c_str(),ios::trunc);

	 int num = 0;
	 for (int i=0; i<usedNodes.size(); i++) {
		 if (usedNodes[i] == 0) {	//singleton
			 outFile2 << g->getCommentOfNode(i) << "\n";
			 outFile2 << g->getSeqOfNode(i) << "\n";
		 } else {
			 num++;
		 }
	 }
	 outFile2.close();

	 ofstream outFile3;
	 outFile3.open(numOfUsedESTsFileName.c_str(),ios::trunc);
	 outFile3 << ">number of used ESTs by EAST\n";
	 outFile3 << num;
	 outFile3.close();
}

/*
 * reconstruct the sequence
 *
 * @return an arraylist: The assembled sequences.
 */
std::vector<std::string> Reconstruction::reconstruct() {
	sPosDebug = vector<int> (g->graphNodes.size());

	vector<string> ret = processLeftEnds();

	//for debug
	std::stringstream out;
	out << printStr << ret.size() << " consensus from above " << leftMostNodes.size() << " left ends:\n";
	printStr = out.str();
	for (int p=0; p<ret.size(); p++) {
		printStr = printStr + ret[p] + "\n";
	}

	//print debug information about the generated consensus.
	cout << "*********************consensus:********************\n";
	cout << printStr << endl;

	return ret;
}

std::vector<std::vector<int> > Reconstruction::genDGraph() {
	/*
	 * Calculate the length of dGraph.
	 */
	int len = 0;
	for (int i=0; i<alignArray.size(); i++) {
		if (alignArray[i]->leftNode != -1) {
			len++;
		}
		if (alignArray[i]->rightNode != -1) {
			len++;
		}
	}

	/*
	 * generate dGraph.
	 */
	vector<vector<int> > dGraph = vector<vector<int> > (len, std::vector<int>(4));
	int indexOfDGraph = 0;
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode != -1) {
			dGraph[indexOfDGraph][0] = curTuple->leftNode;
			dGraph[indexOfDGraph][1] = curTuple->curNode;
			dGraph[indexOfDGraph][2] = abs(curTuple->lDis);	//distance
			dGraph[indexOfDGraph][3] = abs(curTuple->lOvlLen);	//overlap length
			indexOfDGraph++;
		}

		if (curTuple->rightNode != -1) {
			dGraph[indexOfDGraph][0] = curTuple->curNode;
			dGraph[indexOfDGraph][1] = curTuple->rightNode;
			dGraph[indexOfDGraph][2] = curTuple->rDis;	//distance
			dGraph[indexOfDGraph][3] = curTuple->rOvlLen;	//overlap length
			indexOfDGraph++;
		}
	}

	return dGraph;
}


/*
 *  If flag = 0, return number of used nodes for this left end.
 *  If flag = 1, return consensus from this left end.
 *
 * 1. Generate dGraph.
 * 2. For each left-end node, starting from it to calculate positions for each node.
 * In order to get starting positions, it constructs a MST. The weight of the Minimum Spanning tree is
 * the overlap distance instead of overlap length.
 * Then reconstruct the sequence from the set of ESTs.
 *
 *  	Because the Prim algorithm starts from index 0 to generate MST, we have to
 *  		put left-end node to index 0 in order to get the MST we want. If Prim does
 *  		not start from the left-end node, the directed tree will be unconnected.
 *
 */
string Reconstruction::getInfoOfLeftEnd(int leftEnd, vector<vector<int> >& dGraph, int flag) {
	string ret = "";
	sPos = vector<int> (g->graphNodes.size());	//store starting positions of all the nodes


	// Calculate starting positions using minimum spanning tree starting from this left-end node.
	DefGraph primMST = constructMinTree(g->graphNodes.size(), dGraph, leftEnd); //the first param is the total number of ESTs.

	//get starting positions for the nodes in primMST
	getStartPos(leftEnd, primMST, dGraph);


	vector<StartPos> tmpArray;
	for (int j=0; j<sPos.size(); j++) {
		if ((j == leftEnd) || (sPos[j] != 0)) {
			tmpArray.push_back(StartPos(sPos[j], j));
		}
	}


	if (flag == 0) { //get number of used nodes
		std::stringstream out;
		out << getNumUsedNodes(tmpArray);
		ret = out.str();
	} else if (flag == 1) { //get consensus
		printStr = printStr + printLeftEndInfo(leftEnd); //get information of this left end which is used to start reconstruction.
		vector<string> tStr = reconstructSeq(tmpArray);
		if (tStr.size() != 0) {
			//this is an estimated value, not an exact one because some nodes may repeat (both from nodes2 in inclusion list and also exist in the ordinary list)
			//the value in numOfUsedEsts file is an exact value.
			printStr = printStr + tStr[1] + " nodes are used to reconstruct the sequence.(estimated. real number<= the number)\n";
			printStr = printStr + tStr[0] + "\n\n";
			ret = tStr[0];
		}
	}

	return ret;
}

/*
* This method is designed to group left ends into independent sets. All the dependent left ends (include one by one) are put
* into one group.
* In each group, we find the left end with the longest length and the one which will use the most number of ESTs to reconstruct,
* and combine them together to form a consensus by calling the function "processLeftEndsWithInclusion".

* @return all the consensus sequences.
*/
vector<string>Reconstruction::processLeftEnds() {
	vector<string> allOutputContigs;	//store all the generated sequences

	vector<vector<int> > dGraph = genDGraph();
	int sizeOfs = leftMostNodes.size();
	if (sizeOfs == 0) {
		cout << "zero left ends!" << endl;
		exit(1);
	}
	vector<LeftEnd> resultArray(sizeOfs); //store the starting positions of ests
	for (int i=0; i<sizeOfs; i++) { //start for
		int idx = leftMostNodes[i]->curNode;
		string num = getInfoOfLeftEnd(idx, dGraph, 0);
		resultArray[i] = LeftEnd(idx, atoi(num.c_str()), g->getSeqOfNode(idx));
	} //end for
	sort(resultArray.begin(), resultArray.end());

	vector<LeftEnd> allLeftEnds;
	for (int i=sizeOfs-1; i>=0; i--) {
		allLeftEnds.push_back(resultArray[i]);
	}

	while (true) {
		string s1 = allLeftEnds[0].seq;
		vector<LeftEnd> includedEnds;
		includedEnds.push_back(allLeftEnds[0]);
		vector<LeftEnd> excludedEnds;
		for (int i=1; i<allLeftEnds.size(); i++) {
			bool b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1); //if resultArray[i].firstEst is included in s1
			if (b) {
				includedEnds.push_back(allLeftEnds[i]);
			} else {
				excludedEnds.push_back(allLeftEnds[i]);
			}
		}

		string contig = processLeftEndsWithInclusion(includedEnds);
		if ((contig.size() != 0) && (contig.erase(contig.find_last_not_of(" \n\r\t")+1).size() != 0)) { //trim
			allOutputContigs.push_back(contig);
		}
		if (excludedEnds.size() == 0) {
			break;
		} else {
			allLeftEnds = excludedEnds;
		}
	}
	return allOutputContigs;
 }

/*
 * This method is called by "processLeftEnds".
 * This method is used to process those left ends that include each other.
 *
 * We combine the one with the longest first EST and the one using the most number of ESTs together to
 * form a new one and return it.
 */
string Reconstruction::processLeftEndsWithInclusion(vector<LeftEnd>& includeStrs) {
	 if (includeStrs.size() == 0) {
		 return "";
	 } else if (includeStrs.size() == 1) {
		 vector<vector<int> > dGraph = genDGraph();
		 string s = getInfoOfLeftEnd(includeStrs[0].index, dGraph, 1);
		 return s;
	 } else {
		 int maxLen = 0; //the maximal length of the first EST.
		 int idxMaxLen = 0;
		 int maxNumNodes = 0;
		 int idxMaxNumNodes = 0;

		 for (int i=0; i<includeStrs.size(); i++) {
			 int tLen = includeStrs[i].lenOfSeq;
			 if (tLen > maxLen) {
				 maxLen = tLen;
				 idxMaxLen = i;
			 }

			 int num = includeStrs[i].numOfUsedNodes;
			 if (num > maxNumNodes) {
				 maxNumNodes = num;
				 idxMaxNumNodes = i;
			 }
		 }

		 string s1 = includeStrs[idxMaxLen].seq;
		 usedNodes[includeStrs[idxMaxLen].index] = 1;	//mark that the node is used.
		 vector<vector<int> > dGraph = genDGraph();
		 string s2 = getInfoOfLeftEnd(includeStrs[idxMaxNumNodes].index, dGraph, 1);

		 AlignResult strs = alignment.getLocalAlignment(s1, s2);
		 int offset = s1.find(replace(strs.str1, "-", ""));
		 return (s1.substr(0, offset) + s2);
	 }
}

/*
 * reconstruct a sequence which starts from a left end.
 * return: ret[0]-the consensus, ret[1]-the total number of nodes used for reconstruction.
 */
vector<string> Reconstruction::reconstructSeq(vector<StartPos>& a) {
	vector<string> ret(2);
	int sizeOfa = a.size();
	if ((sizeOfa == 0) || (sizeOfa == 1)){ //size=1, this node should be a singleton if it is not used in other place. it shouldn't appear as a consensus sequence.
		return ret;
	}

	sort(a.begin(), a.end());

	vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction.
	std::stringstream out;
	out << addedNodes.size();
	ret[1] = out.str();

	vector<StartPos> resultArray(addedNodes.size());
	for (int i=0; i<addedNodes.size(); i++) {
		UsedNode tmpNode = addedNodes[i];
		resultArray[i] = StartPos(tmpNode.pos, tmpNode.index);
	}
	sort(resultArray.begin(), resultArray.end());

	vector<SingleBase*> bases;
	string tConsensus = g->getSeqOfNode(resultArray[0].index);
	usedNodes[resultArray[0].index] = 1;	//mark that the node is used.
	string curSeq = "";
	int len = resultArray.size() - 1;
	for (int i=1; i<=len; i++) {
		curSeq = g->getSeqOfNode(resultArray[i].index);
		usedNodes[resultArray[i].index] = 1;	//mark that the node is used.
		AlignResult strs = alignment.getLocalAlignment(tConsensus, curSeq);
		tConsensus = replace(tConsensus, replace(strs.str1, "-", ""), strs.str1);
		int offset = (int)tConsensus.find(strs.str1);

		string tSeq = replace(curSeq, replace(strs.str2, "-", ""), strs.str2);
		curSeq = tSeq.substr((int)tSeq.find(strs.str2));
		if (i == 1) {
			int len1 = tConsensus.size();
			int len2 = curSeq.size();
			int end = max(offset+len2, len1);
			for (int j=0; j<end; j++) {
				if ((j < len1) && (j-offset >= 0) && (j-offset < len2)) { //overlap part
					bases.push_back(new SingleBase(tConsensus[j], curSeq[j-offset]));
				} else if ((j-offset < 0) || (j-offset >= len2)) {
					bases.push_back(new SingleBase(tConsensus[j]));
				} else if (j >= len1) {
					bases.push_back(new SingleBase(curSeq[j-offset]));
				}
			}
		} else {
			int len1 = tConsensus.size();
			int len2 = curSeq.size();
			int end = max(offset+len2, len1);
			for (int j=offset; j<end; j++) {
				if ((j < len1) && (j-offset < len2)) { //overlap part
					char c1 = tConsensus[j];
					char c2 = curSeq[j-offset];

					if (c1 != '-') {
						bases[j]->addOneBase(c2);
					} else {
						bases.insert(bases.begin()+j, new SingleBase(c1, c2));
					}
				} else if (j >= len1) {
					bases.push_back(new SingleBase(curSeq[j-offset]));
				}
			}
		}

		tConsensus = getCurConsensus(bases);
	}

	ret[0] = replace(tConsensus, "P", "");

	//release bases
	for (int i=0; i<bases.size(); i++) {
		delete bases[i];
	}
	return ret;
}

int Reconstruction::getNumUsedNodes(vector<StartPos>& a) {
	int sizeOfa = a.size();
	if (sizeOfa == 0) {
		return 0;
	} else if (sizeOfa == 1) {
		return 1;
	}

	sort(a.begin(), a.end());

	vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction.
	return addedNodes.size();
}


string Reconstruction::getCurConsensus(vector<SingleBase*> bases) {
	int len = bases.size();
	string tStr;
	for (int i=0; i<len; i++) {
		tStr.push_back(bases[i]->getCurBase());
	}
	return tStr;
}

/*
 * in string 'str', replace 'old' with 'new'.
 */
string Reconstruction::replace(string& str, const string& old, const string& newstr) {
	size_t found = str.find(old);
	while (found != string::npos) { //there is "old" in the sequence
		str.replace(found, old.size(), newstr);
		found = str.find(old, found+old.size());
	}
	return str;
}


/*
 * Add all the inclusion nodes into the input arraylist.
 * For each element in the arraylist, put its corresponding node just after it.
 */
vector<UsedNode> Reconstruction::addInclusionNodes(vector<StartPos>& input) {
	vector<UsedNode> retList;
	int size = input.size();

	for (int i=0; i<size; i++) {
		int curIdx = input[i].index;
		int pos = input[i].pos;
		retList.push_back(UsedNode(curIdx, pos));

		vector<int> chdIdx = incNodes->containPNode(curIdx); //inclusion children index of the curIdx if exist.
		for (int j=0; j<chdIdx.size(); j++) {
			retList.push_back(UsedNode(chdIdx[j], pos+1));
		}
	}
	sort(retList.begin(), retList.end());
	return retList;
}

/*
 * Construct a directed Miminum spanning tree.
 *
 *  @param nOfNodes number of nodes
 *  @param g a directed graph, the second dimension has three elements:
 *  	index of starting node, index of ending node, weight between them.
 */
DefGraph Reconstruction::constructMinTree(int nOfNodes, vector<vector<int> >& input, int source) {
	// Make a directed graph.
	DefGraph dGraph(nOfNodes);
	for (int j=0; j<input.size(); j++) {
		if (input[j][3] != 0) {	//there is an edge between the nodes
			dGraph[input[j][0]].push_back(Edge(input[j][1], input[j][2]));
		}
	}
	DefGraph mst = Prim(dGraph, source, true);
	return mst;
}

/*
 * Calculate starting positions for each node.
 */
void Reconstruction::getStartPos(int parentNode, DefGraph& tree, vector<vector<int> >& d) {
	for (EdgeIterator it=tree[parentNode].begin(); it!=tree[parentNode].end(); it++) {
		int index = (*it).node;

		int overlapLen = 0;
		for (int i=0; i<d.size(); i++) {
			if ((d[i][0] == parentNode) && (d[i][1] == index)) {
				overlapLen = d[i][3];
				break;
			}
		}

		sPos[index] = sPos[parentNode] + g->getLenOfNode(parentNode) - overlapLen;
		getStartPos(index, tree, d);
	}
}


string Reconstruction::printLeftEndInfo(int leftEnd) {
	string ret = "";
	std::stringstream out;
	out << "Node " << leftEnd << " starts from " << g->getNameOfNode(leftEnd) << endl;
	ret = out.str();
	out.str("");
	out.clear();

	int sp = atoi(g->getNameOfNode(leftEnd).c_str());	//actual starting position of the node
	int ln = g->getLenOfNode(leftEnd);
	int flag = 0;
	for (int k = 0; k < g->getSizeofGraph(); k++) {
		int tmpSp = atoi(g->getNameOfNode(k).c_str( ));
		int tmpLn = g->getLenOfNode(k);
		if ((sp > tmpSp) && (sp < (tmpSp + tmpLn - 1)) && ((sp + ln - 1)
				> (tmpSp + tmpLn - 1))) {
			//if overlap length is less than windowsize, we consider they're not overlapping.
			//if not, we see them as overlapping,so this is not a real left end.
			if ((tmpSp + tmpLn - sp) >= ((g->ovl).d2.getWindowSize())) {
				out << ret << "Node " << leftEnd << " is not a real left-most node. \n";
				ret = out.str();
				out.str("");
				out.clear();
				out << ret << "Overlap length with node " << k << " is " << (tmpSp + tmpLn - sp) << "\n"; //(sp+ln-tmpSp));
				ret = out.str();
				out.str("");
				out.clear();
				flag = 1;
			}
		}
	}

	if (flag == 0) {
		out << ret << "Node " << leftEnd << " is a real left-most node.\n";
		ret = out.str();
	}
	return ret;
}
