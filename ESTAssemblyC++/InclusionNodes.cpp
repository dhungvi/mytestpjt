#include "InclusionNodes.h"
using namespace std;

void InclusionNodes::addNode(int chd, int parent) {
	this->nodes.insert(pair<int, int> (chd, parent));
	this->nodes_nodup.insert(pair<int, int> (chd, parent));
}

int InclusionNodes::getSize() {
	return this->nodes_nodup.size();
}

bool InclusionNodes::containInclusionNode(int idx) {
	map<int, int>::iterator ret;
	ret = nodes_nodup.find(idx);
	if (ret == nodes_nodup.end()) {
		return false;
	} else {
		return true;
	}
}

vector<int> InclusionNodes::containPNode(int pIdx) {
	vector<int> chdIdx;
	multimap<int, int>::iterator it;

	for (it = nodes.begin(); it != nodes.end(); it++) {
		if ((*it).second == pIdx) {
			chdIdx.push_back((*it).first);
		}
	}
	return chdIdx;
}

vector<int> InclusionNodes::getAllChdNodes() {
	vector<int> chdIdx;
	multimap<int, int>::iterator it;
	int tmpValue = -1; //check duplicate children index.

	for (it = nodes_nodup.begin(); it != nodes_nodup.end(); it++) {
		int idx = (*it).first;
		if (idx != tmpValue) {
			chdIdx.push_back(idx);
		}
		tmpValue = idx;
	}
	return chdIdx;
}
