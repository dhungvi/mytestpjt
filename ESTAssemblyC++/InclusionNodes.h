#ifndef ZY_InclusionNodes_HH
#define ZY_InclusionNodes_HH

#include <limits.h>
#include <string>
#include <vector>
#include <map>
#include "Param.h"

class InclusionNodes {
public:
	std::multimap<int, int> nodes;
	std::map<int, int> nodes_nodup;

	void addNode(int chd, int parent);
	int getSize();
	bool containInclusionNode(int idx);
	std::vector<int> containPNode(int pIdx);
	std::vector<int> getAllChdNodes();
};

#endif
