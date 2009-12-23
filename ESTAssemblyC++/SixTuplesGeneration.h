#ifndef ZY_SixTuplesGeneration_HH
#define ZY_SixTuplesGeneration_HH
// get six-tuple for each node.

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include "SixTuple.h"
#include "InclusionNodes.h"
#include "Graph.h"
#include "Param.h"

class SixTuplesGeneration {
public:
	Graph* g;
	InclusionNodes* incNodes;
	std::vector<SixTuple*> leftMostNodes;
	std::vector<SixTuple*> alignArray;
	SixTuplesGeneration() {g=NULL; incNodes=NULL;};
	SixTuplesGeneration(Graph* graph, InclusionNodes* inc);
	inline std::vector<SixTuple*> getAlignArray() {
		return alignArray;
	}

	inline std::vector<SixTuple*> getLeftEnds() {
		return leftMostNodes;
	}

private:
	void init();
	void createAlignArray();
	void processAlignArray();
};

#endif
