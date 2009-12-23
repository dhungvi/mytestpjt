#ifndef ZY_ESTAssembly_HH
#define ZY_ESTAssembly_HH

#include <ctype.h>
#include "Reconstruction.h"
#include "Graph.h"
#include "SixTuplesGeneration.h"
#include "InclusionNodes.h"
#include "Prim.h"

class ESTAssembly {
public:
	Graph* g;	//graph to store all the ests. It is generated in 'readEstFile' function.
	SixTuplesGeneration* gen;
	Reconstruction* rec;
	InclusionNodes* incNodes;

	ESTAssembly();
	void assemble();
	~ESTAssembly();

//private:
	void readEstFile(const std::string& inFileName);
	void readMST(const std::string& inFileName);
	std::vector<std::string> split(const std::string& str, char delimit);
	std::string toUpperCase(const std::string& str);
};

#endif
