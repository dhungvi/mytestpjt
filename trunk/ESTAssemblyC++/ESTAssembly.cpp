#include "ESTAssembly.h"
#include <cstdlib>
using namespace std;


ESTAssembly::ESTAssembly(std::string& estF, std::string& mstF) {
	incNodes = new InclusionNodes();
	g = new Graph(incNodes);
    gen = NULL;
    rec = NULL;
    estFileName = estF;
    mstFileName = mstF;
}

void ESTAssembly::assemble(std::string& con, std::string& sing, std::string& numF) {
	readEstFile(estFileName);
	readMST(mstFileName);
	gen = new SixTuplesGeneration(g, incNodes);
	rec = new Reconstruction(g, gen->getAlignArray(), gen->getLeftEnds(), incNodes, con, sing, numF);
	rec->getConsensus();
}

ESTAssembly::~ESTAssembly() {
	delete incNodes;
	delete g;
	delete gen;
	delete rec;
}
/*
 * Read ests from the input file;
 * Generate a Graph object, all the ests are considered to be one node in the graph;
 * No edge in the graph. Edges will be added in "createAlignArray" function.
 *
 * FASTA format:
 * A sequence in FASTA format begins with a single-line description, followed by lines of
 * sequence data. The description line is distinguished from the sequence data by a greater-than
 * (">") symbol in the first column. The word following the ">" symbol is the identifier of the
 * sequence, and the rest of the line is the description (both are optional). There should be no
 * space between the ">" and the first letter of the identifier. It is recommended that all lines
 * of text be shorter than 80 characters. The sequence ends if another line starting with a ">"
 * appears; this indicates the start of another sequence.
 *
 * EST file comment format:
 * >g001_001169_001679: first one is number of gene, second is index of tarting position, third is
 * 						index of ending position. Index starts from 0.
 */
void ESTAssembly::readEstFile(const string& inFileName) {
	vector<string> ests; //store all the ests.

	ifstream in;
	in.open(inFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);
		while (str.size() != 0) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			// first line is comment line which begins from '>'
			if (str[0] == '>') {	//comment line begins from '>'
				vector<string> paras = split(str, '_');
				ests.push_back(paras[1]);
				ests.push_back(str.substr(1)); //comment

				//get est in the next lines
				getline(in, str);
				string estStr;
				while (str.size() != 0) {
					str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
					if (str.size() != 0)	{
						if (str[0] != '>') {
							estStr.append(str);
						} else  {
							ests.push_back(estStr);
							break;
						}
					}
					getline(in, str);
				}
				if (str.size() == 0) {
					ests.push_back(estStr);
				}
			}
		}
	}
	in.close();

	//generate a graph from the input file
	int i=0;
	while (i<ests.size()) {
		//the sequence is upper-case
		g->addNode(Node(ests[i], ests[i+1], toUpperCase(ests[i+2])));
		i = i+3;
	}
}

/*
 * read a minimum spanning tree from the input MST file
 */
void ESTAssembly::readMST(const string& inFileName) {
	int nOfNodes = g->graphNodes.size();
	if (nOfNodes == 0) {
		cout << "zero nodes!" << endl;
		exit(1);
	}

	vector<vector<int> > nodes(nOfNodes-1, vector<int>(3));	//store edges in MST, there are n-1 edges, n is number of nodes.

	//read mst from the input file
	ifstream in;
	in.open(inFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);

		int curIndex = 0;
		while (str.size() != 0) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str[0] != '#') {	//comment line begins from '#'
				vector<string> paras = split(str, ',');
				if (atoi(paras[0].c_str()) != -1) {	//-1 means root of MST
					int i0 = atoi(paras[0].c_str());
					int i1 = atoi(paras[1].c_str());
					int i2 = atoi(paras[2].c_str());
					nodes[curIndex][0] = i0;
					nodes[curIndex][1] = i1;
					nodes[curIndex][2] = i2;
					curIndex++;
				}
			}
			getline(in, str);
		}
	}
	in.close();

	// Make a undirected MST.
	DefGraph mst(nOfNodes);
	for (int j=0; j<nodes.size(); j++) {
		addEdge(mst, nodes[j][0], nodes[j][1], nodes[j][2], false);
	}

	g->setMst(mst);
}


vector<string> ESTAssembly::split(const string& str, char delimit) {
	vector<string> retStr;
	size_t found = str.find(delimit);
	if (found != string::npos) {
		int start = 0;
		while (found != string::npos) {
			retStr.push_back(str.substr(start, (int) found - start));
			start = found + 1;
			found = str.find(delimit, found + 1);
		}
		retStr.push_back(str.substr(start));
	}
	return retStr;
}

string ESTAssembly::toUpperCase(const string& str) {
	string retStr;
	for (int i=0; i<str.size(); i++) {
		retStr.push_back(toupper(str[i]));
	}
	return retStr;
}
