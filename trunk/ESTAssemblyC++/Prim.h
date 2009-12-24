#include <vector>

#ifndef JK_PRIM_H
#define JK_PRIM_H

class Edge {
public:
  Edge(int node, int weight) {this->node = node; this->weight = weight;}
  int node;
  int weight;
};

class NodeEntry {
public:
  NodeEntry(int node, int source, int distance) {
    this->node = node;
    this->source = source;
    this->distance = distance;
  }
  int node;       // The number of the node held by the heap element.
  int source;     // The number of the closest node that is part of the MST.
  int distance;   // The distance to the closest node that is part of the MST (source)
};

class NodeHeap {
public:
  NodeHeap(int max_size);

  void push(int node, int source, int distance);
  NodeEntry pop();
  void update(int node, int source, int newDistance);
  bool contains(int node);

  int size() {return V.size();}
  bool empty() {return V.size()==0;}

  // For debugging
  void printHash();
  void printElements();
  void printM();

  std::vector<NodeEntry> V;
  std::vector<int> M;

private:
//  std::vector<NodeEntry> V;
//  std::vector<int> M;

  void swapPositions(int p1, int p2);
  void bubbleUp(int currentNode);
  void bubbleDown(int currentNode);
};



// Graph: we will use the adjacency matrix representation
//    G[i] will be the adjacency list for node i
//    G[i][j] is the jth edge touching node i.  (Order is arbitrary.)
//    G[i][j].node is the node on the other end of the edge.
//    G[i][j].distance is the weight of the edge.
typedef std::vector<std::vector<Edge> > DefGraph;
typedef std::vector<Edge>::iterator EdgeIterator;

DefGraph Prim(DefGraph G, int source, bool directed = true, bool forrest=false);

void addNode(DefGraph& G);
void addEdge(DefGraph& G, int s, int t, int weight, bool directed = true);
int degree(DefGraph G, int n); // "Out degree" if directed

#endif
