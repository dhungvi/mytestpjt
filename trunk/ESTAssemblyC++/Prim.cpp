#include "Prim.h"
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <assert.h>
#include <limits.h>

using namespace std;

NodeHeap::NodeHeap(int max_size) {
  M = vector<int>(max_size);
}

void NodeHeap::swapPositions(int p1, int p2) {
  M[V[p1].node] = p2;
  M[V[p2].node] = p1;

  NodeEntry tmp = V[p1];
  V[p1] = V[p2];
  V[p2] = tmp;
}


void NodeHeap::bubbleUp(int currentNode) {
  while (currentNode > 0) {
    int parentIndex = (currentNode-1)/2;
    if (V[currentNode].distance < V[parentIndex].distance) {
      swapPositions(currentNode, parentIndex);
      currentNode = parentIndex;
    }
    else
      break;
  }
}

void NodeHeap::bubbleDown(int currentNode) {
  while (1) {
    int leftChild = 2*currentNode + 1;
    if (leftChild >= (int)V.size())
      break;

    int rightChild = leftChild+1;

    int bestIndex = V[currentNode].distance < V[leftChild].distance ? currentNode : leftChild;
    if (rightChild < (int)V.size() && V[rightChild].distance < V[bestIndex].distance)
      bestIndex = rightChild;

    if (bestIndex == currentNode)
      break;

    swapPositions(bestIndex, currentNode);
    currentNode = bestIndex;
  }
}

bool NodeHeap::contains(int node) {
  return M[node] != -1;
}

void NodeHeap::push(int node, int source, int distance) {
  assert(node < (int)M.size());
  NodeEntry ne(node, source, distance);
  M[node] = V.size();

  int currentNode = V.size();
  V.push_back(ne);
  bubbleUp(currentNode);
}

NodeEntry NodeHeap::pop() {
  NodeEntry ne = V[0];
  M[ne.node] = -1;
  V[0] = V[V.size()-1];
  V.pop_back();
  bubbleDown(0);

  return ne;
}

void NodeHeap::update(int node, int source, int newDistance) {
  int currentNode = M[node];
  assert(currentNode != -1);
  if (newDistance < V[currentNode].distance) {
    V[currentNode].distance = newDistance;
    V[currentNode].source = source;
    bubbleUp(currentNode);
  }
}

DefGraph Prim(DefGraph G, int source, bool directed, bool forrest) {
  DefGraph T(G.size());
  NodeHeap H(G.size());
  for (int i=0; i < (int)G.size(); i++)
    H.push(i, -1, INT_MAX);

  H.update(source,-1,0);
  while (!H.empty()) {
   NodeEntry ne = H.pop();
   if (!forrest && (ne.distance == INT_MAX))
     break;

   if (ne.source != -1) {
     T[ne.source].push_back(Edge(ne.node, ne.distance));
     if (!directed)
       T[ne.node].push_back(Edge(ne.source, ne.distance)); // Comment out for undirected graph
   }

   for (EdgeIterator i = G[ne.node].begin(); i != G[ne.node].end(); i++)
     if (H.contains(i->node)) {
       H.update(i->node, ne.node, i->weight);
     }

  }
  return T;
}

void NodeHeap::printHash() {
  while (!empty()) {
    NodeEntry p = pop();
    cout << p.node << " " << p.source << " " << p.distance << "\n";
  }
}

void NodeHeap::printElements() {
  for (vector<NodeEntry>::iterator i=V.begin(); i!=V.end(); i++) {
    cout << i->node << " " << i->source << " " << i->distance << "\n";
  }
}

void NodeHeap::printM() {
  for (vector<int>::iterator i=M.begin(); i!=M.end(); i++) {
    cout << (*i) << endl;
  }
}

void addNode(DefGraph& G) {
  G.push_back(vector<Edge>());
}

void addEdge(DefGraph& G, int s, int t, int weight, bool directed) {
  assert(s < (int)G.size() && t < (int)G.size());
  G[s].push_back(Edge(t,weight));
  if (!directed)
    G[t].push_back(Edge(s,weight));
}

int degree(DefGraph G, int n) {
  assert(n < (int)G.size());
  return G[n].size();
}

