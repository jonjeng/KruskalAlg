/*****************************************
 Project 5 - CIS 22C
 
 * Contributors:
 * Evan Finnigan
 * Forest Finnigan
 * Jonathan Jeng
 * Abhishek Rajbahndri
 *****************************************/

#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

#ifndef nullptr
#define nullptr 0
#endif

#include <iostream>
#include <limits>
#include "Exception.h"
#include "Dynamic_queue.h"

// My current idea is to use an adjacency matrix and this is potentially more costly, but allows for unknowns, such as the number of edges. Can also represent with an edge list but searching this list is more costly

class Weighted_graph {
private:
    static const double INF;
    int numVertices;
    int numEdges;
    int* graphVertices;
    double** weights;
    
    // Do not implement these functions!
    // By making these private and not implementing them, any attempt
    // to make copies or assignments will result in errors
    Weighted_graph( Weighted_graph const & );
    Weighted_graph &operator=( Weighted_graph );
    
    // your choice
    
public:
    Weighted_graph( int = 10 );
    ~Weighted_graph();
    
    int degree( int ) const;
    int edge_count() const;
    std::pair<double, int> minimum_spanning_tree() const;
    
    bool insert_edge( int, int, double );
    bool erase_edge( int, int );
    void clear_edges();
    
    // Friends
    
    friend std::ostream &operator<<( std::ostream &, Weighted_graph const & );
};

const double Weighted_graph::INF = std::numeric_limits<double>::infinity();

Weighted_graph::Weighted_graph( int n ) :
numVertices(n), numEdges(0) {
    // Assume a graph with vertices numbered from 0 through n − 1. The vertex i is given an ?initial? priority of i.

    weights = new double*[numVertices];
    for (int i = 0; i < numVertices; ++i)
        weights[i] = new double[numVertices];
    for (int i = 0; i < numVertices; i++)
        for (int j = 0; j < numVertices; j++)
            weights[i][j] = INF;
}

Weighted_graph::~Weighted_graph() {
    // The destructor frees up any memory allocated by the constructor.
    delete graphVertices;
    for (int i = 0; i < numVertices; ++i)
        delete [] weights[i];
    delete [] weights;
}

int Weighted_graph::degree(int vertex) const {
    // If the vertex is not in the graph throw an illegal argument exception.
    if (vertex < 0 || vertex >= numVertices)
        throw illegal_argument();
    
    // Count the number of edges incident upon the given vertex
    int degree = 0;
    for (int i = 0; i < numVertices; i++) {
        if (weights[vertex][i] != INF)
            degree++;
    }
    
    return degree;
}

int Weighted_graph::edge_count() const {
    // Returns the number of edges in the graph.
    int count = 0;
    for (int i = 0; i < numVertices; i++)
        for (int j = 0; j < numVertices; j++)
            if (weights[i][j] != INF) count++;
    return (count/2);
}

bool Weighted_graph::insert_edge( int i, int j, double w ) {
    // If i or j are outside the range of vertices or if the weight w is not positive, the argument(s) is/are illegal
    if ((i<0) || (i>numVertices) || (j<0) || (j>numVertices) || (w <= 0))
        throw illegal_argument();
    // If i equals j and are in the graph, return false (self edges are inconsequential)
    if (i==j && i >= 0 && i < numVertices && j >= 0 && j < numVertices)
        return false;
    // If a parallel edge of inconsequential (equal or heavier) weight would be formed, return false
    if (weights[i][j] != INF && w >= weights[i][j])
        return false;
    // Otherwise, the new edge is safe to insert, as there was either no edge there yet, or an edge of greater weight (parallel edge - but shouldn't we only update the weight if it is lower than the preexisting weight?).
    weights[i][j] = w;
    weights[j][i] = w;
    numEdges++;
    return true;
}

bool Weighted_graph::erase_edge(int i, int j) {
    // If i and j are outside the legal range, they are invalid - throw an exception
    if (i < 0 || i > numVertices || j < 0 || j > numVertices) throw illegal_argument();
    
    // If a specified edge exists, remove it (set its weight to INF) and return true.
    if (weights[i][j] != INF || i == j) {
        weights[i][j] = INF;
        weights[j][i] = INF;
        numEdges--;
        return true;
    }
    
    // Otherwise, the erasure cannot be completed so return false
    return false;
}

void Weighted_graph::clear_edges() {
    // Removes all the edges from the graph.
    if (numEdges > 0) {
        numEdges = 0;
        for (int i = 0; i < numVertices; i++)
            for (int j = 0; j < numVertices; j++)
                weights[i][j] = INF;
    }
}


std::pair<double, int> Weighted_graph::minimum_spanning_tree() const {
    // If there aren't at least n + 1 vertices, where n is the number of edges, there is no spanning tree
    if (numEdges >= numVertices - 1) {
        // Use Kruskal's algorithm to find the minimum spanning tree. You will return the weight of the minimum spanning tree and the number of edges that were tested for insertion into the minimum spanning tree.
        int edgesTested = 0, nodesLeft = numVertices, edgeNode1[numEdges], edgeNode2[numEdges], edgeTop = 0, disjointSetEdges[numVertices];
        double MSTWeight = 0.0, edgeWeights[numEdges];
        
        for (int i = 0; i < numVertices; i++)
            disjointSetEdges[i] = -1;
        
        // Store edges separately
        double** weightsCopy = new double*[numVertices];
        for (int i = 0; i < numVertices; i++)
            weightsCopy[i] = new double[numVertices];
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                weightsCopy[i][j] = weights[i][j];
            }
        }
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                if (weightsCopy[i][j] != INF) {
                    edgeNode1[edgeTop] = i;
                    edgeNode2[edgeTop] = j;
                    edgeWeights[edgeTop] = weightsCopy[i][j];
                    edgeTop++;
                    weightsCopy[j][i] = INF;
                }
            }
        }
        
        while (nodesLeft > 0) {
            // While a spanning tree has not been discovered:
            // Find the edge with the least weight (for such small numbers, selection sort is fine, so brute force searching for the minimum is adequate
            double currWeight = edgeWeights[0];
            int currEdgeNode1 = edgeNode1[0], currEdgeNode2 = edgeNode2[0];
            for (int i = 0; i < nodesLeft; i++) {
                if (edgeWeights[i] > currWeight) {
                    currWeight = edgeWeights[i];
                    currEdgeNode1 = edgeNode1[i]; currEdgeNode2 = edgeNode2[i];
                }
            }
            
            
            // Test to see if the edge creates a cycle (such is only possible if the edge to be formed is between two vertices upon which are incident edges that have already been selected), using the union-find algorithm
            int node1 = edgeNode1[edgeTop-1], node2 = edgeNode2[edgeTop-1], root1 = node1, root2 = node2;
            edgesTested++;
            if (edgesTested > 3 && disjointSetEdges[node1] != -1 && disjointSetEdges[node2] != -1) {
                // If the nodes of the new edge have the same root, adding them forms a cycle
                while (root1 >= 0) root1 = disjointSetEdges[root1];
                while (root2 >= 0) root2 = disjointSetEdges[root2];
                if (root1 != root2 || (root1 == -1 && root2 == -1)) {
                    // If the root of node1 has a larger tree (more negative number in disjointSetEdges) than the root of node2,
                    if (root2 < root1) {
                        disjointSetEdges[root2] += disjointSetEdges[root1];
                        disjointSetEdges[root1] = root2;
                    }
                    else {
                        disjointSetEdges[root1] += disjointSetEdges[root2];
                        disjointSetEdges[root2] = root1;
                    }
                }
                MSTWeight += currWeight;
                nodesLeft--;
            }
            edgeTop--;
        }
        return std::pair<double, int>( MSTWeight, edgesTested );
    }
    else throw exception();
}

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
    // Your implementation
    
    return out;
}

#endif
