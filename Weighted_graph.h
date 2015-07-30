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
    numEdges++;
    return true;
}

bool Weighted_graph::erase_edge(int i, int j) {
    // If i and j are outside the legal range, they are invalid - throw an exception
    if (i < 0 || i > numVertices || j < 0 || j > numVertices) throw illegal_argument();
    
    // If a specified edge exists, remove it (set its weight to INF) and return true.
    if (weights[i][j] != INF || i == j) {
        weights[i][j] = INF;
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
    // Use Kruskal's algorithm to find the minimum spanning tree. You will return the weight of the minimum spanning tree and the number of edges that were tested for insertion into the minimum spanning tree.
    
    // If there are no edges or only one edge in the graph, insertion will not create a cycle and so is legal

    // DFS from i
    Dynamic_queue<int> BFSqueue(numEdges + 1);
    int DFS_stack[numEdges+1];
    int top_of_stack = 0;
    int startingNode = -1;
    bool loopBreak = false;
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            if (weights[i][j] != INF) {
                startingNode = i;
                loopBreak = true;
                break;
            }
        }
        if (loopBreak) break;
    }

    return std::pair<double, int>( 0.0, 0 );
}

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
    // Your implementation
    
    return out;
}

#endif
