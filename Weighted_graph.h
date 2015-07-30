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
        int edgesTested = 0, nodesLeft = numVertices;
        double MSTWeight = 0.0;
        
        double** weightsCopy = new double*[numVertices];
        for (int i = 0; i < numVertices; i++)
            weightsCopy[i] = new double[numVertices];
        for (int i = 0; i < numVertices; i++)
            for (int j = 0; j < numVertices; j++)
                weightsCopy[i][j] = weights[i][j];
        
        while (nodesLeft > 0) {
            // Find the smallest weighted edge using the weight matrix
            double minWeight = 0.0;
            int v1 = -1, v2 = -2;
            for (int i = 0; i < numVertices; i++) {
                for (int j = 0; j < numVertices; j++) {
                    if (weightsCopy[i][j] != INF) {
                        weightsCopy[j][i] = INF;
                        if (minWeight == 0.0) {
                            minWeight = weights[i][j];
                            v1 = i; v2 = j;
                        }
                        else if (weights[i][j] < minWeight) {
                            minWeight = weights[i][j];
                            v1 = i; v2 = j;
                        }
                        weightsCopy[i][j] = INF;
                    }
                }
            }
            
            
            
            // Test if the edge creates a cycle
            bool done = false;
            int DFS_stack1[numEdges], DFS_stack2[numEdges];
            DFS_stack1[0] = v1; DFS_stack2[0] = v2;
            int top_of_stack = 0;
            
            while (!done && nodesLeft > 0) {
                int curr = DFS_stack[top_of_stack];
                int vertices[numVertices];
                int vtop = 0;
                for (int i = 0; i < numVertices; i++)
                    vertices[i] = -1;
                double vweights[numVertices];
                for (int i = 0; i < numVertices; i++)
                    vweights[i] = 0.0;
                
                
                // Find adjacent nodes
                for (int i = 0; i < numVertices; i++) {
                    if (weights[curr][i] != INF) {
                        vertices[vtop] = i;
                        vweights[vtop] = weights[curr][i];
                        vtop++;
                    }
                }
                
                // Find incident edge with least weight
                
                
                
            }
            
            // Increment edgesTested and test the obtained minimum edge to see if it creates a cycle
            edgesTested++;
            bool createsCycle = true;
            
            // If the edge does not create a cycle, add the edge to the MST and decrement nodesLeft
            if (!createsCycle) {
                
                nodesLeft--;
            }
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
