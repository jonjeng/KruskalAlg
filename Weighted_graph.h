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

class Weighted_graph {
private:
    static const double INF;
    int numVertices;
    int currSize;
    int numEdges;
    int* graphVertices;
    double** weights;
    
    // Do not implement these functions!
    // By making these private and not implementing them, any attempt
    // to make copies or assignments will result in errors
    Weighted_graph( Weighted_graph const & );
    Weighted_graph &operator=( Weighted_graph );
    
    
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
    graphVertices = new int[numVertices];
    for (int i = 0; i < numVertices; i++) {
        graphVertices[i] = -1;
    }
    weights = new double*[numVertices];
    for (int i = 0; i < numVertices; ++i)
        weights[i] = new double[numVertices];
    for (int i = 0; i < numVertices; i++)
        for (int j = 0; j < numVertices; j++)
            weights[i][j] = INF;
}

Weighted_graph::~Weighted_graph() {
    // The destructor frees up any memory allocated by the constructor.
    delete [] graphVertices;
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
    if ((i<0) || (i>=numVertices) || (j<0) || (j>=numVertices) || (w <= 0))
        throw illegal_argument();
    // If i equals j and bothare in the graph, return false (self edges are inconsequential)
    if (i==j && i >= 0 && i < numVertices && j >= 0 && j < numVertices)
        return false;
    // Otherwise, the new edge is safe to insert, as there was either no edge there yet, or an edge of greater weight (parallel edge - but shouldn't we only update the weight if it is lower than the preexisting weight?).
    int top = 0;
    if (graphVertices[numVertices - 1] != -1)
        top = numVertices - 1;
    else {
        for (int n = 0; n < numVertices; n++) {
            if (graphVertices[n] == -1) {
                top = n;
                break;
            }
        }
    }
    bool isNew_i = true, isNew_j = true;
    for (int n = 0; n <= top; n++) {
        if (graphVertices[n] == i) {
            isNew_i = false;
        }
        if (graphVertices[n] == j) {
            isNew_j = false;
        }
    }
    if (isNew_i) {
        currSize++;
        graphVertices[top++] = i;
    }
    if (isNew_j) {
        currSize++;
        graphVertices[top++] = j;
    }
    if (weights[i][j] == INF) numEdges++;
    weights[i][j] = w;
    weights[j][i] = w;
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
    // Use Kruskal's algorithm to find the minimum spanning tree. You will return the weight of the minimum spanning tree and the number of edges that were tested for insertion into the minimum spanning tree.
    int edgesTested = 0, nodesLeft = currSize, edgeNode1[numEdges], edgeNode2[numEdges], edgeTop = 0, disjointSetEdges[currSize];
    double MSTWeight = 0.0, edgeWeights[numEdges];
    
    for (int i = 0; i < numVertices; i++)
        disjointSetEdges[i] = -1;
    
    // Store edges separately
    double** weightsCopy = new double*[currSize];
    for (int i = 0; i < currSize; i++)
        weightsCopy[i] = new double[currSize];
    for (int i = 0; i < currSize; i++) {
        for (int j = 0; j < currSize; j++) {
            weightsCopy[i][j] = weights[i][j];
        }
    }
    for (int i = 0; i < currSize; i++) {
        for (int j = 0; j < currSize; j++) {
            if (weightsCopy[i][j] != INF) {
                edgeNode1[edgeTop] = i;
                edgeNode2[edgeTop] = j;
                edgeWeights[edgeTop] = weightsCopy[i][j];
                edgeTop++;
                weightsCopy[j][i] = INF;
            }
        }
    }
    edgeTop--;
    bool done = false;
    while (!done) {
        // While a spanning tree has not been discovered:
        // Find the edge with the least weight (for such small numbers, selection sort is fine, so brute force searching for the minimum is adequate
        double currWeight = edgeWeights[0];
        int currEdgeNode1 = edgeNode1[0], currEdgeNode2 = edgeNode2[0], currEdge = 0;
        for (int i = 0; i < numEdges; i++) {
            if (edgeWeights[i] < currWeight) {
                currEdge = i;
                currWeight = edgeWeights[i];
                currEdgeNode1 = edgeNode1[i]; currEdgeNode2 = edgeNode2[i];
            }
        }
        edgeWeights[currEdge] = INF;
        
        
        // Test to see if the edge creates a cycle, using the union-find algorithm
        int /*node1 = edgeNode1[edgeTop], node2 = edgeNode2[edgeTop], */ind_root1 = currEdgeNode1, root1 = disjointSetEdges[ind_root1], ind_root2 = currEdgeNode2, root2 = disjointSetEdges[ind_root2];
        edgesTested++;
        edgeTop--;
        if (disjointSetEdges[currEdgeNode1] == -1) {
            nodesLeft--;
        }
        if (disjointSetEdges[currEdgeNode2] == -1) {
            nodesLeft--;
        }
        
        while (root1 >= 0) {
            ind_root1 = root1;
            root1 = disjointSetEdges[ind_root1];
        }
        while (root2 >= 0) {
            ind_root2 = root2;
            root2 = disjointSetEdges[ind_root2];
        }
        // If the nodes of the edge have different roots, adding them does not form a cycle
        if (root1 != root2 || (ind_root1 != ind_root2 && root1 < 0 && root2 < 0)) {
            // If the root of node1 has a larger tree (more negative number in disjointSetEdges) than the root of node2,
            if (root2 < root1) {
                disjointSetEdges[ind_root2] += disjointSetEdges[ind_root1];
                disjointSetEdges[ind_root1] = ind_root2;
            }
            else {
                disjointSetEdges[ind_root1] += disjointSetEdges[ind_root2];
                disjointSetEdges[ind_root2] = ind_root1;
            }
            
            MSTWeight += currWeight;
        }
        int count = 0;
        for (int z = 0; z < currSize; z++) {
            if (disjointSetEdges[z] < 0) {
                count++;
            }
        }
        if (edgesTested == numEdges || count == 1)
            done = true;
    }
    for (int i = 0; i < currSize; i++)
        delete [] weightsCopy[i];
    delete [] weightsCopy;
    return std::pair<double, int>( MSTWeight, edgesTested );
}

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
    std::cout << "The following is an adjacency matrix representing the graph so far\n (INF indicates that there is no connection between the nodes; a finite value indicates the weight \nof the connection between the nodes):\n";
    for (int i = 0; i < graph.numVertices; i++) {
        for (int j = 0; j < graph.numVertices; j++) {
            std::cout << graph.weights[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    return out;
}

#endif
