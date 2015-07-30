/*****************************************
 * Instructions
 *  - Replace 'uwuserid' with your uWaterloo User ID
 *  - Select the current calendar term and enter the year
 *  - List students with whom you had discussions and who helped you
 *
 * uWaterloo User ID:  uwuserid @uwaterloo.ca
 * Submitted for ECE 250
 * Department of Electrical and Computer Engineering
 * University of Waterloo
 * Calender Term of Submission:  (Winter|Spring|Fall) 201N
 *
 * By submitting this file, I affirm that
 * I am the author of all modifications to
 * the provided code.
 *
 * The following is a list of uWaterloo User IDs of those students
 * I had discussions with in preparing this project:
 *    -
 *
 * The following is a list of uWaterloo User IDs of those students
 * who helped me with this project (describe their help; e.g., debugging):
 *    -
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


class Weighted_graph {
private:
    static const double INF;
    int numVertices;
    int* graphVertices;
    double** weights;
    char* edges;
    
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
numVertices(n) {
    // The constructor takes an argument n which defines a graph with vertices numbered from 0 through n − 1. The vertex i is given an initial priority of i. The default value of n is 10.

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
            if (weights[i][j] != INF) count++;  // But this counts i-j and j-i - how to account for this extraneous count?
    return (numVertices+1);
}

bool Weighted_graph::insert_edge( int i, int j, double w ) {
    // If i or j are outside the range of vertices or if the weight w is not positive, the argument(s) is/are illegal
    if ((i<0) || (i>numVertices) || (j<0) || (j>numVertices) || (w <= 0))
        throw illegal_argument();
    // If i equals j and are in the graph, return false (self edges are inconsequential)
    if (i==j && i >= 0 && i < numVertices && j >= 0 && j < numVertices)
        return false;
    // If an edge between i and j would form a cycle (if there is already a path from i to j, return false
    // BFS from i
    if ()
        return false;
    // Otherwise, either insert a new edge from vertex i to vertex j or, if the edge already exists, update the weight and return true (parallel edge - but shouldn't we only update the weight if it is lower than the preexisting weight?).
    weights[i][j] = w;
    return true;
}

bool Weighted_graph::erase_edge(int i, int j) {
    // If an edge between nodes i and j exists, remove the edge. In this case or if i equals j return true. Otherwise, if no edge exists, return false. If i or j are outside the range 0, ..., n − 1, throw an illegal argument exception.
    return false;
}

void Weighted_graph::clear_edges() {
    // Removes all the edges from the graph.
}

std::pair<double, int> Weighted_graph::minimum_spanning_tree() const {
    // Use Kruskal's algorithm to find the minimum spanning tree. You will return the weight of the minimum spanning tree and the number of edges that were tested for insertion into the minimum spanning tree.
    return std::pair<double, int>( 0.0, 0 );
}

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
    // Your implementation
    
    return out;
}

#endif
