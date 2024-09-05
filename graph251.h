#pragma once

#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

using namespace std;

/// @brief Simple directed graph using an adjacency list.
/// @tparam VertexT vertex type
/// @tparam WeightT edge weight type
template <typename VertexT, typename WeightT>
class graph {
   private:
    size_t numEdges;
    unordered_map<VertexT, map<VertexT, WeightT>> adjacencyList;

   public:
    /// Default constructor
    graph() {
        numEdges = 0;
    }

    /// @brief Add the vertex `v` to the graph, must run in at most O(log |V|).
    /// @param v
    /// @return true if successfully added; false if it existed already
    bool addVertex(VertexT v) {

        // if already exists
        if (adjacencyList.count(v) != 0) { return false; }

        // add vertex, and return true
        adjacencyList.emplace(v, map<VertexT, WeightT>{});
        return true; 
    
    }

    /// @brief Add or overwrite directed edge in the graph, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight edge weight / label
    /// @return true if successfully added or overwritten;
    ///         false if either vertices isn't in graph
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        
        // at least one of the vertices are not in the graph
        if (adjacencyList.count(from) == 0 || adjacencyList.count(to) == 0 ) { return false; }

        // add edge and weight to the map for <from>
        auto retPair = adjacencyList.at(from).emplace(to, weight);

        if (retPair.second == true) { // new edge
            numEdges++;
        }
        else { // overwrite
            adjacencyList.at(from).at(to) = weight; 
        }
        
        return true; // sucessfully added or overwritten
    }

    /// @brief Maybe get the weight associated with a given edge, must run in at most O(log |V|).
    /// @param from starting vertex
    /// @param to ending vertex
    /// @param weight output parameter
    /// @return true if the edge exists, and `weight` is set;
    ///         false if the edge does not exist
    bool getWeight(VertexT from, VertexT to, WeightT& weight) const {

        // vertex doesn't exist
        if (adjacencyList.count(from) == 0 || adjacencyList.count(to) == 0 ) { return false; }

        // check if <to> exists as a neighbor for <from>, then get weight
        if (adjacencyList.count(from) > 0 && adjacencyList.at(from).count(to) > 0) {
            weight = adjacencyList.at(from).at(to);
            return true;
        }
        else return false;
    }

    /// @brief Get the out-neighbors of `v`. Must run in at most O(|V|).
    /// @param v
    /// @return vertices that v has an edge to
    set<VertexT> neighbors(VertexT v) const {
        set<VertexT> S;
        
        // for each vertex in the list, 
        // emplace all adjacent vertices in the set
        for (const auto& e : adjacencyList.at(v)) {
            S.emplace(e.first);
        }

        return S;
    }

    /// @brief Return a vector containing all vertices in the graph
    vector<VertexT> getVertices() const {
        vector<VertexT> verts; 

        // for each vertex in the list, add it to the vector
        for (auto const& e : adjacencyList) {
            verts.push_back(e.first);
        }

        return verts;
    }

    /// @brief Get the number of vertices in the graph. Runs in O(1).
    size_t NumVertices() const { return adjacencyList.size(); }

    /// @brief Get the number of directed edges in the graph. Runs in at most O(|V|).
    size_t NumEdges() const { return numEdges; }
};