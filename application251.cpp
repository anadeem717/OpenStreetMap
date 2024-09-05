#include "application.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm> // for dijkstra

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// for priority queue
class prioritize {
   public:
    bool operator()(const pair<long long, double>& p1,
                    const pair<long long, double>& p2) const {
        return p1.second > p2.second;
    }
};


double INF = numeric_limits<double>::max();


// Builds a graph using Nodes (vertices)
// and computes edges by checking <Footways>
// and footways close to buildings
graph<long long, double> buildGraph(
    const map<long long, Coordinates>& Nodes,
    const vector<FootwayInfo>& Footways,
    const vector<BuildingInfo>& Buildings) {
    graph<long long, double> G;

    // add all the vertices (nodes)
    for (auto const& e : Nodes) {
        G.addVertex(e.first);
    }

    // add Footway edges to the graph
    for (auto const& e : Footways) {
        for (size_t i = 0; i < e.Nodes.size() - 1; i++) {

            auto from = e.Nodes.at(i);   // start node
            auto to = e.Nodes.at(i + 1); // end node

            // find the weight/distance between from and to
            auto weight = distBetween2Points(Nodes.at(from).Lat, Nodes.at(from).Lon,
                                             Nodes.at(to).Lat, Nodes.at(to).Lon);
            
            // add the edges (both ways) to the graph
            G.addEdge(from, to, weight); 
            G.addEdge(to, from, weight);
        }
        
    }

    // next add vertices and edges for buildings
    for (auto const& e : Buildings) {
        G.addVertex(e.Coords.ID);

        // iterate through the Nodes to see any possible
        // ways connecting to building
        for (const auto& i : Nodes) {
            auto distance = distBetween2Points(e.Coords.Lat, e.Coords.Lon,
                                             i.second.Lat, i.second.Lon);

            // Edge created must be within 0.041 miles of a footway
            if (i.second.OnFootway && distance <= 0.041) {
                G.addEdge(e.Coords.ID, i.second.ID, distance);
                G.addEdge(i.second.ID, e.Coords.ID, distance); 
            }
        }
    }

    return G;
}


// Computes the shortest possible path from
// start to target, using the dijkstra algorithm
vector<long long> dijkstra(
    const graph<long long, double>& G,
    long long start,
    long long target,
    const set<long long>& ignoreNodes) {
    vector<long long> path;

    priority_queue<pair<long long, double>, 
                   vector<pair<long long, double>>, 
                   prioritize> 
        worklist;

    // kinda like a copy of the adjList of the graph
    // but with INF for weights
    unordered_map<long long, double> vertexList;

    // map of a vertex->vertex
    // stores predecessors to vertices
    unordered_map<long long, long long> predecessor; 

    // edge case, already at target
    if (start == target) { 
        path.push_back(start);
        return path;
    }

    // emplace each vertex in both maps
    // and add each to priority queue
    for (auto const& v : G.getVertices()) {
        vertexList.emplace(v, INF);
        predecessor.emplace(v,-1);
        worklist.push(make_pair(v, INF));
    }

    // first "visit" of algorithm
    vertexList.at(start) = 0;
    worklist.push(make_pair(start, 0));

    // generate the minimum spanning tree
    // after this, we can use it to find shortest path
    while (!worklist.empty()) {
        auto currV = worklist.top();
        worklist.pop(); 

        // if currV in ingnoreNodes, continue 
        // only if it is not target/start
        if ((ignoreNodes.count(currV.first) > 0) 
             && (currV.first != target) 
             && (currV.first != start)) {
            continue;
        }

        // iterate through neighbors and find alternative paths
        for (auto const& neighbor : G.neighbors(currV.first)) {
            double edgeWeight;
            if (!(G.getWeight(currV.first, neighbor, edgeWeight))) {
                continue;
            }

            double altPathDistance = currV.second + edgeWeight; 

            // if new path is shorter, change the values for
            // <predecessor> and <vertexList> maps
            // add to priority queue
            if (altPathDistance < vertexList.at(neighbor)) {
                vertexList.at(neighbor) = altPathDistance;
                predecessor.at(neighbor) = currV.first;
                worklist.push(make_pair(neighbor, altPathDistance));
            }
        }
    }

    auto pred = target; // start from target, work backwards

    // use the map of predecessors to build path
    while (pred != -1) {
        path.push_back(pred);
        pred = predecessor.at(pred);
    }

    // right now the path is target->start,
    // so we need to flip it
    reverse(path.begin(), path.end());

    // if target unreachable from start 
    if (path.back() != target || path.front() != start) {
        path.clear();
    }

    return path;
}

double pathLength(const graph<long long, double>& G, const vector<long long>& path) {
    double length = 0.0;
    double weight;
    for (size_t i = 0; i + 1 < path.size(); i++) {
        bool res = G.getWeight(path.at(i), path.at(i + 1), weight);
        assert(res);
        length += weight;
    }
    return length;
}

void outputPath(const vector<long long>& path) {
    for (size_t i = 0; i < path.size(); i++) {
        cout << path.at(i);
        if (i != path.size() - 1) {
            cout << "->";
        }
    }
    cout << endl;
}

void application(
    const vector<BuildingInfo>& Buildings,
    const graph<long long, double>& G) {
    string person1Building, person2Building;

    set<long long> buildingNodes;
    for (const auto& building : Buildings) {
        buildingNodes.insert(building.Coords.ID);
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
        cout << "Enter person 2's building (partial name or abbreviation)> ";
        getline(cin, person2Building);

        //
        // find the building coordinates
        //
        bool foundP1 = false;
        bool foundP2 = false;
        Coordinates P1Coords, P2Coords;
        string P1Name, P2Name;

        for (const BuildingInfo& building : Buildings) {
            if (building.Abbrev == person1Building) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (building.Abbrev == person2Building) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        for (const BuildingInfo& building : Buildings) {
            if (!foundP1 &&
                building.Fullname.find(person1Building) != string::npos) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (!foundP2 && building.Fullname.find(person2Building) != string::npos) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        if (!foundP1) {
            cout << "Person 1's building not found" << endl;
        } else if (!foundP2) {
            cout << "Person 2's building not found" << endl;
        } else {
            cout << endl;
            cout << "Person 1's point:" << endl;
            cout << " " << P1Name << endl;
            cout << " (" << P1Coords.Lat << ", " << P1Coords.Lon << ")" << endl;
            cout << "Person 2's point:" << endl;
            cout << " " << P2Name << endl;
            cout << " (" << P2Coords.Lat << ", " << P2Coords.Lon << ")" << endl;

            string destName;
            Coordinates destCoords;

            Coordinates centerCoords = centerBetween2Points(
                P1Coords.Lat, P1Coords.Lon, P2Coords.Lat, P2Coords.Lon);

            double minDestDist = numeric_limits<double>::max();

            for (const BuildingInfo& building : Buildings) {
                double dist = distBetween2Points(
                    centerCoords.Lat, centerCoords.Lon,
                    building.Coords.Lat, building.Coords.Lon);
                if (dist < minDestDist) {
                    minDestDist = dist;
                    destCoords = building.Coords;
                    destName = building.Fullname;
                }
            }

            cout << "Destination Building:" << endl;
            cout << " " << destName << endl;
            cout << " (" << destCoords.Lat << ", " << destCoords.Lon << ")" << endl;

            vector<long long> P1Path = dijkstra(G, P1Coords.ID, destCoords.ID, buildingNodes);
            vector<long long> P2Path = dijkstra(G, P2Coords.ID, destCoords.ID, buildingNodes);

            // This should NEVER happen with how the graph is built
            if (P1Path.empty() || P2Path.empty()) {
                cout << endl;
                cout << "At least one person was unable to reach the destination building. Is an edge missing?" << endl;
                cout << endl;
            } else {
                cout << endl;
                cout << "Person 1's distance to dest: " << pathLength(G, P1Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P1Path);
                cout << endl;
                cout << "Person 2's distance to dest: " << pathLength(G, P2Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P2Path);
            }
        }

        //
        // another navigation?
        //
        cout << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or #> ";
        getline(cin, person1Building);
    }
}