#ifndef GRAPH_ALGS
#define GRAPH_ALGS

#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <algorithm>
#include <limits>
#include "weighted_graph.hpp"
#include "easy_weighted_graph_algorithms.cpp"

//Returns true if the graph is connected, false otherwise.
template <typename vertex>
bool is_connected(const weighted_graph<vertex>& g){
	
	std:: vector<vertex> bft_vertices;
	bool is_connected = true;
	
	if(g.num_vertices () > 0){
		bft_vertices = breadth_first(g, *g.begin());
		//if no. of vertices in bfs == g, then true
		if (bft_vertices.size() == g.num_vertices() ){
			is_connected = true;
			}
		else {
			is_connected = false;
		}
	}
	return is_connected;
}

//Returns a vector of weighted graphs, where each weighted graph is a connected component of the input graph.
template <typename vertex>
std::vector<weighted_graph<vertex>> connected_components(const weighted_graph<vertex>& g){
	
	std::vector<weighted_graph<vertex>> connected_components;
	std::map<vertex, bool> vis;
	std::vector<vertex> dft;
	
	for (auto itr = g.begin(); itr != g.end(); itr++){																			
		//initialise all vis entties as false
				vis[*itr] = false;
	}
	
	//for all vertices of the given graph
	for(auto itr = g.begin(); itr != g.end(); itr++ ){
		//process only if vertex is marked as unvisited
		if (!vis[*itr]) {																								
			weighted_graph<vertex> temp;
			//bft on given graph with any vertex
			dft = depth_first(g, *itr );
			for(int j = 0; j < dft.size(); j++){
				//add dft vertices to temp component graph(temp), push vertex edges to temp component graph
				temp.add_vertex(bft[j]);	
				vis[dft[j]] = true;	     																																	//mark all vertices in bft as visited
				for (auto n_it = g.cneighbours_begin(dft[j]); n_it != g.cneighbours_end(dft[j]); ++n_it){ 					
					temp.add_edge(dft[j],n_it->first , n_it->second);
				}		
			}
			connected_components.push_back(temp);																		
		}	
	}
	return connected_components;	
}

//Returns a map of the vertices of the weighted graph g and their distances from
//the given starting vertex v.
template <typename vertex> 
std::map<vertex, int> dijkstras(const weighted_graph<vertex>& g, const vertex& v){
	
	std::map<vertex, int> dijkstras;
	std::map<vertex, int> dist;
	std::map<vertex, bool> processed;
		
	//initialise dist as infinite and p as false
	for(auto itr = g.begin(); itr != g.end(); itr++ ){
		dist[*itr] = __INT_MAX__; 
		processed[*itr] = false;
	}
			
	//set distance of source vertex as 0
	dist[v] = 0;
	//find shortest path of all vertices
	for(auto itr = g.begin(); itr != g.end(); itr++){
		//pick minimum dist vertex
		int min = __INT_MAX__; 
		vertex u;
  		for (auto itr = g.begin(); itr != g.end(); itr++ ){
    		if (processed[*itr] == false && dist[*itr] <= min){
         		min = dist[*itr], u = *itr;
			}
		}
		//mark the picked vertex as processed
		processed[u] = true;
		//update dist value of adjacent vertices of picked vertex
		for(auto itr = g.begin(); itr != g.end(); itr++ ){
			if(!processed[*itr] && g.are_adjacent(u,*itr) && dist[u] != __INT_MAX__ && dist[u]+g.get_edge_weight(u,*itr) < dist[*itr]){
				dist[*itr] = dist[u] + g.get_edge_weight(u,*itr);
			}
		}
	}
	for(auto itr = g.begin(); itr != g.end(); itr++ ){
		dijkstras.insert(std::pair<vertex,int>(*itr,dist[*itr]));
	}
	return dijkstras;	
}


//Returns a vector containing all the articulation points of the
//input weighted graph g.
template <typename vertex>
std::vector<vertex> articulation_points(const weighted_graph<vertex>& g){

	std::vector<vertex> ap;
	
	//loop thorugh all vertices, remove current from temp graph. test temp graph to see if connected. 
	//If !is_connected(temp), removed vertex is an articulation point.
	for(auto itr = g.begin(); itr != g.end(); itr++ ){
		auto temp = g;
		temp.remove_vertex(*itr);
		if (!is_connected(temp) ){
			ap.push_back(*itr);
		}
	}
	return ap;
}

#endif
