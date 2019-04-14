#include <iostream>
#include <vector>

#include "weighted_graph.hpp"

int main() {

	weighted_graph<int> g;
g.add_vertex(0);
g.add_vertex(1);
g.add_vertex(2);
g.add_vertex(3);
g.add_vertex(4);
g.add_vertex(5);
g.add_vertex(6);
g.add_vertex(7);

g.add_edge(0,1,5);
g.add_edge(0,2,3);
g.add_edge(1,4,2);
g.add_edge(2,5,8);
g.add_edge(3,6,1);
g.add_edge(5,3,9);
g.add_edge(2,3,4);
g.add_edge(1,6,9);
g.add_edge(3,4,2);
g.add_edge(5,7,4);
g.add_edge(3,7,3);
g.add_edge(4,7,4);

		weighted_graph<int> t = g.mst();
	

		//TS_ASSERT_EQUALS(t.num_edges(), t.num_vertices() - 1);
		//TS_ASSERT_EQUALS(t.total_weight(), t.num_vertices() - 1);
		//TS_ASSERT(is_connected(t));
	
	
	
}
