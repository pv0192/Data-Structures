#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>

template <typename vertex>
class weighted_graph {

	private:
	
	std::vector<std::vector<int>> matrix;
	std::vector<vertex> mapping;
	size_t size;
	size_t edges;
	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.

	//The graph_iterator class provides an iterator
	//over the vertices of the graph.
	//This is one of the harder parts, so if you're
	//not too comfortable with C++ leave this for last.
	//If you are, there are many ways of doing this,
	//as long as it passes the tests, it's okay.
	class graph_iterator {

		private:

		int position;
		weighted_graph<vertex> owner;
		
		public:
			graph_iterator(const weighted_graph &);
			graph_iterator(const weighted_graph &, size_t);
			~graph_iterator();
			graph_iterator operator=(const graph_iterator&);
			bool operator==(const graph_iterator&) const;
			bool operator!=(const graph_iterator&) const;
			graph_iterator operator++();
			graph_iterator operator++(int);
			const vertex operator*();
			const vertex* operator->();
	};

	//The neighbour_iterator class provides an iterator
	//over the neighbours of a given vertex. This is
	//probably harder (conceptually) than the graph_iterator.
	//Unless you know how iterators work.
	class neighbour_iterator {

		private:

		int position;
		weighted_graph<vertex> owner;
		vertex ver;

		public:
			neighbour_iterator(const neighbour_iterator&);
			neighbour_iterator(const weighted_graph &, const vertex&);
			neighbour_iterator(const weighted_graph &, const vertex&, size_t);
			~neighbour_iterator();
			neighbour_iterator operator=(const neighbour_iterator& it);
			bool operator==(const neighbour_iterator&) const;
			bool operator!=(const neighbour_iterator&) const;
			neighbour_iterator operator++();
			neighbour_iterator operator++(int);			
			const std::pair<vertex, int> operator*();
			const std::pair<const vertex, int>* operator->();
	};

	public:


	weighted_graph(); //A constructor for weighted_graph. It should start empty. DONE
	~weighted_graph(); //A destructor. Depending on how you do things, this may
					   //not be necessary. DONE
	int get_index(const vertex&) const; //Returns index of vertex passed. DONE
										  
	bool are_adjacent(const vertex&, const vertex&) const; //Returns true if the two vertices are
														   //adjacent, false otherwise. DONE
	bool has_vertex(const vertex&) const; //Returns true if the passed in vertex is 
										  //a vertex of the graph, false otherwise. DONE

	void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges). DONE
	void add_edge(const vertex&, const vertex&, const int&); //Adds an edge between the two vertices
															 //with the given weight (as an int). DONE

	void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges. DONE
	void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists. DONE
	void set_edge_weight(const vertex&, const vertex&, const int&); //Changes the edge weight between the two
																	//vertices to the new weight (the int). DONE

	int get_edge_weight(const vertex&, const vertex&) const; //Returns the weight on the edge between the two vertices. DONE
	int degree(const vertex&) const; //Returns the degree of the vertex. DONE
	int weighted_degree(const vertex&); //Returns the sum of the weights on all the edges incident to the vertex. DONE
	int num_vertices() const; //Returns the total number of vertices in the graph. DONE
	int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight). DONE
	int total_weight(); //Returns the sum of all the edge weights in the graph. DONE

	std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices. DONE
	std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex. DONE
	void test_function(); //Returns index of vertex passed. DONE

	graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
	graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

	neighbour_iterator neighbours_begin(const vertex&); //Returns a neighbour_iterator pointing to the start
														//of the neighbour set for the given vertex.
	neighbour_iterator neighbours_end(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end
													  //of the neighbour set for the given vertex.

	std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they
													//are visited in by a depth-first traversal starting at
													//the given vertex.
	std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they
													  //are visisted in by a breadth-first traversal starting
													  //at the given vertex.

	weighted_graph<vertex> mst(); //Returns a minimum spanning tree of the graph.

};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

// initialise graph iterator
template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g){
	owner = g; 
	position = 0;
}

//initialise graph iterator with given positon
template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos){
	owner = g;
	position = start_pos;
}

//iterator deconstructor
template <typename vertex> weighted_graph<vertex>::graph_iterator::~graph_iterator(){
	
}

//iterator copy constructor, assigns position of given iterator to new iterator
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator=(const graph_iterator& it){ 
	position = it.position;
	return *this;
}

// returns true if two iterators are at the same postion
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator==(const graph_iterator& it) const { 
	return position = it.position; 
}

//returns false if two iterators are not in the same position
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator!=(const graph_iterator& it) const { 
	return position != it.position; 
}

// reduces iterator position by 1 count
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(){ 
	++position;
	return *this;
}

//increments interator position by 1 count
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(int){ 
	position++;
	return *this;
}

 //dereferencing object iterator is pointing to
template <typename vertex> const vertex weighted_graph<vertex>::graph_iterator::operator*(){ 
	return owner.mapping[position]; 
}
 //dereferencing object iterator is pointing to
template <typename vertex> const vertex* weighted_graph<vertex>::graph_iterator::operator->(){ 
	return &owner.mapping[position]; 
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u) {
	owner =g;
	ver = u;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u, size_t start_pos) {
	owner = g;
	ver = u;
	position = start_pos;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::~neighbour_iterator() {}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator=(const neighbour_iterator& it) { 
	position = it.position;
	ver = it.ver;
	return *this; 
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator==(const neighbour_iterator& it) const { 
	return position = it.position; 
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator!=(const neighbour_iterator& it) const { 
	return false;
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++() { auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++(int){ auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }			

template <typename vertex> const std::pair<vertex, int> weighted_graph<vertex>::neighbour_iterator::operator*(){ auto p = std::pair<vertex,int>(); return p; }

template <typename vertex> const std::pair<const vertex, int>* weighted_graph<vertex>::neighbour_iterator::operator->(){ return nullptr; }

// returns iterator pointing to the beginning of vertex vector
template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::begin() {
	return graph_iterator(*this);
}
 //Returns iterator pointing to the end of the vertex vector
template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::end() {
	return graph_iterator(*this, mapping.size());
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_begin(const vertex& u) {
	return neighbour_iterator(*this, vertex());
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_end(const vertex& u) {
	return neighbour_iterator(weighted_graph<vertex>(), vertex());
}

//initialise weighted graph
template <typename vertex> weighted_graph<vertex>::weighted_graph(){
	size = 0;
	edges = 0;
}

template <typename vertex> weighted_graph<vertex>::~weighted_graph(){

}

//function to retrieve index of given vertex
template <typename vertex> int weighted_graph<vertex>::get_index(const vertex& u) const {
	for(int ii = 0; ii<mapping.size(); ii++){	//loop through mapping vector, if vector is found, return vector index at which vertex is found.
		if(mapping.at(ii)==u){
			return ii;
		}
	}
	return -1;
}

//returns true if vertex esists in graph
template <typename vertex> bool weighted_graph<vertex>::has_vertex(const vertex& u) const {
	for(int i = 0; i<mapping.size(); i++){		//loop through mapping vector, if vertex found, return true
		if(mapping.at(i)==u){
			return true;
		}
	}
	return false;
}
//function to return true if given vertices are adjacent
template <typename vertex> bool weighted_graph<vertex>::are_adjacent(const vertex& u, const vertex& v) const {
	bool result = false;						//set initial value to false
	int i = get_index(u);       				//get index of given vertex
	int j = get_index(v);
	if(matrix[i][j] && matrix[j][i] > 0){		//check to see if both directions of edges are >0, if ture, vertices are adjacent.
		result = true;
	}
	return result;
}

//Add new vertex to weighted graph
template <typename vertex> void weighted_graph<vertex>::add_vertex(const vertex& v) {
	mapping.push_back(v);
	matrix.resize(mapping.size());
	for(int i = 0; i<mapping.size(); i++){
		matrix[i].resize(mapping.size());
	}

}

// add edge between existing vertices
template <typename vertex> void weighted_graph<vertex>::add_edge(const vertex& u, const vertex& v, const int& weight) {
	int i = get_index(u); 						//get index of vertex u
	int j = get_index(v);						//get index of vertex v
	matrix[i][j] = weight; 						//set edge weight
	matrix[j][i] = weight;						//as graph is non directional, edge weights need to be set in both directions
	
}

template <typename vertex> void weighted_graph<vertex>::test_function(){
	
}
	
//remove vertex from weighted graph
template <typename vertex> void weighted_graph<vertex>::remove_vertex(const vertex& u) {
	int pos = get_index(u); 					// get index of vertex in mapping vector
	mapping.erase(mapping.begin()+pos); 		//erase vector from mapping vector
	matrix.erase(matrix.begin()+pos); 			//erase edge records from adjacency matrix 
}

//remove edge between vertices
template <typename vertex> void weighted_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	int i = get_index(u); 						//get index
	int j = get_index(v); 						//get index
	matrix[i][j] = 0; 							//set edge to 0
	matrix[j][i] = 0; 							// as graph is non directional, edge weights need to be set to 0 in both directions
}

//set edge weight between vertex u and v with weight.
template <typename vertex> void weighted_graph<vertex>::set_edge_weight(const vertex& u, const vertex& v, const int& weight) {
	int i = get_index(u); 						//get index
	int j = get_index(v); 						//get index
	matrix[i][j] = weight;						// setting edge weight in both directions
	matrix[j][i] = weight; 						// setting edge weight in both directions
}

//retrieve edge weight between vertex u and v.
template <typename vertex> int weighted_graph<vertex>::get_edge_weight(const vertex& u, const vertex& v) const {
	int i = get_index(u); 						//get index
	int j = get_index(v); 						//get index
	return matrix[i][j]; 						// as graph is non directional, returning edge weight of one direction is enough
}

//Find degree of given vertex = number of edges sprouting from given vertex
template <typename vertex> int weighted_graph<vertex>::degree(const vertex& u) const {
	int index = get_index(u);  					// get index of vertex
	int vertices = num_vertices(); 				// find total number of vertices
	int count = 0;
	//loop thorugh all edges of vertex u , if edge > 0, this implies edge exists and counter increments
	for(int i=0;i<vertices;i++){
		if(matrix[index][i] > 0){
			count++;
		}
	}
	return count; //return final count
}

//caluclate weighted degree of a vertex. Add up weight of each edge sprouting from given vertex u. 
template <typename vertex> int weighted_graph<vertex>::weighted_degree(const vertex& u) {
	int index = get_index(u);  					// get index of vertex
	int vertices = num_vertices(); 				// total no. of vertices in graph
	int weight = 0;                             //starting weight 0
												//loop through all edges of give vertex u, if edge weight >0, add to weight degree count.
	for(int i=0 ;i <= vertices; i++){
		if(matrix[index][i] > 0){
			weight += get_edge_weight(index,i);
		}
	}
	return weight;	
}

//Returns total number of verties in graph
template <typename vertex> int weighted_graph<vertex>::num_vertices() const {
	return mapping.size(); 						// return size of vector tracking vertices
}

//Returns total number of edges in adjacency matrix 
template <typename vertex> int weighted_graph<vertex>::num_edges() const {
	int number_edges = 0; 						//start count at 0
												//loop through all edge weight instances, if >0, edge exists and increment counter.
	for (int i = 0; i<mapping.size() ; i++){
		for(int j=0; j<mapping.size(); j++){
			if(matrix[i][j] > 0){
				number_edges++;
			}
		}
	}
	
	return number_edges/2;         				// since graph is non directional, no of edges/2 gives correct number of edges.
}

//Returns total weigh of graph
template <typename vertex> int weighted_graph<vertex>::total_weight() {
	int vertices = num_vertices(); //total number of vertices
	int count = 0;                 //count starts at 0
	                                            
	for(int i = 0; i< vertices; i++){			//loop through all edges, if edge weight >0, edge exists and add to total weight.
		for(int j = 0; j < vertices;j++){
			count += get_edge_weight(i,j); 
		}
	}
	return count/2;             				// since graph is non directional, total weight/2 gives correct total weight.
}

//return vector of vertices
template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_vertices() {
	return mapping;             				//simply return vector that stores all vertices, in our case napped mapping.
}

												//function that returns neighbours of a given vertex u
template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_neighbours(const vertex& u) {
	std::vector<vertex> neighbours; 			// create vector of vertices to store neighbours
	int ver = get_index(u); 					//get index
	int total = mapping.size(); 				// total number of vertices
												 
	for (int i = 0; i < total; i++ ){			//loop thorugh all edges of u, if weight > 0, vertex is a neighbour, add to neighbour vector.
		if (get_edge_weight(ver,i) > 0){
			neighbours.push_back(i);
		}
	}
	return neighbours; //return vector of neighbours
}

//depth first traversal
//This algorithm works in two stages – in the first stage the visited vertices are pushed onto the stack and later on when there is no vertex further to visit those are popped-off.
template <typename vertex> std::vector<vertex> weighted_graph<vertex>::depth_first(const vertex& start_vertex){
	
	int start = get_index(start_vertex); //get index of given starting vertex
	bool visited[mapping.size()];        //create array of booleans to track if vertex has been visited in the size of total number of vertex
	//set all values to false
	for (unsigned i = 0; i < mapping.size(); i++){
		visited[i] = false;
	}
	
	std::stack<int> unprocessed;         //stack of unprocessed vertices, as we need FILO
	unprocessed.push(start);             //add first given vertex to unprocessed
	
	std::vector<int> ordered;            //create vector for ordered traversal
	
	while (!unprocessed.empty()){        //while unprocessed stack in not empty,continue loop
		int n = unprocessed.top();       //assign top most vertex to n
		unprocessed.pop();               //remove top most element from stack
		if (!visited[n]){				 // if vertex n is not visited, add to ordered vector
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = mapping.size(); i != 0; i--){			//add rest of vertices to unprocessed to be processed by loop in reverse order. 
				if (matrix[n][i-1]){
					unprocessed.push(i-1);
				}
			}
		}
	}
		
	return ordered; 
}
	

//breadth first traversal
template <typename vertex> std::vector<vertex> weighted_graph<vertex>::breadth_first(const vertex& start_vertex){
	
	bool visited[mapping.size()];
	int start = get_index(start_vertex);
	for (unsigned i = 0; i < mapping.size(); i++){
		visited[i] = false;
	}
	
	std::queue<int> unprocessed; 			// queue of non visited vertices, as we need FIFO.
	unprocessed.push(start);
	
	std::vector<int> ordered; 				// vector where the bft will be stored
	
	//loop runs as long as unprocessed ques is not empty
	while (!unprocessed.empty()){
		int n = unprocessed.front(); 		// assign first in queue in unprocessed to n
		unprocessed.pop(); 					//remove that element from the queue
		if (!visited[n]){               	//check to see if vertex n has been visited, if not mark as visited, add to ordered vector
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = 0; i < mapping.size(); i++){            // add the rest of vertices to unprocessed queue.
				if (matrix[n][i]){
					unprocessed.push(i);
				}
			}
		}
	}
		
	
	return ordered;	
}

template <typename vertex>	weighted_graph<vertex> weighted_graph<vertex>::mst() {
	weighted_graph<vertex> mst;
	/*mst.add_vertex(mapping[0]);
	int position = 0;
	std::vector<vertex> neighbours = get_neighbours(mapping[position]);
	int graph_size = mapping.size();
	
	
     //while(mst.num_vertices() < graph_size){
	//for(int j=0; j<graph_size;j++){
		for(int i =0; i < neighbours.size();i++){
			std::vector<vertex> neighbours = get_neighbours(mapping[position]);
			if(!mst.has_vertex(neighbours[i])){
				mst.add_vertex(neighbours[i]);
				mst.add_edge(mapping[position],neighbours[i],get_edge_weight(mapping[position],neighbours[i]));
			}
		}
		
		position ++;
		std::cout<< "current mapping index" << position;
	 
	
	
	std::cout<< mst.mapping[4];*/
	return mst;
	//return weighted_graph<vertex>();
}


#endif
