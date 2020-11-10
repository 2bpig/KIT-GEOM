#ifndef ALGO_MAXMATCHING_H
#define ALGO_MAXMATCHING_H
// geom-output and data
#include "View/View_Geom.hpp"
#include "View/View_File.hpp"
// std
// boost
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

// declare constants

//using namespace boost;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> my_graph;

class MaxMatching {
private:
	GeomOutput& go;
	Polyhedron& polyhedron;
    // parameters
    // matching
    void extract_vertices(std::vector<Vertex_handle>& vertices);
    void id2vertex(std::vector<Vertex_handle>& vertices, std::size_t id, Vertex_handle& vertex);
    std::size_t vertex2id(std::vector<Vertex_handle>& vertices, Vertex_handle& vertex);
    void extract_graph(std::vector<Vertex_handle>& vertices, my_graph& g);
    bool v2e(Vertex_handle& v1, Vertex_handle& v2, Halfedge_handle& e);
    void pairs2edges(my_graph& g, std::vector<Vertex_handle>& vertices,
      std::vector<boost::graph_traits<my_graph>::vertex_descriptor>& mate,
      std::vector<Halfedge_handle>& M,
      std::vector<Halfedge_handle>& N
//    std::vector<Halfedge_handle>& O
    );
    void max_match_graph_EDMONDS(my_graph& g,
      std::vector<boost::graph_traits<my_graph>::vertex_descriptor>& mate
//    std::vector<Halfedge_handle>& M, 
//    std::vector<Halfedge_handle>& N, 
//    std::vector<Halfedge_handle>& O
    );
    void max_match_graph_RING(
//    std::vector<Halfedge_handle>& M, 
//    std::vector<Halfedge_handle>& N, 
//    std::vector<Halfedge_handle>& O
    );
public:
    MaxMatching(GeomOutput& geomoutput, Polyhedron& polyhedron);
    void match(
      std::vector<Halfedge_handle>& M, 
      std::vector<Halfedge_handle>& N 
    );
};

#endif
