#include "Algo_MaxMatching.hpp"

////////////////////////////////////////////////
//    constructor:                            //
////////////////////////////////////////////////

MaxMatching::MaxMatching(GeomOutput& geomoutput, Polyhedron& polyhedron) :
  go(geomoutput), polyhedron(polyhedron){
  //
}


////////////////////////////////////////////////
//    private:                                //
////////////////////////////////////////////////

void MaxMatching::extract_vertices(std::vector<Vertex_handle>& vertices) {
    Polyhedron_cgal& p = polyhedron.get_p();
    Vertex_iterator v_i = p.vertices_begin();
    Vertex_iterator v_last = p.vertices_end();
    do {
        vertices.push_back(v_i);
    } while ( ++v_i != v_last );
}

void MaxMatching::id2vertex(std::vector<Vertex_handle>& vertices, std::size_t id, Vertex_handle& vertex) {
    vertex = vertices.at(id);
}

std::size_t MaxMatching::vertex2id(std::vector<Vertex_handle>& vertices, Vertex_handle& vertex) {
    std::size_t id = 0;
    for(Vertex_handle v: vertices){
        if(v == vertex) break;
        id++;
    }
    return id;
}

void MaxMatching::extract_graph(std::vector<Vertex_handle>& vertices, my_graph& g) {
    
    Polyhedron_cgal& p = polyhedron.get_p();
    Edge_iterator e_i = p.edges_begin();
    Edge_iterator e_last = p.edges_end();
    Vertex_handle v_1;
    Vertex_handle v_2;
    std::size_t id_1;
    std::size_t id_2;
    do {
        v_1 = e_i->opposite()->vertex();
        v_2 = e_i->vertex();
        id_1 = vertex2id(vertices, v_1);
        id_2 = vertex2id(vertices, v_2);
        boost::add_edge(id_1,id_2,g);
    } while ( ++e_i != e_last );
}

bool MaxMatching::v2e(Vertex_handle& v1, Vertex_handle& v2, Halfedge_handle& e) {
    Polyhedron_cgal& p = polyhedron.get_p();
	HV_circulator h_v2 = v2->vertex_begin();
	HV_circulator h_v2_begin = v2->vertex_begin();
	Halfedge_handle v1_v2;
    do {
        if(h_v2->opposite()->vertex() == v1) {
            e = h_v2;
            return true; 
        }
    } while ( ++h_v2 != h_v2_begin);
    return false;
}

void MaxMatching::pairs2edges(my_graph& g, std::vector<Vertex_handle>& vertices,
  std::vector<boost::graph_traits<my_graph>::vertex_descriptor>& mate,
  std::vector<Halfedge_handle>& M,
  std::vector<Halfedge_handle>& N
//  std::vector<Halfedge_handle>& O
  ) {
//    Vertex_handle v_1;
//    Vertex_handle v_2;
//    Halfedge_handle e;
    
//    boost::graph_traits<my_graph>::vertex_iterator vi, vi_end;
//    for(boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
//        if (mate[*vi] != boost::graph_traits<my_graph>::null_vertex() && *vi < mate[*vi]) {
//            std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;
//            id2vertex(vertices, *vi, v_1);
//            id2vertex(vertices, mate[*vi], v_2);
//            v2e(v_1, v_2, e);
//            M.push_back(e);
//        }
//    }
    
    Polyhedron_cgal& p = polyhedron.get_p();
    Edge_iterator e_i = p.edges_begin();
    Edge_iterator e_last = p.edges_end();
    Vertex_handle v_1;
    Vertex_handle v_2;
    std::size_t id_1;
    std::size_t id_2;
    do {
        v_1 = e_i->opposite()->vertex();
        v_2 = e_i->vertex();
        id_1 = vertex2id(vertices, v_1);
        id_2 = vertex2id(vertices, v_2);
        if(mate[id_1] == id_2) {
            M.push_back(e_i);
        } else {
            N.push_back(e_i);
        }
    } while ( ++e_i != e_last );
}

void MaxMatching::max_match_graph_EDMONDS(my_graph& g,
  std::vector<boost::graph_traits<my_graph>::vertex_descriptor>& mate
//  std::vector<Halfedge_handle>& M, 
//  std::vector<Halfedge_handle>& N, 
//  std::vector<Halfedge_handle>& O
  ) {
    
//    std::vector<boost::graph_traits<my_graph>::vertex_descriptor> mate(polyhedron.get_edge_num());
    bool success = checked_edmonds_maximum_cardinality_matching(g, &mate[0]);
    assert(success);
    
}

void MaxMatching::max_match_graph_RING(
//  std::vector<Halfedge_handle>& M, 
//  std::vector<Halfedge_handle>& N, 
//  std::vector<Halfedge_handle>& O
  ) {
    
    Polyhedron_cgal& p = polyhedron.get_p();
    Halfedge_handle e = p.edges_begin();
    HF_circulator h_f1 = e->facet_begin();
    HF_circulator h_f1_begin = e->facet_begin();
    // start the traversal from f1
    std::size_t order_e = 0;
	h_f1++;
    do {
        // match the edges around the f1
    	++ order_e;
    } while ( ++h_f1 != h_f1_begin);
    // traversal the external ring one by one
    // split the ring, and forming some tree structure
}

////////////////////////////////////////////////
//    public:                                 //
////////////////////////////////////////////////

void MaxMatching::match(
  std::vector<Halfedge_handle>& M, 
  std::vector<Halfedge_handle>& N
  ) {
    polyhedron.update_nums();
    std::vector<Vertex_handle> vertices;
    extract_vertices(vertices);
    my_graph g(polyhedron.get_edge_num());
    extract_graph(vertices, g);
    
    std::vector<boost::graph_traits<my_graph>::vertex_descriptor> mate(polyhedron.get_edge_num());
    max_match_graph_EDMONDS(g, mate);
    pairs2edges(g, vertices, mate, M, N);
}
