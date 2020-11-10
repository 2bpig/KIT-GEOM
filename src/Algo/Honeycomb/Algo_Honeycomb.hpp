// geom-output and data
#include "View/View_Geom.hpp"
#include "View/View_File.hpp"
// std
#include <list>
#include <map>

// declare constants
const int max_iter = 8;
const double max_alpha_edge = 0.5;
const double max_alpha_facet = 0.2;
const double max_beta_edge = 0.5; // alpha_edge
const double max_beta_vertex = 0.2; // alpha_facet
const double tolerance = 0.1;
const double error_prop = 0.01;
const double angle_prop = 0.01;

class Honeycomb {
private:
	GeomOutput go;
    /*
                offset                  direction                   section
                
    mode 0 :    uniform     min_dist    uniform     old_normal      uniform     bisection
    mode 1 :    uniform     avg_dist    uniform     old_normal      uniform     bisection
    
    mode 2 :    uniform     avg_dist    uniform     avg_normal      uniform     bisection
    mode 3 :    uniform     avg_dist    weighted    avg_normal      uniform     bisection
    
    mode 4 :    uniform     avg_dist    uniform     old_normal      weighted    section
    */
    int mode;
    int iter_num;
    // parameters
    double alpha_edge;
    double alpha_facet;
    double beta_edge; // alpha_edge
    double beta_vertex; // alpha_facet
    void update_alpha_facet_with_iter();
    void update_alpha_edge_with_iter();
    // honeycomb scheme
    void cal_edge_normal(Halfedge_handle& edge,
      std::map<Halfedge_handle, Vector_3, Hedge_cmp>& hedge_to_normal_map,
      std::map<Facet_handle, Point_3, Facet_cmp>& facet_to_center_map);
    Vector_3 cal_normal_for_lowest_plane_honeycomb(Point_3 p1, Point_3 p2, Point_3 min_p, 
      std::vector<Point_3>& p_list, double error_dist);
    bool cal_new_vertices_new_normal_honeycomb(Vector_3& normal_facet, Point_3 center, 
      HV_circulator& h_v, const HV_circulator& h_v_begin, std::size_t max_order, 
      std::vector<Point_3>& new_vertices, std::vector<Vector_3>& normal_edges, 
      std::vector<Plane_3>& plane_edges);
    void cal_new_vertices_honeycomb(Vector_3& normal_facet, 
      HV_circulator& h_v, const HV_circulator& h_v_begin, std::size_t max_order, 
      std::vector<Point_3>& new_vertices, std::vector<Plane_3>& plane_edges);
    bool pull_facet( Polyhedron_cgal& p, Facet_iterator& f, Point_3 f_center, 
      std::map<Halfedge_handle, Vector_3, Hedge_cmp>& hedge_to_normal_map);
    void pull_facets(Polyhedron_cgal& p);
    // kernel scheme
    void subdivid_edge(Polyhedron_cgal& p, Halfedge_handle& edge);
    Point_3 cal_new_vertex_kernel(Point_3& old_v, 
      HV_circulator& h_v, const HV_circulator& h_v_begin, std::size_t neighbor_size,
      std::vector<Point_3>& neighbor_vertices, Vector_3& push_vector);
    void push_vertex( Polyhedron_cgal& p, Vertex_iterator& v);
    void push_vertices(Polyhedron_cgal& p);
    
public:
    Honeycomb(GeomOutput& geomoutput);
    // honeycomb and kernel
    void honeycomb_iter(Polyhedron& polyhedron, const int iter, const int hc_mode, double alpha);
    void kernel_iter(Polyhedron& polyhedron, const int iter, const int hc_mode, double beta);
    void algo(Polyhedron& polyhedron);
};
