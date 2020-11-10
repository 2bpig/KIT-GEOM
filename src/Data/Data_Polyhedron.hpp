#ifndef DATA_POLYHEDRON_H
#define DATA_POLYHEDRON_H

#include "Data/Data_Primitive.hpp"
// cgal
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Modifier_base.h>

typedef CGAL::Polyhedron_traits_with_normals_3<gvKernel> 	normTraits;
typedef CGAL::Polyhedron_3<normTraits>                     	Polyhedron_cgal;
typedef Polyhedron_cgal::Point_3				            Point;
typedef Polyhedron_cgal::Vertex                             Vertex;
typedef Polyhedron_cgal::Vertex_iterator                    Vertex_iterator;
typedef Polyhedron_cgal::Vertex_handle		                Vertex_handle;
typedef Polyhedron_cgal::Halfedge_iterator                  Halfedge_iterator;
typedef Polyhedron_cgal::Halfedge_handle		            Halfedge_handle;
typedef Polyhedron_cgal::Edge_iterator                      Edge_iterator;
typedef Polyhedron_cgal::Facet_iterator                     Facet_iterator;
typedef Polyhedron_cgal::Facet_handle		                Facet_handle;
typedef Polyhedron_cgal::HalfedgeDS                         HalfedgeDS;
typedef Polyhedron_cgal::Halfedge_around_vertex_circulator 	        HV_circulator;
typedef Polyhedron_cgal::Halfedge_around_vertex_const_circulator 	HV_const_circulator;
typedef Polyhedron_cgal::Halfedge_around_facet_circulator           HF_circulator;
typedef Polyhedron_cgal::Halfedge_around_facet_const_circulator     HF_const_circulator;


struct Hedge_cmp{
    bool operator()(Halfedge_handle e1, Halfedge_handle e2) const{
        return &*e1 < &*e2;
    }
};

struct Facet_cmp{
    bool operator()(Facet_handle f1, Facet_handle f2) const{
        return &*f1 < &*f2;
    }
};

struct Normal_vector {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        return CGAL::unit_normal(
          h->next()->vertex()->point(),
          h->next()->next()->vertex()->point(),
		  h->vertex()->point());
    }
};

class Polyhedron {
private:
    Polyhedron_cgal p;
    std::size_t vertex_num;
    std::size_t edge_num;
    std::size_t facet_num;
    
    void trans_to_origin();
    
public:
    Polyhedron();
    // init polyhedron
    void clear();
    Halfedge_handle init_tetrahedron();
    Halfedge_handle init_cube();
    Halfedge_handle init_octahedron();
    Halfedge_handle init_icosahedron();
    Halfedge_handle init_random_polyhedron();
    // parameter of polyhedron
    bool update_nums();
    std::size_t get_vertex_num() { return vertex_num; }
    std::size_t get_edge_num() { return edge_num; }
    std::size_t get_facet_num() { return facet_num; }
    Polyhedron_cgal& get_p() { return p; }
    // modify polyhedron
    void split_facet_in_center();
    void cal_facets_normal();
    //static
    static Point_3 cal_facet_center(Facet_iterator& f);
};

#endif
