#include "Algo_Dualizing.hpp"

////////////////////////////////////////////////
//    constructor:                            //
////////////////////////////////////////////////

Dualizing::Dualizing(GeomOutput& geomoutput, Polyhedron& polyhedron) :
  go(geomoutput), polyhedron(polyhedron){
  //
}


////////////////////////////////////////////////
//    private:                                //
////////////////////////////////////////////////

Point_3 Dualizing::facet2vertex(Facet_handle& f) {
    Vector_3 n = f->plane();
    Point_3 p = f->halfedge()->vertex()->point();
    double f_dist = (p-origin_point)*n;
    if(f_dist <= 0.0) {
        std::cout<< "Unexpected facet through origin ! " <<std::endl;
        return Point_3(inf_double, inf_double, inf_double);
    }
    double v_dist = 110.0/f_dist;
    Vector_3 v = -v_dist*n;
    return Point_3(v.x(), v.y(), v.z());
}

void Dualizing::facets2vertices(std::vector<Point_3>& dual_points) {
    polyhedron.cal_facets_normal();
    Polyhedron_cgal& p = polyhedron.get_p();
    Facet_iterator f_i = p.facets_begin();
    Facet_iterator f_last = p.facets_end();
    do {
	    dual_points.push_back(facet2vertex(f_i));
    } while ( ++f_i != f_last );
}

////////////////////////////////////////////////
//    public:                                 //
////////////////////////////////////////////////

void Dualizing::convexDual(Polyhedron& dual) {
    Polyhedron_cgal& dual_p = dual.get_p();
	CGAL_precondition(dual_p.is_valid());
	std::vector<Point_3> dual_points;
	facets2vertices(dual_points);
	CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), dual_p);
	CGAL_postcondition(dual_p.is_valid());
}
