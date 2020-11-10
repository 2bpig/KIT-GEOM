#include "Algo_EdgecuttingHC.hpp"

EdgecuttingHC::EdgecuttingHC(GeomOutput& geomoutput, Polyhedron& polyhedron) :
  go(geomoutput), polyhedron(polyhedron) {
    mode = 0;
    iter_num = 0;
    alpha = 0.5;
}

void EdgecuttingHC::algo() {
    GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, 0, "initial");
    FileIO::write_off_polyhedron(polyhedron);
    // TODO store the queue of polyhedron, for undo/redo
    
    prime();
    
    polyhedron.update_nums();
    GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, 0, "final");
    FileIO::write_off_polyhedron(polyhedron);
}

void EdgecuttingHC::prime() {
    
    Polyhedron_cgal& p = polyhedron.get_p();
	std::size_t iter = 1;
    do {
        std::vector<Point_3> new_points;
        Vertex_iterator v_i = p.vertices_begin();
        Vertex_iterator v_last = p.vertices_end();
        do {
            // find the center in each neighboring facet
            // fix the new vertex with the given ratio in the corresponding neighboring facet
            // find the plane according to the neighboring edge and the two new vertices besides
            // find the new vertex as intersection of these planes
	        new_points.push_back();
        } while ( ++v_i != v_last );
        
            // build polyhedron convex on the new vertices
        build_convex_hull(Polyhedron_cgal& p, std::vector<Point_3>& new_points);
	    
        polyhedron.update_nums();
        GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "p");
        FileIO::write_off_polyhedron(polyhedron);
        
    } while(iter++ < max_iter);
}

void EdgecuttingHC::build_convex_hull(Polyhedron_cgal& p, std::vector<Point_3>& new_points) {
    polyhedron.clear();
	CGAL_precondition(p.is_valid());
	CGAL::convex_hull_3(new_points.begin(), new_points.end(), p);
	CGAL_postcondition(p.is_valid());
}
