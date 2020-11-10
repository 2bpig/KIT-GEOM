#ifndef ALGO_DUALIZING_H
#define ALGO_DUALIZING_H
// geom-output and data
#include "View/View_Geom.hpp"
#include "View/View_File.hpp"
// std

// declare constants

class Dualizing {
private:
	GeomOutput& go;
	Polyhedron& polyhedron;
    // parameters
    // dualizing
    Point_3 facet2vertex(Facet_handle& f);
    void facets2vertices(std::vector<Point_3>& dual_points);
public:
    Dualizing(GeomOutput& geomoutput, Polyhedron& polyhedron);
    void convexDual(Polyhedron& dual);
};

#endif
