#ifndef VIEW_GEOM_H
#define VIEW_GEOM_H

#include <iostream>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#ifndef CGAL_USE_GEOMVIEW
// not use geomview
class GeomOutput {
private:
    bool cgal_use_geoview;
    int id;
    
public:
    // get usability immediately after constructing new object
    GeomOutput(int id) : cgal_use_geoview(false) {
    this->id = id;
	  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
    }
    bool getUsability() { return cgal_use_geoview; }

};
#else
// use geomview
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include "Data/Data_Polyhedron.hpp"
// std
#include <unistd.h> // for sleep()

// declare constants
const int gray_fac = 220;
const CGAL::Color color_gray = CGAL::Color(gray_fac,gray_fac,gray_fac);
const CGAL::Color color0 = CGAL::Color(255, 255, 255);
const CGAL::Color color1 = CGAL::Color(200, 200, 0);
const CGAL::Color color2 = CGAL::Color(200, 0, 200);
const CGAL::Color color3 = CGAL::Color(100, 100, 100);

class GeomOutput {
private:
    bool cgal_use_geoview;
    // construct the geom_view
    int id;
    CGAL::Color color;
    CGAL::Geomview_stream gv = 
      CGAL::Geomview_stream(CGAL::Bbox_3(-100, -100, -100, 600, 600, 600));
    void init_color();
    void init_geomoutput();
    
public:
    GeomOutput(int id);
    bool getUsability() { return cgal_use_geoview; }
    int getId() { return id; }
    GeomOutput& operator<<(const CGAL::Color& color) { gv<<color; }
    GeomOutput& operator<<(Polyhedron& polyhedron) { gv<<polyhedron.get_p(); }
    GeomOutput& operator<<(const Point_3& point) { gv<<point; }
    GeomOutput& operator<<(const Triangle_3& triangle) { gv<<triangle; }
    void clearView();
    
    static void draw_polyhedron(
      Polyhedron& p0, Polyhedron& p1,
      GeomOutput& go_0, GeomOutput& go_1,
      const int iter, const char* info);
    static void print_polyhedron_info(Polyhedron& p);
};

#endif
#endif
