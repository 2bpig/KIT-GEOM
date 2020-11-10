#ifndef VIEW_FILE_H
#define VIEW_FILE_H

// std
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h> // for sleep()
// cgal
#include "Data/Data_Polyhedron.hpp"

// declare constants
const std::string filename_off("test2.off");
const std::string filename_obj("test2.obj");
const std::string filename_stl("test.stl");
//const int gray_fac = 220;
//const CGAL::Color color_gray = CGAL::Color(gray_fac,gray_fac,gray_fac);
/*
template <class HDS>
class ExtractHDS4STL : public CGAL::Modifier_base<HDS> {
private:
    std::istream& in_stl;
    
public:
    ExtractHDS4STL(std::istream& in) : in_stl(in) {}
    
    void operator()( HDS& hds);
    
     {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( 3, 1, 6);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( 0, 0, 0));
        B.add_vertex( Point( 1, 0, 0));
        B.add_vertex( Point( 0, 1, 0));
        B.begin_facet();
        B.add_vertex_to_facet( 0);
        B.add_vertex_to_facet( 1);
        B.add_vertex_to_facet( 2);
        B.end_facet();
        B.end_surface();
    }
    
};
*/
class FileIO {
private:
    
    static bool read_stl( std::istream& input,
      std::vector< std::array<double,3> >& points,
      std::vector< std::array<int,3> >& facets );
    
public:

    FileIO();
    
    static bool write_off_polyhedron(Polyhedron& p);
    
    static bool read_off_polyhedron(Polyhedron& p);
    
    static void write_obj_polyhedron(Polyhedron& p);
    
    static void read_obj_polyhedron(Polyhedron& p);
    
    static void write_stl_polyhedron(Polyhedron& p);
    
    static void read_stl_polyhedron(Polyhedron& p);
      
    static void print_polyhedron_info(Polyhedron& p);
    
    static bool split_polyhedron_inFile();
};

#endif
