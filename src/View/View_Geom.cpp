#include "View_Geom.hpp"

#ifndef CGAL_USE_GEOMVIEW
#else

////////////////////////////////////////////////
//    constructor:                            //
////////////////////////////////////////////////

GeomOutput::GeomOutput(int id) : cgal_use_geoview(true) {
    this->id = id;
    init_color();
    init_geomoutput();
}

void GeomOutput::init_color() {
    switch(id){
      case 0: color = color0; break;
      case 1: color = color1; break;
      case 2: color = color2; break;
      case 3: color = color3; break;
      default: color = color0;
    }
}

void GeomOutput::init_geomoutput() {
	gv.set_line_width(4);
//  gv.set_trace(true);
    gv.set_bg_color(color);
	gv.set_wired(true);
	gv.set_face_color(color_gray);
//	gv.set_edge_color(color0);
//	gv.set_vertex_color(CGAL::GREEN);
    gv.clear();
}


////////////////////////////////////////////////
//    public:                                 //
////////////////////////////////////////////////

void GeomOutput::clearView() {
    gv.clear();
}


////////////////////////////////////////////////
//    static:                                 //
////////////////////////////////////////////////

void GeomOutput::draw_polyhedron(
  Polyhedron& p0, Polyhedron& p1,
  GeomOutput& go_0, GeomOutput& go_1,
  const int iter, const char* info) {
    // print info
	std::cout << "Drawing P" << iter << " : " << info << ".\n";
    GeomOutput::print_polyhedron_info(p0);
    // draw polyhedron
	go_0 << p0;
//	go_1 << p1;
//  sleep(5);
	std::cout << "Enter any key to clear view !" << std::endl;
	char nextIter;
	std::cin >> nextIter;
	go_0.clearView();
//	go_1.clearView();
}

void GeomOutput::print_polyhedron_info(Polyhedron& p) {
	std::cout << "facet num : " << p.get_facet_num() << std::endl;
	std::cout << "edge num : " << p.get_edge_num() << std::endl;
	std::cout << "vertex num : " << p.get_vertex_num() << std::endl;
}

//void print_p_info(int order_num, Point_3& p_info) {
//  std::cout << order_v+1 << " : " << h_v->opposite()->vertex()->point() << " -> " << h_v->vertex()->point() << std::endl;
//  std::cout << order_v << "new_v : " << new_v << std::endl;
//	std::cout << "facet center : " << center << std::endl;
//	std::cout << "facet norm : " << normal_facet << std::endl;
//	char ch;
//	std::cin >> ch;
//}

#endif
