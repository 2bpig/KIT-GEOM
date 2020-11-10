#include "Data_Polyhedron.hpp"

////////////////////////////////////////////////
//    public static:                          //
////////////////////////////////////////////////

Point_3 Primitive::intersection_planes(Plane_3& plane1, Plane_3& plane2, Plane_3& plane3){
    double a_plane1 = plane1.a();
    double b_plane1 = plane1.b();
    double c_plane1 = plane1.c();
    double d_plane1 = plane1.d();
    double a_plane2 = plane2.a();
    double b_plane2 = plane2.b();
    double c_plane2 = plane2.c();
    double d_plane2 = plane2.d();
    double a_plane3 = plane3.a();
    double b_plane3 = plane3.b();
    double c_plane3 = plane3.c();
    double d_plane3 = plane3.d();
    // chose z axis
    int axis = 0;
    double det_z = std::abs(a_plane1*b_plane2 - a_plane2*b_plane1);
    double det_y = std::abs(c_plane1*a_plane2 - c_plane2*a_plane1);
    double det_x = std::abs(b_plane1*c_plane2 - b_plane2*c_plane1);
    double max_det = std::max(det_z, std::max(det_y,det_x));
    if(det_z>=max_det) axis = 3;
    else if(det_y>=max_det) axis = 2;
	else if(det_x>=max_det) axis = 1;
	// build new plane1/2/3 with new z
	Plane_3 new_plane1;
	Plane_3 new_plane2;
	Plane_3 new_plane3;
	if(axis == 1) {         // b c a
	    new_plane1 = Plane_3(b_plane1, c_plane1, a_plane1, d_plane1);
	    new_plane2 = Plane_3(b_plane2, c_plane2, a_plane2, d_plane2);
	    new_plane3 = Plane_3(b_plane3, c_plane3, a_plane3, d_plane3);
	} else if(axis == 2) {  // c a b
	    new_plane1 = Plane_3(c_plane1, a_plane1, b_plane1, d_plane1);
	    new_plane2 = Plane_3(c_plane2, a_plane2, b_plane2, d_plane2);
	    new_plane3 = Plane_3(c_plane3, a_plane3, b_plane3, d_plane3);
	} else {                // a b c as default
	    new_plane1 = Plane_3(a_plane1, b_plane1, c_plane1, d_plane1);
	    new_plane2 = Plane_3(a_plane2, b_plane2, c_plane2, d_plane2);
	    new_plane3 = Plane_3(a_plane3, b_plane3, c_plane3, d_plane3);
	}
	// intersect 3 new planes
    CGAL::cpp11::result_of<Intersect_3(Plane_3, Plane_3)>::type intersect_result1 =
      CGAL::intersection(new_plane1, new_plane2);
    Line_3 line(origin_point, unit_point);
    if(intersect_result1)
        line = boost::apply_visitor( intersect_line_visitor(), *intersect_result1);
//    if(line == Line_3(origin_point, unit_point))
//        std::cout << "Unexpected: neighbouring F_e doesn't intersect with each other" 
//                  << std::endl;
    CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type intersect_result2 =
      CGAL::intersection(line, new_plane3);
    Point_3 point(inf_double, inf_double, inf_double);
    if(intersect_result2)
        point = boost::apply_visitor( intersect_visitor(), *intersect_result2);
    if(point == inf_point)
        std::cout << "Unexpected: intersect_line of neighbouring F_es doesn't intersect with F_e or F_f" << std::endl;
	// rebuild the intersect point
	if(axis == 1) {         // b c a
	    point = Point_3(point.z(), point.x(), point.y());
	} else if(axis == 2) {  // c a b
	    point = Point_3(point.y(), point.z(), point.x());
	} else {                // a b c as default
	    
	}
    return point;
}

Point_3 Primitive::intersection_planes(Plane_3& plane, Line_3& line){
    double a_line = line.to_vector().x();
    double b_line = line.to_vector().y();
    double c_line = line.to_vector().z();
    double a_plane = plane.a();
    double b_plane = plane.b();
    double c_plane = plane.c();
    double d_plane = plane.d();
    // chose z axis
    int axis = 0;
    double max_coord = std::max(c_line, std::max(b_line, a_line));
    if(c_line>=max_coord) axis = 3;
    else if(b_line>=max_coord) axis = 2;
	else if(a_line>=max_coord) axis = 1;
	// build new plane and line with new z
	Line_3 new_line;
	Plane_3 new_plane;
	if(axis == 1) {         // b c a
	    new_line = line.transform(t_bca);
	    new_plane = Plane_3(b_plane, c_plane, a_plane, d_plane);
	} else if(axis == 2) {  // c a b
	    new_line = line.transform(t_cab);
	    new_plane = Plane_3(c_plane, a_plane, b_plane, d_plane);
	} else {                // a b c as default
	    new_line = line;
	    new_plane = Plane_3(a_plane, b_plane, c_plane, d_plane);
	}
	// intersect new plane with new line
    CGAL::cpp11::result_of<Intersect_3(Line_3, Plane_3)>::type intersect_result =
      CGAL::intersection(new_line, new_plane);
    Point_3 point(inf_double, inf_double, inf_double);
    if(intersect_result)
        point = boost::apply_visitor( intersect_visitor(), *intersect_result);
    if(point == inf_point)
        std::cout << "Unexpected: intersect_line of neighbouring F_es doesn't intersect with F_e or F_f" << std::endl;
	// rebuild the intersect point
	if(axis == 1) {         // b c a
	    point = Point_3(point.z(), point.x(), point.y());
	} else if(axis == 2) {  // c a b
	    point = Point_3(point.y(), point.z(), point.x());
	} else {                // a b c as default
	    
	}
    return point;
}

Point_3 Primitive::random_scale_point(Point_3& point, double random, double a, double b){
    // scale in [a, b)
    double scale = a + random*(b-a);
	double new_x = scale*point.x();
	double new_y = scale*point.y();
	double new_z = scale*point.z();
    return Point_3(new_x, new_y, new_z);
}
