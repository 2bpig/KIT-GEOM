#ifndef DATA_PRIMITIVE_H
#define DATA_PRIMITIVE_H
// cgal
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/intersections.h>
#include <CGAL/point_generators_3.h>

// std lib
#include <vector>
#include <array>
#include <list>
#include <map>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>

// boost
#include <boost/variant.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Cartesian<double>     						gvKernel;
typedef gvKernel::Vector_3                               	Vector_3;
typedef gvKernel::Intersect_3                               Intersect_3;
typedef CGAL::Point_3<gvKernel>     						Point_3;
typedef CGAL::Line_3<gvKernel>     							Line_3;
typedef CGAL::Segment_3<gvKernel>     						Segment_3;
typedef CGAL::Ray_3<gvKernel>     							Ray_3;
typedef CGAL::Plane_3<gvKernel>     						Plane_3;
typedef CGAL::Triangle_3<gvKernel>     						Triangle_3;
typedef CGAL::Aff_transformation_3<gvKernel>                Aff_transformation_3;

typedef CGAL::Creator_uniform_3<double, Point_3>                Pt_creator;
typedef CGAL::Random_points_on_sphere_3<Point_3, Pt_creator>    Random_p_sphere;

// declare constants
const double inf_double = std::numeric_limits<double>::infinity();
const Point_3 inf_point(inf_double, inf_double, inf_double);
const Point_3 origin_point(0.0, 0.0, 0.0);
const Point_3 unit_point(1.0, 1.0, 1.0);
const Vector_3 origin_vector = CGAL::ORIGIN - origin_point;

const Aff_transformation_3 t_bca(0,1,0,0,
                                 0,0,1,0,
                                 1,0,0,0,1);
const Aff_transformation_3 t_cab(0,0,1,0,
                                 1,0,0,0,
                                 0,1,0,0,1);

class Primitive {
public:
    //static
    static Point_3 intersection_planes
      (Plane_3& plane1, Plane_3& plane2, Plane_3& plane3);
    static Point_3 intersection_planes
      (Plane_3& plane, Line_3& line);
    
    static Point_3 random_scale_point
      (Point_3& point, double random, double a, double b);
};

class intersect_visitor : public boost::static_visitor<Point_3> {
public:
    Point_3 operator()(Point_3& point) const { return point; }
    Point_3 operator()(Line_3& line) const { return Point_3(inf_double, inf_double, inf_double); }
    Point_3 operator()(Plane_3& plane) const { return Point_3(inf_double, inf_double, inf_double); }
    Point_3 operator()(Ray_3& ray) const { return Point_3(inf_double, inf_double, inf_double); }
};

class intersect_line_visitor : public boost::static_visitor<Line_3> {
public:
    Line_3 operator()(Line_3& line) const { return line; }
    Line_3 operator()(Plane_3& plane) const { return Line_3(origin_point, unit_point); }
    Line_3 operator()(Ray_3& ray) const { return Line_3(origin_point, unit_point); }
};

#endif
