#ifndef ALGO_EDGECUTTING_H
#define ALGO_EDGECUTTING_H
// geom-output and data
#include "View/View_Geom.hpp"
#include "View/View_File.hpp"
// algo
#include "Algo/Algo_Dualizing.hpp"
// std

// declare constants
const double tolerance = 0.1;
const double error_prop = 0.01;
const double angle_prop = 0.01;

class Edgecutting {
private:
	GeomOutput& go;
	Polyhedron& polyhedron;
    // parameters
    double depth;
    bool cal_depth;
	double depth_ratio;
	double max_depth_ratio;
    double min_dist;
    double sphere_r;
    bool cut_sphere;
    bool check_min_dist(double dist);
    // edgecutting scheme
    void find_cutting_normal_MID(Halfedge_handle e, 
      Vector_3& cutting_normal);
    bool find_cutting_pos(Halfedge_handle e,
      Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
      Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
      std::vector<double>& neighborV_cutPlane_dist_list,
      Vector_3& cutting_normal, Point_3& cutting_pos);
    bool find_cutting_pos_SPHERE(Halfedge_handle e,
      Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
      Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
      std::vector<double>& neighborV_cutPlane_dist_list,
      Vector_3& cutting_normal, Point_3& cutting_pos);
    bool find_cutting_plane(Halfedge_handle e,
      Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
      Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
      std::vector<double>& neighborV_cutPlane_dist_list,
      Vector_3& cutting_normal, Plane_3& cutting_plane);
    void intersect_cutting_plane(Plane_3 cutting_plane,
      Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
      Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
      std::vector<Point_3>& new_vertices, 
      std::vector<Halfedge_handle>& end1_halfedges,
      std::vector<Halfedge_handle>& end2_halfedges);
    void construct_cutting_facet(Halfedge_handle& e,
      std::vector<Halfedge_handle>& end1_halfedges,
      std::vector<Halfedge_handle>& end2_halfedges,
      std::vector<Point_3> new_vertices,
      std::vector<Halfedge_handle>& new_edges,
      Vector_3 cutting_normal);
public:
    Edgecutting(GeomOutput& geomoutput, Polyhedron& polyhedron);
    void init_sphere_r();
    void init_sphere_r_dual();
	void update_depth_ratio(double ratio);
	void update_cal_depth(bool cal);
    bool cut_edge(Halfedge_handle& e, std::vector<Halfedge_handle>& new_edges);
    bool interpolate_vertex4edge(Halfedge_handle& e, std::vector<Halfedge_handle>& new_edges);
};

#endif
