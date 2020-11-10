#include "Algo_Edgecutting.hpp"

Edgecutting::Edgecutting(GeomOutput& geomoutput, Polyhedron& polyhedron) : 
  go(geomoutput), polyhedron(polyhedron) {
    depth = 0.0;
    cal_depth = true;
	depth_ratio = 0.4;
	max_depth_ratio = 0.9;
    min_dist = 1e-13;
    sphere_r = 0.0;
    cut_sphere = false;
}

bool Edgecutting::check_min_dist(double dist) {
    if(min_dist < dist){
        min_dist = dist;
	    std::cout << "min min distance : " << min_dist << std::endl;
        return false;
    }
    return true;
}
/*
bool Edgecutting::find_cutting_normal_TILED(Halfedge_handle e, 
  Vector_3& cutting_normal) {
    Vector_3 n1 = e->facet()->plane();
    Vector_3 n2 = e->opposite()->facet()->plane();
    Vector_3 v1end1 = e->vertex()->point() - e->next()->vertex()->point();
    Vector_3 v2end1 = e->vertex()->point() - e->opposite()->prev()->vertex()->point();
    double 
}
*/

void Edgecutting::find_cutting_normal_MID(Halfedge_handle e, 
  Vector_3& cutting_normal) {
    double n1_ratio = 0.5;
    double n2_ratio = 1 - n1_ratio;
    Vector_3 n1 = e->facet()->plane();
    Vector_3 n2 = e->opposite()->facet()->plane();
	cutting_normal = n1_ratio*n1 + n2_ratio*n2;
    cutting_normal = cutting_normal/std::sqrt(cutting_normal.squared_length());
}

bool Edgecutting::find_cutting_pos(Halfedge_handle e,
  Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
  Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
  std::vector<double>& neighborV_cutPlane_dist_list,
  Vector_3& cutting_normal, Point_3& cutting_pos) {
    // find the nearest neighbor in the direction of normal
	double neighborV_cutPlane_dist = 0.0;
	double min_neighborV_cutPlane_dist = inf_double;
    Vector_3 v2end1;
    std::size_t order_v = 0;
	h_end1++;
	do {
    	v2end1 = e_end1->point() - h_end1->opposite()->vertex()->point();
    	neighborV_cutPlane_dist = v2end1*cutting_normal;
    	// TODO re-calculate the cutting normal, or give up cutting this edge
    	if(neighborV_cutPlane_dist < min_dist) return false;
    	if(neighborV_cutPlane_dist < min_neighborV_cutPlane_dist) {
    	    min_neighborV_cutPlane_dist = neighborV_cutPlane_dist;
    	}
        neighborV_cutPlane_dist_list.push_back(neighborV_cutPlane_dist);
    	++ order_v;
    } while ( ++h_end1 != h_end1_begin);
    h_end2++;
    do {
    	v2end1 = e_end1->point() - h_end2->opposite()->vertex()->point();
    	neighborV_cutPlane_dist = v2end1*cutting_normal;
    	// TODO re-calculate the cutting normal, or give up cutting this edge
    	if(neighborV_cutPlane_dist < min_dist) return false;
    	if(neighborV_cutPlane_dist < min_neighborV_cutPlane_dist) {
    	    min_neighborV_cutPlane_dist = neighborV_cutPlane_dist;
    	}
        neighborV_cutPlane_dist_list.push_back(neighborV_cutPlane_dist);
    	++ order_v;
    } while ( ++h_end2 != h_end2_begin);
    CGAL_assertion( neighborV_cutPlane_dist_list.size() == order_v);
    // decide the cutting depth
	double cutting_depth;
    if(cal_depth) {
        cutting_depth = min_neighborV_cutPlane_dist*depth_ratio;
        depth = cutting_depth;
    } else {
        cutting_depth = std::min(depth, min_neighborV_cutPlane_dist*max_depth_ratio);
    }
	cutting_pos = e_end1->point() - cutting_normal*cutting_depth;
	return true;
}

bool Edgecutting::find_cutting_pos_SPHERE(Halfedge_handle e,
  Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
  Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
  std::vector<double>& neighborV_cutPlane_dist_list,
  Vector_3& cutting_normal, Point_3& cutting_pos) {
	cutting_pos = origin_point + cutting_normal*sphere_r;
    Vector_3 pos2end1 = e_end1->point() - cutting_pos;
    double pos_cutPlane_dist = pos2end1*cutting_normal;
    // find the nearest neighbor in the direction of normal
	double neighborV_cutPlane_dist = 0.0;
	double min_neighborV_cutPlane_dist = inf_double;
    Vector_3 v2end1;
    std::size_t order_v = 0;
	h_end1++;
	do {
    	v2end1 = e_end1->point() - h_end1->opposite()->vertex()->point();
    	neighborV_cutPlane_dist = v2end1*cutting_normal;
    	// TODO re-calculate the cutting normal, or give up cutting this edge
    	if(neighborV_cutPlane_dist < min_dist) return false;
    	if(neighborV_cutPlane_dist < min_neighborV_cutPlane_dist) {
    	    min_neighborV_cutPlane_dist = neighborV_cutPlane_dist;
    	}
        neighborV_cutPlane_dist_list.push_back(neighborV_cutPlane_dist);
    	++ order_v;
    } while ( ++h_end1 != h_end1_begin);
    h_end2++;
    do {
    	v2end1 = e_end1->point() - h_end2->opposite()->vertex()->point();
    	neighborV_cutPlane_dist = v2end1*cutting_normal;
    	// TODO re-calculate the cutting normal, or give up cutting this edge
    	if(neighborV_cutPlane_dist < min_dist) return false;
    	if(neighborV_cutPlane_dist < min_neighborV_cutPlane_dist) {
    	    min_neighborV_cutPlane_dist = neighborV_cutPlane_dist;
    	}
        neighborV_cutPlane_dist_list.push_back(neighborV_cutPlane_dist);
    	++ order_v;
    } while ( ++h_end2 != h_end2_begin);
    CGAL_assertion( neighborV_cutPlane_dist_list.size() == order_v);
    if(pos_cutPlane_dist <= 0.0 || pos_cutPlane_dist >= min_neighborV_cutPlane_dist){
        std::cout << "fail to cut sphere ! " << pos_cutPlane_dist << " > " << min_neighborV_cutPlane_dist << std::endl;
        return false;
    }
	return true;
}

bool Edgecutting::find_cutting_plane(Halfedge_handle e,
  Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
  Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
  std::vector<double>& neighborV_cutPlane_dist_list,
  Vector_3& cutting_normal, Plane_3& cutting_plane) {
    // decide the cutting normal
    find_cutting_normal_MID(e, cutting_normal);
    Vector_3 end22end1 = e_end1->point() - e_end2->point();
	check_min_dist(std::abs(end22end1*cutting_normal));
	// decide the cutting pos : find_cutting_pos_SPHERE(...)
	Point_3 cutting_pos;
	if(!find_cutting_pos(e, 
	  e_end1, h_end1_begin, h_end1, e_end2, h_end2_begin, h_end2,
      neighborV_cutPlane_dist_list, cutting_normal, cutting_pos)) return false;
	cutting_plane = Plane_3(cutting_pos, cutting_normal);
	return true;
}

void Edgecutting::intersect_cutting_plane(Plane_3 cutting_plane,
  Vertex_handle e_end1, HV_circulator h_end1_begin, HV_circulator& h_end1,
  Vertex_handle e_end2, HV_circulator h_end2_begin, HV_circulator& h_end2,
  std::vector<Point_3>& new_vertices, 
  std::vector<Halfedge_handle>& end1_halfedges,
  std::vector<Halfedge_handle>& end2_halfedges) {
	Line_3 v2end1_line;
	Line_3 v2end2_line;
	Point_3 new_v;
	std::size_t order_v = 0;
	h_end1++;
    do {
        v2end1_line = Line_3(h_end1->opposite()->vertex()->point(), e_end1->point());
    	new_v = Primitive::intersection_planes(cutting_plane, v2end1_line);
        new_vertices.push_back(new_v);
        end1_halfedges.push_back(h_end1);
    	++ order_v;
    } while ( ++h_end1 != h_end1_begin);
	h_end2++;
    do {
        v2end2_line = Line_3(h_end2->opposite()->vertex()->point(), e_end2->point());
    	new_v = Primitive::intersection_planes(cutting_plane, v2end2_line);
        new_vertices.push_back(new_v);
        end2_halfedges.push_back(h_end2);
    	++ order_v;
    } while ( ++h_end2 != h_end2_begin);
    CGAL_assertion( new_vertices.size() == order_v);
}

void Edgecutting::construct_cutting_facet(Halfedge_handle& e,
  std::vector<Halfedge_handle>& end1_halfedges,
  std::vector<Halfedge_handle>& end2_halfedges,
  std::vector<Point_3> new_vertices,
  std::vector<Halfedge_handle>& new_edges,
  Vector_3 cutting_normal) {
    Polyhedron_cgal& p = polyhedron.get_p();
	CGAL_precondition(p.is_valid());
	
	// split the edges adjacent to both end of the edge
    std::size_t order_v = 0;
    for( Halfedge_handle end1_halfedge : end1_halfedges ) {
        p.split_edge(end1_halfedge->opposite());
        end1_halfedge->vertex()->point() = new_vertices[order_v];
    	++ order_v;
    }
    for( Halfedge_handle end2_halfedge : end2_halfedges ) {
        p.split_edge(end2_halfedge->opposite());
        end2_halfedge->vertex()->point() = new_vertices[order_v];
    	++ order_v;
    }
    // split the facets adjacent to the edge
    std::size_t order_v1 = end1_halfedges.size();
    order_v = 0;
    for( Halfedge_handle end1_halfedge : end1_halfedges ) {
        order_v++;
        if( order_v < order_v1 ) {
            p.split_facet(end1_halfedge, end1_halfedge->next()->next());
        } else {
            p.split_facet(end1_halfedge, end1_halfedge->next()->next()->next());
        }
        new_edges.push_back(end1_halfedge->next());
    }
    std::size_t order_v2 = end2_halfedges.size();
    order_v = 0;
    for( Halfedge_handle end2_halfedge : end2_halfedges ) {
        order_v++;
        if( order_v < order_v2 ) {
            p.split_facet(end2_halfedge, end2_halfedge->next()->next());
        } else {
            p.split_facet(end2_halfedge, end2_halfedge->next()->next()->next());
        }
        new_edges.push_back(end2_halfedge->next());
    }
    // remove the edge and construct the new cutting facet
    Halfedge_handle second_center = e->prev();
    Halfedge_handle e_cutting_facet = second_center->prev();
    p.erase_center_vertex(e);
    p.erase_center_vertex(second_center);
    e_cutting_facet->facet()->plane() = cutting_normal;
	CGAL_postcondition(p.is_valid());
}

bool Edgecutting::cut_edge(Halfedge_handle& e, 
  std::vector<Halfedge_handle>& new_edges) {
	// prepare the ends of edge and the circulator of ends
    Vertex_handle e_end1 = e->vertex();
    Vertex_handle e_end2 = e->opposite()->vertex();
	HV_circulator h_end1 = e->vertex_begin();
	HV_circulator h_end1_begin = e->vertex_begin();
	HV_circulator h_end2 = e->opposite()->vertex_begin();
	HV_circulator h_end2_begin = e->opposite()->vertex_begin();
    // find the cutting plane with normal and depth
	std::vector<double> neighborV_cutPlane_dist_list;
	Vector_3 cutting_normal;
	Plane_3 cutting_plane;
    if(!find_cutting_plane(e,
      e_end1, h_end1_begin, h_end1, e_end2, h_end2_begin, h_end2,
      neighborV_cutPlane_dist_list,
      cutting_normal, cutting_plane)) {
        std::cout << "fail to cut current edge !" << std::endl;
        return false;
    }
    // intersect the cutting plane with all adjacent edges
    std::vector<Point_3> new_vertices;
	std::vector<Halfedge_handle> end1_halfedges;
	std::vector<Halfedge_handle> end2_halfedges;
    intersect_cutting_plane(cutting_plane,
      e_end1, h_end1_begin, h_end1, e_end2, h_end2_begin, h_end2,
      new_vertices, end1_halfedges, end2_halfedges);
    // split the edges adjacent to both ends, and construct the new cutting facet
    construct_cutting_facet(e, 
      end1_halfedges, end2_halfedges,
      new_vertices, new_edges, cutting_normal);
    return true;
}

void Edgecutting::init_sphere_r() {
    Polyhedron_cgal& p = polyhedron.get_p();
    Halfedge_handle e = p.edges_begin();
    Vector_3 n = e->facet()->plane();
    Vector_3 v = e->vertex()->point() - origin_point;
    sphere_r = n*v;
    cut_sphere = true;
}

void Edgecutting::update_depth_ratio(double ratio) {
    depth_ratio = ratio;
}

void Edgecutting::update_cal_depth(bool cal) {
    cal_depth = cal;
}


////////////////////////////////////////////////
//    DUAL:                                   //
////////////////////////////////////////////////

bool Edgecutting::interpolate_vertex4edge(Halfedge_handle& e, 
  std::vector<Halfedge_handle>& new_edges) {
	// prepare the ends of edge and the circulator of ends
    Vertex_handle e_end1 = e->vertex();
    Vertex_handle e_end2 = e->opposite()->vertex();
    Vertex_handle f1_v = e->next()->vertex();
    Vertex_handle f2_v = e->opposite()->next()->vertex();
    // find the pos of new vertex on the sphere
    Vector_3 end2_end1 = e_end1->point() - e_end2->point();
    Vector_3 origin_end2 = e_end2->point() - origin_point;
    Vector_3 mid_v = end2_end1*0.5 + origin_end2;
    mid_v = sphere_r*mid_v/std::sqrt(mid_v.squared_length());
    
    if(!true) {
        std::cout << "fail to cut current edge !" << std::endl;
        return false;
    }
    // construct new edges around new vertex
    Polyhedron_cgal& p = polyhedron.get_p();
	CGAL_precondition(p.is_valid());
	
    p.split_edge(e->opposite());
    e->vertex()->point() = Point_3(mid_v.x(), mid_v.y(), mid_v.z());
    Halfedge_handle next_e = e->next();
    p.split_facet(e, next_e->next());
    p.split_facet(next_e->opposite(), e->opposite()->next());
    new_edges.push_back(e);
    new_edges.push_back(e->next());
    new_edges.push_back(next_e);
    new_edges.push_back(e->opposite()->prev());
	
	CGAL_postcondition(p.is_valid());
    
    return true;
}

void Edgecutting::init_sphere_r_dual() {
    Polyhedron_cgal& p = polyhedron.get_p();
    Halfedge_handle e = p.edges_begin();
    Vector_3 v = e->vertex()->point() - origin_point;
    sphere_r = std::sqrt(v.squared_length());
    cut_sphere = true;
}
