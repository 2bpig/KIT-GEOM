#include "Algo_Honeycomb.hpp"

Honeycomb::Honeycomb(GeomOutput& geomoutput) : go(geomoutput) {
    mode = 0;
    iter_num = 0;
}

void Honeycomb::cal_edge_normal(Halfedge_handle& edge, 
  std::map<Halfedge_handle, Vector_3, Hedge_cmp>& hedge_to_normal_map,
  std::map<Facet_handle, Point_3, Facet_cmp>& facet_to_center_map) {
    // first, check in hedge -> normal map, for both hedge of current edge
    Vector_3 normal_edge;
    std::map<Halfedge_handle, Vector_3, Hedge_cmp>::iterator hedge_map_iter;
    hedge_map_iter = hedge_to_normal_map.find(edge);
    if(hedge_map_iter != hedge_to_normal_map.end()) return;
    Halfedge_handle edge_opposite = edge->opposite();
    hedge_map_iter = hedge_to_normal_map.find(edge_opposite);
    if(hedge_map_iter != hedge_to_normal_map.end()) {
        normal_edge = hedge_map_iter->second;
        hedge_to_normal_map.insert(std::pair<Halfedge_handle, Vector_3>(edge, normal_edge));
        return;
    }
    // not exists, cal normal of edge and insert in hedge -> normal map
    Vector_3 normal_left = edge->facet()->plane();
    Vector_3 normal_right = edge_opposite->facet()->plane();
//    if(mode != 4) {
//        normal_edge = alpha_edge*normal_left + (1.0-alpha_edge)*normal_right;
//    }else{
        // edge center
        Point_3 edge_to = edge->vertex()->point();
        Point_3 edge_from = edge_opposite->vertex()->point();
        Point_3 edge_center = edge_from + 0.5*(edge_to - edge_from);
        // left facet center
        std::map<Facet_handle, Point_3, Facet_cmp>::iterator facet_map_iter;
        Point_3 center_left;
        facet_map_iter = facet_to_center_map.find(edge->facet());
        if(facet_map_iter == facet_to_center_map.end()) {
            HF_circulator h_f_left = edge->facet()->facet_begin();
            Vector_3 vec_left( 0.0, 0.0, 0.0 );
            std::size_t order_f_left = 0;
            do {
                vec_left = vec_left + ( h_f_left->vertex()->point() - CGAL::ORIGIN);
                ++ order_f_left;
            } while ( ++h_f_left != edge->facet()->facet_begin());
            CGAL_assertion(order_f_left >= 3);
            center_left =  CGAL::ORIGIN + (vec_left / static_cast<double>(order_f_left));
            facet_to_center_map.insert(std::pair<Facet_handle, Point_3>(edge->facet(), center_left));
        } else {
            center_left = facet_map_iter->second;
        }
        double left_dist = std::sqrt((center_left - edge_center).squared_length());
        // right facet center
        Point_3 center_right;
        facet_map_iter = facet_to_center_map.find(edge_opposite->facet());
        if(facet_map_iter == facet_to_center_map.end()) {
            HF_circulator h_f_right = edge_opposite->facet()->facet_begin();
            Vector_3 vec_right( 0.0, 0.0, 0.0 );
            std::size_t order_f_right = 0;
            do {
                vec_right = vec_right + ( h_f_right->vertex()->point() - CGAL::ORIGIN);
                ++ order_f_right;
            } while ( ++h_f_right != edge_opposite->facet()->facet_begin());
            CGAL_assertion(order_f_right >= 3);
            center_right =  CGAL::ORIGIN + (vec_right / static_cast<double>(order_f_right));
            facet_to_center_map.insert(std::pair<Facet_handle, Point_3>(edge_opposite->facet(), center_right));
        } else {
            center_right = facet_map_iter->second;
        }
        double right_dist = std::sqrt((center_right - edge_center).squared_length());
    // edge normal
    if(mode != 4) {
        normal_edge = alpha_edge*normal_left + (1.0-alpha_edge)*normal_right;
    }else{
        double left_ratio = std::min((1.0-tolerance), std::max((left_dist/(left_dist + right_dist)), tolerance));
        normal_edge = left_ratio*normal_left + (1.0-left_ratio)*normal_right;
    }
    hedge_to_normal_map.insert(std::pair<Halfedge_handle, Vector_3>(edge, normal_edge));
    hedge_to_normal_map.insert(std::pair<Halfedge_handle, Vector_3>(edge_opposite, normal_edge));
}

Vector_3 Honeycomb::cal_normal_for_lowest_plane_honeycomb(
  Point_3 p1, Point_3 p2, Point_3 min_p, std::vector<Point_3>& p_list, double error_dist) {
    std::size_t p_num = p_list.size();
    Point_3 p0 = min_p;
    Vector_3 p0_normal = CGAL::unit_normal(min_p, p1, p2);
    std::size_t p_order = 0;
    do {
        Point_3 p = p_list[p_order];
        Vector_3 p_vector = p - p1;
        double dist = p_vector*p0_normal;
        if(dist<(-0.1*error_dist)) {
            p0 = p;
            p0_normal = CGAL::unit_normal(p, p1, p2);
        }
    } while (++p_order < p_num);
    return p0_normal;
}


bool Honeycomb::cal_new_vertices_new_normal_honeycomb(
  Vector_3& normal_facet, Point_3 center, 
  HV_circulator& h_v, const HV_circulator& h_v_begin, std::size_t max_order, 
  std::vector<Point_3>& new_vertices, std::vector<Vector_3>& normal_edges, 
  std::vector<Plane_3>& plane_edges) {
    // min and tolerant distance of all the new intersection points to the facet
    std::vector<Point_3> facet_v;
    Point_3 facet_v_local;
    std::vector<Point_3> intersect_v;
    Point_3 intersect_v_local;
    Point_3 min_intersect_v;
    double avg_dist = 0.0;
    double min_dist = inf_double;
    double min_dist_local = inf_double;
    std::vector<double> edge_length_list;
    double total_length = 0.0;
    double min_length = inf_double;
    std::size_t order_v = 0;
    std::size_t sub_order_v = 0;
    do {
        Point_3 facet_vertex = h_v->prev()->vertex()->point();
        Point_3 prev_vertex = h_v->opposite()->next()->vertex()->point();
        double edge_length = std::sqrt((facet_vertex - prev_vertex).squared_length());
        total_length += edge_length;
        edge_length_list.push_back(edge_length);
        if(edge_length<min_length) min_length = edge_length;
        do {
        	Point_3 new_v = Primitive::intersection_planes(
        	  plane_edges[order_v],
        	  plane_edges[(order_v+1)%max_order],
              plane_edges[(order_v+2+sub_order_v)%max_order]);
            if(new_v!=inf_point){
            	Vector_3 new_vector = new_v - center;
            	double dist = new_vector*normal_facet;
            	if((dist>=0.0) && (dist<min_dist_local)){
            		min_dist_local = dist;
            		intersect_v_local = new_v;
                    if(dist<min_dist){
                        min_dist = dist;
                        min_intersect_v = new_v;
                    }
        	    }
            }
        } while(sub_order_v++ < (max_order - 3));
        avg_dist += min_dist_local;
        min_dist_local = inf_double;
        intersect_v.push_back(intersect_v_local);
        facet_v.push_back(facet_vertex);
        sub_order_v = 0;
    	++ order_v;
    } while ( ++h_v != h_v_begin );
    CGAL_assertion( intersect_v.size() == max_order);
    CGAL_assertion( facet_v.size() == max_order);
    // check ratio between min_dist and min edge of facet
//    double ratio = min_dist/min_length;
//    if(ratio_sq < error_prop){
        // do not rotate the facet to create new vertices
        // TODO
//        std::cout << "min_edge : " << min_length << std::endl;
//        return false;
//    }
    avg_dist = avg_dist/static_cast<double>(order_v);
    avg_dist = std::min(alpha_facet*avg_dist, (1.0-tolerance)*min_dist);
    // adjust facet_v and intersection_v, avg intersection_v
    std::vector<Point_3> adj_facet_v;
    std::vector<Point_3> adj_intersect_v;
    std::vector<Point_3> avg_new_v;
    Point_3 avg_v_local;
    double tolerance_dist = tolerance*min_dist;
    double error_dist = error_prop*min_dist;
    order_v = 0;    
    do {
        intersect_v_local = intersect_v[order_v];
        facet_v_local = facet_v[order_v];
        Vector_3 facet_intersect_vector = intersect_v_local - facet_v_local;
        double normal_dist = facet_intersect_vector*normal_facet;
        // adjust along the line of facet_v and intersection_v
        // adjust along the normal_facet
        Vector_3 tolerance_vector = tolerance_dist*normal_facet; //facet_intersect_vector*tolerance_dist/normal_dist;
        Vector_3 avg_vector = facet_intersect_vector*avg_dist/normal_dist;
        avg_v_local = facet_v_local + avg_vector;
        intersect_v_local = intersect_v_local - tolerance_vector;
        facet_v_local = facet_v_local + tolerance_vector;
        adj_intersect_v.push_back(intersect_v_local);
        adj_facet_v.push_back(facet_v_local);
        avg_new_v.push_back(avg_v_local);
    	++ order_v;
    } while ( ++h_v != h_v_begin);
    CGAL_assertion( adj_intersect_v.size() == max_order);
    CGAL_assertion( adj_facet_v.size() == max_order);
    // avg normal of all the plane of edges
	Vector_3 avg_normal = origin_vector;
	double normal_weight = 0.0;
    order_v = 0;
    do {
        if(mode == 3) normal_weight = edge_length_list[order_v];
        else normal_weight = 1.0;
        avg_normal = avg_normal + normal_edges[order_v]*normal_weight;
        ++ order_v;
    } while ( ++h_v != h_v_begin );
    avg_normal = avg_normal - origin_vector;
    avg_normal = avg_normal/std::sqrt(avg_normal.squared_length());
    // check difference between avg-normal and normal-facet
    double diff_sq = (avg_normal-normal_facet).squared_length();
//    if(diff_sq < angle_prop*angle_prop){
    if(diff_sq < angle_prop*angle_prop){
        // do not rotate the facet to create new vertices
        order_v = 0;    
        do {
            new_vertices.push_back(avg_new_v[order_v]);
            ++ order_v;
        } while ( ++h_v != h_v_begin);
        return true;
    }
    // min and avg distance of all the new intersection points to avg-normal plane
    double avg_dist_avg_normal = 0.0;
    double min_dist_avg_normal = inf_double;
    order_v = 0;
    do {
        intersect_v_local = adj_intersect_v[order_v];
        Vector_3 v_vector = intersect_v_local - center;
        double dist = v_vector*avg_normal;
        if(dist<min_dist_avg_normal) min_dist_avg_normal = dist;
        avg_dist_avg_normal += dist;
    	++ order_v;
    } while ( ++h_v != h_v_begin );
    avg_dist_avg_normal = avg_dist_avg_normal/static_cast<double>(order_v);
    avg_dist_avg_normal = std::min(alpha_facet*avg_dist_avg_normal, min_dist_avg_normal);
    // max distance of all the vertices of facet to avg-normal plane
    // lowest avg-normal plane according to all the vertices of facet
    double max_dist_avg_normal = 0.0;
    std::vector<std::size_t> max_vertex_orders;
    std::size_t max_vertex_num = 0;
    order_v = 0;
    do {
        std::size_t vertex_order = order_v;
        facet_v_local = adj_facet_v[order_v];
        Vector_3 v_vector = facet_v_local - center;
        double dist = v_vector*avg_normal;
        if(dist>max_dist_avg_normal+error_dist){
            max_vertex_num = 1;
            max_dist_avg_normal = dist;
            max_vertex_orders.clear();
            max_vertex_orders.push_back(vertex_order);
        } else if(dist>=max_dist_avg_normal-error_dist){
            max_vertex_num++;
            max_dist_avg_normal = std::max(dist, max_dist_avg_normal);
            max_vertex_orders.push_back(vertex_order);
        }
        ++ order_v;
    } while ( ++h_v != h_v_begin );
    CGAL_assertion( max_vertex_orders.size() == max_vertex_num);
    CGAL_assertion( max_vertex_num <= max_order);
    // new plane with avg-normal or adjusted avg-normal
	Vector_3 new_normal;
	Point_3 new_center;
    if(max_dist_avg_normal<avg_dist_avg_normal) {          // avg avg-normal plane
        // TODO std::cout << "avg" << std::endl;
        new_normal = avg_normal;
        new_center = center + avg_normal*avg_dist_avg_normal;
    } else if(max_dist_avg_normal<=min_dist_avg_normal) {   // lowest avg-normal plane
        // TODO std::cout << "lowest" << std::endl;
        new_normal = avg_normal;
        new_center = center + avg_normal*max_dist_avg_normal;
    } else {                                                // adjusted avg-normal plane
        if(max_vertex_num==1) {
            // TODO std::cout << "adj_1" << std::endl;
	        // significant vertex with max distance to the avg-normal plane
	        std::size_t max_vertex_order = max_vertex_orders[0];
	        Point_3 max_vertex = adj_facet_v[max_vertex_order];
	        Point_3 prev_vertex = adj_facet_v[(max_vertex_order-1)%max_order];
	        Point_3 next_vertex = adj_facet_v[(max_vertex_order+1)%max_order];
	        // two planes neighboring to the significant vertex
	        Vector_3 prev_normal = cal_normal_for_lowest_plane_honeycomb(max_vertex, prev_vertex, 
	          min_intersect_v, adj_intersect_v, error_dist);
	        Vector_3 next_normal = cal_normal_for_lowest_plane_honeycomb(next_vertex, max_vertex, 
	          min_intersect_v, adj_intersect_v, error_dist);
	        // find nearest point from avg normal to the segment between two normals of the planes
            Vector_3 prev_next_vector = next_normal - prev_normal;
            Vector_3 prev_avg_vector = avg_normal - prev_normal;
            double length_ratio = (prev_avg_vector*prev_next_vector)/prev_next_vector.squared_length();
	        std::cout << "ratio : " << length_ratio << std::endl;
	        length_ratio = std::min(std::max(length_ratio, 0.0), 1.0);
            new_normal = prev_normal + prev_next_vector*length_ratio;
	        new_center = max_vertex;
        } else if(max_vertex_num>=2 && max_vertex_num<max_order) {
            // TODO
	        std::cout << "adj_2" << std::endl;
//	        std::cout << "first : " << max_vertex_orders[0] << std::endl;
//	        std::cout << "last : " << max_vertex_orders[max_vertex_num-1] << std::endl;
//	        std::cout << "i_v num : " << adj_intersect_v.size() << std::endl;
	        // two significant vertices with max distance to the avg-normal plane
	        Point_3 max_vertex_1 = adj_facet_v[max_vertex_orders[0]];
	        Point_3 max_vertex_2 = adj_facet_v[max_vertex_orders[max_vertex_num-1]];
	        // build new plane according to two significant vertices
	        new_normal = cal_normal_for_lowest_plane_honeycomb(max_vertex_2, max_vertex_1, 
	          min_intersect_v, adj_intersect_v, error_dist);
	        new_center = max_vertex_1;
        } else {
            // TODO std::cout << "adj_m : " << max_vertex_num << " : " << max_order << std::endl;
            // do not rotate the facet to create new vertices
            return false;
        }
    }
	// build new plane and find new vertices around the new facet
	Plane_3 new_plane_facet(new_center, new_normal);
    order_v = 0;
    do {
    	Point_3 new_v = Primitive::intersection_planes(
    	  plane_edges[order_v],
    	  plane_edges[(order_v+1)%max_order],
    	  new_plane_facet);
        new_vertices.push_back(new_v);
    	++ order_v;
    } while ( ++h_v != h_v_begin);
    CGAL_assertion( new_vertices.size() == max_order);
    return true;
}

void Honeycomb::cal_new_vertices_honeycomb(Vector_3& normal_facet, 
  HV_circulator& h_v, const HV_circulator& h_v_begin, std::size_t max_order, 
  std::vector<Point_3>& new_vertices, std::vector<Plane_3>& plane_edges) {
    double avg_dist = 0.0;
    double min_dist = inf_double;
    double min_dist_local = inf_double;
    Point_3 min_v(inf_double, inf_double, inf_double);
    Point_3 min_v_local;
    Vector_3 min_vector;
    std::size_t order_v = 0;
    std::size_t sub_order_v = 0;
    do {
        do {
        	Point_3 new_v = Primitive::intersection_planes(
        	  plane_edges[order_v],
        	  plane_edges[(order_v+1)%max_order],
              plane_edges[(order_v+2+sub_order_v)%max_order]);
            if(new_v!=inf_point){
            	Vector_3 new_vector = new_v - h_v->prev()->vertex()->point();
            	double dist = new_vector*normal_facet;
            	if((dist>=0.0) && (dist<min_dist_local)){
            		min_dist_local = dist;
            		min_v_local = new_v;
                    if(dist<min_dist){
                        min_dist = dist;
            		    min_v = h_v->prev()->vertex()->point();
                        min_vector = new_vector;
                    }
        	    }
            }
        } while(sub_order_v++ < (max_order - 3));
        avg_dist += min_dist_local;
        min_dist_local = inf_double;
        sub_order_v = 0;
    	++ order_v;
    } while ( ++h_v != h_v_begin);
    if(mode == 0) {
        avg_dist = alpha_facet*min_dist;
    } else if(mode == 1||mode==4) {
        avg_dist = std::min(alpha_facet*avg_dist/static_cast<double>(max_order), (1-tolerance)*min_dist);
    }
    min_v = min_v + min_vector*avg_dist/min_dist;
	Plane_3 plane_facet(min_v, normal_facet);
    order_v = 0;
    do {
    	Point_3 new_v = Primitive::intersection_planes(plane_edges[order_v],plane_edges[(order_v+1)%max_order],plane_facet);
        new_vertices.push_back(new_v);
        
        	// TODO
            if(iter_num == 12) {
                //std::cout << "vertex -> new_point : " << Vector_3(h_v->opposite()->vertex()->point(), new_vertices[order_v]) << std::endl;
                Point_3 prev_v = h_v->opposite()->next()->vertex()->point();
                Point_3 curr_v = h_v->opposite()->vertex()->point();
                Point_3 next_v = h_v->prev()->prev()->vertex()->point();
                go<<CGAL::Color(255, 0, 0);
                go<<new_v;
                go<<curr_v;
                go<<CGAL::Color(0, 255, 0);
	            //std::cout << "normal 0 : " << plane_edges[order_v].orthogonal_vector() << std::endl;
                go<<prev_v;
                Point_3 new_on_plane0 = plane_edges[order_v].projection(new_v);
                Point_3 prev_on_plane0 = plane_edges[order_v].projection(prev_v);
                Point_3 curr_on_plane0 = plane_edges[order_v].projection(curr_v);
                Point_3 p1_plane0 = plane_edges[order_v].point();
                Point_3 p2_plane0 = p1_plane0 + 100.0*plane_edges[order_v].base1();
                Point_3 p3_plane0 = p1_plane0 + 100.0*plane_edges[order_v].base2();
//                go<<Triangle_3(new_on_plane0,prev_on_plane0,curr_on_plane0);
                go<<Triangle_3(p1_plane0,p2_plane0,p3_plane0);
                
                go<<CGAL::Color(0, 0, 255);
	            //std::cout << "normal 1 : " << plane_edges[(order_v+1)%max_order].orthogonal_vector() << std::endl;
                go<<next_v;
                Point_3 new_on_plane1 = plane_edges[(order_v+1)%max_order].projection(new_v);
                Point_3 curr_on_plane1 = plane_edges[(order_v+1)%max_order].projection(curr_v);
                Point_3 next_on_plane1 = plane_edges[(order_v+1)%max_order].projection(next_v);
                Point_3 p1_plane1 = plane_edges[(order_v+1)%max_order].projection(p1_plane0);
                Point_3 p2_plane1 = plane_edges[(order_v+1)%max_order].projection(p2_plane0);
                Point_3 p3_plane1 = plane_edges[(order_v+1)%max_order].projection(p3_plane0);
//                go<<Triangle_3(new_on_plane1,curr_on_plane1,next_on_plane1);
                go<<Triangle_3(p1_plane1,p2_plane1,p3_plane1);
                
                go<<CGAL::Color(220, 220, 220);
	            std::cout << "Enter a key to next vertex" << std::endl;
	            char ch1;
	            std::cin >> ch1;
        	}
    	++ order_v;
    } while ( ++h_v != h_v_begin);
}

bool Honeycomb::pull_facet( Polyhedron_cgal& p, Facet_iterator& f, Point_3 f_center, 
  std::map<Halfedge_handle, Vector_3, Hedge_cmp>& hedge_to_normal_map) {
	Vector_3 normal_facet = f->plane();
//	HF_circulator h_f = f->facet_begin();
    // the normal for the origin facet remains in new facets
    Halfedge_handle new_center = p.create_center_vertex( f->halfedge()); 
    new_center->vertex()->point() = f_center; // TODO can be removed
    // calculate normal for edges, HF/HV in clock-wise rotation
    std::vector<Plane_3> plane_edges;
    std::vector<Vector_3> normal_edges;
    std::vector<Point_3> new_vertices;
	HV_circulator h_v = new_center->vertex()->vertex_begin();
    std::size_t order_v = 0;
    do {
    	Halfedge_handle edge = h_v->opposite()->next();
    	Vector_3 normal_edge = hedge_to_normal_map.at(edge);
    	normal_edges.push_back(normal_edge);
    	Plane_3 plane_edge(edge->vertex()->point(), normal_edge);
    	plane_edges.push_back(plane_edge);
        ++ order_v;
    } while ( ++h_v != new_center->vertex()->vertex_begin());
    // calculate new vertices
    bool is_new_vertices = true;
    if(mode == 2 || mode == 3) {
        is_new_vertices = cal_new_vertices_new_normal_honeycomb(
          normal_facet, f_center, h_v, new_center->vertex()->vertex_begin(), 
          order_v, new_vertices, normal_edges, plane_edges);
    } else if(mode <=1 || mode == 4) {
        cal_new_vertices_honeycomb(normal_facet, h_v, new_center->vertex()->vertex_begin(), 
          order_v, new_vertices, plane_edges);
    }
    if(is_new_vertices) {
        // create new vertices
        order_v = 0;
        do {
        	p.split_edge(h_v);
        	h_v->opposite()->vertex()->point() = new_vertices[order_v];
            ++ order_v;
        } while ( ++h_v != new_center->vertex()->vertex_begin());
        // create new facets
        order_v = 0;
        do {
        	p.split_facet(h_v->opposite(), h_v->opposite()->next()->next()->next());
            ++ order_v;
        } while ( ++h_v != new_center->vertex()->vertex_begin());
    }
    p.erase_center_vertex(h_v);
    return is_new_vertices;
}

void Honeycomb::pull_facets(Polyhedron_cgal& p) {
	CGAL_precondition(p.is_valid());
	std::size_t facets_num = p.size_of_facets();
    if ( facets_num <= 0) return;
    // new halfedges/facets are appended at the end.
    Edge_iterator last_e = p.edges_end();
    --last_e;
    Facet_iterator last_f = p.facets_end();
    --last_f;
    // cal normal of edges around facet
	std::map<Facet_handle, Point_3, Facet_cmp> facet_to_center_map;
	std::map<Halfedge_handle, Vector_3, Hedge_cmp> hedge_to_normal_map;
    Facet_iterator f = p.facets_begin();
    do {
        HF_circulator h_f = f->facet_begin();
        do {
            Halfedge_handle e = h_f;
            cal_edge_normal(e, hedge_to_normal_map, facet_to_center_map);
        } while ( ++h_f != f->facet_begin());
    } while ( f++ != last_f);
    // pull facets
    f = p.facets_begin();
    bool is_pull = true;
    int unpulled_facets_num = 0;
    int facet_index = 0;
//    go << p;
    do {
        is_pull = pull_facet( p, f, facet_to_center_map.at(f), hedge_to_normal_map);
        if(!is_pull) ++unpulled_facets_num;
    } while ( f++ != last_f);
//    	go.clearView();
    double unpulled = static_cast<double>(unpulled_facets_num)/static_cast<double>(facets_num);
    // TODO
    std::cout << "unpulled facets : " << unpulled << std::endl;
	// remove old edges
    Edge_iterator e = p.edges_begin();
    ++last_e;
    while ( e != last_e) {
        Halfedge_handle old_edge = e;
        ++e; // incr. before join since join destroys current edge
        // test the neighboring facet
        /*
        Vector_3 normal_0 = CGAL::unit_normal(
          old_edge->next()->vertex()->point(),
          old_edge->next()->next()->vertex()->point(),
		  old_edge->vertex()->point());
        Vector_3 normal_1 = CGAL::unit_normal(
          old_edge->opposite()->next()->vertex()->point(),
          old_edge->opposite()->next()->next()->vertex()->point(),
		  old_edge->opposite()->vertex()->point());
        double diff_sq = (normal_0-normal_1).squared_length();
        if(diff_sq < angle_prop*angle_prop){
            p.join_facet( old_edge);
        }
        */
        p.join_facet( old_edge);
    };
	CGAL_postcondition(p.is_valid());
}

void Honeycomb::honeycomb_iter(Polyhedron& polyhedron, const int iter, const int hc_mode, double alpha) {
    mode = hc_mode;
    iter_num = std::max(iter,1);
    update_alpha_facet_with_iter();
    update_alpha_edge_with_iter();
    alpha_facet = alpha;
    polyhedron.cal_facets_normal();
    pull_facets(polyhedron.get_p());
}

void Honeycomb::subdivid_edge(Polyhedron_cgal& p, Halfedge_handle& edge) {
    Vector_3 edge_vector = edge->vertex()->point() - edge->opposite()->vertex()->point();
    Point_3 subdivid_point = edge->opposite()->vertex()->point() + edge_vector*beta_edge;
    p.split_edge(edge);
    edge->opposite()->vertex()->point() = subdivid_point;
}

Point_3 Honeycomb::cal_new_vertex_kernel(Point_3& old_v, 
  HV_circulator& h_v, const HV_circulator& h_v_begin, 
  std::size_t neighbor_size, std::vector<Point_3>& neighbor_vertices, Vector_3& push_vector) {
    Ray_3 v_push_ray(old_v, old_v + push_vector);
    std::size_t order_v = 0;
    std::size_t sub_order_v = 0;
    Vector_3 min_vector(0.0, 0.0, 0.0);
    double min_dist_2 = inf_double;
    double min_dist = inf_double;
    double min_dist_local = inf_double;
    double avg_dist = 0.0;
    do {
        do {
        	boost::optional< boost::variant< Point_3, Ray_3 > > intersect_result = 
              CGAL::intersection(Plane_3(neighbor_vertices[order_v], 
                neighbor_vertices[(order_v+1)%neighbor_size], 
                neighbor_vertices[(order_v+2+sub_order_v)%neighbor_size]), v_push_ray);
            Point_3 intersect_point(inf_double, inf_double, inf_double);
            if(intersect_result) 
                intersect_point = boost::apply_visitor( intersect_visitor(), *intersect_result);
            if(intersect_point!=inf_point){
            	Vector_3 new_vector = intersect_point - old_v;
            	double dist_2 = new_vector.squared_length();
            	double dist = std::sqrt(dist_2);
            	if((dist>=0.0) && (dist<min_dist_local)){
            	    min_dist_local = dist;
            	    if(dist<min_dist){
                		min_dist = dist;
                		min_vector = new_vector;
            	    }
        	    }
            }
        } while(sub_order_v++ < (neighbor_size - 3));
        avg_dist += min_dist_local;
        min_dist_local = inf_double;
        sub_order_v = 0;
    	++ order_v;
    } while ( ++h_v != h_v_begin);
    if(mode == 0) {
        avg_dist = beta_vertex*min_dist;
    } else if(mode >= 1) {
        avg_dist = std::min(beta_vertex*avg_dist/static_cast<double>(neighbor_size), (1-tolerance)*min_dist);
    }
    return old_v + min_vector*avg_dist/min_dist;
}

void Honeycomb::push_vertex( Polyhedron_cgal& p, Vertex_iterator& v) {
    // store the vertices around v
    Vector_3 v_vector = v->point() - origin_point;
//  double v_origin_dist = std::sqrt(v_vector.squared_length());
	HV_circulator h_v = v->vertex_begin();
    std::vector<Point_3> neighbor_vertices;
    std::vector<Vector_3> neighbor_normals;
//	Vector_3 avg_vertex_vector = origin_vector;
	Vector_3 avg_normal_vector = origin_vector;
	double normal_weight = 0.0;
    std::size_t order_v = 0;
    do {
        Point_3 neighbor_vertex = h_v->opposite()->vertex()->point();
    	neighbor_vertices.push_back(neighbor_vertex);
//    	avg_vertex_vector = avg_vertex_vector + (neighbor_vertex-origin_point);
        Point_3 prev_neighbor_vertex = h_v->opposite()->next()->vertex()->point();
//        double edge_length = std::sqrt((neighbor_vertex - prev_neighbor_vertex).squared_length());
        double edge_length = (neighbor_vertex - prev_neighbor_vertex).squared_length();
    	Vector_3 neighbor_normal = h_v->opposite()->facet()->plane();
    	neighbor_normals.push_back(neighbor_normal);
    	normal_weight = mode == 3? edge_length*edge_length : 1.0;
    	avg_normal_vector = avg_normal_vector + neighbor_normal*normal_weight;
        ++ order_v;
    } while ( ++h_v != v->vertex_begin());
    std::size_t neighbor_size = neighbor_vertices.size();
    CGAL_assertion( order_v == neighbor_size);
//    avg_vertex_vector = avg_vertex_vector - origin_vector;
//    avg_vertex_vector = avg_vertex_vector/std::sqrt(avg_vertex_vector.squared_length());
    avg_normal_vector = avg_normal_vector - origin_vector;
    Point_3 old_vertex;
    Vector_3 push_vector;
    if(mode <= 1) {         // push old vertex towards the origin
        old_vertex = v->point();
        push_vector = -1*v_vector;
    }else if(mode == 2 || mode == 3) {
        // push old vertex along the avg_normal according to the surrounding facets
        old_vertex = v->point();
        push_vector = -1*avg_normal_vector;
        // push rotated vertex along the direction of avg_vertex according to the surrounding vertices
//        old_vertex = origin_point + avg_vertex_vector*v_origin_dist;
//        push_vector = -1*avg_vertex_vector;
    }
    // calculate the closest point cutting the origin line among all the planes for neighbor vertices
    Point_3 new_v = cal_new_vertex_kernel(old_vertex, h_v, v->vertex_begin(), neighbor_size, 
      neighbor_vertices, push_vector);
    // split facets around v if necessary
    order_v = 0;
    do {
        if(h_v->opposite()->next()->vertex()->point() != neighbor_vertices[(order_v-1)%neighbor_size]) {
    	    p.split_facet(h_v->opposite(), h_v->opposite()->prev()->prev());
        }
        ++ order_v;
    } while ( ++h_v != v->vertex_begin());
    // push vertex to the new position
    v->point() = new_v;
}

void Honeycomb::push_vertices(Polyhedron_cgal& p) {
	CGAL_precondition(p.is_valid());
    if ( p.size_of_vertices() <= 0) return;
    // new vertices/halfedges are appended at the end.
    Vertex_iterator last_v = p.vertices_end();
    --last_v;
    Edge_iterator last_e = p.edges_end();
    --last_e;
    // sub-divid edges
    Edge_iterator e = p.edges_begin();
    bool is_last_e = false;
    do {
        Halfedge_handle edge = e;
        // incr. before subdivid since subdivid destroys current edge
        if(e++ == last_e) is_last_e = true; 
        subdivid_edge(p, edge);
    } while (!is_last_e);
    // push vertices
    Vertex_iterator v = p.vertices_begin();
    do {
        push_vertex( p, v);
    } while ( v++ != last_v);
	CGAL_postcondition(p.is_valid());
}

void Honeycomb::kernel_iter(Polyhedron& polyhedron, const int iter, const int hc_mode, double beta) {
    mode = hc_mode;
    iter_num = std::max(iter,1) ;
    update_alpha_facet_with_iter();
    update_alpha_edge_with_iter();
    beta_vertex = beta;
    polyhedron.cal_facets_normal();
    push_vertices(polyhedron.get_p());
}

void Honeycomb::update_alpha_facet_with_iter() {
    alpha_facet = max_alpha_facet;
    beta_vertex = alpha_facet;
}

void Honeycomb::update_alpha_edge_with_iter() {
    alpha_edge = max_alpha_edge;
    beta_edge = max_beta_edge;
}

void Honeycomb::algo(Polyhedron& polyhedron) {
    GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, 0, "initial");
    FileIO::write_off_polyhedron(polyhedron);
    // TODO store the queue of polyhedron, for undo/redo
    
	int iter = 1;
    do {
        // TODO
	    std::cout << "Enter the mode ( h / k ) for next iteration!" << std::endl;
	    char ch;
	    std::cin >> ch;
	    if(ch=='h'){
            int honeycomb_iter_num = 1;
            while(honeycomb_iter_num-- > 0) {
//                double offset_fac = 0.5*(pow(1.1, static_cast<double>(iter-1)));
                honeycomb_iter(polyhedron, iter, 0, 0.8);
                polyhedron.update_nums();
                GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "honeycomb");
            }
        } else if(ch=='k'){
            int kernel_iter_num = 1;
            while(kernel_iter_num-- > 0) {
//                double offset_fac = 0.4*(pow(1.1, static_cast<double>(iter-1)));
                kernel_iter(polyhedron, iter, 0, 0.9);
                polyhedron.update_nums();
                GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "kernel");
            }
        }
        FileIO::write_off_polyhedron(polyhedron);
    } while(iter++ < max_iter);
}

