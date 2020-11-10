#include "Data_Polyhedron.hpp"

////////////////////////////////////////////////
//    private:                                //
////////////////////////////////////////////////

void Polyhedron::trans_to_origin() {
    Vertex_iterator vertex_iter = p.vertices_begin();
    Vector_3 vec( 0.0, 0.0, 0.0);
    std::size_t vec_num = 0;
    do {
        vec = vec + ( vertex_iter->point() - CGAL::ORIGIN);
        ++ vec_num;
    } while ( ++vertex_iter != p.vertices_end());
    if(vec_num!=vertex_num){
        std::cout << "unexpected update of vertex num" << std::endl;
        update_nums();
    }
    vec = (CGAL::ORIGIN + (vec / static_cast<double>(vec_num))) - origin_point;
    Vertex_iterator trans_vertex_iter = p.vertices_begin();
    do {
        trans_vertex_iter->point() = trans_vertex_iter->point() - vec;
    } while ( ++trans_vertex_iter != p.vertices_end());
}

////////////////////////////////////////////////
//    constructor:                            //
////////////////////////////////////////////////

Polyhedron::Polyhedron() {
    vertex_num = 0;
    edge_num = 0;
    facet_num = 0;
}

////////////////////////////////////////////////
//    public:                                 //
////////////////////////////////////////////////

void Polyhedron::clear() {
    p.clear();
    vertex_num = 0;
    edge_num = 0;
    facet_num = 0;
}

Halfedge_handle Polyhedron::init_tetrahedron() {
	CGAL_precondition(p.is_valid());
	
	double x = 50.0;
	double y = 43.3;
	double z = 40.8;
    Point v0(  0.0,  0.0,  0.0 ); // 0 0 0
    Point v1(    x,  0.0,  0.0 ); // x 0 0
    Point v2(x*0.5,    y,  0.0 ); // 0 y 0
    Point v3(x*0.5,y/3.0,    z ); // 0 0 z
	Halfedge_handle h0 = p.make_tetrahedron(v0, v1, v3, v2)->next();
    vertex_num = 4;
    edge_num = 6;
    facet_num = 4;
	
    trans_to_origin();
	CGAL_postcondition(p.is_valid());
	return h0;
}

Halfedge_handle Polyhedron::init_cube() {
	CGAL_precondition(p.is_valid());
	
	double x = 200.0;
	double y = 100.0;
	double z = 50.0;
    Point v0(  0.0,  0.0,  0.0 ); // 0 0 0
    Point v1(    x,  0.0,  0.0 ); // x 0 0
    Point v2(  0.0,    y,  0.0 ); // 0 y 0
    Point v3(  0.0,  0.0,    z ); // 0 0 z
    Point v4(    x,  0.0,    z ); // x 0 z
    Point v5(  0.0,    y,    z ); // 0 y z
    Point v6(    x,    y,  0.0 ); // x y 0
    Point v7(    x,    y,    z ); // x y z
	Halfedge_handle h0 = p.make_tetrahedron(v0, v1, v3, v2)->next();
	Halfedge_handle h1 = h0->next()->opposite()->next();
	p.split_edge(h0->next());
	p.split_edge(h1->next());
	p.split_edge(h1);
	h0->next()->vertex()->point() = v4;
	h1->next()->vertex()->point() = v5;
	h1->opposite()->vertex()->point() = v6;
	Halfedge_handle h2 = p.split_facet(h1->next(), h1->next()->next()->next());
	Halfedge_handle h3 = p.split_edge(h2);
	h3->vertex()->point() = v7;
	p.split_facet(h3, h2->next()->next());
    vertex_num = 8;
    edge_num = 12;
    facet_num = 6;
	
    trans_to_origin();
	CGAL_postcondition(p.is_valid());
	return h0;
}

Halfedge_handle Polyhedron::init_octahedron() {
	CGAL_precondition(p.is_valid());
	
	double x = 50.0;
	double y = 50.0;
	double z = 35.4;
    Point v0(  0.0,  0.0,  0.0 ); // 0 0 0
    Point v1(    x,  0.0,  0.0 ); // x 0 0
    Point v2(  0.0,    y,  0.0 ); // 0 y 0
    Point v3(    x,    y,  0.0 ); // x y 0
    Point v4(x*0.5,y*0.5,    z ); // x y z
    Point v5(x*0.5,y*0.5,   -z ); // x y-z
	Halfedge_handle h0 = p.make_tetrahedron(v0, v1, v4, v2);
	Halfedge_handle h1 = h0->prev()->opposite()->next();
	p.split_edge(h1);
	h1->opposite()->vertex()->point() = v3;
	p.split_facet(h1->next(), h1->prev());
    Halfedge_handle h2 = p.create_center_vertex( h1->opposite()); 
    h2->vertex()->point() = v5;
    vertex_num = 6;
    edge_num = 12;
    facet_num = 8;
    
    trans_to_origin();
	CGAL_postcondition(p.is_valid());
	return h0;
}

Halfedge_handle Polyhedron::init_icosahedron() {
	CGAL_precondition(p.is_valid());
	
	std::vector<Point_3> points;
	double factor = 50.0;
	double golden_ratio = (1.0 + std::sqrt(5.0))/2.0;
    Point v0(  0.0+factor*golden_ratio,  factor+factor*golden_ratio,  factor*golden_ratio ); //
    Point v1(  0.0+factor*golden_ratio, -factor+factor*golden_ratio,  factor*golden_ratio ); //
    Point v2(  0.0+factor*golden_ratio,  factor+factor*golden_ratio, -factor*golden_ratio ); //
    Point v3(  0.0+factor*golden_ratio, -factor+factor*golden_ratio, -factor*golden_ratio ); //
    
    Point v4(  factor*golden_ratio+factor*golden_ratio,  0.0+factor*golden_ratio,  factor ); //
    Point v5(  factor*golden_ratio+factor*golden_ratio,  0.0+factor*golden_ratio, -factor ); //
    Point v6( -factor*golden_ratio+factor*golden_ratio,  0.0+factor*golden_ratio,  factor ); //
    Point v7( -factor*golden_ratio+factor*golden_ratio,  0.0+factor*golden_ratio, -factor ); //
    
    Point v8(  factor+factor*golden_ratio,  factor*golden_ratio+factor*golden_ratio,  0.0 ); //
    Point v9( -factor+factor*golden_ratio,  factor*golden_ratio+factor*golden_ratio,  0.0 ); //
    Point v10( factor+factor*golden_ratio, -factor*golden_ratio+factor*golden_ratio,  0.0 ); //
    Point v11(-factor+factor*golden_ratio, -factor*golden_ratio+factor*golden_ratio,  0.0 ); //
	
	points.push_back(v0);
	points.push_back(v1);
	points.push_back(v2);
	points.push_back(v3);
	points.push_back(v4);
	points.push_back(v5);
	points.push_back(v6);
	points.push_back(v7);
	points.push_back(v8);
	points.push_back(v9);
	points.push_back(v10);
	points.push_back(v11);
	
	CGAL::convex_hull_3(points.begin(), points.end(),p);
	Halfedge_handle h0 = p.halfedges_begin();
    vertex_num = 12;
    edge_num = 30;
    facet_num = 20;
	
    trans_to_origin();
	CGAL_postcondition(p.is_valid());
	return h0;
}

Halfedge_handle Polyhedron::init_random_polyhedron() {
    //interactive getting the num of vertex
	std::cout << "Enter the number of points !" << std::endl;
	int num;
	std::cin >> num;
	//init process
	CGAL_precondition(p.is_valid());
	
    std::srand(std::time(0)); // use current time as seed for random generator
	std::vector<Point_3> points;
	Random_p_sphere g(200.0);
	CGAL::cpp11::copy_n(g, num, std::back_inserter(points));
	double a = 1.0;
	double b = 1.0;
	for(std::size_t i = 0; i<points.size(); ++i) {
	    double random = std::rand() /(RAND_MAX+1u);
	    points[i] = Primitive::random_scale_point(points[i], random, a, b);
	}
	CGAL::convex_hull_3(points.begin(), points.end(),p);
	Halfedge_handle h0 = p.halfedges_begin();
	update_nums();
	
    trans_to_origin();
	CGAL_postcondition(p.is_valid());
	return h0;
}

bool Polyhedron::update_nums() {
    std::size_t v_num = p.size_of_vertices();
    std::size_t e_num = p.size_of_halfedges()/2;
    std::size_t f_num = p.size_of_facets();
    bool is_changed = false;
    if(vertex_num!=v_num){
        is_changed = true;
        vertex_num = v_num;
    }
    if(edge_num!=e_num){
        is_changed = true;
        edge_num = e_num;
    }
    if(facet_num!=f_num){
        is_changed = true;
        facet_num = f_num;
    }
    return is_changed;
}

void Polyhedron::split_facet_in_center() {
    CGAL_precondition(p.is_valid());
    if ( p.size_of_vertices() <= 0) return;
    // new vertices/halfedges are appended at the end.
    Facet_iterator last_f = p.facets_end();
    --last_f;
    Facet_iterator f = p.facets_begin();
    Point_3 f_center;
    Halfedge_handle new_center;
    do {
        f_center = cal_facet_center(f);
        new_center = p.create_center_vertex( f->halfedge()); 
        new_center->vertex()->point() = f_center;
    } while ( f++ != last_f);
	CGAL_postcondition(p.is_valid());
}

void Polyhedron::cal_facets_normal() {
    std::transform( p.facets_begin(), p.facets_end(), p.planes_begin(), Normal_vector());
}

////////////////////////////////////////////////
//    static:                                 //
////////////////////////////////////////////////

Point_3 Polyhedron::cal_facet_center(Facet_iterator& f) {
    Vector_3 sum_vertices( 0.0, 0.0, 0.0 );
    std::size_t valence_f = 0;
    HF_circulator h_f = f->facet_begin();
    do {
        sum_vertices = sum_vertices + ( h_f->vertex()->point() - CGAL::ORIGIN);
        ++ valence_f;
    } while ( ++h_f != f->facet_begin());
    CGAL_assertion(valence_f >= 3);
    return CGAL::ORIGIN + (sum_vertices / static_cast<double>(valence_f));
}

