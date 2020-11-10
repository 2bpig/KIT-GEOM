#include "Algo_Edgecutting48.hpp"

Edgecutting48::Edgecutting48(GeomOutput& geomoutput, Polyhedron& polyhedron) :
  go(geomoutput), polyhedron(polyhedron) {
    mode = 0;
    iter_num = 0;
}

void Edgecutting48::algo() {
    GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, 0, "initial");
    FileIO::write_off_polyhedron(polyhedron);
    // TODO store the queue of polyhedron, for undo/redo
    
    dual();//prime();//
    
    polyhedron.update_nums();
    GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, 0, "final");
    FileIO::write_off_polyhedron(polyhedron);
}

void Edgecutting48::prime() {
    polyhedron.cal_facets_normal();
	Edgecutting ec(go, polyhedron);
	ec.init_sphere_r();
    ec.update_depth_ratio(0.3);
    
	MaxMatching mm(go, polyhedron);
    std::vector<Halfedge_handle> M;
    std::vector<Halfedge_handle> N;
    std::vector<Halfedge_handle> new_M;
    std::vector<Halfedge_handle> new_N;
    mm.match(new_M, new_N);
	std::size_t iter = 1;
    do {
        M.swap(new_M);
        N.swap(new_N);
        new_M.clear();
        new_N.clear();
        ec.update_cal_depth(true);
        ec.update_depth_ratio(0.25+0.2*iter/10.0);//std::sqrt(iter/10.0)
        for(Halfedge_handle e_M : M) {
            if(!ec.cut_edge(e_M, new_M)) return;
            ec.update_cal_depth(false);
//        polyhedron.update_nums();
//        GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "e");
        }
        ec.update_cal_depth(true);
        for(Halfedge_handle e_N : N) {
            if(!ec.cut_edge(e_N, new_N)) return;
            ec.update_cal_depth(false);
//        polyhedron.update_nums();
//        GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "e");
        }
        polyhedron.update_nums();
        GeomOutput::draw_polyhedron(polyhedron, polyhedron, go, go, iter, "p");
        FileIO::write_off_polyhedron(polyhedron);
        
    } while(iter++ < max_iter);
}

void Edgecutting48::dual() {
    // dualize polyhedron to the dual space
    Dualizing dualizing(go, polyhedron);
    Polyhedron dual_polyhedron;
    dualizing.convexDual(dual_polyhedron);
    dual_polyhedron.update_nums();
    GeomOutput::draw_polyhedron(dual_polyhedron, dual_polyhedron, go, go, 0, "dual");
    // 48 in dual space
    dual_polyhedron.cal_facets_normal();
	Edgecutting ec(go, dual_polyhedron);
	ec.init_sphere_r_dual();
    
	MaxMatching mm(go, dual_polyhedron);
    std::vector<Halfedge_handle> M;
    std::vector<Halfedge_handle> N;
    std::vector<Halfedge_handle> new_M;
    std::vector<Halfedge_handle> new_N;
    mm.match(new_M, new_N);
	std::size_t iter = 1;
    do {
        M.swap(new_M);
        N.swap(new_N);
        new_M.clear();
        new_N.clear();
        for(Halfedge_handle e_M : M) {
            if(!ec.interpolate_vertex4edge(e_M, new_M)) goto end_48_dual;
        }
        for(Halfedge_handle e_N : N) {
            if(!ec.interpolate_vertex4edge(e_N, new_N)) goto end_48_dual;
        }
        dual_polyhedron.update_nums();
        GeomOutput::draw_polyhedron(dual_polyhedron, dual_polyhedron, go, go, iter, "e");
        FileIO::write_off_polyhedron(dual_polyhedron);
    } while(iter++ < max_iter);
end_48_dual:
    // dualize the polyhedron back to the prime space
    Dualizing dualizing_dual(go, dual_polyhedron);
    polyhedron.clear();
    dualizing_dual.convexDual(polyhedron);
}
