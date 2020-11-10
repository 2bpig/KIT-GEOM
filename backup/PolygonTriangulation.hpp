// Copyright 2018
//   Pawel Herman - Karlsruhe Institute of Technology - pherman@ira.uka.de
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include "happah/geometry/TriangleMesh.hpp"

namespace happah {

//DECLARATIONS

template<class Vertex>
class PolygonGraph;

namespace pog {

struct HalfEdge;

struct PolygonEnumerator;

struct SpokesWalker;
//TODO remove
template<class Vertex>
struct Domain;

}//namespace pog

pog::PolygonEnumerator make_polygon_enumerator(const std::vector<pog::HalfEdge>& edges, const Indices& offsets);

template<class Vertex>
PolygonGraph<Vertex> make_polygon_graph(std::vector<Vertex> vertices, std::vector<pog::HalfEdge> edges, Indices offsets);

template<class Vertex>
TriangleMesh<Vertex> make_triangle_mesh(const PolygonGraph<Vertex>& graph);

namespace detail {
//TODO move
using Point = typename Vertex::SPACE::POINT;
//TODO refactory parameters
template<class Vertex>
float eval_domain(pog::Domain<Vertex>& domain);

}//namespace detail

template<class Vertex>
hpuint size(const PolygonGraph<Vertex>& graph);

//DEFINITIONS

namespace pog {

struct HalfEdge {
     HalfEdge(int8_t offset, int8_t next, int8_t previous, hpindex adjacentFace, int8_t opposite, hpindex oppositeFace, hpindex vertex)
          : m_adjacentFace(adjacentFace), m_next(next), m_offset(offset), m_opposite(opposite), m_previous(previous), m_oppositeFace(oppositeFace), m_vertex(vertex) {}

     auto getId() const { return m_offset + m_adjacentFace; }

     auto getNext() const { return m_next + m_adjacentFace; }

     auto getOpposite() const { return m_opposite + m_oppositeFace; }

     auto getPrevious() const { return m_previous + m_adjacentFace; }

     auto& getVertex() const { return m_vertex; }

private:
     hpindex m_adjacentFace;//index in m_edges to first edge of adjacent face
     int8_t m_next;//offset to next edge from m_adjacentFace
     int8_t m_offset;//offset to current edge from m_adjacentFace
     int8_t m_opposite;//offset to opposite edge from m_oppositeFace
     int8_t m_previous;//offset to previous edge from m_adjacentFace
     hpindex m_oppositeFace;//index in m_edges to first edge of opposite face
     hpindex m_vertex;//index in m_vertices to vertex to which edge points

};//HalfEdge

}//namespace pog

//NOTE: Edges are arranged counter-clockwise and edges of the same face are stored together sequentially.
template<class Vertex>
struct PolygonGraph {
     PolygonGraph(std::vector<Vertex> vertices, std::vector<pog::HalfEdge> edges, Indices offsets)
          : m_edges(std::move(edges)), m_offsets(std::move(offsets)), m_vertices(std::move(vertices)) {}

     auto& getEdges() const { return m_edges; }

     auto getNumberOfPolygons() const { return hpuint(m_offsets.size() - 1); }

     auto& getOffsets() const { return m_offsets; }

     auto& getVertex(hpuint v) const { return m_vertices[v]; }

     auto& getVertices() const { return m_vertices; }

private:
     std::vector<pog::HalfEdge> m_edges;
     Indices m_offsets;
     std::vector<Vertex> m_vertices;

};//PolygonGraph

namespace pog {

class PolygonEnumerator {
public:
     PolygonEnumerator(const std::vector<HalfEdge>& edges, const Indices& offsets)
          : m_edges(edges), m_i(std::begin(offsets)), m_end(std::end(offsets) - 1) {}

     explicit operator bool() const { return m_i != m_end; }

     auto operator*() const {
          auto temp = [](auto i) { return (*i).getVertex(); };

          return std::make_tuple(transform<HalfEdge>(std::begin(m_edges) + m_i[0], temp), transform<HalfEdge>(std::begin(m_edges) + m_i[1], temp));
     }

     auto& operator++() {
          ++m_i;
          return *this;
     }

private:
     using Iterator = typename Indices::const_iterator;

     const std::vector<HalfEdge>& m_edges;
     Iterator m_end;
     Iterator m_i;

};//PolygonEnumerator

struct SpokesWalker {
     SpokesWalker(const std::vector<HalfEdge>& edges, hpindex e)
          : m_e(e), m_edges(edges) {}

     auto operator==(const SpokesWalker& walker) const { return m_e == walker.m_e; }
     
     auto operator!=(const SpokesWalker& walker) const { return !(*this == walker); }
     
     auto operator*() const { return m_e; }

     auto& operator++() {
          m_e = m_edges[m_edges[m_e].getPrevious()].getOpposite();
          return *this;
     }

     auto& operator--() {
          m_e = m_edges[m_edges[m_e].getOpposite()].getNext();
          return *this;
     }

     auto& operator+=(hpuint n) {
          while(n--) ++(*this);
          return *this;
     }

     auto& operator-=(hpuint n) {
          while(n--) --(*this);
          return *this;
     }

     auto operator+(hpuint n) const {
          auto copy = *this;
          return copy += n;
     }

     auto operator-(hpuint n) const {
          auto copy = *this;
          return copy -= n;
     }

private:
     hpindex m_e;
     const std::vector<HalfEdge>& m_edges;

};//SpokesWalker

//TODO
template<class Vertex>
struct Domain {
    Domain(std::vector<Vertex>& d_vertices, int d_key)
        : v_list(d_vertices), key(d_key) {}

    Domain(std::vector<Vertex>& d_vertices, int begin_vid, int end_vid)
        : v_list(d_vertices) {
        key = hash_vids(begin_vid, end_vid);
    }

    //TODO
    int hash_vids(int begin_vid, int end_vid) {
        //Cantor pairing function
        //hash two integers(0-100) to one(0-40300)
        return (begin_vid+end_vid)*(begin_vid+end_vid+1)/2 + begin_vid;
        //Enter Szudzik's function
        //hash two integers(0-100) to one(0-10200)
        //return begin_vid<=end_vid ? end_vid*end_vid+end_vid+begin_vid : end_vid+begin_vid*begin_vid;
    }

    auto& getVList() const { return v_list; }

    std::size_t getSize() const { return v_list.size(); }

    bool isEmpty() const { return getSize()<=2; }

    int getKey() const { return key; }

    void setWeight(float d_w) { weight = d_w; }

    float getWeight() { return weight; }

private:
    std::vector<Vertex> v_list;
    int key;
    float weight;

};//Domain

}//namespace pog

pog::PolygonEnumerator make_polygon_enumerator(const std::vector<pog::HalfEdge>& edges, const Indices& offsets) { return { edges, offsets }; }

template<class Vertex>
PolygonGraph<Vertex> make_polygon_graph(std::vector<Vertex> vertices, std::vector<pog::HalfEdge> edges, Indices offsets) { return { std::move(vertices), std::move(edges), std::move(offsets) }; }

template<class Vertex, class Weighter0, class Weighter1>
TriangleMesh<Vertex> make_triangle_mesh(const PolygonGraph<Vertex>& graph, Weighter0&& weighter0, Weighter1&& weighter1) {
     auto indices = Triples<hpindex>();

    auto eval_triangle = [&](auto& point0, auto& point1, auto& splitPoint){
        //TODO
        return (hpreal)1.0;
    };

    auto eval_bitriangle = [&](auto& point0, auto& point1, auto& splitPoint, auto& domainPoint){
        //TODO
        return (hpreal)0.0;
    };

     visit(make_polygon_enumerator(graph.getEdges(), graph.getOffsets()), [&](auto begin, auto end) {
         //TODO:
    	 //compute triangle's vertex indices for face [begin, end];
    	 //*begin returns index to vertex;
    	 //use begin, end just as any random access iterator
         
        //init the indices of vertices within the face [begin, end]
        std::vector<hpindex> v_indices;
        while(begin!=end) {
            v_indices.push_back(*begin);
            begin++;
        }
        hpindex begin_index = 0;
        hpindex end_index = v_indices.size()-1;
        //init root domain
        //pog::Domain<Vertex> root_domain(d_vertices, begin_vid, end_vid);
        //std::unordered_map<int, pog::Domain<Vertex>> domain_map;
        //domain_map.insert({root_domain.getKey(), root_domain});
        //init the hashmap <<d_begin_id, d_end_id>, <d_weight, d_split_id>>
        auto begin_end_2_w_sp_map = make_map<std::pair<hpreal, hpindex>>(0);

        //evaluate top-down from root-domain
        eval_domain(graph, begin_index, end_index, null, v_indices.begin(), eval_triangle, eval_bitriangle, begin_end_2_w_sp_map, (hpreal)std::numeric_limits<float>::infinity());
     });

     return make_triangle_mesh(graph.getVertices(), std::move(indices));
}

namespace detail {

using Map = std::unordered_map<std::pair<hpuint, hpuint>, std::pair<hpreal, hpindex>>;

template<class Vertex, class Weighter0, class Weighter1>
hpreal eval_domain(const PolygonGraph<Vertex>& graph, const hpindex begin_index, const hpindex end_index, const hpindex dv_index, std::vector<hpindex>::iterator& vid_begin, Weighter0&& weighter0, Weighter1&& weighter1, Map& w_sp_map, hpreal& min_global_eval) {
    hpindex d_size = end_index - begin_index;
    if(d_size <= 2) return (hpreal)0.0;
    using Point = typename Vertex::SPACE::POINT;
    auto getP_withId = [&](auto index) {
        auto& vertex = graph->getVertex(*(vid_begin+index));
        return vertex->position;
    };
    Point p0 = getP_withId(end_index);
    Point p1 = getP_withId(begin_index);
    Point p_d = getP_withId(dv_index);

    hpreal w_lo = 0.0;
    hpreal w_hi = 0.0;
    hpreal w_tri = 0.0;
    hpreal w_bitri = 0.0;
    hpreal w_d = 0.0;
    hpreal min_w_d = (hpreal)std::numeric_limits<float>::infinity();
    hpindex min_split_index = 1;
    for( hpindex split_index = begin_index+1; split_index < end_index; ++split_index) {
        //evaluate low domain
        auto key_lo = std::make_pair(begin_index, split_index);
        auto w_sp_lo = w_sp_map.find(key_lo);
        if(w_sp_lo == w_sp_map.end()) {
            w_lo = eval_domain(graph, begin_index, split_index, end_index, vid_begin, weighter0, weighter1, w_sp_map, min_global_eval);
        }else{
            w_lo = w_sp_lo->second.first;
        }
        //evaluate high domain
        auto key_hi = std::make_pair(split_index, end_index);
        auto w_sp_hi = w_sp_map.find(key_hi);
        if(w_sp_hi == w_sp_map.end()) {
            w_hi = eval_domain(graph, split_index, end_index, begin_index, vid_begin, weighter0, weighter1, w_sp_map, min_global_eval);
        }else{
            w_hi = w_sp_hi->second.first;
        }
        //evaluate domain
        Point p_split = getP_withId(split_index);
        w_tri = weighter0(p0, p1, p_split);//eval_triangle(tri_vertices);
        w_bitri = weighter1(p0, p1, p_split, p_d);//eval_bitriangle(bitri_vertices);
        w_d = w_tri + w_bitri + w_lo + w_hi;
        //update the min state
        if(w_d<min_w_d) {
            min_w_d=w_d;
            min_split_index=split_index;
        }
    }
    auto key = std::make_pair(begin_index, end_index);
    begin_end_2_w_sp_map[key] = std::make_pair(min_w_d, min_split_index);
    return min_w_d;
}

}//namespace detail

template<class Vertex>
hpuint size(const PolygonGraph<Vertex>& graph) { return graph.getNumberOfPolygons(); }

}//namespace happah

