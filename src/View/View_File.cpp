#include "View_File.hpp"

    /*
void ExtractHDS4STL::operator()( HDS& hds);
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
    */

////////////////////////////////////////////////
//    constructor:                            //
////////////////////////////////////////////////

FileIO::FileIO() {

}

////////////////////////////////////////////////
//    private:                                //
////////////////////////////////////////////////
bool FileIO::read_stl( std::istream& input,
  std::vector< std::array<double,3> >& points,
  std::vector< std::array<int,3> >& facets) {
  /*
    bool verbose = false;
    std::string s, solid("solid");
    std::map<std::array<double,3>, int> pmap;
    int index = 0;
    std::array<int,3> ijk;
    std::array<double,3> p;

    char line[80];
    for(int i=0;i < 80; i++){
      boost::uint8_t c;
      input.read(reinterpret_cast<char*>(&c), sizeof(c));
      line[i]=c;
      if(i==5){
        s = std::string(line,5);
        if(s == solid){
          break;
        }
      }
    }

    if(s!= solid){
      boost::uint32_t N32;
      input.read(reinterpret_cast<char*>(&N32), sizeof(N32));
      unsigned int N = N32;

      for(unsigned int i=0; i < N; i++){
        float normal[3];
        input.read(reinterpret_cast<char*>(&normal[0]), sizeof(normal[0]));
        input.read(reinterpret_cast<char*>(&normal[1]), sizeof(normal[1]));
        input.read(reinterpret_cast<char*>(&normal[2]), sizeof(normal[2]));

        for(int j=0; j < 3; j++){
          float x,y,z;
          input.read(reinterpret_cast<char*>(&x), sizeof(x));
          input.read(reinterpret_cast<char*>(&y), sizeof(y));
          input.read(reinterpret_cast<char*>(&z), sizeof(z));
          p[0]=x; p[1]=y; p[2]=z;
          std::map<std::array<double,3>, int>::iterator iti=
            pmap.insert(std::make_pair(p,-1)).first;
          if(iti->second==-1){
            ijk[j] = index;
            iti->second = index++;
            points.push_back(p);
          } else {
            ijk[j] = iti->second;
          }
        }
        if((ijk[0] != ijk[1]) && 
           (ijk[0] != ijk[2]) &&
           (ijk[1] != ijk[2])){
          facets.push_back(ijk);
        }else{
          if(verbose){
            std::cerr << "ignore degenerate face" << std::endl;
          }
        }
        char c;
        input.read(reinterpret_cast<char*>(&c), sizeof(c));
        input.read(reinterpret_cast<char*>(&c), sizeof(c));
      }
      return true;
    } else {
      std::string facet("facet"),
        outer("outer"),
        loop("loop"),
        vertex("vertex"),
        endloop("endloop"),
        endsolid("endsolid");

      while(input >> s){
        if(s == endsolid){
          //std::cerr << "found endsolid" << std::endl;
        } else if(s == facet){
          //std::cerr << "found facet" << std::endl;
          std::getline(input, s); // ignore the normal
          input >> s;
          if(s != outer){
            if (verbose)
              std::cerr << "Expect 'outer' and got " << s << std::endl;
            return false;
          }
          input >> s;
          if(s != loop){
            if (verbose)
              std::cerr << "Expect 'loop' and got " << s << std::endl;
            return false;
          }
          int count = 0;
          do {
            input >> s;
            if(s == vertex){
              //      std::cerr << "found vertex" << std::endl;
              if(count < 3){
                input >> p[0] >> p[1] >> p[2];
                std::map<std::array<double,3>, int>::iterator iti=
                  pmap.insert(std::make_pair(p,-1)).first;
                if(iti->second==-1){
                  ijk[count] = index;
                  iti->second = index++;
                  points.push_back(p);
                } else {
                  ijk[count] = iti->second;
                }
                ++count;
              } else {
                if (verbose)
                  std::cerr << "We can only read triangulated surfaces" << std::endl;
                return false;
              }
            }
          }while(s != endloop);
          
          facets.push_back(ijk);
        }
      }
      return true;
    }
    */
}

////////////////////////////////////////////////
//    public:                                 //
////////////////////////////////////////////////

////////////////////////////////////////////////
//    static:                                 //
////////////////////////////////////////////////

bool FileIO::split_polyhedron_inFile() {
    Polyhedron p;
    if(!read_off_polyhedron(p)) return false;
    p.split_facet_in_center();
    if(!write_off_polyhedron(p)) return false;
    return true;
}

bool FileIO::write_off_polyhedron(Polyhedron& p) {
	std::cout << "Enter the output filename.off :" << std::endl;
	char path_off[50];
	std::cin >> path_off;
	if(strlen(path_off)>4) {
        std::ofstream outoff(path_off);
//      std::ofstream outoff(filename_off);
        outoff << p.get_p();
        outoff.close();
        return true;
    }
    return false;
}
      
bool FileIO::read_off_polyhedron(Polyhedron& p) {
    std::cout << "Enter the input filename.off :" << std::endl;
	char path_off[50];
	std::cin >> path_off;
	if(strlen(path_off)>4) {
        std::ifstream inoff(path_off);
//      std::ifstream inoff(filename_off);
        inoff >> p.get_p();
        inoff.close();
        p.update_nums();
        return true;
    }
    return false;
}

void FileIO::write_obj_polyhedron(Polyhedron& p) {
    std::ofstream outobj(filename_obj);
    outobj << p.get_p();
    outobj.close();
}
      
void FileIO::read_obj_polyhedron(Polyhedron& p) {
    std::ifstream inobj(filename_obj);
    inobj >> p.get_p();
    inobj.close();
}

void FileIO::write_stl_polyhedron(Polyhedron& p) {
    CGAL::Polygon_mesh_processing::triangulate_faces(p.get_p());
    std::ofstream outstl(filename_stl);
    
    outstl.close();
}

void FileIO::read_stl_polyhedron(Polyhedron& p) {

    typedef typename std::array<double,3> p_point;
    typedef typename std::array<int,3> p_facet;
    
    std::ifstream instl(filename_stl);
    
    std::vector<p_point> p_points;
    std::vector<p_facet> p_facets;
    
//    read_stl(instl, p_points, p_facets);
    instl.close();
    
	std::cout << "vertex num : " << p_points.size() << std::endl;
	std::cout << "facet num : " << p_facets.size() << std::endl;
//    print_polyhedron_info(p);
}

void FileIO::print_polyhedron_info(Polyhedron& p) {
	std::cout << "facet num : " << p.get_facet_num() << std::endl;
	std::cout << "edge num : " << p.get_edge_num() << std::endl;
	std::cout << "vertex num : " << p.get_vertex_num() << std::endl;
}

/*
    // Write polyhedron in Object File Format (OFF).
    CGAL::set_ascii_mode( std::cout);
    std::cout << "OFF" << std::endl << P.size_of_vertices() << ' '
              << P.size_of_facets() << " 0" << std::endl;
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<Point_3>( std::cout, "\n"));
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        std::cout << CGAL::circulator_size(j) << ' ';
        do {
            std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        std::cout << std::endl;
    }
*/
