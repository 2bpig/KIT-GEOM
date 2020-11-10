// algo
#include "Algo_Edgecutting.hpp"
#include "Algo_MaxMatching.hpp"

// declare constants
const std::size_t max_iter = 3;//8;

class Edgecutting48 {
private:
	GeomOutput& go;
	Polyhedron& polyhedron;
    std::size_t mode;
    std::size_t iter_num;
    // parameters
    // edgecutting 4-8 scheme
    void prime();
    void dual();
public:
    Edgecutting48(GeomOutput& geomoutput, Polyhedron& polyhedron);
    void algo();
};
