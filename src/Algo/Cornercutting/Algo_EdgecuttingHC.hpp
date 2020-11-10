// algo

// declare constants
const std::size_t max_iter = 7;

class EdgecuttingHC {
private:
	GeomOutput& go;
	Polyhedron& polyhedron;
    std::size_t mode;
    std::size_t iter_num;
    // parameters
    double alpha;
    // edgecutting honeycomb scheme
    void prime();
public:
    EdgecuttingHC(GeomOutput& geomoutput, Polyhedron& polyhedron);
    void algo();
};
