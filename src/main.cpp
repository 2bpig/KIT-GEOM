//#include <config.h>

// algo
//#include "Algo/Honeycomb/Algo_Honeycomb.hpp"
#include "Algo/Cornercutting/Algo_Edgecutting48.hpp"



int main()
{
	// setup geomview
	GeomOutput geomoutput0(0);
	if(!geomoutput0.getUsability()) return 0;
	
	// build polyhedron
    Polyhedron p0;
    p0.init_icosahedron();//init_tetrahedron();//init_cube();//
	// algo
    Edgecutting48 ec48(geomoutput0, p0);
    ec48.algo();
        
    // clear geomview and exit
    geomoutput0.clearView();
    
	return 0;
}
