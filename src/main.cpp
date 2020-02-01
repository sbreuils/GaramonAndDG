#include <iostream>
#include <cstdlib>
#include <list>
#include <c2ga/Mvec.hpp>

#include "c2gaTools.hpp"
#include "Geogebra_c2ga.hpp"

#include "Entry.hpp"

using multivector = c2ga::Mvec<double>;

using namespace c2ga;

/*As suggested, this function merely visualizes
four points defined in Conformal Geometric Algebra */

void vizuFourPoints(){
  Viewer_c2ga viewer; /* visualized used*/
  
  multivector p1,p2,p3,p4;
  
  p1 = point(1.0,0.0);
  p2 = point(0.0,1.0);
  p3 = point(0.0,0.0);
  p4 = point(1.0,1.0);
  
  std::vector<multivector> vectorOfPoints = {p1,p2,p3,p4};
  
  for(multivector point: vectorOfPoints) /*same thing as for point in objectsToBeVisualized*/
    viewer.push(point,"",255,0,0); /* arguments are: multivector, description, R, V, B */
  
  viewer.display();
  viewer.render("output.html"); /* resulting file in build/output.html*/
}




//void


int main(){
    
  //vizuFourPoints() ;

  //discretization() ;

  //rayThroughCircle() ;

  //snail() ;

  vizuFourPoints();
  
  return 0;
}
