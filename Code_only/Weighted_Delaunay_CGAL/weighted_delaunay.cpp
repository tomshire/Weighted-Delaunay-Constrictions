#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <cassert>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>  Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Vertex_iterator                                 Vertex_iterator;
typedef Rt::Vertex_handle                                   Vertex_handle;
typedef Rt::Cell_handle                                     Cell_handle;

//~ Need the following additional lines
#include <iostream>
#include <fstream>

using std::cin;
using std::cout;
using std::ios;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::vector;

int main()
{
  /*~ Read in the data points one-by-one from a formatted text file
    called "Input_Data.txt". This text file contains a varying number
    of rows, each corresponding to a data point, and four columns.
    Columns 1-3 are the x, y and z coordinates and column 4 contains
    the radii (i.e., the weights).*/

  ifstream inputfile ("Input_Data.txt", ios::in);
  vector<Weighted_point> P;
  int number_of_points = 0;
  inputfile.precision(100); //~ An insanely high precision!
  double minrad = 1e100; //~ The minimum particle radius
  double maxrad = -1.0; //~ and the maximum radius

  if (!inputfile) {
    cerr << "Error: unable to open Input_Data.txt.\n";
    exit(1);
  } else {
    while(!inputfile.eof()) {
      //~ Find the smallest and largest particles in the system

      double x,y,z,r;
      inputfile >> x >> y >> z >> r;
      if (r < minrad) minrad = r;
      if (r > maxrad) maxrad = r;
    }

    inputfile.clear();
    inputfile.seekg( 0, std::ios::beg );

    while(!inputfile.eof()) {
      /*~ Note that this method reads in the last line twice;
	however this does not matter as the point is counted
	only once in the tessellation*/
      
      double x,y,z,r;
      inputfile >> x >> y >> z >> r;
      Point p(x, y, z); //~ The point is simply the coordinates
      Weight w = r*r; //~ and the weight is the radius squared
      
      P.push_back(Weighted_point(p, w));
      ++number_of_points;
    }
  }
  
  //~ Decrement the number of points so it is correct
  number_of_points--;

  //~ Find the vertices of the regular triangulation and store in T
  Rt T;
  T.insert (P.begin(), P.end());

  //~ Ensure that the output is valid and reasonable
  assert( T.is_valid() );
  assert( T.dimension() == 3 );

  //~ Delete the cells containing infinite vertices
  vector<Cell_handle> infinite_cells;
  T.incident_cells(T.infinite_vertex(),std::back_inserter(infinite_cells));

  for(int ii = 0; ii < infinite_cells.size(); ii++)
    {
      Cell_handle c = infinite_cells[ii];
      Cell_handle e = c->neighbor(c->index(T.infinite_vertex()));

      T.tds().delete_cell(c);
    }

  //~ Ensure that the output is still valid and reasonable
  assert( T.is_valid() );
  assert( T.dimension() == 3 );

  //~ Display the numbers of points, vertices and cells
  cout << "Number of data points:\t\t\t" << number_of_points << endl;
  cout << "Number of vertices:\t\t\t" << T.number_of_vertices() << endl;
  cout << "Number of infinite cells deleted:\t" << infinite_cells.size() << endl;
  cout << "Number of cells remaining:\t\t" << T.number_of_cells() << endl;

  /*~ "The information stored in the iostream is: the dimension,
    the number of vertices, the number of cells, the indices of
    the vertices of each cell, then the indices of the neighbors
    of each cell, where the index corresponds to the preceding 
    list of cells." - from the "TriangulationDataStructure_3"
    page of the manual.

    Write the iostream to a formatted text file (imaginatively
    titled "Output_Data.txt").

    Ordering of "Output_Data.txt" file:
    - Dimension (always 3)
    - Number of vertices
    - List of vertices
    - Number of cells
    - Indices of the vertices of each cell
    - Indices of the neighbours of each cell
  */
  ofstream outputfile ("Output_Data.txt", ios::out);
  outputfile.setf(ios::scientific,ios::floatfield);
  outputfile.precision(16);
  // outputfile << minrad << endl;
  outputfile << T;
 
  return 0;
}
