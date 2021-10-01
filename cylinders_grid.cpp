#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//Generates cartesian grid and transforms it to represent a void space between two cylinders of radius R1 and R2

std::string mesh_name = "Cylinders-grid";

using namespace INMOST;
const Storage::real pi = 3.1415926535897932384626433832795;

#define V_ID(x, y, z) ((x)*(ny+1)*(nz+1) + (y)*(nz+1) + (z))

int main(int argc, char *argv[]) 
{
  if (argc < 6)
  {
    std::cout << "Usage: " << argv[0] << " grid nx ny nz R1 R2 [shiftx] [shifty] [rotatexy (degree)]" << std::endl;
    return -1;
  }
  double shift[2] = { 0,0 }, rotate = 0;
  std::string fname = std::string(argv[1]);
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[4]);
  double R1 = atof(argv[5]); //bottom left corner cylinder radius
  double R2 = atof(argv[6]); //top right corner cylinder radius
  if (argc > 7) shift[0] = atof(argv[7]);
  if (argc > 8) shift[1] = atof(argv[8]);
  if (argc > 9) rotate = atof(argv[9])/180.0*pi;
  double cr = cos(rotate), sr = sin(rotate);
  Mesh * m = new Mesh();

  if (R1 < 0 || R1 > 1) std::cout << "R1 is " << R1 << " should be 0 <= R1 <= 1 " << std::endl;
  if (R2 < 0 || R2 > 1) std::cout << "R2 is " << R2 << " should be 0 <= R2 <= 1 " << std::endl;
  
  // Create i-j-k structure of nodes
  ElementArray<Node> newverts(m);
  newverts.reserve((nx+1)*(ny+1)*(nz+1));
  double alpha, beta, sign = 1;
  //double Rc1 = 0.5 + (R1 - R2) / std::sqrt(8.0), a1 = (2.0 * Rc1 - 1.0) / (4.0 * Rc1 * (Rc1 - 1)), b1 = -0.5 - 2 * a1;
  //double Rc2 = 0.5 + (R2 - R1) / std::sqrt(8.0), a2 = (2.0 * Rc2 - 1.0) / (4.0 * Rc2 * (Rc2 - 1)), b2 = -0.5 - 2 * a2;
  for (int i = 0; i <= nx; i++)
	  for (int j = 0; j <= ny; j++)
		  for (int k = 0; k <= nz; k++)
		  {
			  Storage::real xyz[3];
			  xyz[0] = i * 1.0 / nx;
			  xyz[1] = j * 1.0 / ny;
			  xyz[2] = k * 1.0 / nz;
			  //projection of nodes
			  if (fabs(R1) && fabs(R2))
			  {
				  if (xyz[1] <= 0.5)
				  {
					  alpha = xyz[0] / 0.5;
					  beta = xyz[1] / 0.5;
					  xyz[0] = (1.0 - beta) * R1 * sin(alpha * pi / 4.0) + beta * (alpha / 2.0 + 0.5 * sign * (R1 - R2) * sin(alpha * pi / 2.0));
					  xyz[1] = (1.0 - beta) * R1 * cos(alpha * pi / 4.0) + beta * (1.0 - alpha / 2.0 + 0.5 * sign * (R1 - R2) * sin(alpha * pi / 2.0));
				  }
				  else if (xyz[1] >= 0.5)
				  {
					  alpha = (1 - xyz[0]) / 0.5;
					  beta = (1 - xyz[1]) / 0.5;
					  xyz[0] = 1.0 - ((1.0 - beta) * R2 * sin(alpha * pi / 4.0) + beta * (alpha / 2.0 - 0.5 * sign * (R1 - R2) * sin(alpha * pi / 2.0)));
					  xyz[1] = 1.0 - ((1.0 - beta) * R2 * cos(alpha * pi / 4.0) + beta * (1.0 - alpha / 2.0 - 0.5 * sign * (R1 - R2) * sin(alpha * pi / 2.0)));
				  }
			  }
			  else if (fabs(R1))
			  {
				  alpha = xyz[1];
				  beta = xyz[0];
				  if (xyz[1] <= 0.5)
				  {
					  xyz[0] = (1.0 - beta) * R1 * cos(alpha * pi / 2.0) + beta * 1.0;
					  xyz[1] = (1.0 - beta) * R1 * sin(alpha * pi / 2.0) + beta * alpha * 2.0;
				  }
				  else
				  {
					  xyz[0] = (1.0 - beta) * R1 * cos(alpha * pi / 2.0) + beta * (2.0 - alpha * 2.0);
					  xyz[1] = (1.0 - beta) * R1 * sin(alpha * pi / 2.0) + beta * 1.0;
				  }
			  }
			  else if (fabs(R2))
			  {
				  alpha = xyz[1];
				  beta = xyz[0];
				  if (xyz[1] <= 0.5)
				  {
					  xyz[0] = 1.0 - ((1.0 - beta) * R2 * cos(alpha * pi / 2.0) + beta * 1.0);
					  xyz[1] = 1.0 - ((1.0 - beta) * R2 * sin(alpha * pi / 2.0) + beta * alpha * 2.0);
				  }
				  else
				  {
					  xyz[0] = 1.0 - ((1.0 - beta) * R2 * cos(alpha * pi / 2.0) + beta * (2.0 - alpha * 2.0));
					  xyz[1] = 1.0 - ((1.0 - beta) * R2 * sin(alpha * pi / 2.0) + beta * 1.0);
				  }
			  }
			  if (rotate)
			  {
				  double x0 = xyz[0] - 0.5, y0 = xyz[1] - 0.5;
				  xyz[0] = x0 * cr + y0 * sr + 0.5;
				  xyz[1] = -x0 * sr + y0 * cr + 0.5;
			  }
			  xyz[0] += shift[0];
			  xyz[1] += shift[1];
			  newverts.push_back(m->CreateNode(xyz)); // Create node in the mesh with index V_ID(i,j,k)
		  }

  const INMOST_DATA_INTEGER_TYPE face_nodes[24] = { 0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6 };
  const INMOST_DATA_INTEGER_TYPE num_nodes[6] = { 4,       4,       4,       4,       4,       4 };

  // Create i-j-k structure of elements
  for (int i = 1; i <= nx; i++)
	  for (int j = 1; j <= ny; j++)
		  for (int k = 1; k <= nz; k++)
		  {
			  // Create local array of eight nodes                            /*      (4)*-------*(6)  */
			  // using representation on the right figure                     /*        /|      /|     */
			  ElementArray<Node> verts(m);                                    /*       /       / |     */
			  verts.push_back(newverts[V_ID(i - 1, j - 1, k - 1)]); // 0      /*      /  |    /  |     */
			  verts.push_back(newverts[V_ID(i - 0, j - 1, k - 1)]); // 1      /*  (5)*-------*(7)|     */
			  verts.push_back(newverts[V_ID(i - 1, j - 0, k - 1)]); // 2      /*     |   |   |   |     */
			  verts.push_back(newverts[V_ID(i - 0, j - 0, k - 1)]); // 3      /*     |       |   |     */
			  verts.push_back(newverts[V_ID(i - 1, j - 1, k - 0)]); // 4      /*     |   |   |   |     */
			  verts.push_back(newverts[V_ID(i - 0, j - 1, k - 0)]); // 5      /*     |(0)*- -|- -*(2)  */
			  verts.push_back(newverts[V_ID(i - 1, j - 0, k - 0)]); // 6      /*     |  /    |  /      */
			  verts.push_back(newverts[V_ID(i - 0, j - 0, k - 0)]); // 7      /*     |       | /       */
																			  /*     |/      |/        */
			  m->CreateCell(verts, face_nodes, num_nodes, 6); // Create the cubic cell in the mesh
		  }

	

  Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("GRIDNAME",DATA_BULK,MESH,NONE));
  std::stringstream str;
  str << mesh_name << "_" << nx << "x" << ny << "x" << nz;
  mesh_name = str.str();
  name.replace(name.begin(),name.end(),mesh_name.begin(),mesh_name.end());

  
  m->Save(fname.c_str());
  
  std::cout << "File " << fname << " written " << std::endl;
  


  delete m;
  return 0;
}
