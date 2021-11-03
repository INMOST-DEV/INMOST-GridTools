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
  if (argc < 5)
  {
    std::cout << "Usage: " << argv[0] << " grid nx ny nz [R=1] [L=10] [rotatexy=0(degree)] [rotatexz=0(degree)" << std::endl;
    return -1;
  }
  double L = 10, rotatexy = 0, rotatexz = 0;
  std::string fname = std::string(argv[1]);
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[4]);
  double R = 1;
  if (argc > 5) R = atof(argv[5]); //bottom left corner cylinder radius
  if (argc > 6) L = atof(argv[6]);
  if (argc > 7) rotatexy = atof(argv[7]) / 180.0 * pi;
  if (argc > 8) rotatexz = atof(argv[8]) / 180.0 * pi;
  double crxy = cos(rotatexy), srxy = sin(rotatexy);
  double crxz = cos(rotatexz), srxz = sin(rotatexz);
  Mesh * m = new Mesh();

  if (R < 0) std::cout << "R is " << R << std::endl;
  // Create i-j-k structure of nodes
  ElementArray<Node> newverts(m);
  newverts.reserve((nx+1)*(ny+1)*(nz+1));
  double alpha, beta, d, xr, yr, zr, dr;
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
			  alpha = k * 1.0 / nz, beta = j * 1.0 / ny, d = i * 1.0 / nx;
			  xr = -R / sqrt(3.0);
			  yr = R * (2.0 * beta - 1.0) / sqrt(3.0);
			  zr = R * (2.0 * alpha - 1.0) / sqrt(3.0);
			  dr = sqrt(xr * xr + yr * yr + zr * zr);
			  xr = R * xr / dr;
			  yr = R * yr / dr;
			  zr = R * zr / dr;
			  xyz[0] = (1 - d) * (-L) + d * xr;
			  xyz[1] = (1 - d) * (-L + 2 * L * beta) + d * yr;
			  xyz[2] = (1 - d) * (-L + 2 * L * alpha) + d * zr;
			  
			  //xyz[0] = (1 - d) * (-L) + d * (-R * cos(pi / 2.0 * (0.5 - alpha)) );
			  //xyz[1] = (1 - d) * (-L + 2 * L * beta) + d * (-R * sin(pi / 2.0 * (0.5 - alpha)) * sin(pi / 2.0 * (0.5 - beta)));
			  //xyz[2] = (1 - d) * (-L + 2 * L * alpha) + d * (-R * sin(pi / 2.0 * (0.5 - alpha)) * cos(pi / 2.0 * (0.5 - beta)));

			  if (rotatexy)
			  {
				  double x0 = xyz[0], y0 = xyz[1];
				  xyz[0] =  x0 * crxy + y0 * srxy;
				  xyz[1] = -x0 * srxy + y0 * crxy;
			  }
			  if (rotatexz)
			  {
				  double x0 = xyz[0], z0 = xyz[2];
				  xyz[0] =  x0 * crxz + z0 * srxz;
				  xyz[2] = -x0 * srxz + z0 * crxz;
			  }
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
