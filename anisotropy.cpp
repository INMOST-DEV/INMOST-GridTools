#include "inmost.h"

using namespace INMOST;

//todo: want to separate all faces into those that share edges
//see todo: inside text


int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] << " mesh [tag=PERM]" << std::endl;
		return -1;
	}
	
	Mesh::Initialize(&argc,&argv);
	
	std::string tag_name = "PERM";
	if( argc > 2) tag_name = std::string(argv[2]);
	
	
	std::cout << "Load: " << argv[1] << std::endl;
	Mesh m;
	m.SetFileOption("VERBOSITY","2");
	m.Load(argv[1]);
	
	if( !m.HaveTag(tag_name) )
	{
		std::cout << "mesh does not have tag named " << tag_name << std::endl;
		exit(-1);
	}
	
	TagRealArray tag_K = m.GetTag(tag_name);
	TagReal tag_smax = m.CreateTag("smax", DATA_REAL, CELL, NONE, 1);
	TagReal tag_smin = m.CreateTag("smin", DATA_REAL, CELL, NONE, 1);

	std::cout << "Start:" << std::endl;
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	std::cout << "Edges: " << m.NumberOfEdges() << std::endl;

	
	double ani = 1;
#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
		rMatrix U(3, 3), S(3, 3), V(3, 3), K(3, 3);
		double smax, smin, ani_loc = 1;
		Cell n;
#if defined(USE_OMP)
#pragma omp for
#endif
		for (int q = 0; q < m.CellLastLocalID(); ++q) if (m.isValidCell(q))
		{
			n = m.CellByLocalID(q);
			K = rMatrix::FromTensor(tag_K[n].data(), tag_K[n].size(), 3);
			K.SVD(U, S, V);
			smax = smin = S(0, 0);
			for (int j = 0; j < 3; ++j)
			{
				smax = std::max(smax, S(j, j));
				smin = std::min(smin, S(j, j));
			}
			tag_smax[n] = smax;
			tag_smin[n] = smin;

			ani_loc = std::max(ani_loc, smax / smin);
		}
#if defined(USE_OMP)
#pragma omp critical
#endif
		ani = std::max(ani, ani_loc);
	}
	std::cout << "Local cell anisotropy: " << ani << std::endl;
	ani = 1;

#if defined(USE_OMP)
#pragma omp parallel
#endif
	{
		double ani_loc = 1, smax, smin;
		Cell n;
#if defined(USE_OMP)
#pragma omp for
#endif
		for (int q = 0; q < m.CellLastLocalID(); ++q) if (m.isValidCell(q))
		{
			n = m.CellByLocalID(q);
			smax = tag_smax[n];
			smin = tag_smin[n];
			ElementArray<Cell> adj = n.NeighbouringCells();
			for (ElementArray<Cell>::iterator jt = adj.begin(); jt != adj.end(); ++jt)
			{
				smax = std::max(smax, tag_smax[*jt]);
				smin = std::min(smin, tag_smin[*jt]);
			}
			ani_loc = std::max(ani_loc, std::max(smax / tag_smin[n], tag_smax[n] / smin));
		}
#if defined(USE_OMP)
#pragma omp critical
#endif
		ani = std::max(ani, ani_loc);
	}

	std::cout << "Adjacent cell anisotropy: " << ani << std::endl;
	
	Mesh::Finalize();
	return 0;
}
