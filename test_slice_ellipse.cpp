#include "inmost.h"
#include "slice_func.h"
using namespace INMOST;


//ellipse (x - ec)^T eA (x - ec) = 1
double eA[6] =
{
	4.0, 0.1, 0.0,
	     6.0, 0.1,
		      24.0
};
double ec[3] = { 0.5,0.5,0.5 };

class SliceEllipse : public Slice
{
	double A[6], c[3];
	Storage::real func(Storage::real x, Storage::real y, Storage::real z) const
	{
		Storage::real xc[3], Axc[3];
		xc[0] = x - c[0];
		xc[1] = y - c[1];
		xc[2] = z - c[2];
		Axc[0] = A[0] * xc[0] + A[1] * xc[1] + A[2] * xc[2];
		Axc[1] = A[1] * xc[0] + A[3] * xc[1] + A[4] * xc[2];
		Axc[2] = A[2] * xc[0] + A[4] * xc[1] + A[5] * xc[2];
		
		return 1.0 - (xc[0] * Axc[0] + xc[1] * Axc[1] + xc[2] * Axc[2]);
	}
public:
	SliceEllipse(double Ain[6], double cin[3]) :Slice() 
	{
		std::copy(Ain, Ain + 6, A);
		std::copy(cin, cin + 3, c);
	}
	SliceEllipse(const SliceEllipse &b) :Slice(b) 
	{
		std::copy(b.A, b.A + 6, A);
		std::copy(b.c, b.c + 3, c);
	}
	SliceEllipse & operator =(SliceEllipse const & b) 
	{ 
		Slice::operator =(b);
		if (this != &b)
		{
			std::copy(b.A, b.A + 6, A);
			std::copy(b.c, b.c + 3, c);
		}
		return *this;
	}
	Storage::real LevelFunction(Storage::real p[3]) const {return func(p[0],p[1],p[2]);}
};



int main(int argc, char ** argv)
{
	if( argc < 2 )
	{
		std::cout << "Usage: " << argv[0] << " mesh  [mesh_out=grid.pmf]" << std::endl;
		return -1;
	}
	
	std::string grid_out = "grid.pmf";

	if( argc > 2 ) grid_out = std::string(argv[2]);

	
	

	Mesh m;
	m.Load(argv[1]);
	//m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER | MARK_ON_ERROR | ADJACENT_DUPLICATE);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	
	SliceEllipse(eA, ec).SliceMesh(m, true);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	    
	m.Save(grid_out);
	return 0;
}
