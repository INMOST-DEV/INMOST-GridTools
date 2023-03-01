#include "inmost.h"
#include "slice_func.h"
using namespace INMOST;


//ellipse (x - c)^T A (x - c) = 1
class Ellipse
{
	double A[6], c[3];
public:
	Ellipse()
	{
		std::fill(A, A + 6, 0.0);
		std::fill(c, c + 3, 0.0);
	}
	Ellipse(double Axx, double Axy, double Axz, double Ayy, double Ayz, double Azz, double cx, double cy, double cz)
	{
		A[0] = Axx;
		A[1] = Axy;
		A[2] = Axz;
		A[3] = Ayy;
		A[4] = Ayz;
		A[5] = Azz;
		c[0] = cx;
		c[1] = cy;
		c[2] = cz;
	}
	Ellipse(double Ain[6], double cin[3])
	{
		std::copy(Ain, Ain + 6, A);
		std::copy(cin, cin + 3, c);
	}
	Ellipse(const Ellipse& b)
	{
		std::copy(b.A, b.A + 6, A);
		std::copy(b.c, b.c + 3, c);
	}
	Ellipse& operator =(Ellipse const& b)
	{
		if (this != &b)
		{
			std::copy(b.A, b.A + 6, A);
			std::copy(b.c, b.c + 3, c);
		}
		return *this;
	}
	double func(double x, double y, double z) const
	{
		double xc[3], Axc[3];
		xc[0] = x - c[0];
		xc[1] = y - c[1];
		xc[2] = z - c[2];
		Axc[0] = A[0] * xc[0] + A[1] * xc[1] + A[2] * xc[2];
		Axc[1] = A[1] * xc[0] + A[3] * xc[1] + A[4] * xc[2];
		Axc[2] = A[2] * xc[0] + A[4] * xc[1] + A[5] * xc[2];
		return 1.0 - (xc[0] * Axc[0] + xc[1] * Axc[1] + xc[2] * Axc[2]);
	}
};

class SliceEllipse : public Slice
{
	std::vector<Ellipse> ellipses;
	double func(double x, double y, double z) const
	{
		double fmax = -1.0e+20;
		for (size_t k = 0; k < ellipses.size(); ++k)
			fmax = std::max(fmax, ellipses[k].func(x, y, z));
		return fmax;
	}
public:
	SliceEllipse() :Slice() {}
	SliceEllipse(const SliceEllipse &b) :Slice(b), ellipses(b.ellipses) {}
	void AddEllipse(const Ellipse& E) { ellipses.push_back(E); }
	SliceEllipse & operator =(SliceEllipse const & b) 
	{ 
		Slice::operator =(b);
		if (this != &b)
			ellipses = b.ellipses;
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
	
	std::string grid_out = "grid.vtk";

	if( argc > 2 ) grid_out = std::string(argv[2]);

	
	

	Mesh m;
	m.Load(argv[1]);
	//m.SetTopologyCheck(NEED_TEST_CLOSURE|PROHIBIT_MULTILINE|PROHIBIT_MULTIPOLYGON|GRID_CONFORMITY|DEGENERATE_EDGE|DEGENERATE_FACE|DEGENERATE_CELL | FACE_EDGES_ORDER | MARK_ON_ERROR | ADJACENT_DUPLICATE);
	//m.RemTopologyCheck(THROW_EXCEPTION);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;

	
	SliceEllipse slice;
	slice.AddEllipse(Ellipse(4.0, 0.1, 0.0, 6.0, 0.1, 24.0, 0.5, 0.5, 0.5));
	slice.AddEllipse(Ellipse(6.0, 1.0, -10.0, 4.0, 1.0, 32.0, 0.7, 0.5, 0.8));
	slice.SliceMesh(m, true);
	
	std::cout << "Cells: " << m.NumberOfCells() << std::endl;
	std::cout << "Faces: " << m.NumberOfFaces() << std::endl;
	    
	m.Save(grid_out);
	return 0;
}
