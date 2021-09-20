#include "inmost.h"


using namespace INMOST;

const double pi = 3.1415926535897932384626433832795;


int main(int argc, char *argv[]) 
{
	
	if( argc < 2 )
	{
		printf("Usage: %s meshin [meshout=out.pmf]\n",argv[0]);
		return -1;
	}
	std::string infile(argv[1]);
	std::string outfile = "out.pmf";
	if( argc > 2 ) outfile = std::string(argv[2]);  
	
	Mesh mesh; 
	mesh.Load(infile);
	
	for (Mesh::iteratorNode it = mesh.BeginNode(); it != mesh.EndNode(); ++it)
	{
		Storage::real_array c = it->Coords();
		std::swap(c[1], c[2]);
	}
	  
	std::cout << "Save file " << outfile << std::endl;
	mesh.Save(outfile);
	return 0;
}
