#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

std::string mesh_name = "CubeWithChannels";

using namespace INMOST;

const bool DEBUG_VALUES = false;

#define V_ID(x, y, z) ((x-localstart[0])*(localsize[1]+1)*(localsize[2]+1) + (y-localstart[1])*(localsize[2]+1) + (z-localstart[2]))

constexpr double pi() { return std::atan(1)*4.; }

// Create single channel
// n         - is the number of points in the generated curve
// pt_start  - starting point of the curve
// pt_end    - ending point of the curve
// m, ay, az - coefficients in the curve definition: (t, y0+(y1-y0)*t+ay*sin(pi*n*t), z0+(z1-z0)*t+az*sin(pi*n*t))
// curve     - vector of curve points (size=3*n)
void create_curve(int n,
		  const std::vector<double>& pt_start,
		  const std::vector<double>& pt_end,
		  int m, double ay, double az,
		  std::vector<double>& curve)
{
    assert(n > 0);
    assert(pt_start.size() == 3);
    assert(pt_end  .size() == 3);

    curve.resize(3*n);

    if( n == 1 )
    {
	curve = pt_start;
	return;
    }

    double Dx = pt_end[0] - pt_start[0];
    double Dy = pt_end[1] - pt_start[1];
    double Dz = pt_end[2] - pt_start[2];
    
    double h = 1.0 / (n-1.);
    for(int i=0; i < n; ++i)
    {
	double t = i*h; // in [0,1]
	curve[3*i  ] = pt_start[0] + Dx*t;
	curve[3*i+1] = pt_start[1] + Dy*t + ay*std::sin(pi()*m*t);
	curve[3*i+2] = pt_start[2] + Dz*t + az*std::sin(pi()*m*t);
    }
}

// Add a pipe of channels
// npy, npz  - the size of the pipe (npy x npz curves)
// n         - is the number of points in the generated curves
// p_start   - starting point of the pipe
// p_end     - ending point of the pipe
// m, ay, az - coefficients in the curve definition: (t, y0+(y1-y0)*t+ay*sin(pi*n*t), z0+(z1-z0)*t+az*sin(pi*n*t))
// pipe      - vector of curves
void add_pipe(int npy, int npz,
	      int n,
	      const std::vector<double>& p_start,
	      const std::vector<double>& p_end,
	      int m, double ay, double az,
	      std::vector<std::vector<double>>& pipe)
{
    assert(npy > 0);
    assert(npz > 0);
    assert(n > 0);
    assert(p_start.size() == 3);
    assert(p_end  .size() == 3);

    // By the current design the pipe will be originated from points near pt_start and end near pt_end
    // according to the following rule:
    // y0i = pt_start[1] + a*(-0.5+i/(npy-1)), i=0,npy-1
    // z0i = pt_start[2] + a*(-0.5+i/(npz-1)), i=0,npz-1
    double a = 0.1;
    for(int k=0; k < npz; ++k)
    {
	double z0 = p_start[2], z1 = p_end[2];
	if( npz > 1 )
	{
	    z0 += a*(-0.5 + k/(npz-1.));
	    z1 += a*(-0.5 + k/(npz-1.));
	}
	for(int  j=0; j < npy; ++j)
	{
	    double y0 = p_start[1], y1 = p_end[1];
	    if( npy > 1 )
	    {
		y0 += a*(-0.5 + j/(npy-1.));
		y1 += a*(-0.5 + j/(npy-1.));
	    }
	    std::vector<double> pt_start({p_start[0], y0, z0});
	    std::vector<double> pt_end  ({p_end  [0], y1, z1});
	    std::vector<double> curve;
	    create_curve(n, pt_start, pt_end, m, ay, az, curve);
	    pipe.push_back(curve);
	}
    }
}

// Create channels
// n         - is the number of points in the generated curves
// alpha     - factor for wibling
// channels  - vector of curves
void create_channels(int n,
		     double alpha,
		     std::vector<std::vector<double>>& channels)
{
    assert(n > 0);

    //add_pipe(5, 5, n, {0,0.25,0.25}, {1.,0.25,0.15}, 5, 0.10*alpha, 0.05*alpha, channels);
    //add_pipe(5, 5, n, {0,0.75,0.25}, {1.,0.80,0.15}, 4, 0.05*alpha, 0.10*alpha, channels);
    //add_pipe(5, 5, n, {0,0.25,0.75}, {1.,0.25,0.80}, 5, 0.10*alpha, 0.05*alpha, channels);
    //add_pipe(5, 5, n, {0,0.75,0.75}, {1.,0.80,0.85}, 6, 0.05*alpha, 0.10*alpha, channels);

    //add_pipe(5, 5, n, {0,0.25,0.25}, {1.,0.25,0.80}, 5, 0.10*alpha, 0.05*alpha, channels);
    add_pipe(20, 20, n, {0.05,0.25,0.25}, {0.95,0.25,0.25}, 5, 0.10*alpha, 0.05*alpha, channels);
    add_pipe(20, 20, n, {0.05,0.75,0.25}, {0.95,0.75,0.25}, 4, 0.05*alpha, 0.10*alpha, channels);
    add_pipe(20, 20, n, {0.05,0.25,0.75}, {0.95,0.25,0.75}, 5, 0.10*alpha, 0.05*alpha, channels);
    add_pipe(20, 20, n, {0.05,0.75,0.75}, {0.95,0.75,0.75}, 6, 0.05*alpha, 0.10*alpha, channels);
    // add_pipe(10, 10, n, {0.0,0.25,0.25}, {1.,0.25,0.25}, 5, 0.10*alpha, 0.05*alpha, channels);
    // add_pipe(10, 10, n, {0.0,0.75,0.25}, {1.,0.75,0.25}, 4, 0.05*alpha, 0.10*alpha, channels);
    // add_pipe(10, 10, n, {0.0,0.25,0.75}, {1.,0.25,0.75}, 5, 0.10*alpha, 0.05*alpha, channels);
    // add_pipe(10, 10, n, {0.0,0.75,0.75}, {1.,0.75,0.75}, 6, 0.05*alpha, 0.10*alpha, channels);
}


// Generate cubic grid in parallel
Mesh *ParallelGenerator(INMOST_MPI_Comm comm, int nx, int ny, int nz)
{
    double hx = 1.0 / static_cast<double>(nx);
    double hy = 1.0 / static_cast<double>(ny);
    double hz = 1.0 / static_cast<double>(nz);
    
    int procs_per_axis[3] = {1, 1, 1};
    int sizes[3] = {nx, ny, nz};
    Mesh *m = new Mesh(); // Create a mesh to be constructed

    m->SetCommunicator(comm); // Set the MPI communicator, usually MPI_COMM_WORLD

#if defined(USE_MPI)
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
#endif

    int rank = m->GetProcessorRank();    // Get the rank of the current process
    int size = m->GetProcessorsNumber(); // Get the number of processors used in communicator comm

    // Compute the configuration of processes connection
    {
	int divsize = size;
	std::vector<int> divs;
	while( divsize > 1 )
	{
	    for(int k = 2; k <= divsize; ++k)
		if( divsize % k == 0 )
		{
		    divs.push_back(k);
		    divsize /= k;
		    break;
		}
	}

	int elements_per_procs[3] = {sizes[0], sizes[1], sizes[2]};
	for(std::vector<int>::reverse_iterator it = divs.rbegin(); it != divs.rend(); ++it)
	{
	    int *max = std::max_element(elements_per_procs+0, elements_per_procs+3);
	    procs_per_axis[max-elements_per_procs] *= *it;
	    (*max) /= *it;
	}
    }

    int proc_coords[3] = {rank % procs_per_axis[0],
			  rank / procs_per_axis[0] % procs_per_axis[1],
			  rank / (procs_per_axis[0]*procs_per_axis[1])};

    int localsize[3], localstart[3], localend[3];
    int avgsize[3] = {
	(int)ceil((double)sizes[0] / procs_per_axis[0]),
	(int)ceil((double)sizes[1] / procs_per_axis[1]),
	(int)ceil((double)sizes[2] / procs_per_axis[2])
    };

    for(int j = 0; j < 3; ++j)
    {
	localstart[j] = avgsize[j] * proc_coords[j];
	if( proc_coords[j] == procs_per_axis[j] - 1 )
	    localsize[j] = sizes[j] - avgsize[j] * (procs_per_axis[j]-1);
	else
	    localsize[j] = avgsize[j];
	localend[j] = localstart[j] + localsize[j];
    }

    // Create i-j-k structure of nodes
    ElementArray<Node> newverts(m);
    newverts.reserve(localsize[0]*localsize[1]*localsize[2]);

    if( DEBUG_VALUES)
    {
	if( size == 1 )
	{
	    std::cout << "ixs="<<localstart[0]<<" ixe="<<localend[0]<<std::endl;
	    std::cout << "jys="<<localstart[1]<<" jye="<<localend[1]<<std::endl;
	    std::cout << "kzs="<<localstart[2]<<" kze="<<localend[2]<<std::endl;
	} else
	{
	    std::cout << "Proc " << rank << ": ixs="<<localstart[0]<<" ixe="<<localend[0]<<std::endl;
	    std::cout << "Proc " << rank << ": jys="<<localstart[1]<<" jye="<<localend[1]<<std::endl;
	    std::cout << "Proc " << rank << ": kzs="<<localstart[2]<<" kze="<<localend[2]<<std::endl;
	}
    }
    
    for(int i = localstart[0]; i <= localend[0]; ++i)
	for(int j = localstart[1]; j <= localend[1]; ++j)
	    for(int k = localstart[2]; k <= localend[2]; ++k)
	    {
		Storage::real xyz[3];
		xyz[0] = i * hx;
		xyz[1] = j * hy;
		xyz[2] = k * hz;
		newverts.push_back(m->CreateNode(xyz)); // Create node in the mesh with index V_ID(i,j,k)
	    }

    // Create i-j-k structure of elements
    // Define six cube faces assuming verts are numerated in the way presented above
    const INMOST_DATA_INTEGER_TYPE face_nodes[24] = {0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6};
    const INMOST_DATA_INTEGER_TYPE num_nodes[6]   = {4,       4,       4,       4,       4,       4};
    for(int i = localstart[0]+1; i <= localend[0]; ++i)
	for(int j = localstart[1]+1; j <= localend[1]; ++j)
	    for(int k = localstart[2]+1; k <= localend[2]; ++k)
	    {
		// Create local array of eight nodes                      /*      (4)*-------*(6)  */
		// using representation on the right figure               /*        /|      /|     */
		ElementArray<Node> verts(m);                              /*       /       / |     */
		verts.push_back(newverts[V_ID(i-1, j-1, k-1)]); // 0      /*      /  |    /  |     */
		verts.push_back(newverts[V_ID(i-0, j-1, k-1)]); // 1      /*  (5)*-------*(7)|     */
		verts.push_back(newverts[V_ID(i-1, j-0, k-1)]); // 2      /*     |   |   |   |     */
		verts.push_back(newverts[V_ID(i-0, j-0, k-1)]); // 3      /*     |       |   |     */
		verts.push_back(newverts[V_ID(i-1, j-1, k-0)]); // 4      /*     |   |   |   |     */
		verts.push_back(newverts[V_ID(i-0, j-1, k-0)]); // 5      /*     |(0)*- -|- -*(2)  */
		verts.push_back(newverts[V_ID(i-1, j-0, k-0)]); // 6      /*     |  /    |  /      */
		verts.push_back(newverts[V_ID(i-0, j-0, k-0)]); // 7      /*     |       | /       */
                                                                          /*     |/      |/        */
		// Create the cubic cell in the mesh                      /*  (1)*-------*(3)      */
		m->CreateCell(verts, face_nodes, num_nodes, 6);
	    }

    m->ResolveShared(); // Resolve duplicate nodes
    return m;
}

void define_perms(Mesh *mesh,
		  const std::vector<std::vector<double>>& channels,
		  double k_base, double k_channel)
{
    Tag mat = mesh->CreateTag("MATERIAL", DATA_INTEGER, CELL, NONE, 1);
    Storage::real cnt[3];
    Tag perm = mesh->CreateTag("K", DATA_REAL, CELL, NONE, 1);
    for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
    {
	ElementArray<Node> nodes = it->getNodes();
	// Find the bounding box for the cell
	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = xmax = nodes[0].Coords()[0];
	ymin = ymax = nodes[0].Coords()[1];
	zmin = zmax = nodes[0].Coords()[2];
	for(int i=1; i < 6; ++i) // Since we have a parallelepiped 5 point should be enough
	{
	    xmin = std::min(xmin, nodes[i].Coords()[0]);
	    xmax = std::max(xmax, nodes[i].Coords()[0]);
	    ymin = std::min(ymin, nodes[i].Coords()[1]);
	    ymax = std::max(ymax, nodes[i].Coords()[1]);
	    zmin = std::min(zmin, nodes[i].Coords()[2]);
	    zmax = std::max(zmax, nodes[i].Coords()[2]);
	}
	if( DEBUG_VALUES )
	    std::cout << "Cell bounding box: ["<<xmin<<","<<xmax<<"]x["<<ymin<<","<<ymax<<"]x["<<zmin<<","<<zmax<<"]\n";

	// Go over all polylines and check if any vertex is contained in the box
	bool in_channel = false;
	for(int ic=0; ic < channels.size(); ++ic)
	{
	    if( in_channel ) break;
	    // Use the knowledge that x-coordinates of each curve are monotonically increasing
	    double hc = channels[ic][3] - channels[ic][0];
	    int i_start = std::max(0, static_cast<int>(std::floor(xmin/hc)-1));
	    for(int i=0; i < channels[ic].size()/3; ++i)
	    {
		double x = channels[ic][3*i+0];
		if( x < xmin ) continue;
		if( x > xmax ) break;
		double y = channels[ic][3*i+1];
		if( y < ymin || y > ymax ) continue;
		double z = channels[ic][3*i+2];
		if( z < zmin || z > zmax ) continue;
		in_channel = true;
		break;
	    }
	}

	if( in_channel )
	{
	    it->RealDF(perm) = k_channel;
	    it->IntegerDF(mat) = 1;
	} else {
	    it->RealDF(perm) = k_base;
	    it->IntegerDF(mat) = 0;
	}
    }
}


int main(int argc, char *argv[]) 
{
    if (argc < 2)
    {
	std::cout << "Usage: " << argv[0] << " nx|ny [nz = 1] [k_base = 0.001] [k_channel=10.] [alpha=0.0]" << std::endl;
	return -1;
    }

    Mesh::Initialize(&argc, &argv);

    // Get command line options
    double alpha = 0.0;       // Amplitude of vibration for channels
    double k_channel = 10.;   // Permeability in channels
    double k_base = 0.001;    // Permeability in a matrix
    int nz = 1;               // Default number of layers in z direction
    int nx = 8;               // Default number of layers in x- and y- directions

    if (argc > 1 )  nx = atoi(argv[1]);
    if( argc > 2 )  nz = atoi(argv[2]);

    if (argc > 3 )  k_base    = atof(argv[3]);
    if (argc > 4 )  k_channel = atof(argv[4]);
    if (argc > 5 )  alpha     = atof(argv[5]);


    // Construct a mesh in parallel
    double t_loc = Timer();
    Mesh *mesh = ParallelGenerator(INMOST_MPI_COMM_WORLD, nx, nx, nz);
    t_loc = Timer() - t_loc;
    double t_tot = mesh->Integrate(t_loc);
    t_loc = t_tot / static_cast<double>(mesh->GetProcessorsNumber());

    int nodes_tot = mesh->TotalNumberOf(NODE);
    int edges_tot = mesh->TotalNumberOf(EDGE);
    int faces_tot = mesh->TotalNumberOf(FACE);
    int cells_tot = mesh->TotalNumberOf(CELL);
    int rank = mesh->GetProcessorRank();    // Get the rank of the current process
    if( rank == 0 )
    {
	std::cout << "Grid: " << nx << " x " << nx << " x " << nz << std::endl;
	std::cout << "Processors: " << mesh->GetProcessorsNumber() << std::endl;
	std::cout << "Mesh generator time: total = " << t_tot << "  average = " << t_loc << std::endl;
	std::cout << "Nodes: " << nodes_tot << std::endl
		  << "Edges: " << edges_tot << std::endl 
		  << "Faces: " << faces_tot << std::endl 
		  << "Cells: " << cells_tot << std::endl;
    }

#if defined(USE_MPI)
    MPI_Barrier(mesh->GetCommunicator());
#endif

    // Define curves
    std::vector<std::vector<double>> channels;

    create_channels(std::max(100,2*nx), alpha, channels);

    // Save channels
    if( rank == 0 )
    {
	std::ofstream of("channels.txt");
	of << "# of channels: " << channels.size() << std::endl;
	for(int i=0; i < channels.size(); ++i)
	{
	    std::vector<double>& ch = channels[i];
	    assert( ch.size()%3 == 0);
	    for(int k=0; k < ch.size()/3; ++k)
		of << ch[3*k] << "  " << ch[3*k+1] << "  " << ch[3*k+2] << std::endl;
	    of << std::endl << std::endl;
	}
	of.close();
    }
  
#if defined(USE_MPI)
    MPI_Barrier(mesh->GetCommunicator());
#endif

    t_loc = Timer();
    if( rank == 0 )
	std::cout << "Define permeability" << std::endl;
    define_perms(mesh, channels, k_base, k_channel);

//===============================
  
#if defined(USE_MPI)
    mesh->ExchangeGhost(1, NODE); // Construct Ghost cells in 1 layer connected via nodes
    MPI_Barrier(mesh->GetCommunicator());
#endif

    if( mesh->GetProcessorRank() == 0 )
	std::cout << "I'm ready!\n";

    std::string filename("grid");
    t_loc = Timer();
    if( mesh->GetProcessorsNumber() == 1 )
    {
        mesh->Save(filename + ".vtk");
	mesh->Save(filename + ".pmf");
	mesh->Save(filename + ".gmv");
    }
    else
    {
        mesh->Save(filename + ".pvtk");
    }

#if defined(USE_MPI)
    MPI_Barrier(mesh->GetCommunicator());
#endif

    t_loc = Timer() - t_loc;
    if( mesh->GetProcessorRank() == 0 )
	std::cout << "Saving to file \"" << filename << "\", time: " << t_loc << std::endl;

    delete mesh;
    Mesh::Finalize();
    return 0;
}
