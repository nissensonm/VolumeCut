/*		MRC CLASS
 *
 *		-- Dong Si, last updated on Oct 2014
 *
 *
 *		Old Dominion University
 *		Department of Computer Science
 *		Engineering & Computational Sciences Bldg,
 *		4700 Elkhorn Ave, Suite 3318, Norfolk, VA
 *
 */

#ifndef MRC_H
#define MRC_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <math.h>
#include <limits.h>
#include "protein.h"
#include "geometry.h"


using namespace std;

/*
 *		This is an improtant aspect....you should change this according to the mode of the map file you are working on
 *		the default is 2 (float) but you can change it to any type
 */
typedef float vxlDataType;

/*
 *	The Maximum Size of the grid (map) we handle
 */
#define MAXLEN	2000



//data structure for deleting the small voxel groups in map
struct Position
{
    int x;
    int y;
    int z;

    //initializer
    Position() : x(0), y(0), z(0) {}
};
struct Node
{
    Position pos;
    double density;
    bool traveled;
    bool isHlxEnd;

    //initializer
    Node() : density(0.0), traveled(false), isHlxEnd(false){}
};




//local peaks counter in the density map
//for each cell store the coordinate and its local peak counter
struct peakCell{
	double peak;
	Coordinate pnt;
};
// Determine priority.
bool operator<(const peakCell &a, const peakCell &b)
{
	return a.peak  < b.peak ;
}



/*
 *		MRC HEADER SRTUCTURE
 *		source http://www2.mrc-lmb.cam.ac.uk/image2000.html
 */

/*
Map/Image Header Format

Length = 1024 bytes, organized as 56 LONG words followed
by space for 10 80 byte text labels.

1		NX number of columns (fastest changing in map)
2		NY number of rows
3		NZ number of sections (slowest changing in map)
4		MODE data type :	0 image : signed 8-bit bytes range -128 to 127
							1 image : 16-bit halfwords
							2 image : 32-bit reals
							3 transform : complex 16-bit integers
							4 transform : complex 32-bit reals

5		NXSTART number of first column in map (Default = 0)
6		NYSTART number of first row in map
7		NZSTART number of first section in map
8		MX number of intervals along X
9		MY number of intervals along Y
10		MZ number of intervals along Z
11-13	CELLA cell dimensions in angstroms
14-16	CELLB cell angles in degrees
17		MAPC axis corresp to cols (1,2,3 for X,Y,Z)
18		MAPR axis corresp to rows (1,2,3 for X,Y,Z)
19		MAPS axis corresp to sections (1,2,3 for X,Y,Z)
20		DMIN minimum density value
21		DMAX maximum density value
22		DMEAN mean density value
23		ISPG space group number 0 or 1 (default=0)
24		NSYMBT number of bytes used for symmetry data (0 or 80)
25-49	EXTRA extra space used for anything - 0 by default
50-52	ORIGIN origin in X,Y,Z used for transforms
53		MAP character string 'MAP ' to identify file type
54		MACHST machine stamp
55		RMS rms deviation of map from mean density
56		NLABL number of labels being used
57-256	LABEL(20,10) 10 80-character text labels
*/
struct MRC_HEADER {


  int    nx;            /* # of columns ( fastest changing in the map    */
  int    ny;            /* # of rows                                     */
  int    nz;            /* # of sections (slowest changing in the map    */

  int    mode;          /* data type
                              0 = image data in bytes
                              1 = image data in short integer
                              2 = image data in floats
                              3 = complex data in complex short integers
                              4 = complex data in complex reals          */

  int    nxstart;       /* number of first column in map (default = 0)   */
  int    nystart;       /* number of first row in map (default = 0)      */
  int    nzstart;       /* number of first ssection in map (default = 0) */

  int    mx;            /* number of intervals along X                   */
  int    my;            /* number of intervals along Y                   */
  int    mz;            /* number of intervals along Z                   */

  float  xlength;       /* cell dimensions in X (angstrom)               */
  float  ylength;       /* cell dimensions in Y (angstrom)               */
  float  zlength;       /* cell dimensions in Z (angstrom)               */

  float  alpha;         /* cell angles between Y and Z                   */
  float  beta;          /* cell angles between X and Z                   */
  float  gamma;         /* cell angles between X and Y                   */

  int    mapc;          /* number of axis corresponding to columns (X)   */
  int    mapr;          /* number of axis corresponding to rows (Y)      */
  int    maps;          /* number of axis corresponding to sections (Z)  */

  float  amin;          /* minimum density value                         */
  float  amax;          /* maximum density value                         */
  float  amean;         /* mean density value                            */

  int    ispg;          /* space group number (0 for images)             */
  int    nsymbt;        /* # of bytes for symmetry operators             */

  int    extra[25];     /* user defined storage space                    */

  float  xorigin;       /* X phase origin                                */
  float  yorigin;       /* Y phase origin                                */
  float	 zorigin;		/* Z phase origin								 */

  char map[4];			/* constant string "MAP "						*/

  int machineStamp;		/* machine stamp in ccp4 convention: big endian: 0x11110000, little endian 0x44440000 */

  float rms;			/* rms deviation of map from mean density		*/


  int    nlabl;         /* # of labels being used in the MRC header      */

  char   label[10][80]; /* actual text labels                            */

/*
 * NOTE: In some cases what follows the MRC header are symmetry records stored as
 *       text as in the International Tables operators. These are separated by
 *       a * and grouped into 'lines' of 80 characters.
 */

  //initializer
  MRC_HEADER() : nx(0), ny(0), nz(0), mode(2), nxstart(0), nystart(0), nzstart(0), mx(0), my(0), mz(0), xlength(0), \
				ylength(0), zlength(0), alpha(0), beta(0), gamma(0), mapc(0), mapr(0), maps(0), amin(0), amax(0),	\
				amean(0), ispg(0), nsymbt(0), xorigin(0.0), yorigin(0.0), zorigin(0.0), nlabl(0) {}

};

/*
 *		END of MRC HEADER STRUCTURE
 */

struct Gradient
{
  double dx;
  double dy;
  double dz;  //
  double da;  //

  //initializer
  Gradient() : dx(0.0), dy(0.0), dz(0.0), da(0.0) {}
};

struct Tensor
{
  vector<vector<double> > Hmatrix;      //Hessian matrix
  vector<double> Evalue;                //Eigenvalue
  vector<vector<double> > Evector;      //Eigenvector
};

struct Thickness
{
  double t1;
  double t2;
  double t3;

  //initializer
  Thickness() : t1(0.0), t2(0.0), t3(0.0) {}
};

/*
 *		DENSITY MAP : CLASS Definitioan
 */

class Map
{
public:

	MRC_HEADER	hdr;

	/*
	 *		actual data : The Body Of The Map
	 *
	 *		the cube where we store actual values of voxels
	 *		The occurence of the points on the cube
	 *		nx			the fastest point changes	....	but this would be correspondant to x coordinate from pdb
	 *		ny																  correspondant to y coordinate from pdb
	 *		nz			the slice...the slowest point changes				  correspondant to z coordinate from pdb
	 *
	 *					 nx	ny nz	 nx ny nz	 nx ny nz		nx ny nz	nx ny nz	 nx ny nz	 nx ny nz	nx ny nz	ny ny nz
	 *		i.e			[0, 0, 0]	[1, 0, 0]	[2, 0, 0].......[0, 1, 0]	[1, 1, 0]	[2, 1, 0]....[0, 0, 1] [1, 0, 1]	[2, 0, 1]
	 *
	 *		The map is structured in a way the first slice (nz) written first then the second slice....
	 *		the size of the cube will be nx * ny * nz (rows X cols X depth)
	 *		nx : Rows
	 *		ny : Cols
	 *		nz : Depths (slices)
	 *
	 */
    double apixX;								//Angstrom per pixel ratio for X direction
	double apixY;
	double apixZ;
	vector<vector<vector<vxlDataType> > >	cube;   // voxel density value
	vector<vector<vector<Gradient> > >	 grad;      // voxel gradient
	vector<vector<vector<Tensor> > >   tens;        // voxel tenser
	vector<vector<vector<Thickness> > >   thick;    // voxel thickness
	vector<vector<vector<double> > >   dt;           // Distance Transform
	vector<vector<vector<double> > >   dr;           // DT value of the Distance Ridge / Medial Axis

    vector<vector<vector<Node> > >	node;   //vector for all the voxels in map, for filerting small groups


    void setApix();								//set Apix values
	void read(string);							//read the density map ... given the name of the density map.
	void write(string);							//write back the density on a given file
	void createCube(short, short, short);		//create the grid of the size by given dimensions	(rows, cols, slices)
	short numRows();							//returns number of rows in grid3D
	short numCols();							//returns number of cols in grid3D
	short numSlcs();							//returns number of slices in grid3D (depth)




private:
	short slcLen;		//how many cell in each slice
	short sizeOfVxl;	//the size of the voxel,	depends on (mode) if bytes, short, floats,....
};
/*
 *		DENSITY MAP : END of CLASS Definitioan
 */




/*
 *		DENSITY MAP : CLASS Implementation
 */
void Map::read (string mrcFname)
{
	ifstream inMapF;	//map file
	int iRow,
		iCol,
		iDepth;				//regular counters

	//open given mrc file
	inMapF.open (mrcFname.c_str (), ios::binary);

	if (!inMapF.is_open ()){
		cout<<"============================== in MRC::read (string) =========================="<<endl;
		cout<<"Can't open given Map file ( "<<mrcFname<<" ). "<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}
	/*
	 *		Read Map Header
	 */
	inMapF.read ((char*)(&hdr), sizeof(MRC_HEADER));

    if ( hdr.nx <= 0 || hdr.nx >= MAXLEN ||
		hdr.ny <= 0 || hdr.ny >= MAXLEN ||
		hdr.nz <= 0 || hdr.nz >= MAXLEN )
	{
		cout<<"============================== in MRC::read (string) =========================="<<endl;
		cout<<"Strange header of the file. One of (nx,ny,nz) exceeds MAXLEN ( "<<MAXLEN<<" )."<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}

	/*
	 *		Set Apix ratios
	 */
	setApix();

	createCube(hdr.nx, hdr.ny , hdr.nz);		//create the grid (rows , Cols, Depth) --> nx X ny X nz

	/*
	 *		Read map data (Voxels)
	 */

	for (iDepth=0; iDepth<hdr.nz; iDepth++)
		for (iCol=0; iCol<hdr.ny; iCol++)
			for (iRow=0; iRow<hdr.nx; iRow++)
				// Read one cell at a time
				inMapF.read ((char *) &cube[iRow][iCol][iDepth], sizeOfVxl);
}
/////////////////////////////////////////////////////////////////////////////////
void Map::write (string outFileName)
{

	ofstream outMapF;

	//open out file
	outMapF.open (outFileName.c_str (), ios::binary);

	if (!outMapF.is_open ()){
		cout<<"============================= in MRC::write (string) =========================="<<endl;
		cout<<"Can't open given Map file ( "<<outFileName<<" ). "<<endl;
		cout<<"==============================================================================="<<endl;
		exit(1);
	}
	/*
	 *		write the header first
	 */
	outMapF.write ((char *) &hdr, sizeof(MRC_HEADER));

	/*
	 *		write the grid
	 */
	int iRow,
		iCol,
		iDepth;		//regular counters

	for (iDepth = 0; iDepth <numSlcs(); iDepth++)
		for (iCol=0; iCol<numCols(); iCol++)
			for (iRow=0; iRow<numRows(); iRow++)
				//write cell at a time
				outMapF.write ((char *) &cube[iRow][iCol][iDepth], sizeOfVxl);

	outMapF.close();
}

////////////////////////////////////////////////////////////////////////////////////
void Map::setApix ()
{
	apixX = hdr.xlength / (double) hdr.mx;
	apixY = hdr.ylength / (double) hdr.my;
	apixZ = hdr.zlength / (double) hdr.mz;
}
///////////////////////////////////////////////////////////////////////////////////
void Map::createCube(short nX, short nY, short nZ)
{

	// set slice length
	slcLen = nX * nY;

	cube.resize(nX);			//resize rows according to nx

	for(int iRow=0; iRow<nX; iRow++){
		cube[iRow].resize(nY);		// resize cols according to ny ..  generate ny cols for each row

		for (int iCol=0; iCol<nY; iCol++)
			cube[iRow][iCol].resize(nZ);	//resize depth according to nz .. generate nz slices for each cell
	}

	//set the size of voxel (the size of data type used)
	sizeOfVxl = sizeof(cube[0][0][0]);
}

////////////////////////////////////////////////////////////////////////////////////
short Map::numRows ()
{
	return 	cube.size();
}
////////////////////////////////////////////////////////////////////////////////////
short Map::numCols()
{
	if (cube.size ())
		return cube[0].size();
	else
		return 0;
}
////////////////////////////////////////////////////////////////////////////////////
short Map::numSlcs()
{
	if (cube.size ())
		return cube[0][0].size();
	else
		return 0;
}

/*
 *		DENSITY MAP : End of CLASS Implementation
 */
#endif


