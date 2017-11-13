#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>

#include "include/protein.h"
#include "include/mrc.h"

using namespace std;

// Generates coordinates along a backbone chain. See function for details. 
void cutVolumeAroundChain(Protein pdb, Map mrcF, Map & outMap, vector<Coordinate> axis, float distance, float apixX, float apixY, float apixZ);

// Display chains to users. Also does some error checking if MRC or PDB file are not formatted correctly. 
string getChainToCut(string pdbString);

// Generate the trace to define spheres and cut volume from.
// Axis is the coordinates backbone is stored into. 
void generateTraceToCutAlong(vector<Coordinate> &axis, Protein &pdb, int option);

// Arguments passed in should be PDB file, MRC file, Radius
int main(int argc, char *argv[])
{	
	if (argc != 5){
		cout << "Must invoke Volume Cut with the following parameters:\n";
		cout << "PDBFileName MRCFileName Radius MRCResolution\n";
		return -1;
	}
	
	// Declaring variables at the outset. 
	string pdbStr, mrcStr, pdbOUTStr, mrcOUTStr, chainID;
	Protein pdb;		
	Map inMRC, outMRC;
    float apixX=1, apixY=1, apixZ=1;
	float radius = 1;
	float resolution;
	
	// Axis holds the coordinates along the backbone chain that volume is cut around. 
	vector<Coordinate> axis;
	
	// Error handling for PDB and MRC parameters done when the file is read in later.
	pdbStr = argv[1];
	mrcStr = argv[2];
	
	// Check if radius is a valid number
	istringstream ss(argv[3]);
	if (!(ss >> radius) || radius < 0)
	{
		cout << argv[3] << " is not a valid number. \n";
		cout << "Must enter a valid, positive number for the radius. Setting radius to the default value of 1. \n";
		radius = 1;
	}
	
	// Check if resolution is a valid number
	istringstream ssRes(argv[4]);
	if (!(ssRes >> resolution))
	{
		cout << argv[4] << " is not a valid resolution. \n";
		cout << "Must enter a valid resolution value. Setting method to cut to high resolution. \n";
		resolution = 1;
	}
	
	// Read in the PDB file.
	cout<<"Reading in PDB File to get Chain Information."<<endl;
	chainID = getChainToCut(pdbStr);
	if (chainID == "BADSELECTION")
	{
		cout << "Something went wrong when trying to select a chain to cut\n or the PDB file is not valid or missing." << endl;
		return -1;
	}
	
	// Read in the PDB file with the matching chain information retrieved from getChainToCut().
	cout<<"Reading in PDB File w/ Selected Chain."<<endl;
	pdb.read(pdbStr + ".pdb", chainID);
	
	// Read in the MRC file.
	cout<<"Reading in MRC File."<<endl;
	inMRC.read(mrcStr + ".mrc");

	//calculate apix
	apixX = inMRC.hdr.xlength / inMRC.hdr.mx;
	apixY = inMRC.hdr.ylength / inMRC.hdr.my;
	apixZ = inMRC.hdr.zlength / inMRC.hdr.mz;
	
	// Write pdb cut chain
	pdbOUTStr = pdbStr + "_" + chainID + "_cut.pdb";
	mrcOUTStr = mrcStr + "_" + chainID + "_cut.mrc";
	cout<<"Writing PDB Cut file as " << pdbOUTStr << endl;
	pdb.writePDB(pdbOUTStr, 1, pdb.numOfAA());
	
	// Create MRC that will be written.
	outMRC.hdr = inMRC.hdr;
	outMRC.createCube(inMRC.numRows(), inMRC.numCols(), inMRC.numSlcs());

	// Select the points along the backbone to cut. 
	// If the resolution is considered medium, 
	//	we cut down the middle of the helicies.
	// Otherwise, trace point to point for helicies.
	generateTraceToCutAlong(axis, pdb, resolution);
		
	// Write all points along the axis. 
	cutVolumeAroundChain(pdb, inMRC, outMRC, axis, radius, apixX, apixY, apixZ);
	
	cout<< "Writing MRC Cut file as " << mrcOUTStr << endl;
    outMRC.write(mrcOUTStr);
	
	return 0;
}

// Generate the backbone trace. 
// Note that for volume surrounding helices,
// medium resolution data is cut differently than high resolution data.
// High ( < 5 A ) -> Cut Volume point to point.
// Medium ( > 5 A ) -> Cut Volume from center of a helix.
void generateTraceToCutAlong(vector<Coordinate> &axis, Protein &pdb, int resolution)
{
	// If high resolution ( resolution is less than 5 A ) 
	if (resolution < 5)
	{
		for (int iAA=0; iAA<pdb.numOfAA() - 1; iAA++)
		{
			// Walk along an interpolated line from one alpha carbon to the next.
			for (float distance = 0; 
					distance <= getDistance(pdb.AAs[iAA + 1].atoms[pdb.getAtomIndx(iAA + 1, " CA ")].coord, 
						pdb.AAs[iAA].atoms[pdb.getAtomIndx(iAA, " CA ")].coord); 
					distance += 0.1)
			{
				// Set the next point along the interpolated line.
				Coordinate pointOnLineBtwACarbons = pointOnLine(pdb.AAs[iAA + 1].atoms[pdb.getAtomIndx(iAA + 1, " CA ")].coord, 
					pdb.AAs[iAA].atoms[pdb.getAtomIndx(iAA, " CA ")].coord, distance);
				
				axis.push_back(pointOnLineBtwACarbons);
			}
		}
	}
	else // If medium resolution ( resolution is greater than 5 A )
	{
		// First store all points along the center of helicies
		for (int i=0; i<pdb.hlces.size(); i++){
			//set the axis segments along the center of the helicies.
			pdb.setAxisSegments(pdb.hlces[i].startIndx, pdb.hlces[i].endIndx, 0.1, axis);
		}
		
		// Find all the indicies for loops and sheets.
		pair<short, short> loopIndces;
		// Loops (startIndx, endIndx)
		vector<pair<short, short> > loops;		
		// Find all indicies that are NOT part of a helix and store them
		if (pdb.hlces.size()){
			int i = 0;
			// Find + store starting loops indicies
			loopIndces.first = 0;
			loopIndces.second = pdb.hlces [0].startIndx -1;
			if (loopIndces.second >= loopIndces.first)
				loops.push_back (loopIndces);
			// Check + store middle loops
			for (i=0; i<pdb.hlces.size ()-1; i++){
				loopIndces.first = pdb.hlces [i].endIndx +1;
				loopIndces.second = pdb.hlces [i+1].startIndx -1;

				if (loopIndces.second >= loopIndces.first)
					loops.push_back (loopIndces);
			}
			// Check + store last loop
			loopIndces.first = pdb.hlces [i].endIndx+1;
			loopIndces.second = pdb.numOfAA()-1;

			if (loopIndces.second >= loopIndces.first)
				loops.push_back (loopIndces);
		}
		else {
			// Treat all AAs as one loop if there is no SS
			loopIndces.first = 0;                       
			loopIndces.second = pdb.numOfAA()-1;
			loops.push_back(loopIndces);
		}

		// With all indicies for helicies, loops, and sheets found, interpolate lines between 
		//	each coordinate.
		for (int i=0; i<loops.size(); i++){
			for (int j = loops[i].first; j < loops[i].second; j++)
			{
				// Walk along an interpolated line from one alpha carbon to the next.
				for (float distance = 0; 
						distance <= getDistance(pdb.AAs[j + 1].atoms[pdb.getAtomIndx(j + 1, " CA ")].coord, 
							pdb.AAs[j].atoms[pdb.getAtomIndx(j, " CA ")].coord); 
						distance += 0.1)
				{
					// Set the next point along the interpolated line.
					Coordinate pointOnLineBtwACarbons = pointOnLine(pdb.AAs[j + 1].atoms[pdb.getAtomIndx(j + 1, " CA ")].coord, 
						pdb.AAs[j].atoms[pdb.getAtomIndx(j, " CA ")].coord, distance);
					
					axis.push_back(pointOnLineBtwACarbons);
				}
			}
		}
	}
}

// Cuts volume around the points along the axis that was passed in.
// Does so by defining a sphere around each coordinate on the axis, then 
//	only copies voxels inside each sphere to the output MRC file 
void cutVolumeAroundChain(Protein pdb, Map mrcF, Map & outMap, vector<Coordinate> axis, float radius, float apixX, float apixY, float apixZ)
{
	Coordinate p;
	
	// Set up a flag for each possible location in the MRC file. 
	// Used so we only need to check distances if the location has already been written to.
	vector<vector<vector<char> > > mrcOutWroteFlag(mrcF.numRows(), vector<vector<char> >(mrcF.numCols(), vector<char>(mrcF.numSlcs())));

	// Used to hold min and max possible coordinates in x, y, and z
	// that we need to check. 
	float minCoordinate[3];
	float maxCoordinate[3];
		
	for (int j = 0; j < axis.size() - 1; j++)
	{
		// Pick the upper and lower bounds in the MRC file to examine for x, y, and z coordinates.
		// This allows us to only have to examine a very narrow box of points surrounding the individual point,
		//	rather than comparing distances across the entire MRC file.
		// This effectively inscribes the sphere inside a sphere so we only need to check the minimum number of voxels. 
		minCoordinate[0] = axis[j].x - radius;
		maxCoordinate[0] = axis[j].x + radius;
		
		minCoordinate[1] = axis[j].y - radius;
		maxCoordinate[1] = axis[j].y + radius;
		
		minCoordinate[2] = axis[j].z - radius;
		maxCoordinate[2] = axis[j].z + radius;
		
		// Error handling; makes sure we don't check values that are outside the bounds of the MRC File 
		// Never check below the lower bounds of the MRC file.
		if ((minCoordinate[0] - mrcF.hdr.xorigin ) / apixX < 0)
			minCoordinate[0] = mrcF.hdr.xorigin;
		if ((minCoordinate[1] - mrcF.hdr.yorigin) / apixY < 0)
			minCoordinate[1] = mrcF.hdr.xorigin;
		if ((minCoordinate[2] - mrcF.hdr.zorigin) / apixZ < 0)
			minCoordinate[2] = mrcF.hdr.xorigin;
		
		// Never check past the upper bounds of the MRC file.
		if ((maxCoordinate[0] - mrcF.hdr.xorigin) / apixX > mrcF.numRows())
			maxCoordinate[0] = mrcF.numRows() * apixZ + mrcF.hdr.xorigin;
		if ((maxCoordinate[1] - mrcF.hdr.yorigin) / apixY > mrcF.numCols())
			maxCoordinate[1] = mrcF.numCols() * apixY + mrcF.hdr.yorigin;
		if ((maxCoordinate[2] - mrcF.hdr.zorigin) / apixZ > mrcF.numSlcs())
			maxCoordinate[2] = mrcF.numSlcs() * apixZ + mrcF.hdr.zorigin;
		
		for (float x=minCoordinate[0]; x<maxCoordinate[0]; x++){
			p.x = x;
			for (float y=minCoordinate[1]; y<maxCoordinate[1]; y++){
				p.y = y;
				for (float z=minCoordinate[2]; z<maxCoordinate[2]; z++){
					p.z = z;
				
					// Doing the calculations here saves us having to re-calculate it 4+ times below.
					// Using the formula: coordinate = origin + index * spacing
					// We have the coordinate, origin, and spacing, we want to solve for the index:
					int xIndx = round((x - mrcF.hdr.xorigin) / apixX);
					int yIndx = round((y - mrcF.hdr.yorigin) / apixY);
					int zIndx = round((z - mrcF.hdr.zorigin) / apixZ);
					
					//First check if the location in the MRC file was written already (never check distance again).
					//	If it was NOT written to, then check if it's within range.
					if (mrcOutWroteFlag[xIndx][yIndx][zIndx] == 0 
						&& getDistance(p, axis[j])<= radius)
					{		
						//store the points included in calculations
						outMap.cube[xIndx][yIndx][zIndx] 
							= mrcF.cube[xIndx][yIndx][zIndx];	
							
						// Set the written flag so this point is never checked again.
						mrcOutWroteFlag[xIndx][yIndx][zIndx] = 1;
					}
				}
			}
		}
	}	
}

// Finds all chains that exist in a PDB file and finds if the chain is unique or a duplicate.
// Also finds the size of each chain. 
string getChainToCut(string pdbStr)
{
	string chainID;
	// First item of each vector is an original chain, subsequent items are duplicates
	vector<vector<string> > chains; 
	
	// Used to store the # of atoms and AA. 
	vector<int> sizeOfChainsAtoms;
	vector<int> sizeOfChainsAA;
	
	Protein pdb;
	pdb.read(pdbStr + ".pdb");
	for (int i = 0; i < pdb.header.size(); i++)
	{
		// Find all parts of the header that have chains listed. 
		int indexStartOfChain = pdb.header[i].find("CHAIN:"); 		
		if (indexStartOfChain != string::npos)
		{
			// Move index forward to first chain.
			indexStartOfChain += 7;
			string firstChain = pdb.header[i].substr(indexStartOfChain, 1);
			
			// Create a temp vector, store in in chains, then add the non-duplicate 
			//	chain to it first.
			vector<string> temp;
			chains.push_back(temp);
			// Added to size() - 1, because that was the last created.
			chains[chains.size() - 1].push_back(firstChain);
			
			// Now check if there are any duplicate chains declared.
			if (pdb.header[i].substr(indexStartOfChain + 1, 1) == ",")
			{
				// Continue to add until we hit an empty space, at which point 
				//	there are no more duplicate chains to add. 
				for (int x = indexStartOfChain + 3 ; ; x += 3)
				{
					string dupChain = pdb.header[i].substr(x, 1);
					if (dupChain == " ")
						break;
					chains[chains.size() - 1].push_back(dupChain);
				}
			}
		}
	}
	
	// Now read in data for each of our non-duplicate chains,
	//	such as # of atoms and # of AAs.
	for (int currChain = 0; currChain < chains.size(); currChain++)
	{
		// Read in the chain, add the number of AA and create an int
		// 	in the atoms vector so it can be used to sum the total # of atoms
		//	in the chain. .
		pdb.read(pdbStr + ".pdb", chains[currChain][0]);
		sizeOfChainsAtoms.push_back(0);
		sizeOfChainsAA.push_back(pdb.numOfAA());
		
		// Iterate over each AA to get total # of atoms. 
		for (int iAA=0; iAA<pdb.numOfAA(); iAA++)
		{
			sizeOfChainsAtoms[currChain] += pdb.numOfAtoms(iAA);
		}
	}
	
	// Display all chains and related information.
	cout << "Chain\t# Atoms\t# Amino Acid\tDuplicate Chains " << endl;
	for (int i = 0; i < chains.size(); i++)
	{
		cout << chains[i][0] << "\t" << sizeOfChainsAtoms[i] <<"\t"<< sizeOfChainsAA[i] << "\t\t";
		
		if (chains[i].size() == 1)
		{
			cout << "No Duplicate Chains";
		}
		
		for (int x = 1; x < chains[i].size(); x++)
		{
			cout << chains[i][x] << " ";

		}
		cout << "\n";
	}
	
	cout<<"Please enter the ID of chain you want to cut the volume around: ";
    cin>>chainID;
	
	// Do some error handling to make sure chain that was entered correctly.
	for (int i = 0; i < chains.size(); i++)
	{
		// Verify it's a real chain that exists. 
		if (chains[i][0] == chainID)
		{
			return chainID;
		}
		else 
		{
			// Check if it's a duplicate chain someone wants to cut.
			for (int x = 0; x < chains[i].size(); x++)
			{
				if (chains[i][x] == chainID)
				{
					return chainID;
				}
			}
		}
	}
	
	// If we got this far, something went wrong (user input chain did not match
	//	an existing chain). Return following string to notify main so we exit. 
	return "BADSELECTION";
} 
