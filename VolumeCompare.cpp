#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>

#include "include/skeleton_overall.h"
#include "include/MRC.h"

using namespace std;

// Define threshold to filter noise.
float const THRESHOLD_ORIGINAL = 0.25;
float const THRESHOLD_SIM = 0.25; 

// Calculate the Hausdorff distance
float calculateHausdorff(vector<Coordinate> &setA, vector<Coordinate> &setB); 

int getTimeTakenInMS(struct timeval &start, struct timeval &end);

// Arguments passed in should be Cut MRC and Simulated MRC
int main(int argc, char *argv[])
{	
	if (argc != 3)
	{
		cout << "Must invoke Volume Comparison with the following parameters:\n";
		cout << "CutMRC SimulatedMRC\n";
		return -1;
	}

	string mrcOrigStr, mrcSimStr;
	Map originalMRC, simulatedMRC;

	struct timeval start, end, startTotal;
	
	// Error handling for PDB and MRC done when the file is read in.
	mrcOrigStr = argv[1];
	mrcSimStr = argv[2];
	
	gettimeofday(&startTotal, NULL);

	// Read in the MRC file.
	cout << "Reading in MRC Files: " << (mrcOrigStr + ".mrc") << " and " << (mrcSimStr + ".mrc") << endl;
	originalMRC.read(mrcOrigStr + ".mrc");
	
	// Read in the Simulated MRC file.
	simulatedMRC.read(mrcSimStr + ".mrc");
	
	gettimeofday(&start, NULL);
	cout << "Files read in successfully, now Normalizing... "; 
	// Normalize the density data.
	originalMRC.normalize();
	simulatedMRC.normalize();
	
	cout << " Thresholding... ";
	// Set the threshold for both simulated and original data.
	originalMRC.filterize(THRESHOLD_ORIGINAL);
	simulatedMRC.filterize(THRESHOLD_SIM);
		
	int totalVoxels = originalMRC.numRows() * originalMRC.numCols() * originalMRC.numSlcs();	
	
	// Used to keep track of coordinates for calculating Hausdorff distance (done to reduce run time)
	vector<Coordinate> maskSimCoordinates;
	vector<Coordinate> originalCoordinates;
	
	// Estimate by using the simulatedMRC as a mask (anything)
	for (int rows = 0; rows < originalMRC.numRows(); rows++)
	{
		for (int cols = 0; cols < originalMRC.numCols(); cols++)
		{
			for (int slcs = 0; slcs < originalMRC.numSlcs(); slcs++)
			{
				// (Delete me later) Temp debug code to make sure thresholding is working...
				// if (simulatedMRC.cube[rows][cols][slcs] < 0.25 && simulatedMRC.cube[rows][cols][slcs] > 0.25)
				// {
				//	cout << " R / C / S " << rows << " " << cols << " " << slcs << " " << simulatedMRC.cube[rows][cols][slcs] << endl;
				// }
				
				// Write all "original" MRC coordinates.
				if (originalMRC.cube[rows][cols][slcs] != 0)
				{
					Coordinate temp;
					temp.x = rows;
					temp.y = cols;
					temp.z = slcs;
					originalCoordinates.push_back(temp);
				}
				
				// Write all "simulated" MRC coordinates.
				// Creating a "mask", if it is non-0, then we have a match.
				if (simulatedMRC.cube[rows][cols][slcs] != 0)
				{
					Coordinate temp;
					temp.x = rows;
					temp.y = cols;
					temp.z = slcs;
					maskSimCoordinates.push_back(temp);
				}
			}
		}		
	}
	cout << " Calculating the OSHD..." << endl;
	cout << "One Sided Hausdorff Distance Simulated to Original Data: " << calculateHausdorff(maskSimCoordinates, originalCoordinates) << endl;
	gettimeofday(&end, NULL);
	
	cout << "OSHD Time: " << getTimeTakenInMS(start, end) << " total time: " << getTimeTakenInMS(startTotal, end) << endl;
	
} 

// Calculates the one-sided Hausdorff distance from coordinate setA to setB
float calculateHausdorff(vector<Coordinate> &setA, vector<Coordinate> &setB)
{
	if (setA.size() == 0 || setB.size() == 0)
		return -1;

	float maxDist = 0;
	for (int a = 0; a < setA.size(); a++)
	{
		// Set number to max value (to simulate infinity). 
		float minDist = numeric_limits<float>::max();
		// Iterate through setB and compare all points to setB to find minDistance. 
		for (int b = 0; b < setB.size(); b++)
		{
			float distance = getDistance(setA[a], setB[b]);
			if (distance < minDist)
			{
				minDist = distance;
			}
		}
		if (minDist > maxDist)
		{
			maxDist = minDist;
		}
		
	}
	
	return maxDist;
}

// Calculate the time taken, returns ms taken as an int. 
int getTimeTakenInMS(struct timeval &start, struct timeval &end) 
{
	long seconds, useconds;    
	seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    return ((seconds) * 1000 + useconds/1000.0) + 0.5;
}
