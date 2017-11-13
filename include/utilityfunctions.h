#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <cmath>
#include <sstream>
#include <algorithm>		//used mainly to generate permutations
#include <vector>
#include <iostream>

//#include "constants.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32
#include <time.h>					//used to find elapsed time
#else
#include <sys/time.h>
#endif

using namespace std;
#define SHOW_ERRORS 0		//1 : if you want the program to show you errors. 0: if you won't
#define BEEP 007			//the code of the beep sound of the internal speaker of the computer

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List Of Data Structures //////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct OnePermutation
{//used to store one permutation
	vector<int> permutation;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double getRadius (char atomType,char method = 'O');				//given the type of atom and the method ... return the radius of this atom
void errMsg(string className, string methodName, string msg, bool important = false);		//errMsg

template <class T> inline std::string toString (const T& t);				//convert a data type to string
#define MAX(x, y)   ((x) >= y ? (x) : (y))				// return max
#define MIN(x, y)	((x) < y ? (x) : (y))				// return min
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//error msg
// given the name of the class, method, and the error msg, print out the msg on the screen
//if the error msg is important then the program should stop working untill the user hit any key
void errMsg(string className, string methodName, string msg, bool important){

	cout<<endl;
	cout<<"================================ in "<<className<<"."<<methodName<<"() ===================="<<endl;
	cout<< msg<< endl;
	cout<<"==============================================================================="<<endl;
	if (SHOW_ERRORS || important)
	{
		putchar(BEEP);
		cout<<"Press any key..."<<endl;
		getchar();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//method V = VDW;  C= calculated; O= Covalent (default)
double getRadius(char atomType,char method)
{
	/*
	// covalent bond radius (2008 values) /www.webelements.com
	#define RADIUS_C        0.76
	#define RADIUS_N        0.71
	#define RADIUS_O        0.66
	#define RADIUS_H        0.31
	#define RADIUS_S        1.03
	#define RADIUS_P        1.07

	  // calculated radius  //www.webelements.com
	#define RADIUS_C        0.67
	#define RADIUS_N        0.56
	#define RADIUS_O        0.48
	#define RADIUS_H        0.53
	#define RADIUS_S        0.83
	#define RADIUS_P        0.98


	// van der waals force radius
	#define RADIUS_C        1.7
	#define RADIUS_N        1.55
	#define RADIUS_O        1.52
	#define RADIUS_H        1.2
	#define RADIUS_S        1.85
	#define RADIUS_P        1.9

	*/
	switch (method){
		case 'V': switch (atomType){
						/*
						case 'N': return 1.55;
						case 'C': return 1.7;
						case 'O': return 1.52;
						case 'P': return 1.9;
						case 'S': return 1.85;
						case 'H': return 1.2;
						*/
						/*
						case 'N': return 1.5;
						case 'C': return 1.7;
						case 'O': return 1.4;
						case 'P': return 1.9;
						case 'S': return 1.8;
						case 'H': return 1.0;
						*/
						///
						///			as in SCWRL 3.0 Dunbrack paper
						///
						case 'N': return 1.3;
						case 'C': return 1.6;
						case 'O': return 1.3;
						case 'P': return 1.8;
						case 'S': return 1.7;
						case 'H': return 1.0;

					}
		case 'C': switch (atomType){
						/*
						case 'N': return 0.56;
						case 'C': return 0.67;
						case 'O': return 0.48;
						case 'P': return 0.98;
						case 'S': return 0.83;
						case 'H': return 0.53;
						*/
						//Dr. Weitao
						case 'N': return 0.70;
						case 'C': return 0.77;
						case 'O': return 0.66;
						case 'P': return 1.1;
						case 'S': return 1.04;
						case 'H': return 0.32;

				  }
		case 'O': switch (atomType){

						case 'N': return 0.71;
						case 'C': return 0.76;
						case 'O': return 0.66;
						case 'P': return 1.07;
						case 'S': return 1.03;
						case 'H': return 0.31;

			/*
						// Dr. Weitao's values
						case 'N': return 0.73;
						case 'C': return 0.77;
						case 'O': return 0.74;
						case 'P': return 1.1;
						case 'S': return 1.03;
						case 'H': return 0.3;
			*/
				  }
	}

	cout<<"============================= in getRadius(char, char) ==============="<<endl;
	cout<<"The atomtype or method class sent are incorrect"<<endl;
	cout<<"======================================================================"<<endl;
	exit(1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T>
inline std::string toString (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

#endif


