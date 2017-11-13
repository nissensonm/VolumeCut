#ifndef PROTEIN_H
#define PROTEIN_H

#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include <stdio.h>


//#include "constants.h"
#include "geometry.h"
#include "utilityfunctions.h"


using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// Some Constants and Data Structures ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// list of missing AAs if any.
struct missingAA
{
	string chr3;				//AA in 3 letters representation
	char chr1;					//AA in one letter representation
	int num;					//sequence number
	char resInsertion;			//Code for insertions
	char SStype;				//the type of the AA...H helix, S strand, or L loop

	missingAA(): chr3(""), chr1(' '), num(0), resInsertion(' '), SStype('L') {}

};

								//////////////////////////////////////////////////////////////////////////
								//////////////////// Helices Secondary Structure Struct //////////////////
								//////////////////////////////////////////////////////////////////////////
struct HelicesSecondaryStructure
{
	int startIndx;		//The index (AA vector indx) of the first AA in the helix
	int endIndx;		//The index (AA vector indx) of the last AA in the helix
	int serialNum;		//Helix serial Number
	int hlxID;			//Helix identifier
	int nAAcur;			//The length (the missing AA's are not counted) (# of AA)  of the helix
	int nAA;			//The Actual length of the helix...missing AA's are counted
	string comment;		//any comment

	int type;			//Contains the type of the helix

						//The type as follows:
						// 1  Right-handed alpha (default)		6  Left-handed alpha
						// 2  Right-handed omega				7  Left-handed omega
						// 3  Right-handed pi					8  Left-handed gamma
						// 4  Right-handed gamma				9  2/7 ribbon/helix
						// 5  Right-handed 3/10					10  Polyproline

	//Initializer
	HelicesSecondaryStructure() : startIndx(-1), endIndx(-1), serialNum(0), hlxID(0), nAAcur(0), nAA(0), comment(""), type(1) {}
};
								//////////////// END OF Helices SECONDARY STRUCTURES /////////////////////

								//////////////////////////////////////////////////////////////////////////
								//////////////////// Sheets Secondary Structure Struct ///////////////////
								//////////////////////////////////////////////////////////////////////////
struct SheetsSecondaryStructure
{
	//store the information of the strand
	int startIndx;		//The index (AA vector indx) of the first AA in the strand
	int endIndx;		//The index (AA vector indx) of the last AA in the strand
	int strandNum;		//strand number
	string sheetID;		//sheet identifier
	int nStrand;		//number of strands in current sheet

	int nAAcur;			//The length (the missing AA's are not counted) (# of AA)  of the beta strand
	int nAA;			//The Actual length of the strand...missing AA's are counted
	int sense;			//Contains strand sense with respect to previous
						//0 first strand			1 parrallel		-1 antiparallel
	//The following fields identify two atoms involved in a hydrogen bond, the first in the current strand and the second in the previous strand.  \
    //These fields should be blank for strand 1 (the first strand in a sheet).
	string curAtomName;		//the name of the atom in the current strand involved in a hydrogen bond
	string prevAtomName;	//the name of the atom in the previous strand involved in a hydrogen bond
	int AACurIndx;			//the index of the AA in the current strand involved in a hydrogen bond
	int	AAPrevIndx;			//the index of the AA in the previous strand involved in the hydrogen bond


	//Initializer
	SheetsSecondaryStructure() : startIndx(-1), endIndx(-1), strandNum(0), sheetID(""), nStrand(0), nAAcur(0), nAA(0), sense(1), AACurIndx(-1), AAPrevIndx(-1) {}
};
								//////////////// END OF Sheets SECONDARY STRUCTURES /////////////////////

								/////////////////////////////////////////////////////////////////////////
								////////////////////// ATOM STRUCTURE ///////////////////////////////////
								/////////////////////////////////////////////////////////////////////////
struct Atom
{
	string name;								//Atom name ...considered as 4 characters
	char locIndicator;							//Alternate_location_indicator			position (16,1)
	string occupancy;							//Occupancy							//	position (54,6)		right alignment
	string tempFactor;							//Temperature factor				//	position (60,6)		right alignment
	string charge;								//Charge							//	position (78,2)

	Coordinate coord;							// X,Y,Z position of the atom

	char type;									//The type of the atom..(N,C,O....)
	bool isSideChain;							//side chain flag....true if it is a side chain atom...false otherwise

	//Initializer
	Atom() : name("NONE"), locIndicator(' '), occupancy(""), tempFactor(""), charge(""), type(' '), isSideChain(false) {}
};
								////////////////// END OF ATOM STRUCTURE ////////////////////////////////


								/////////////////////////////////////////////////////////////////////////
								/////////////////// AMINO ACID STRUCTURE ////////////////////////////////
								/////////////////////////////////////////////////////////////////////////
struct AminoAcid
{
	vector<Atom> atoms;			//list of atoms from atom struct
	/*
					atoms structure will look like
					N
					N-H atoms if any
					Ca
					Ca-H atoms if any
					CB if any
					CB-H atoms if any
					SideChain atoms followed by their H atoms
					...
					...
					C
					O
					OXT
	*/
	char chr1;					//AA name (1-character)
	string chr3;				//AA name (3-character)
	string chain;				//the current chain
	int num;					//AA serial number...same as the number in PDB file
	char resInsertion;			//Code for insertions of residues	//	position (26,1)
	Coordinate coord;			//could be the center of side-chain or the all-atoms center
	char SStype;				//stores the type of secondary structure this AA represent
								//  H : Helix
								//	S : Sheet
								//  L : Loop or other;
	Coordinate ScEndPoint;		//The end point that represents the end of the side chain.... a vector from CA to this point represent the size and
								// direction of the side chain
	char whtInCoord;			//what is stored in coord variable...could be
								//	S : side-chain center
								//	A : All-atoms center
								//  B : backbone-atoms center
								//	N : Nothing Interesting (initial)
	char atomsIncluded;			//H heavy atoms, A all, N nothing (initial)..included in calculation stored in coord
	double gyration;			//stores the gyration radius of amino acid (for side chain)...if it is equal to 0 then gyCoord is invalid
	Coordinate gyCoord;		    // cotains the Coordinate of the radius of gyration
	Torsion angles;				//data structure to store phi psi angles

	//Initializer
	AminoAcid() : chr1(' '), chr3(""), chain(""), num(0), resInsertion(' '), whtInCoord('N'), atomsIncluded('N'), gyration(0), SStype('L') {}
};
								////////////////// END OF AMINO ACID STRUCTURE /////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// End Of Constants and Data Structures /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// List of Functions ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ////////// End Of The List //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// PROTEIN CLASS PROTOTYPE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Protein
{
public:
	Protein();									//The default constructor
	Protein(string, string = "");				//Another Constructor..given file name (path) and chain specifier
	void initialize();							//initialize variables


	string ID;									//The ID of the protein that the portion came from
	Coordinate centOfCharge;					//center of charge
	vector<AminoAcid> AAs;						//Contains the list of AA's that form the Portion
	vector<missingAA> missingAAs;				// a list of amino acids are missing from 3D structure
	vector<string> header;						//The header information from the PDB file
	vector<HelicesSecondaryStructure> hlces;			//The list of hlces in the portion
	vector<SheetsSecondaryStructure> sheets;			//The list of sheets in the portion
	vector<int> AAsCollide;						//stores the Indeces of first 2 AA's that collide with each other ...
	vector<int> atomsCollide;					//stores the Indeces of atoms that collide with each other within the 2 AA's collide from above variable


	void read(string, string = "");					//given a pdb file path and a target chain....reads the PDB file
	void writePDB(string, int, int, bool = false);	//write a specified range of the portion to a PDB file

	inline int numOfAA();						//returns the number of AA in the Protein
	inline int numOfAtoms(int);					//given the indx of AA....returns the number of atoms in that AA
	int numOfSCAtoms(int, char = 'H');			//get the number of side chain atoms, H: heavy or A: all
	int getAAIndx(int);							//given AA num (as in PDB) returns the indx of this AA in the AA vector .. .-1 if not found
	int getAtomIndx(int, string);				//given AA indx and the name of atom (or substring)...returns the indx of this atom or -1 if not exist

	Coordinate getAtomCoordinate(int, string);	//return the coordinate of a given atom (atom name or a substring of it) in a particular AA



	void fillMissingAAs();						//fill missing AAs if any in missingAAs data structure

	void setSCCenter(int, char = 'H');			//given an index of an AA....this stores the center of sidechain in AAs.coord and set whtInCoord to S

	void setBBCenter(int, char = 'H');			//given AA indx ... copmute the center(mass center) of all backbone atoms
	void setRgyration(int, char ='H');			//set the radius of gyration and the gyration coordinate...H heavy atoms included or A all atoms

	void setAxisSegments(int, int, double, vector<Coordinate> &);		//stor the axis of a segment every "given distance" as a set of points
	void buildSS();								//build secondary structure vector.....

	void setScEndPoint(int, char = 'H');//set the side chain end point for a given AA, H: heavy atoms to be included A: all atoms


private:
	string path;								//the path from where the protein came
	char chr3ToChr1(string);					//converts from 3-letters format to 1-letter format
	char getAtomType(string);					//returns the type of given Atom...N, C, O, or S...
	bool isSideChainAtom(string);				//true if the atom is a sidechain atom

	AminoAcid reOrderAtoms(int);

	bool isHeavyAtom(char);						//check if the atom is heavy or not
	void sortHlces(vector<HelicesSecondaryStructure> &, const unsigned int, unsigned int);	//sort hlces according to the indeces of the first AA

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// END OF PROTEIN CLASS PROTOTYPE /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// PROTEIN Class Implementation /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Protein::Protein()		//constructor
{
	//initialize variables
	initialize();
}
////////////////////////////////////////////////////////////////////////////////////////
Protein::Protein(string filePath, string chain)		//constructor
{
	read(filePath, chain);
}

////////////////////////////////////////////////////////////////////////////////////////////
void Protein::initialize()
{
	path = "Unknown";
	ID = "Unknown";
	AAs.clear();
	missingAAs.clear ();
	header.clear();
	hlces.clear();
	sheets.clear();
	AAsCollide.clear ();
	atomsCollide.clear();

}
/////////////////////////////////////////////////////////////////////////////////////////////

// read PDB file
void Protein::read(string fileName, string chain)
{
	string line,tmpChain;

	ifstream infile;

	bool exitReading = false;

	//initialize variables .... every time u wanna read a protein
	initialize();

	path = fileName;		//the path is where the portion came from

	//get portion ID
	if (fileName.length() > 8)
		if ((atoi(fileName.substr (fileName.length () - 9,1).c_str ()) != 0) &&		// "\" symbol in windows
			(fileName.substr (fileName.length () - 9,1) != "/") &&					// "/" symbol
			(atoi(fileName.substr (fileName.length ()- 9,1).c_str()) != 92 ))		// "\" symbol in linux
			ID = fileName.substr(fileName.length()-9,4);
		else
			//the path does not contain "\" or "/" but it may contain the chain ID
			ID = fileName.substr (fileName.length () - 8,4);
	else
		if (fileName.length() == 8)
			ID = fileName.substr (fileName.length ()-8,4);
		else
			//the ID in this case could be wrong since the name is irrelative to the ID
			ID = fileName.substr (0, fileName.length() - 4);

	infile.open(fileName.c_str());
	if (!infile)
	{
		errMsg("Protein", "read", "Unable to open "+ fileName);

		exit(1);
	}
	else
	{
		AminoAcid aa;
		Atom tmpAtom;
		while((!infile.eof()) && (!exitReading))
		{

			getline(infile, line);													//get a line from pdb file

			if ((line.substr(0,4) == "ATOM") &&										//ATOM word token
				(line[17] != ' ') &&												//not a DNA
				((line.substr(21,1).c_str() == chain) || (chain == "")) &&			//same target chain
				(!exitReading))														//remain atoms to be read
			{
				/*****
						Samples:
						ATOM      1  N   VAL A   3      18.481  17.489  41.966  1.00 33.08           N
						ATOM      6  CG1 VAL A   3      14.652  17.181  41.934  1.00 33.32           C
				*****/

				/*****
						 begins with first chain if no target chain was sent
				*****/
				if ((numOfAA() == 0))
				{
					tmpChain = line.substr(21,1).c_str();
					aa.num = -1000; //atoi(line.substr(22,4).c_str());

				}
				/*****
						the current version of the system does not work with portions have code of residue inserting....
						the portion will be all the part before the first AA has a code of inserting
				*****/
				if (line[26] != ' ')
				{
					cout<<ID<<" has code of residue inserting"<<endl;
					exitReading = true;
				}

				/*****
						saves the first chain, whtever this chain is A,B,or C...etc		terminate if u reach
				*****/
				if (tmpChain != line.substr(21, 1).c_str() ||								//another chain or
					//(atoi(line.substr(22,4).c_str()) < aa.num) ||						//AA with less seq num
					(line[26] != ' ') ||												//Code for inserting residue
					(line.substr(0, 3) == "TER"))										//TER reserved word
				{
					exitReading = true;					
				}
				else
				{
					/*****
							check if we reach a new AA and it is not a dublicate AA (by locIndicator)
					*****/
					bool isNew=true;
					for (int n=0; n<AAs.size(); n++)
					{
					    if (atoi(line.substr(22,4).c_str()) == AAs[n].num)
					    {
					        isNew = false;
					    }
					}

					if (isNew == true)
					{
						aa.chr3			= line.substr(17,3);
						aa.num			= atoi(line.substr(22,4).c_str());
						aa.chain		= line.substr(21,1).c_str();
						aa.chr1			= chr3ToChr1(line.substr(17,3).c_str ());
						aa.resInsertion = line[26];

						/*****
								if it is the only alternative take it...or it has more than one alternative...take alternative A.. or B.. or C
						*****/
						if ((line[16] == ' ') || (line[16] == 'A') || (line[16] == 'B') || (line[16] == 'C'))
							AAs.push_back(aa);							//Push the information of the AA
					}

					/*****
							if the atom is not for a duplicated AA then take it or it is alternative A.. or B.. or C
					*****/
					if ((line[16] == ' ') || (line[16] == 'A') || (line[16] == 'B') || (line[16] == 'C'))
					{
						tmpAtom.name         = line.substr(12,4);
						tmpAtom.locIndicator = line[16];
						tmpAtom.coord.x		 = atof(line.substr(30,8).c_str());
						tmpAtom.coord.y		 = atof(line.substr(38,8).c_str());
						tmpAtom.coord.z		 = atof(line.substr(46,8).c_str());
					//	tmpAtom.type		 = line[77];
						/*****
								if one of atoms has no type then get it
						*****/
					//	if (int (tmpAtom.type) == 0)
							tmpAtom.type = getAtomType(tmpAtom.name);
						/*****
								Extra Information for atom
						*****/
						if (line.length()>54)
							tmpAtom.occupancy    = line.substr(54,6).c_str();
						if (line.length()>60)
							tmpAtom.tempFactor	 = line.substr(60,6).c_str();
						if (line.length()>78)
							tmpAtom.charge       = line.substr(78,2).c_str();
						/*****
								check if it is a side chain atom..set the flag
						*****/
						tmpAtom.isSideChain  = isSideChainAtom(tmpAtom.name);
						/*****
								push the atom to the AAs DataStructure
						*****/
						AAs[AAs.size()-1].atoms.push_back(tmpAtom);

					}
				}
			}
			else
			{

				/*****
						save header information from the beginning of the Pdb file..exclude any chain before the specified one (if any has been specified)
				*****/
				if ((!exitReading) && (numOfAA()==0) && (line.substr(0,4) != "ATOM"))
				{
					//cout<<line<<endl;
					header.push_back(line);
				}
			}
		}

		/*****
				close the file
		*****/
		infile.close();

		/*****
				reOrder Atoms...so every atom is followed by its H atoms in the AAs vector..last three atoms are C, O, and OXT
				and set the Side chain end point (where the line from Ca to this end point represents the direction and the length of the side chin
		*****/
		int i;
		for (i=0;i<numOfAA();i++)
		{
			AAs[i] = reOrderAtoms(i);
			setScEndPoint(i);
		}

		/*****
				build the secondary structure information
		*****/
		buildSS();

		/*****
				set the information of secondary structure for each AA
		*****/
		for (i=0;i<numOfAA(); i++)
		{
			for (int m=0; m<sheets.size (); m++)
				if ((i>= sheets[m].startIndx) && (i <= sheets[m].endIndx ))
					AAs[i].SStype = 'S';
			for (int k=0; k<hlces.size (); k++)
				if ((i>=hlces[k].startIndx ) && (i<= hlces[k].endIndx  ))
					AAs[i].SStype = 'H';
		}

		/*****
					fill missingAAs data structure for any missing residue from the 3D structure
		*****/
		fillMissingAAs();
	}
}
////////////////////////////////////////////////////////////////////////////////
void Protein::writePDB(string outFile,int startAARank, int endAARank, bool wHeader)
// The whole Protein could be written to a file or a portion of it
//startAARank is the rank of the first AA to be written to outfile...for example 4 is the 4th AA
//endAARank is the rank of the last AA would be written to the outfile PDB file
//if you want to print with the header info in the original pdb file..then send wHeader = true;
{
	int tmpNumOfAA = numOfAA();

	//assure from the ranks given
	if ((startAARank>=1) && (startAARank <= tmpNumOfAA) && (endAARank <= tmpNumOfAA))
	{
		int i,
			j;
		ofstream out;
		out.open(outFile.c_str());
		if (!out) {
			errMsg("Protein", "writePDBFile", "Unable to open " + outFile, true);
		}

		/*****
				write header information...if the header information is chosen to be written then no need to
				do the next operation which writes the secondary structures in the specified range if the SS
				has been built using buildSS() before;
		*****/
		if (wHeader)
			for (i=0;i<header.size();i++)
				out<<header[i]<<endl;
		else
		{
			/*****
					In the case of buildSS() has been called, SS existance, and wHeader flag = false....write the list of SS first
			*****/
			int hCounter = 1;	//counter for helices
			int sCounter = 1;	//counter for sheets

			/*****
					write hlces if any
			*****/
			for (i=0;i<hlces.size();i++)
			{
				//The whole hlx is within the range
				if ((startAARank - 1 <= hlces[i].startIndx) && (endAARank - 1 >= hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hlces[i].serialNum<<" "<<setw(3)<<right<<hlces[i].hlxID<<" "<<AAs[hlces[i].startIndx].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "
					   <<setw(4)<<right<<AAs[hlces[i].startIndx].num<<AAs[hlces[i].startIndx].resInsertion<<" "<<AAs[hlces[i].endIndx].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right<<AAs[hlces[i].endIndx].num
					  <<AAs[hlces[i].endIndx].resInsertion<<setw(2)<<right<<hlces[i].type<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right<<hlces[i].nAAcur<<endl;
						++hCounter;
						continue;
				}

				//both ends within the range
				if ((startAARank > hlces [i].startIndx ) && (endAARank < hlces [i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "<<
					//new AA1
					AAs[startAARank-1].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "<<setw(4)<<right<<AAs[startAARank-1].num<<AAs[startAARank-1].resInsertion<<" "<<
					//new AA2
					AAs[endAARank-1].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right<<AAs[endAARank-1].num<<AAs[endAARank-1].resInsertion<<setw(2)<<right<<
					hlces[i].type <<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right<<endAARank - startAARank + 1<<endl;
					++hCounter;
					continue;
				}
				//The lower end is within the range
				if ((startAARank - 1<= hlces[i].startIndx) && (endAARank > hlces[i].startIndx) && (endAARank - 1 < hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "<<AAs[hlces[i].startIndx].chr3<<" "<<
						AAs[hlces[i].startIndx ].chain<<" "<<setw(4)<<right<<AAs[hlces[i].startIndx].num<<AAs[hlces[i].startIndx].resInsertion<<" "
					   //new AA
					   <<AAs[endAARank - 1].chr3<<" "<<AAs[hlces[i].endIndx].chain<<" "<<setw(4)<<right
		 			   //endAARank - 1 is the new end right now
					   <<AAs[endAARank - 1].num<<AAs[endAARank-1].resInsertion<<setw(2)<<right<<hlces[i].type<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right
					   //The new length
					   <<endAARank - hlces[i].startIndx<<endl;
					   ++hCounter;
					   continue;
				}

				//The upper end is within the range
				if ((startAARank -1 > hlces[i].startIndx) && (startAARank -1 <= hlces[i].endIndx) && (endAARank - 1 >= hlces[i].endIndx))
				{
					out<<"HELIX  "<<setw(3)<<right<<hCounter<<" "<<setw(3)<<right<<hCounter<<" "
						//New AA
						<<AAs[startAARank - 1].chr3<<" "<<AAs[hlces[i].startIndx].chain<<" "<<setw(4)<<right
						//startAARank - 1 is the new start right now
						<<AAs[startAARank -1 ].num<<AAs[startAARank-1].resInsertion<<" "<<AAs[hlces[i].endIndx].chr3<<" "<<AAs[hlces[i].endIndx].chain
						<<" "<<setw(4)<<right<<AAs[hlces[i].endIndx].num<<AAs[hlces[i].endIndx ].resInsertion<<setw(2)<<right<<hlces[i].type
						<<setw(30)<<left<<hlces[i].comment<<" "<<setw(5)<<right
					  //new length
					  <<hlces[i].endIndx - startAARank + 2<<endl;
						++hCounter;
						continue;
				}
			}

			/*****
					write sheets if any...the operation of hlces above should be done for sheets as well
			*****/
			for (i=0;i<sheets.size();i++)
			{
				//if ((startAARank - 1 <= sheets[i].startIndx) && (endAARank - 1 >= sheets[i].endIndx))
				//{
					out<<"SHEET  "<<setw(3)<<right<<sheets[i].strandNum<<" "<<sheets[i].sheetID<<setw(2)<<sheets[i].nStrand<<" "
						<<AAs[sheets[i].startIndx ].chr3<<" "<<AAs[sheets[i].startIndx].chain<<setw(4)<<AAs[sheets[i].startIndx].num
						<<AAs[sheets[i].startIndx ].resInsertion<<" "<<AAs[sheets[i].endIndx].chr3<<" "<<AAs[sheets[i].endIndx].chain
						<<setw(4)<<AAs[sheets[i].endIndx].num<<AAs[sheets[i].endIndx].resInsertion<<setw(2)<<sheets[i].sense<<" ";
					//write information for AA's involved in Hydrogen Bond starting from second strand
					if (sheets[i].strandNum != 1)
					{
						out<<sheets[i].curAtomName<<AAs[sheets[i].AACurIndx].chr3<<" "<<AAs[sheets[i].AACurIndx].chain<<setw(4)<<AAs[sheets[i].AACurIndx].num
							<<AAs[sheets[i].AACurIndx ].resInsertion <<" "<<sheets[i].prevAtomName <<AAs[sheets[i].AAPrevIndx ].chr3 <<" "
							<<AAs[sheets[i].AAPrevIndx ].chain <<setw(4)<<AAs[sheets[i].AAPrevIndx ].num<<AAs[sheets[i].AAPrevIndx ].resInsertion;

					}
					out<<endl;

					++sCounter;
				//}
				//the 2 cases left r to be written later
			}
		}

		/*****
				write Atoms
		*****/
		int atomsCounter = 1;
		for(i=startAARank-1;i<endAARank;i++)
		{
			for (j=0;j<numOfAtoms(i);j++)
			{

				out<<"ATOM"<<setw(7)<<atomsCounter<<" "<<setw(4)<<AAs[i].atoms[j].name<<AAs[i].atoms [j].locIndicator<<setw(3)<<right<<AAs[i].chr3
					<<setw(2)<<right<<AAs[i].chain<<setw(4)<<right<<AAs[i].num<<AAs[i].resInsertion<<setw(3)
					<<" "<<setw(8)<<right<<setiosflags(ios::fixed) << setprecision(3)<<AAs[i].atoms[j].coord.x<<setw(8)<<right<<AAs[i].atoms[j].coord.y
					<<setw(8)<<right<<AAs[i].atoms[j].coord.z<<setw(6)<<right<<AAs[i].atoms[j].occupancy<<setw(6)<<right<<AAs[i].atoms[j].tempFactor
					<<setw(10)<<" "<<setw(2)<<right<<AAs[i].atoms[j].type <<setw(2)<<AAs[i].atoms[j].charge<<endl;
				atomsCounter++;
			}
		}

		out.close();
	}
	else
	{
		string eMsg = "The given range (";
		eMsg += toString(startAARank);
		eMsg += ", ";
		eMsg += toString(endAARank);
		eMsg += ") is incorrect..or the portion os empty (no.AA= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "writePDB", eMsg, true);

	}
}


////////////////////////////////////////////////////////////////////////////
//number of AA in the whole portion
inline int Protein::numOfAA()
{
	return AAs.size();
}

///////////////////////////////////////////////////////////////////////////////
//number of atoms for a particular AA
inline int Protein::numOfAtoms(int AAIndx)
{
	if ((AAIndx>=0) && (AAIndx<numOfAA()))
		return AAs[AAIndx].atoms.size();
	else
	{
		string eMsg = "The given index (";
		eMsg += toString (AAIndx);
		eMsg += ") is out of range..or the portion is empty (no.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "numOfAtoms", eMsg);
	}

	return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
//atomsincluded: H heavy atoms, A all
int Protein::numOfSCAtoms(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		int cntr=0, i;
		/*****
				look for Heavy atoms
		*****/
		if (atomsIncluded == 'H')
			for (i=0; i<numOfAtoms(AAIndx); i++)
				if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
					cntr++;
		/*****
				look for All atoms in the side chain (H atom will be included)
		*****/
		if (atomsIncluded == 'A')
			for (i = 0; i<numOfAtoms(AAIndx); i++)
				if (AAs[AAIndx].atoms[i].isSideChain)
					cntr++;
		return cntr;
	}
	else
	{
		string eMsg = "AAindx given (";
		eMsg += toString(AAIndx);
		eMsg += "is out of range..";

		errMsg("Protein", "numOfSCAtoms", eMsg);
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////////
//given AA sequence number...return the indx in AAs
int Protein::getAAIndx(int AANum)
{
	for (int i=0;i<numOfAA();i++)
		if (AAs[i].num == AANum)
			return i;

	//Not Found
	return -1;
}
///////////////////////////////////////////////////////////////////////////////
//given AA indx in AAs and atom name (or a substring of atom name)...return atom indx in AAs.atoms
int Protein::getAtomIndx(int AAIndx, string atomName)
{

	if ((AAIndx>=0) && (AAIndx < numOfAA()))
	{
		for (int i=0;i<numOfAtoms(AAIndx);i++)
		{
			if (AAs[AAIndx].atoms[i].name.find(atomName) != AAs[AAIndx].atoms[i].name.npos)
				return i;		//The indx of the atom
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "getAtomIndx", eMsg);
	}

	//Not found
	return -1;

}


///////////////////////////////////////////////////////////////////////////////////////
//given the indx of AA in AAs and the name of the atom in that AA...return the coordinate of that atom
Coordinate Protein::getAtomCoordinate(int AAIndx, string atomName)
{
	Coordinate coord;

	int indx = getAtomIndx(AAIndx,atomName);
	if ( indx != -1 ) //-1 means not found
		return AAs[AAIndx].atoms[indx].coord;
	else
	{
		string eMsg = "No such atom ( " + atomName;
		eMsg += ") found ";
		if (AAIndx >=0 && AAIndx < numOfAA()){

			eMsg += "in AA (";
			eMsg += AAs[AAIndx].chr3 ;
			eMsg += "-";
			eMsg += toString(AAs[AAIndx].num);
			eMsg += " )";

			errMsg("Protein", "getAtomCoordinate", eMsg);

		}
		else{

			eMsg += "...pad AA indx (";
			eMsg += toString(AAIndx);
			eMsg += ").";

			errMsg("Protein", "getAtomCoordinate", eMsg);
		}

		coord.x = -999.0;
		coord.y = -999.0;
		coord.z = -999.0;

	}

	//Not Found....return the initial values 0 0 0
	return coord;
}
///////////////////////////////////////////////////////////////////////////////////////
void Protein::fillMissingAAs ()
{


	if ((missingAAs.size () == 0) && (numOfAA()))		//fill in the case it is not filled before
	{
		bool stop = false;
		int i = 0;
		size_t found;
		while ((i<header.size()) && (!stop))
		{
			found = header[i].find("MISSING RESIDUES");		//REMARK 465 referes to missing residues
			if ( found != header[i].npos)
				stop = true;
			i++;
		}

		if (stop)		//found some missing residues
		{
			found = header[i].find("RES C SSSEQI");
			while (header[i].find("RES C SSSEQI") == header[i].npos)
				i++;
			i++;

			missingAA tmpmissAA;
			string model;
			for ( int cntr = i; (cntr<header.size() && header[cntr].substr (0, 10) == "REMARK 465"); cntr++)
			{
				if (header[cntr].substr (19,1) == AAs[0].chain)		//should be in the same chain
				{
					if (missingAAs.size () == 0)
						model = header[cntr].substr (13, 1).c_str ();

					if (header[cntr].substr (13,1).c_str () == model)		//if in the same model
					{
						tmpmissAA.chr3 = header[cntr].substr(15, 3).c_str ();
						tmpmissAA.chr1 = chr3ToChr1(tmpmissAA.chr3);
						tmpmissAA.num  = atoi(header[cntr].substr (22, 4).c_str ());
						tmpmissAA.resInsertion = header[cntr][26];


						bool checkStrands = true;

						//check helices
						int SScntr = 0;
						while(SScntr < hlces.size())
						{
							if ((tmpmissAA.num >= AAs[hlces[SScntr].startIndx].num) &&
								(tmpmissAA.num <= AAs[hlces[SScntr].endIndx].num))
							{
								checkStrands = false;
								tmpmissAA.SStype = 'H';
								break;

							}
							SScntr++;
						}
						if (checkStrands)
						{
							SScntr = 0;
							while (SScntr < sheets.size ())
							{
								if ((tmpmissAA.num >= AAs[sheets[SScntr].startIndx].num) &&
									(tmpmissAA.num <= AAs[sheets[SScntr].endIndx].num))
								{
									tmpmissAA.SStype = 'S';
									break;
								}
								SScntr++;
							}
						}
						missingAAs.push_back (tmpmissAA);
					}
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////
//calculate the center of sidechain and store it in AAs.coord
void Protein::setSCCenter(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			double segmaX  = 0,
					segmaY = 0,
					segmaZ = 0;
			int numSC = 0;  //counter for number of side chain atoms included...heavy atoms or all

			int i;
			if (atomsIncluded == 'H')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
					if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numSC++;
					}
			}
			if (atomsIncluded == 'A')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
					if (AAs[AAIndx].atoms[i].isSideChain)
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numSC++;
					}
			}
			//if sidechain atoms are exist
			if (numSC)
			{
				AAs[AAIndx].coord.x = segmaX/numSC;
				AAs[AAIndx].coord.y = segmaY/numSC;
				AAs[AAIndx].coord.z = segmaZ/numSC;
			}
			//if no side chain atoms were found
			else
				//The center of SC is considered to be Ca
				AAs[AAIndx].coord = getAtomCoordinate(AAIndx," CA ");

			AAs[AAIndx].whtInCoord = 'S';		//set the flag to indicate that side chain mass center is stored in coord variable
			AAs[AAIndx].atomsIncluded = atomsIncluded;	//set the type of atoms were included in the calculations
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setSCCenter", eMsg);

	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::setBBCenter (int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'B') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			double segmaX  = 0,
					segmaY = 0,
					segmaZ = 0;

			int numOfAtomsIncluded = 0;
			int i;

			//heavy atoms in calculation
			if (atomsIncluded == 'H')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					//heavy and not a side chain atom
					if ((isHeavyAtom(AAs[AAIndx].atoms[i].type)) && (!AAs[AAIndx].atoms [i].isSideChain))
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numOfAtomsIncluded++;
					}
				}
			}
			//all atoms in calculation
			if (atomsIncluded == 'A')
			{
				for (i=0;i<numOfAtoms(AAIndx);i++)
				{
					//not a side chain atom
					if (!AAs[AAIndx].atoms [i].isSideChain)
					{
						segmaX += AAs[AAIndx].atoms[i].coord.x;
						segmaY += AAs[AAIndx].atoms[i].coord.y;
						segmaZ += AAs[AAIndx].atoms[i].coord.z;
						numOfAtomsIncluded++;
					}
				}
			}

			//if AA has atoms
			if (numOfAtomsIncluded)
			{
				AAs[AAIndx].coord.x = segmaX/numOfAtomsIncluded;
				AAs[AAIndx].coord.y = segmaY/numOfAtomsIncluded;
				AAs[AAIndx].coord.z = segmaZ/numOfAtomsIncluded;
			}

			AAs[AAIndx].whtInCoord = 'B';				//set the flag to indicate that Backbone mass center is stored in the coord variable
			AAs[AAIndx].atomsIncluded = atomsIncluded;  //set the type of atoms included in the calculations
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setBBCenter", eMsg);

	}
}
//////////////////////////////////////////////////////////////////////////////////////
//atomsincluded H heavy A all
void Protein::setRgyration(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if (AAs[AAIndx].gyration == 0)
		{
			int numSC = numOfSCAtoms(AAIndx, atomsIncluded);
			if (numSC == 0)
			{
				AAs[AAIndx].gyration = sqrt(getRadius('C', 'C')); // gyration radius is the same as the radius of Ca	.... calculated radius is used
				AAs[AAIndx].gyCoord .x = AAs[AAIndx].gyCoord .y = AAs[AAIndx].gyCoord .z = AAs[AAIndx].gyration ;
			}
			else
			{
				int numAtoms = numOfAtoms(AAIndx);
				//calculate mass center if it is not calculated before
				if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
					setSCCenter(AAIndx);

				Coordinate massCenter = AAs[AAIndx].coord;
				//cout<<"mass center = "<<massCenter.x <<" "<<massCenter.y <<" "<<massCenter.z<<endl;

				int i;

				double gyrationR = 0;
				//initialize gyration coordinate
				AAs[AAIndx].gyCoord .x = 0;
				AAs[AAIndx].gyCoord .y = 0;
				AAs[AAIndx].gyCoord .z = 0;
				double atomR,
						rx,
						ry,
						rz,
						rr;
				if (atomsIncluded == 'A')
				{
					for (i=0; i<numAtoms; i++)
					{
						if (AAs[AAIndx].atoms[i].isSideChain)
						{
							atomR = getRadius(AAs[AAIndx].atoms [i].type , 'C');
							rx = AAs[AAIndx].atoms[i].coord .x - massCenter.x;
							ry = AAs[AAIndx].atoms[i].coord .y - massCenter.y;
							rz = AAs[AAIndx].atoms[i].coord .z - massCenter.z;

							rr = sqrt(rx*rx + ry*ry + rz*rz) + atomR;
							gyrationR += rr * rr;

							AAs[AAIndx].gyCoord .x += (rx + atomR) * (rx + atomR);
							AAs[AAIndx].gyCoord .y += (ry + atomR) * (ry + atomR);
							AAs[AAIndx].gyCoord .z += (rz + atomR) * (rz + atomR);
						}
					}
				}
				if (atomsIncluded == 'H')
				{
					for (i=0; i<numAtoms; i++)
					{
						if ((AAs[AAIndx].atoms[i].isSideChain) && (isHeavyAtom(AAs[AAIndx].atoms[i].type)))
						{
							atomR = getRadius(AAs[AAIndx].atoms [i].type , 'C');
							rx = AAs[AAIndx].atoms[i].coord .x - massCenter.x;
							ry = AAs[AAIndx].atoms[i].coord .y - massCenter.y;
							rz = AAs[AAIndx].atoms[i].coord .z - massCenter.z;

							rr = sqrt(rx*rx + ry*ry + rz*rz) + atomR;
							gyrationR += rr * rr;

							AAs[AAIndx].gyCoord .x += (rx + atomR) * (rx + atomR);
							AAs[AAIndx].gyCoord .y += (ry + atomR) * (ry + atomR);
							AAs[AAIndx].gyCoord .z += (rz + atomR) * (rz + atomR);


							//cout<<rx<<" "<<ry<<" "<<rz<<" "<<gyrationR<<endl;
						}
					}
				}

				AAs[AAIndx].gyCoord .x = sqrt(AAs[AAIndx].gyCoord .x / numSC);
				AAs[AAIndx].gyCoord .y = sqrt(AAs[AAIndx].gyCoord .y / numSC);
				AAs[AAIndx].gyCoord .z = sqrt(AAs[AAIndx].gyCoord .z / numSC);

				gyrationR = gyrationR/numSC;
				AAs[AAIndx].gyration = sqrt(gyrationR);
			}
		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setRgyration", eMsg);

	}
}

///////////////////////////////////////////////////////////////////////////////////////
//the axis will be represented by set of points, the distance b/w the points will be given by "dist"
//except the last portion of the axis which could be less than this dist
void Protein::setAxisSegments (int startIndx, int endIndx, double dist, vector<Coordinate> &axis){

	string eMsg;
	//find the length of protein segment (in terms of # of AAs)
	short segLength = endIndx - startIndx + 1;
	//axis.clear();		//clear the data structure

	int i;

	if (segLength>=4){

		//set Backbone center for all AAs
		for (i=startIndx; i<=endIndx; i++)
			setBBCenter(i, 'H');		//work on heavy atoms only

		vector<Coordinate> triangles;

		//find center of triangles for each 3 consecutive AAs
		for (i=startIndx; i<=endIndx-2; i++)
			triangles.push_back (triangleCenter(AAs[i].coord, AAs[i+1].coord , AAs[i+2].coord));

		//find N terminal...if missing Ca ..error otherwise
		Coordinate terminalP,			//terminal point
					newP;				//the point where the line intersects the point
		short indx = getAtomIndx(startIndx, " N  ");
		if (indx != -1)
			terminalP = AAs[startIndx].atoms [indx].coord;
		else{
			indx = getAtomIndx(startIndx, " CA ");
			if (indx != -1)
				terminalP = AAs[startIndx].atoms [indx].coord;
			else{

				eMsg  = "First AA (";
				eMsg += AAs[startIndx].chr3 ;
				eMsg += "-";
				eMsg += toString(AAs[startIndx].num);
				eMsg += ") has no N or Ca atoms..";

				errMsg("Protein", "setAxisSegments", eMsg);

			}
		}
		//find the projection of this point with the first two triangles
		linePointIntersection(triangles[0], triangles[1], terminalP, newP);
		triangles.insert (triangles.begin (), newP);

		//find C terminal..if missing Ca....error otherwise
		indx = getAtomIndx(endIndx, " C  ");
		if (indx != -1)
			terminalP = AAs[endIndx].atoms [indx].coord;
		else{
			indx = getAtomIndx(endIndx, " CA ");
			if (indx != -1)
				terminalP = AAs[endIndx].atoms [indx].coord;
			else{

				eMsg  = "Last AA ("+AAs[endIndx].chr3 ;
				eMsg += "-";
				eMsg += toString(AAs[endIndx].num) ;
				eMsg += ") has no C or Ca atoms..";

				errMsg("Protein", "setAxisSegments", eMsg);
			}
		}
		//find the projection of this point with last two triangles
		linePointIntersection(triangles[triangles.size ()-2], triangles[triangles.size ()-1], terminalP, newP);
		triangles.push_back (newP);



		short 	sIndx	= 0,		//start indx
				eIndx	= triangles.size ()-1;		//end indx


		Coordinate curP = triangles[sIndx];	//current point we are working on...first point in the set of triangles

		axis.push_back(curP);


		int lastPointIndx = -1;
		i= sIndx + 1;
		while (i <= eIndx){
			if (getDistance(curP, triangles[i]) > dist){
				axis.push_back (pointOnLine(curP, triangles[i], dist));
				curP = axis[axis.size()-1];
				lastPointIndx = i;
			}
			else
				i +=1;
		}

		// add the last point if needed
		if (lastPointIndx <= eIndx)
			if (getDistance(curP, triangles[triangles.size ()-1]) > 0.0)
				axis.push_back (triangles[triangles.size ()-1]);


		/*
		cout<<"# of points in triangles = "<<triangles.size ()<<endl;
		double totalDist=0;
		for (i=0; i<triangles.size()-1; i++){
			cout<<i+1<<" distance = "<<getDistance(triangles[i], triangles[i+1])<<endl;
			totalDist += getDistance(triangles[i], triangles[i+1]);
		}

		cout<<totalDist<<endl;
		cout<<"=============="<<endl;
		for (i=0; i<axis.size()-1; i++)
			cout<<i+1<<"  distance = "<<getDistance(axis[i], axis[i+1])<<endl;
		*/

	}
	else{

		errMsg("Protein", "setAxisSegments", "Be sure the segment sent is longer than 3 AAs..no segments ");
	}

}
/////////////////////////////////////////////////////////////////////////////////////
void Protein::buildSS()
{

	if (numOfAA())
	{
		//string line;
		int i, j;

		//clear both vectors
		hlces.clear();
		sheets.clear();

		for (i=0; i<header.size(); i++)
		{

			HelicesSecondaryStructure tmpHlx;
			SheetsSecondaryStructure  tmpStrand;

			// if encounter HELIX
			if ((header[i].substr(0,5) == "HELIX") && (header[i].substr(19,1).c_str() == AAs[0].chain))
			{
				tmpHlx.serialNum		= atoi(header[i].substr(7, 3).c_str ());
				tmpHlx.hlxID			= atoi(header[i].substr(11, 3).c_str ());
				tmpHlx.startIndx		= getAAIndx(atoi(header[i].substr(21,4).c_str()));
				tmpHlx.endIndx			= getAAIndx(atoi(header[i].substr(33,4).c_str()));
				tmpHlx.type				= atoi(header[i].substr(38,2).c_str());
				tmpHlx.comment			= header[i].substr(40, 30);

				if ((tmpHlx.startIndx != -1) &&									//does the AA exist
					(tmpHlx.endIndx != -1) &&									//does the 2nd AA exist
					(header[i].substr(31,1) == header[i].substr(19,1)) )	//check the other chain ... normally the hlx should be in the same chain
				{

					tmpHlx.nAA				= abs(AAs[tmpHlx.endIndx ].num - AAs[tmpHlx.startIndx].num) + 1;
					tmpHlx.nAAcur		= tmpHlx.endIndx - tmpHlx.startIndx + 1;
					hlces.push_back(tmpHlx);
				}

			}

			//build sheets vector
			if ((header[i].substr(0,5) == "SHEET") && (header[i].substr(21,1).c_str() == AAs[0].chain))
			{

				tmpStrand.strandNum = atoi(header[i].substr (7, 3).c_str ());
				tmpStrand.sheetID	= header[i].substr (11, 3).c_str ();
				tmpStrand.nStrand	= atoi(header[i].substr (14, 2).c_str ());
				tmpStrand.startIndx = getAAIndx(atoi(header[i].substr(22,4).c_str()));
				tmpStrand.endIndx	= getAAIndx(atoi(header[i].substr(33,4).c_str()));

				tmpStrand.sense		= atoi(header[i].substr(38,2).c_str());

				if (header[i].length () > 44){

					tmpStrand.curAtomName = header[i].substr (41, 4);
					tmpStrand.AACurIndx		= getAAIndx(atoi(header[i].substr(50, 4).c_str()));
					tmpStrand.prevAtomName = header[i].substr (56, 4);
					tmpStrand.AAPrevIndx	= getAAIndx(atoi(header[i].substr(65, 4).c_str()));
				}

				if ((tmpStrand.startIndx != -1) &&
					(tmpStrand.endIndx != -1))
				{
					tmpStrand.nAAcur		=  tmpStrand.endIndx - tmpStrand.startIndx + 1;
					tmpStrand.nAA			= abs(AAs[tmpStrand.endIndx ].num - AAs[tmpStrand.startIndx].num) + 1;
					sheets.push_back(tmpStrand);
				}
			}
		}

		//sort hlces and strands according to the index of the first AA assigned to them
		if (hlces.size () > 0)
		{
			sortHlces(hlces, 0, hlces.size()-1);
			//remove duplicates
			i=0;
			while (i <hlces.size())
			{
				j = i+1;
				while (j< hlces.size())
				{
					if (hlces[i].startIndx == hlces[j].startIndx ){
						hlces.erase (hlces.begin () + j);		//delete and then decrement counter
						j--;
					}
					j++;
				}
				i++;
			}

		}
		if (sheets.size ())
		{
			//sortStrands(sheets, 0, sheets.size ()-1);

			//remove duplicates
			i=0;

			while (i < sheets.size())
			{
				j = i+1;
				while (j < sheets.size())
				{
					if (sheets[i].startIndx  == sheets[j].startIndx){
						sheets.erase (sheets.begin () + j);		//delete and then decrement counter
						j--;
					}

					j++;
				}
				i++;
			}

		}
	}
	else
	{
		string eMsg = "Be sure the portion (" + path;
		eMsg += ") contains AAs or you have read the pdb file..";

		errMsg("Protein", "buildSS", eMsg);

	}


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::setScEndPoint(int AAIndx, char atomsIncluded)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		if ((AAs[AAIndx].whtInCoord != 'S') || (AAs[AAIndx].atomsIncluded != atomsIncluded))
		{
			int CaIndx = getAtomIndx(AAIndx, "CA");

			if (CaIndx != -1)
			{
				setSCCenter(AAIndx, atomsIncluded);		//set the center of side chain according to the atom included heavy or all
				setRgyration(AAIndx, atomsIncluded);	//set the radius of gyration for the side chain

				Coordinate CaCoordinate = AAs[AAIndx].atoms [CaIndx].coord;

				Vectors v(AAs[AAIndx].coord, CaCoordinate);				//the vector from Ca atom to the massCenter
				double factorNum = v.length() / AAs[AAIndx].gyration;	// the number is used to divide the vector v by to get the exact vector to be added to the vector v
																		// to add the portion of the vector after the massCenter toward the end of the sidechain. this
																		//addition will be equal to the Rg
				Vectors vToBeAdded = v;
				vToBeAdded.divide (factorNum);							//get the vector of length equal to Rg to be added to the current vector (from CA to massCenter)
				v += vToBeAdded;

				AAs[AAIndx].ScEndPoint.x = v.getX() + CaCoordinate.x;
				AAs[AAIndx].ScEndPoint.y = v.getY() + CaCoordinate.y;
				AAs[AAIndx].ScEndPoint.z = v.getZ() + CaCoordinate.z;
			}
			else
			{
				string eMsg = "Ca atom is not found (" + AAs[AAIndx].chr3;
				eMsg += "-";
				eMsg += toString(AAs[AAIndx].num);
				eMsg += ") ....no EndPoint has been computed";

				errMsg("Protein", "setScEndPoint", eMsg);
			}


		}
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "setScEndPoint", eMsg);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
char Protein::chr3ToChr1(string chr3)
{

	if ( !chr3.compare("ALA"))
	return 'A';
	else if ( !chr3.compare("ASX"))
	return 'B';
	else if ( !chr3.compare("CYS"))
	return 'C';
	else if ( !chr3.compare("ASP"))
	return 'D';
	else if ( !chr3.compare("GLU") )
	return 'E';
	else if ( !chr3.compare("PHE") )
	return 'F';
	else if ( !chr3.compare("GLY") )
	return 'G';
	else if ( !chr3.compare("HIS") )
	return 'H';
	else if ( !chr3.compare("ILE") )
	return 'I';
	else if ( !chr3.compare("XLE"))
	return 'J';
	else if ( !chr3.compare("LYS") )
	return 'K';
	else if ( !chr3.compare("LEU") )
	return 'L';
	else if ( !chr3.compare("MET") )
	return 'M';
	else if ( !chr3.compare("ASN") )
	return 'N';
	else if ( !chr3.compare("PYL"))
	return 'O';
	else if ( !chr3.compare("PRO") )
	return 'P';
	else if ( !chr3.compare("GLN") )
	return 'Q';
	else if ( !chr3.compare("ARG") )
	return 'R';
	else if ( !chr3.compare("SER") )
	return 'S';
	else if ( !chr3.compare("THR") )
	return 'T';
	else if ( !chr3.compare("SEC"))
	return 'U';
	else if ( !chr3.compare("VAL") )
	return 'V';
	else if ( !chr3.compare("TRP") )
	return 'W';
	else if ( !chr3.compare("TYR") )
	return 'Y';
	else if ( !chr3.compare("GLX"))
	return 'Z';
	else
	{
		string eMsg = "Unknown character (X) has been returned for AA (" + chr3;
		eMsg += " )..";

		errMsg("Protein", "Chr3toChr1", eMsg);

		return 'X';		// Unknown Characters;
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////
//returns the type of the atom
char Protein::getAtomType(string atomName)
{
	size_t found;

	found = atomName.find("N");
	if (found!= atomName.npos)
		return 'N';

	found = atomName.find("C");
	if (found !=atomName.npos)
		return 'C';

	found = atomName.find("O");
	if (found != atomName.npos)
		return 'O';

	found = atomName.find("P");
	if (found != atomName.npos)
		return 'P';

	found = atomName.find("S");
	if (found != atomName.npos)
		return 'S';

	found = atomName.find("H");
	if (found != atomName.npos)
		return 'H';

	errMsg("Protein", "getAtomType", "Unknown atom type.....");

	exit(1);

}
//////////////////////////////////////////////////////////////////////////////////////////
bool Protein::isSideChainAtom(string atomName)
{
	if ((atomName == " N  ") || (atomName == " O  ") || (atomName == " C  ") || (atomName == " CA "))
		return false;

	if ((atomName == " H  ") || (atomName == " H1 ") || (atomName == " H2 ") || (atomName == " H3 "))
		return false;

	if ((atomName == " HA ") || (atomName == " HA1") || (atomName == " HA2") || (atomName == " HA3"))
		return false;
	if (atomName == " OXT")
		return false;

	return true;
}
////////////////////////////////////////////////////////////////////////////////////////////
AminoAcid Protein::reOrderAtoms (int AAIndx)
{
	if ((AAIndx >= 0) && (AAIndx < numOfAA()))
	{
		int	tmpNumOfAtoms = numOfAtoms(AAIndx);

		AminoAcid tmpAA;
		tmpAA = AAs[AAIndx];
		tmpAA.atoms .clear();
		for (int i=0;i<tmpNumOfAtoms;i++)
		{
			if ((AAs[AAIndx].atoms [i].name != " C  ") && (AAs[AAIndx].atoms [i].name != " O  ") && (AAs[AAIndx].atoms [i].name != " OXT"))
			{
				if (AAs[AAIndx].atoms [i].type != 'H')
				{
					size_t found;


					tmpAA.atoms .push_back (AAs[AAIndx].atoms [i]);		//push the atom
					string suffix = AAs[AAIndx].atoms [i].name .substr (2,2).c_str();
					//if the fourth letter of the atom name is not empty
					if (suffix.substr(suffix.length()-1, 1) != " ")
					{
						for (int j=i+1; j<tmpNumOfAtoms; j++)
						{
							found = AAs[AAIndx].atoms [j].name.find(suffix);

							if ((found != AAs[AAIndx].atoms[j].name.npos) &&
								(AAs[AAIndx].atoms [j].name != " C  ") &&
								(AAs[AAIndx].atoms [j].name != " O  ") &&
								(AAs[AAIndx].atoms [j].name != " OXT"))
							{
								tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);
							}
						}
					}
					else
						//if the third letter of the atom name is not empty ...for Ca and Cb or other 2 letters atom
						if (suffix.substr(suffix.length()-2,1) != " ")
						{
							for (int j=i+1; j<tmpNumOfAtoms; j++)
							{
								found = AAs[AAIndx].atoms [j].name.substr(2,2).find(suffix.substr(0,1).c_str());

								if ((found != AAs[AAIndx].atoms[j].name.npos) &&
									(AAs[AAIndx].atoms [j].name != " C  ") &&
									(AAs[AAIndx].atoms [j].name != " O  ") &&
									(AAs[AAIndx].atoms [j].name != " OXT"))
								{
									tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);
								}
							}
						}
						else
							//for N...push all H atoms
							if (suffix == "  ")
							{
								for (int j=i+1; j<tmpNumOfAtoms; j++)
								{
									if ((AAs[AAIndx].atoms [j].name == " H  ") || (AAs[AAIndx].atoms [j].name == " H1 ")
										|| (AAs[AAIndx].atoms [j].name == " H2 ") || (AAs[AAIndx].atoms [j].name == " H3 "))
										tmpAA.atoms .push_back (AAs[AAIndx].atoms [j]);

								}
							}
							else
							{
								errMsg("Protein", "reOrderAtoms", "Error occured in the name of the AA " + AAs[AAIndx].chr3);
								exit(1);

							}

				}
			}
		}

		//Always...the last 3 atoms are C, O, and OXT
		int cAtomIndx = getAtomIndx(AAIndx," C  ");
		if (cAtomIndx != -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[cAtomIndx]);		//push C atom
		int oAtomIndx = getAtomIndx(AAIndx," O  ");
		if (oAtomIndx != -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[oAtomIndx]);		//push O atom
		int oxtIndx = getAtomIndx(AAIndx," OXT");
		if ( oxtIndx!= -1)
			tmpAA.atoms.push_back(AAs[AAIndx].atoms[oxtIndx]);						//push OXT if exist

		return tmpAA;
	}
	else
	{
		string eMsg = "The given index (";
		eMsg += toString(AAIndx);
		eMsg += ") is out of range..or the portion is empty (No.AAs= ";
		eMsg += toString(numOfAA());
		eMsg += ")";

		errMsg("Protein", "reOrderAtoms", eMsg);

		exit(1);

	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
bool Protein::isHeavyAtom(char type)
{
	if ((type == 'C') || (type == 'N') || (type == 'O') || (type == 'S'))
		return true;

	return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Protein::sortHlces(vector<HelicesSecondaryStructure> & a, const unsigned int leftArg, unsigned int rightArg)
{
  if (leftArg < rightArg)
  {

    int pivotvalue = a[leftArg].startIndx;
    int left = leftArg - 1;
    int right = rightArg + 1;

  for(;;)
  {

    while (a[--right].startIndx > pivotvalue);
    while (a[++left].startIndx < pivotvalue);

    if (left >= right) break;

    HelicesSecondaryStructure temp = a[right];
    a[right] = a[left];
    a[left] = temp;
  }

  int pivot = right;
  sortHlces(a, leftArg, pivot);
  sortHlces(a, pivot + 1, rightArg);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// END OF PROTEIN CLASS IMPLEMENTATION ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif
