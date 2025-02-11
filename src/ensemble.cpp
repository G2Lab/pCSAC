/*
   Generate chromatin ensemble of reference state.
   Author: Yun Xu
   Email: yxu7@uic.edu
   Data: Sep 20, 2011
*/
/*
   We treat the chromosome structure as a polymer chain.
   1. Bead-stick model consisting of n nodes.
   2. Each node is represented as a sphere of diameter 30 nm to take into 
      account the exclude volume effect.
   3. The sphere is connected by rigid bonds of the persistence length 
      211 nm.
   4. The state space is discrete. uniformely distributed on the sphere.

 The reference state models the generic interaction that expected from 
   random polymer chain and not due to biological interactions for studying
   the chromatin structure. For our study, we randomly generate chromatin 
   fiber ensemble constrained in the nucleus sphere. Our model is obtained 
   by:
   1. Grow polymer chain that models the chromatin fiber in the sphere of 
      diameter 6 micronmeter that models the nucleus.
   2. The location where the chromatin fiber start at random location inside
      the sphere.
   3. Growing a self-avoid chain to its full length with varying bending and
      torsion angles.
   Program unit is amstrong, 10e-10m
*/

#include <iostream>
#include <ctime>
#include "XYMath.h"
#include "XYFile.h"
#include "XYUtility.h"
#include "tree.hh"
#include "ConfigFile.h"
#include "XYEnsemble.h"
#include "XYSO3Sequence.h"

void usage(){
	printf ("Usage: ensemble -conf ConfigFile -prefix AA -samplesize 1 \n");
	printf ("       -conf          ConfigFile\n");
	printf ("       -prefix        PrefixLetters(options, default AA)\n");
	printf ("       -samplesize    Size of ensemble\n");
	exit(0);
}

int main (int argc, char *argv[])
{
	const size_t NEWSIZE=1024;
	char cConfigFile[NEWSIZE] = "";
//	char cOutptsFile[NEWSIZE] = "";
//	char cOutlenFile[NEWSIZE] = "";
	
	// the position file including start and end for each node
	char cStartEndFile[NEWSIZE] = "";
	char cPrefix[NEWSIZE] = "";

	int iSampleSize = 1;
	/* read argv
	* http://stackoverflow.com/questions/441547/most-efficient-way-to-process-arguments-from-the-command-line-in-c  
	* Loop over command-line args
	* (Actually I usually use an ordinary integer loop variable and compare
	* args[i] instead of *i -- don't tell anyone! ;)
	*/ 

	vector<string> args(argv + 1, argv + argc);
	for (vector<string>::iterator i = args.begin(); i != args.end(); ++i) {
		if (*i == "-h" || *i == "--help" ) {
			usage();
			return 0;
		} else if (*i == "-conf") {
			strcpy(cConfigFile, (*++i).c_str());
			cout << "Reading ConfigFile: " << cConfigFile << endl ;
		} else if (*i == "-prefix"){
			strcpy(cPrefix, (*++i).c_str());
		} else if (*i == "-samplesize"){
			iSampleSize = atoi((*++i).c_str());
		}
	}

	if (strcmp(cConfigFile, "") == 0) {
		usage();
	}
	
	if (strcmp(cPrefix ,"") == 0){
		strcpy(cPrefix, "AA");
	}

	CXYFile kFile;
	if (!kFile.FileExists(cConfigFile) ) {
		cout << "File " << cConfigFile << " do not exist" << endl ;
		exit(0);
	}

	// read configuration file
	ConfigFile config(cConfigFile);
//	float fBindingAngleBeg, fBindingAngleEnd, fTorsionAngleBeg, fTorsionAngleEnd ;
//	int iNumBindingAngle, iNumTorsionAngle;
	float fPersistenLength, fCollisionLength;
	float fPackingDensity;
	float fNucleusSphereLength;
	char cOutPath[NEWSIZE];
	int iNumNodes;
	int iNumSamplePoints;
	string sSegLenFile, sOutPath;
	char cSegLenFile[NEWSIZE];
	string sStartEndFile;
	char cContIndFile[NEWSIZE];
	string sContIndFile;

	config.readInto(sOutPath, "out_path"); strcpy(cOutPath, sOutPath.c_str());
//	if (!strcmp(cOutptsFile, "")) {
//		cout << "Out point file: " << cOutPath << "/" << cOutptsFile << endl;
//	}
//	if (!strcmp(cOutlenFile, "")) {
//		cout << "Out length file: " << cOutPath << "/" << cOutlenFile << endl;
//	}


	config.readInto( fPersistenLength, "persistence_length", 2110.0f);
	config.readInto( fCollisionLength, "collision_length", 300.0f);
	config.readInto( fPackingDensity, "packing_density", 0.07f);
	config.readInto( fNucleusSphereLength, "nucleus_sphere_diameter", 60000.0f);
	
	config.readInto( iNumNodes, "number_nodes", 100);
	config.readInto( iNumSamplePoints, "number_sample_points", 72);
	
	config.readInto( sStartEndFile, "start_end_file");
	strcpy(cStartEndFile, sStartEndFile.c_str());

	config.readInto( sSegLenFile, "seg_len_file");
	strcpy(cSegLenFile, sSegLenFile.c_str());
  
	config.readInto( sContIndFile,  "cont_index_file"); 
	strcpy(cContIndFile, sContIndFile.c_str());

	
	CXYEnsemble ensemble(
		cOutPath,
		fPersistenLength,
		fCollisionLength,
		fPackingDensity,
		fNucleusSphereLength,
		iNumNodes,
		iNumSamplePoints,
		cStartEndFile,
    cSegLenFile,
		cContIndFile);

	char buff[NEWSIZE];
	char cFN[NEWSIZE];
	
	
	int i = 0;
  int iCount = 0;
	while (i<iSampleSize) {
		sprintf(buff, "%s_%04d", cPrefix, i);
		cout<< buff << endl;
		if (ensemble.GrowOneChain()) {
			sprintf(cFN,"%s.pts",buff);
			ensemble.WriteChain(cFN);
      sprintf(cFN, "%s.dst",buff);
      ensemble.WriteContactDistance(cFN);
//			sprintf(cFN,"%s.wgt",buff);
//      ensemble.WriteWeight(cFN);
      i++;
      iCount = 0;
		} else {
      iCount ++;
    }
    
    if (iCount >100) {
      cout << "Can not finish generating!\n";
      break;
    }


	}

	
//	CXYSO3Sequence sO3sequence(iNumSamplePoints);
//	sO3sequence.SetSO3Sequence();
//	sO3sequence.printSO3Sequence();
	
  return 0;
}

