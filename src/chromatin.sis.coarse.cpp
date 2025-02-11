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
#include "XYSIS.h"
#include "XYSO3Sequence.h"

void usage(){
//	printf ("Usage: ensemble -configfile ConfigFile -outpts ptsfile -outlen lenfile\n");
	printf ("Usage: chromatin.sis.coarse -conf ConfigFile -prefix AA\n");
	printf ("       -conf          ConfigFile\n");
	printf ("       -prefix        PrefixLetters(options, default AA)\n");
	exit(0);
}

int main (int argc, char *argv[])
{	
	const size_t NEWSIZE=1024;
	char cConfigFile[NEWSIZE] = "";
	
	// the position file including start and end for each node
	char cStartEndFile[NEWSIZE] = "";
	char cPrefix[NEWSIZE] = "";

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
		cout << "ERROR TYPE 1" << endl;
		exit(0);
	}

	// read configuration file
	ConfigFile config(cConfigFile);
//	float fBindingAngleBeg, fBindingAngleEnd, fTorsionAngleBeg, fTorsionAngleEnd ;
//	int iNumBindingAngle, iNumTorsionAngle;
	float fPersistenLength, fCollisionLength;
	float fPackingDensity;
	float fNucleusSphereLength;
	//		char cSegLenFile[NEWSIZE];
	char cOutPath[NEWSIZE];
	int iNumNodes;
	int iNumSamplePoints;
	float fthreshold_A, fthreshold_B;
	string sSegLenFile, sOutPath;
	string sStartEndFile;
	int m_interaction_distance;
	float fmc_coeff;

	config.readInto(sOutPath, "out_path"); //strcpy(cOutPath, sOutPath.c_str());
	sprintf(cOutPath, "%s/%s",sOutPath.c_str(),cPrefix);

	config.readInto( fPersistenLength, "persistence_length", 2110.0f);
	config.readInto( fCollisionLength, "collision_length", 300.0f);
	config.readInto( fPackingDensity, "packing_density", 0.07f);
	config.readInto( fNucleusSphereLength, "nucleus_sphere_diameter", 60000.0f);
	config.readInto( fthreshold_A, "threshold_A", 0.4f); //collision tolerance
	config.readInto( fthreshold_B, "threshold_B", 0.3f); //wrong interactions tolerance
	config.readInto( m_interaction_distance, "interaction_distance", 800);	
	config.readInto( fmc_coeff, "monte_carlo_coefficient", 0.2f);
	config.readInto( iNumNodes, "number_nodes", 100);
	config.readInto( iNumSamplePoints, "number_sample_points", 72);
	
	config.readInto( sStartEndFile, "start_end_file");
	strcpy(cStartEndFile, sStartEndFile.c_str());


	// for sis 
	int iMmax;
	float fRho_1, fRho_2, fRho_3, fTau_t, fAdjust;
//	char cDistFile[NEWSIZE];
//	string sDistFile ;
	config.readInto( iMmax, "M_max", 10);
	
	config.readInto( fRho_1, "rho_1", 0.1f);
	config.readInto( fRho_2, "rho_2", 1.f);
	config.readInto( fRho_3, "rho_3", 0.1f);
	config.readInto( fTau_t, "tau_t", 1.f);
	config.readInto( fAdjust, "adjust", 1e-6f);
//	config.readInto( sDistFile,  "dist_file"); strcpy(cDistFile, sDistFile.c_str());
  
    float fCutoff;
    config.readInto( fCutoff, "cutoff", 75.0f);

    char cPvalFile[NEWSIZE]; string sPvalFile;
    config.readInto( sPvalFile, "pval_file"); 
    strcpy(cPvalFile, sPvalFile.c_str());

    char cContIndFile[NEWSIZE]; string sContIndFile;  
    config.readInto( sContIndFile,  "cont_index_file");
    strcpy(cContIndFile, sContIndFile.c_str());


	CXYSIS SIS(
        cOutPath,
        fPersistenLength,
        fCollisionLength,
        fCutoff,
        fPackingDensity,
        fNucleusSphereLength,
        iNumNodes,
        iNumSamplePoints,
        cStartEndFile,
        iMmax,
        fRho_1,
        fRho_2,
        fRho_3,
        fTau_t,
        fAdjust,
        cPvalFile,
        cContIndFile,
	fthreshold_A,
	fthreshold_B,
	m_interaction_distance,
	fmc_coeff
        );

//	SIS.SIS_Algorithm();
	SIS.GrowChainSegment_SISAlgorithm_2();
//	SIS.GrowChainSegment_SISAlgorithm_2(1);
//	SIS.GrowChainSegment_SISAlgorithm_2(-1);
		

	
  return 0;
}	

