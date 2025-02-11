#include "XYFile.h"

#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;
//----------------------------------------------------------------------------
CXYFile::CXYFile(void)
{
}
//----------------------------------------------------------------------------
CXYFile::~CXYFile(void)
{
}
//----------------------------------------------------------------------------
CXYMatrixf CXYFile::ReadMatrix(const char* acfname)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
             << endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);

	vector<vector<string>*>::iterator p = data.begin();	
	int iRow = data.size();
	int iCol = (*p)->size();
	CXYMatrixf kM(iRow,iCol);
	for (int i = 0; i < iRow; i++)
	{
		for (int j = 0; j < iCol; j++)
		{
			kM[i][j] = (float)atof(data.at(i)->at(j).c_str());
		}
	}

	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;                                  // Be sure to delete
	} 

	return kM;
}
//----------------------------------------------------------------------------
CXYVector<float> CXYFile::ReadVector(const char* acfname)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
             << endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);

	vector<vector<string>*>::iterator p = data.begin();	
	int iRow = data.size();
	CXYVector<float> kV(iRow);
	for (int i = 0; i < iRow; i++)
	{
		kV[i] = (float)atof((data.at(i))->at(0).c_str());
	}

	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;                                  // Be sure to delete
	} 

	return kV;
}
//----------------------------------------------------------------------------
CXYVector<int> CXYFile::ReadVectorInt(const char* acfname)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
		<< endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);
	
	vector<vector<string>*>::iterator p = data.begin();	
	int iRow = data.size();
	CXYVector<int> kV(iRow);
	for (int i = 0; i < iRow; i++)
	{
		kV[i] = atoi((data.at(i))->at(0).c_str());
	}
	
	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;                                  // Be sure to delete
	} 
	
	return kV;
}
//----------------------------------------------------------------------------
void CXYFile::WriteMatrix(const char* acfname, const char* _mode, const CXYMatrixf& rkM)
{
	FILE * pFile;

	int iRow, iCol;
	rkM.GetSize(iRow,iCol);
	pFile = fopen (acfname,_mode);
	for (int i=0 ; i<iRow; i++)
	{
		for (int j=0; j<iCol; j++)
		{
			fprintf (pFile, "% .3E\t",rkM[i][j]);
		}
		fprintf(pFile, "\n");
	}
	fclose (pFile);
}
//----------------------------------------------------------------------------
void CXYFile::WriteVector(const char* acfname, const char* _mode, const CXYVectorf& rkV)
{
	FILE * pFile;

	pFile = fopen (acfname,_mode);
	for (int i=0 ; i<rkV.GetSize(); i++)
	{
		fprintf (pFile, "% .3E\n",rkV[i]);
	}
	fclose (pFile);

}
//----------------------------------------------------------------------------
void CXYFile::ReadCSV(istream& in, vector<vector<string>*>& data) {

	vector<string>* p = NULL;
	string tmp;

	while (!in.eof( )) {
		getline(in, tmp, '\n');                     // Grab the next line
    //the following line trims white space from the beginning of the string
    tmp.erase(tmp.begin(), find_if(tmp.begin(), tmp.end(), not1(ptr_fun<int, int>(isspace)))); 
    if (tmp.find("#") == 0) {                   // Skip "#" comment line
      continue;
    }
    
		CXYStringTokenizer st(tmp);
    
		p = new vector<string>( );
		while (st.hasMoreTokens( )) {
			st.nextToken(tmp);
			p->push_back(tmp);
		}
		if (p->size() > 0)	// skip empty line
		{
			data.push_back(p);
		}
	}
}
//----------------------------------------------------------------------------
/*
CXYPDB CXYFile::ReadAtomsFromPDB(const char* acfname)
{
	CXYPDB kP;
	filebuf BuffFile;
	BuffFile.open(acfname,ios::in);
	istream f(&BuffFile);

	if (! f)
	{
		cerr << "can't open input file \"" << acfname << "\""
             << endl;
		exit(EXIT_FAILURE);
	}

	f.seekg(0);
	char acBuffLine[81];
	char acRecName[] ="ATOM  ";
	
	while (!f.eof())
	{
		f.getline(acBuffLine,81);
		if (! memcmp(acRecName,acBuffLine,6))
		{
			CXYPDBAtom Atom;
			Atom.ReadLine(acBuffLine);
			(kP.GetPDBAtoms())->push_back(Atom);
		}
	}
	return kP;
}
*/
//----------------------------------------------------------------------------
/*
void CXYFile::WriteAtomsToPDB(const char* acfname, ios_base::openmode _mode, CXYPDB& rkP)
{
	vector<CXYPDBAtom >::iterator p;
	CXYPDBAtom Atom;

	ofstream OutF;
	OutF.open(acfname, _mode | ios::out);
	if (! OutF)
	{
		cerr << "can't open output file \"" << acfname << "\""
             << endl;
		exit(EXIT_FAILURE);
	}

	for (p = rkP.GetPDBAtoms()->begin(); p != rkP.GetPDBAtoms()->end(); p++)
	{
		Atom = *p;

		OutF.setf(ios_base::fixed);
		OutF.width(6);	OutF << left  << Atom.GetRecName() ;

		OutF.width(5);	OutF << right << Atom.GetSerial();
		
		OutF.width(1);	OutF << "";

		OutF.width(2);	OutF << right << Atom.GetAtomicSym();
		OutF.width(1);	OutF << Atom.GetRemoteInd();
		OutF.width(1);	OutF << Atom.GetBranchDes();
		
		OutF.width(1);	OutF << Atom.GetAltLoc();
		OutF.width(3);	OutF << Atom.GetResName();
		OutF.width(1);	OutF << "";
		OutF.width(1);	OutF << Atom.GetChainID();
		OutF.width(4);	OutF << right << Atom.GetResSeq();
		OutF.width(1);	OutF << Atom.GetICode();
		OutF.width(3);	OutF << "";
		OutF.precision(3);
		OutF.width(8);	OutF << Atom.GetX();
		OutF.width(8);	OutF << Atom.GetY();
		OutF.width(8);	OutF << Atom.GetZ();
		OutF.precision(2);
		OutF.width(6);	OutF << Atom.GetOccupancy();
		OutF.width(6);	OutF << Atom.GetTempFactor();
		OutF.width(6);	OutF << "";
		OutF.width(4);	OutF << left  << Atom.GetSegID();
		OutF.width(2);	OutF << right << Atom.GetElement();
		OutF.width(2);	OutF << right << Atom.GetCharge();

		OutF << endl;
	}

	OutF.close();

}*/
void CXYFile::WriteAtomsToPDB(const char* acfname, const char* _mode, CXYPDB& rkP)
{
	vector<CXYPDBAtom >::iterator p;
	CXYPDBAtom Atom;

	FILE * pFile;

	pFile = fopen (acfname, _mode);
	for (p = rkP.GetPDBAtoms()->begin(); p != rkP.GetPDBAtoms()->end(); p++)
	{
		Atom = *p;
		fprintf (
			pFile, 
			"%-6s%5d%1s%2s%1s%1s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n",
			Atom.GetRecName(),
			Atom.GetSerial(),
			" ",
			Atom.GetAtomicSym(),
			Atom.GetRemoteInd(),
			Atom.GetBranchDes(),
			Atom.GetAltLoc(),
			Atom.GetResName(),
			" ",
			Atom.GetChainID(),
			Atom.GetResSeq(),
			Atom.GetICode(),
			" ",
			Atom.GetX(),
			Atom.GetY(),
			Atom.GetZ(),
			Atom.GetOccupancy(),
			Atom.GetTempFactor(),
			" ",
			Atom.GetSegID(),
			Atom.GetElement(),
			Atom.GetCharge());
	}
	fclose (pFile);
}

void CXYFile::WriteAtomsToPDB( const char* acfname, const char* _mode, vector<CXYPDB>& rkvPDBs )
{
	FILE * pFile;
	pFile = fopen (acfname, _mode);
	CXYPDBAtom Atom;

	int iCount = rkvPDBs.size();
	for (int i=0; i<iCount; i++)
	{
		fprintf(pFile,"MODEL %4d\n",i);
		CXYPDB kPDB = rkvPDBs.at(i);
		for (vector<CXYPDBAtom>::iterator p = kPDB.GetPDBAtoms()->begin(); p != kPDB.GetPDBAtoms()->end(); p++)
		{
			Atom = *p;
			fprintf (
							 pFile, 
							 "%-6s%5d%1s%2s%1s%1s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n",
							 Atom.GetRecName(),
							 Atom.GetSerial(),
							 " ",
							 Atom.GetAtomicSym(),
							 Atom.GetRemoteInd(),
							 Atom.GetBranchDes(),
							 Atom.GetAltLoc(),
							 Atom.GetResName(),
							 " ",
							 Atom.GetChainID(),
							 Atom.GetResSeq(),
							 Atom.GetICode(),
							 " ",
							 Atom.GetX(),
							 Atom.GetY(),
							 Atom.GetZ(),
							 Atom.GetOccupancy(),
							 Atom.GetTempFactor(),
							 " ",
							 Atom.GetSegID(),
							 Atom.GetElement(),
							 Atom.GetCharge());
		}
		fprintf(pFile,"TER   %5d      %3s %1s%4d%1s\n",
						Atom.GetSerial()+1,
						Atom.GetResName(), 
						Atom.GetChainID(), 
						Atom.GetResSeq(), 
						Atom.GetICode());
		fprintf(pFile,"ENDMDL\n");
	}
	for (int i=iCount-1; i>=0; i--)
	{
		fprintf(pFile,"MODEL %4d\n",i);
		CXYPDB kPDB = rkvPDBs.at(i);
		for (vector<CXYPDBAtom>::iterator p = kPDB.GetPDBAtoms()->begin(); p != kPDB.GetPDBAtoms()->end(); p++)
		{
			Atom = *p;
			fprintf (
							 pFile, 
							 "%-6s%5d%1s%2s%1s%1s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s\n",
							 Atom.GetRecName(),
							 Atom.GetSerial(),
							 " ",
							 Atom.GetAtomicSym(),
							 Atom.GetRemoteInd(),
							 Atom.GetBranchDes(),
							 Atom.GetAltLoc(),
							 Atom.GetResName(),
							 " ",
							 Atom.GetChainID(),
							 Atom.GetResSeq(),
							 Atom.GetICode(),
							 " ",
							 Atom.GetX(),
							 Atom.GetY(),
							 Atom.GetZ(),
							 Atom.GetOccupancy(),
							 Atom.GetTempFactor(),
							 " ",
							 Atom.GetSegID(),
							 Atom.GetElement(),
							 Atom.GetCharge());
		}
		fprintf(pFile,"TER   %5d      %3s %1s%4d%1s\n",
						Atom.GetSerial()+1,
						Atom.GetResName(), 
						Atom.GetChainID(), 
						Atom.GetResSeq(), 
						Atom.GetICode());
		fprintf(pFile,"ENDMDL\n");
	}
	fclose (pFile);
}
//----------------------------------------------------------------------------
void CXYFile::WriteComment(const char* acfname, const char* _mode, const char* acComment)
{
	FILE * pFile;

	pFile = fopen (acfname, _mode);
	fprintf (pFile, "%s",acComment);
	fclose (pFile);

}
//----------------------------------------------------------------------------
void CXYFile::WriteVector(const char* acfname, const char* _mode, const CXYVector<CXYPoint3D<float> >& rkV)
{
   FILE * pFile;

   pFile = fopen (acfname, _mode);
   for (int i=0 ; i<rkV.GetSize(); i++)
   {
     fprintf (pFile, "% 4.4f\t\t% 4.4f\t\t% 4.4f\n",rkV[i].X(),rkV[i].Y(),rkV[i].Z());
   }
   fclose (pFile);

}
//----------------------------------------------------------------------------
CXYPDB CXYFile::ReadAtomsFromPDB(const char* acfname)
{

	FILE * pFile;
	char buff [81];
	char acRecName[] ="ATOM  ";
	CXYPDB kP;

	pFile = fopen (acfname , "r");
	if (pFile == NULL) perror ("Error opening file");
	else {
		while(fgets(buff,81,pFile))
		{
			if (! memcmp(acRecName,buff,6))
			{
				CXYPDBAtom Atom;
				Atom.ReadLine(buff);
				(kP.GetPDBAtoms())->push_back(Atom);
			}
		}
		fclose (pFile);
	}
	return kP;
}
//----------------------------------------------------------------------------
//void CXYFile::MakeDirectory(char* sPath, int imode)
//{
//#ifdef WIN32
//	//		#include <direct.h>
//	_mkdir(sPath);
//#else 
//	//		#include <sys/stat.h>
//	mkdir (sPath, imode);
//#endif
//}

//----------------------------------------------------------------------------
int CXYFile::MakeDirectory(char* sPath, int imode)
{
#ifdef WIN32
	//		#include <direct.h>
	return(_mkdir(sPath));
#else 
	//		#include <sys/stat.h>
  //		#include <errno.h>
  struct stat st;
  if(stat(sPath,&st) != 0){  // check if the directory exist
    if (mkdir ( sPath, imode) == -1){
      cerr << "Error: " << strerror(errno) << endl;
      return (EXIT_FAILURE);
    }
  }
	return 0;
#endif
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
CXYMatrixf CXYFile::ReadMatrix(const char* acfname, int iNumRowSkip, int iNumColSkip)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
		<< endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);
	
	vector<vector<string>*>::iterator p = data.begin() + iNumRowSkip;	
	int iRow = data.size() - iNumRowSkip;
	int iCol = (*p)->size() - iNumColSkip;
	CXYMatrixf kM(iRow,iCol);
	for (int i = 0; i < iRow; i++)
	{
		for (int j = 0; j < iCol; j++)
		{
			kM[i][j] = (float)atof(data.at(i+iNumRowSkip)->at(j+iNumColSkip).c_str());
		}
	}
	
	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;                                  // Be sure to delete
	} 
	
	return kM;
}
//----------------------------------------------------------------------------
bool CXYFile::FileExists(const char* cFileName){
	FILE* fp = NULL;
	
	//will not work if you do not have read permissions
	//to the file, but if you don't have read, it
	//may as well not exist to begin with.
	
	fp = fopen( cFileName, "rb" );
	if( fp != NULL )
	{
		fclose( fp );
		return true;
	}
	
	return false;
}
//----------------------------------------------------------------------------
bool CXYFile::DirExists(const char* cDirName){
#ifdef WIN32

	struct _stat statBuffer;
	return (_tstat(szPath, &statBuffer) >= 0 &&	 // make sure it exists
					statBuffer.st_mode & S_IFDIR);	 // and it's not a file
#else
  struct stat st;
	return (stat(cDirName,&st) == 0) ? true: false;
#endif
}


//----------------------------------------------------------------------------
// Sparse Matrix format
//----------------------------------------------------------------------------
// RowSize ColSize
// RowIndex_i ColIndex_j Element_ij
//----------------------------------------------------------------------------
CXYMatrixf CXYFile::ReadSparseMatrix(const char* acfname)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
		<< endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);
	
	vector<vector<string>*>::iterator p = data.begin();
	
	int iNumRow = atoi(data.at(0)->at(0).c_str());
	int iNumCol = atoi(data.at(0)->at(1).c_str());
	
	CXYMatrixf kM(iNumRow, iNumCol);
	for (unsigned int i=1; i < data.size(); i++) {
		int ind_i = atoi(data.at(i)->at(0).c_str());
		int ind_j = atoi(data.at(i)->at(1).c_str());
		kM[ind_i][ind_j] = (float) atof(data.at(i)->at(2).c_str());
	}
	
	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;     // Be sure to delete
	} 

	return kM;
}

//----------------------------------------------------------------------------
void CXYFile::WriteSparseMatrix(const char* acfname, const char* _mode, const CXYMatrixf& rkM)
{
	FILE * pFile;
	
	int iRow, iCol;
	rkM.GetSize(iRow,iCol);
	pFile = fopen (acfname,_mode);
	fprintf(pFile, "%d\t%d\n",iRow,iCol);
	for (int i=0 ; i<iRow; i++)
	{
		for (int j=0; j<iCol; j++)
		{
			if(rkM[i][j] > CXYMath<float>::ZERO_TOLERANCE)
			{
				fprintf (pFile, "%d\t%d\t%.3f\n",i,j,rkM[i][j]);
			}
		}
	}
	fclose (pFile);
}

//----------------------------------------------------------------------------
// RowSize ColSize
// RowIndex_i ColIndex_j Element_ij
//----------------------------------------------------------------------------

vector<CXYTriple<int, int, float> > CXYFile::ReadSparseMatrixToTriple(const char* acfname)
{
	ifstream inF(acfname);
	if (!inF)
	{
		cerr << "can't open output file \"" << acfname << "\""
		<< endl;
		exit(EXIT_FAILURE);
	}
	vector<vector<string> * > data;
	ReadCSV(inF, data);
	
	vector<vector<string>*>::iterator p = data.begin();

	vector<CXYTriple<int,int,float> > kTriple;
	for (unsigned int i=1; i < data.size(); i++) {
		int ind_i = atoi(data.at(i)->at(0).c_str());
		int ind_j = atoi(data.at(i)->at(1).c_str());
		float value = (float) atof(data.at(i)->at(2).c_str());
		kTriple.push_back(CXYTriple<int, int, float>(ind_i, ind_j, value));
	}
	
	for (p = data.begin( ); p != data.end( ); ++p) 
	{
		delete *p;     // Be sure to delete
	} 
	
	return kTriple;
}

