#ifndef XYFILE_H
#define XYFILE_H
#include <string>
#include <iostream>
#include <vector>
#include "XYMatrix.h"
#include "XYVector.h"
#include "XYPoint3D.h"
#include "XYPDB.h"
#include "XYUtility.h"

#ifdef WIN32
#include <direct.h>
#else 
#include <sys/stat.h>
#include <errno.h>
#endif


#pragma warning (disable: 4267)
#pragma warning (disable: 4996)

using namespace std;

class CXYFile
{
public:
	CXYFile(void);
	~CXYFile(void);

	static void ReadCSV(istream& in, vector<vector<string>*>& data);
	static CXYMatrixf ReadMatrix(const char* acfname);
	static CXYMatrixf ReadMatrix(const char* acfname, int iNumRowSkip, int iNumColSkip);
	static CXYMatrixf ReadSparseMatrix(const char* acfname);
	static void WriteMatrix(const char* acfname, const char* _mode, const CXYMatrixf& rkM);
	static void WriteSparseMatrix(const char* acfname, const char* _mode, const CXYMatrixf& rkM);
  vector<CXYTriple<int, int, float> > ReadSparseMatrixToTriple(const char* acfname);


  
	static CXYVector<float> ReadVector(const char* acfname);
	static CXYVector<int> ReadVectorInt(const char* acfname);
	static void WriteVector(const char* acfname, const char* _mode, const CXYVectorf& rkV);
	static void WriteVector(const char* acfname, const char* _mode, const CXYVector<CXYPoint3D<float> >& rkV);

	static CXYPDB ReadAtomsFromPDB(const char* acfname);
	//static void WriteAtomsToPDB(const char* acfname, ios_base::openmode _mode, CXYPDB& rkP);
	static void WriteAtomsToPDB(const char* acfname, const char* _mode, CXYPDB& rkP);
	static void WriteAtomsToPDB(const char* acfname, const char* _mode, vector<CXYPDB>& rkvPDBs);

	static void WriteComment(const char* acfname, const char* _mode, const char* acComment);
	static int MakeDirectory(char* sPath, int imode);
	static bool FileExists(const char* cFileName);
	static bool DirExists(const char* cDirName);
};

// String tokenizer class.
class CXYStringTokenizer {

public:

   CXYStringTokenizer(const string& s, const char* delim = NULL) :
      str_(s), count_(-1), begin_(0), end_(0) {

      if (!delim)
         delim_ = " \f\n\r\t\v";  //default to whitespace
      else
         delim_ = delim;

      // Point to the first token
      begin_ = str_.find_first_not_of(delim_);
      end_ = str_.find_first_of(delim_, begin_);
   }

   size_t countTokens( ) {
     if (count_ >= 0) // return if we've already counted
       return(count_);

     string::size_type n = 0;
     string::size_type i = 0;

     for (;;) {
        // advance to the first token
        if ((i = str_.find_first_not_of(delim_, i)) == string::npos)
           break;
        // advance to the next delimiter
        i = str_.find_first_of(delim_, i+1);
        n++;
        if (i == string::npos)
          break;
     }
     return (count_ = n);
   }
   bool hasMoreTokens( ) {return(begin_ != end_);}
   void nextToken(string& s) {
     if ( begin_ != (int) string::npos &&  end_ != (int) string::npos) {
        s = str_.substr(begin_, end_-begin_);
        begin_ = str_.find_first_not_of(delim_, end_);
        end_ = str_.find_first_of(delim_, begin_);
     }
     else if (begin_ != (int) string::npos && end_ == (int) string::npos)
     {
        s = str_.substr(begin_, str_.length( )-begin_);
        begin_ = str_.find_first_not_of(delim_, end_);
     }

   }

private:
   CXYStringTokenizer( ) {};
   string delim_;
   string str_;
   int count_;
   int begin_;
   int end_;
};

#endif
