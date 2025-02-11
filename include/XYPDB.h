#ifndef XYPDB_H
#define	XYPDB_H

#include <vector>
#include <algorithm>

#include "XYPDBAtom.h"
#include "XYPoint3D.h"

using namespace std;
//class CXYPDBAtom{};
class CXYPDB 
{
public:
	CXYPDB(void);
	~CXYPDB(void);

	const CXYPDBAtom* operator[] (int iInd) const;
	CXYPDBAtom* operator[] (int iInd);

	vector<CXYPDBAtom >* GetPDBAtoms();
	const vector<CXYPDBAtom >* GetPDBAtoms() const;
	int Size() const;
	//void ReadAtomsFromPDB(const char* fname);
	//void WriteAtomsToPDB(const char* fname);

	// assignment
	CXYPDB& operator= (const CXYPDB& rkPDB);
	
	// Side Chain of CA Model
	void ReduceToSCOfCA(void);
	void PostDealPDB(char* sLevel);
    // this function need ResSeq is uniq
    int GetPDBIndex(int iResSeq);
private:
	vector<CXYPDBAtom> m_kPDBAtoms;

};

#endif
