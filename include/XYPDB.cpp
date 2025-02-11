#include "XYPDB.h"

#include <iostream>
#include <fstream>
#include <cstring>

//#include "XYPDBAtom.h"
//----------------------------------------------------------------------------
CXYPDB::CXYPDB(void)
{

}
//----------------------------------------------------------------------------
CXYPDB::~CXYPDB(void)
{
	
}
//----------------------------------------------------------------------------
int CXYPDB::Size() const
{
	return (int)m_kPDBAtoms.size();
}
//----------------------------------------------------------------------------
void CXYPDB::ReduceToSCOfCA(void)
{
	vector<CXYPDBAtom >::iterator p;

	vector<CXYPDBAtom > lsAtoms;
	for(p = m_kPDBAtoms.begin(); p != m_kPDBAtoms.end() ; p++)
	{
		CXYPDBAtom Atom = (*p);
		if (memcmp ( Atom.GetName(), "CA", 2*sizeof(char)) == 0)
		{
			lsAtoms.push_back(Atom);
		}
	}
	m_kPDBAtoms = lsAtoms;
}
//----------------------------------------------------------------------------
const CXYPDBAtom* CXYPDB::operator[] (int iInd) const
{
	assert(0 <= iInd && iInd < Size());
	return &m_kPDBAtoms[iInd];
}
//----------------------------------------------------------------------------
CXYPDBAtom* CXYPDB::operator[] (int iInd)
{
	assert(0 <= iInd && iInd < Size());
	return &m_kPDBAtoms[iInd];
}
//----------------------------------------------------------------------------
vector<CXYPDBAtom >* CXYPDB::GetPDBAtoms()
{
	return &m_kPDBAtoms;
}
//----------------------------------------------------------------------------
const vector<CXYPDBAtom >* CXYPDB::GetPDBAtoms() const
{
	return &m_kPDBAtoms;
}
//----------------------------------------------------------------------------
CXYPDB& CXYPDB::operator= (const CXYPDB& rkPDB)
{
	m_kPDBAtoms.clear();

	//vector<CXYPDBAtom >::iterator p;
	//vector<CXYPDBAtom >* kvPDBAtoms = rkPDB.GetPDBAtoms();

	//vector<CXYPDBAtom > lsAtoms;
	//for(p = (*kvPDBAtoms).begin(); p != (*kvPDBAtoms).end(); p++)
	//{
	//	m_kPDBAtoms.push_back(*p);
	//}
	for (int i=0; i<rkPDB.Size(); i++)
	{
		m_kPDBAtoms.push_back((*rkPDB[i]));
	}
	return *this;
}
//----------------------------------------------------------------------------
int CXYPDB::GetPDBIndex(int iResSeq)
{
  int iPDBInd = 0;
  for (int i=0; i<Size(); i++)
  {
    if ( m_kPDBAtoms[i].GetResSeq() == iResSeq )
    {
      iPDBInd = i;
      break;
    }
  }
  return iPDBInd;
}
//----------------------------------------------------------------------------
// deal with different level, e.g. residue level or atom level
// default is residue level
void CXYPDB::PostDealPDB(char* sLevel)
{
	int level_flg;
	if (strcmp(sLevel,"residue") == 0)
	{
		level_flg = 1;
	}
	if (strcmp(sLevel,"atom") == 0)
	{
		level_flg = 2;
	}

	switch (level_flg)
	{
	case 1:
		// reduce to CA atom for each residue
		ReduceToSCOfCA();
		break;
	case 2:
		//do nothing
		// need clean the pdb file, do not include H atom
		break;
	default:
		// default is CA atom
		ReduceToSCOfCA();
		break;
	}
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
