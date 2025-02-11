#include "XYPDBAtom.h"

#include <iostream>
#include <fstream>
#include <cstring>

#include "XYUtility.h"

using namespace std;
//----------------------------------------------------------------------------
CXYPDBAtom::CXYPDBAtom(void):
	m_iSerial(0),
	m_iResSeq(0),
	m_fX(0.0), m_fY(0.0), m_fZ(0.0),
	m_fOccupancy(0.0),
	m_fTempFactor(0.0),
	m_point(0.0,0.0,0.0)
{
	Initialize();
}

//----------------------------------------------------------------------------
CXYPDBAtom::~CXYPDBAtom(void)
{
}

//----------------------------------------------------------------------------
void CXYPDBAtom::Initialize(void)
{
	memset(m_acRecName,'\0',7*sizeof(char));
		
	memset(m_acAtom,   '\0',5*sizeof(char));
	memset(m_acAtomicSym, '\0', 3*sizeof(char));
	memset(m_acRemoteInd, '\0', 2*sizeof(char));
	memset(m_acBranchDes, '\0', 2*sizeof(char));

	memset(m_acAltLoc, '\0',  2*sizeof(char));
	memset(m_acResName,'\0',4*sizeof(char));
	memset(m_acChainID,'\0',2*sizeof(char));
	memset(m_acICode,  '\0',2*sizeof(char));
	memset(m_acSegID,  '\0',5*sizeof(char));
	memset(m_acElement,'\0',3*sizeof(char));
	memset(m_acCharge ,'\0',3*sizeof(char));
}

//----------------------------------------------------------------------------
/*
void CXYPDBAtom::ReadLine(const char *acBuffLine)
{
	char acSerial[6];
	char acResSeq[5];
	char acReal8[9];
	char acReal6[7];


	memset(acSerial, '\0', 6*sizeof(char));
	memset(acResSeq, '\0', 5*sizeof(char));
	memset(acReal8,	'\0', 9*sizeof(char));
	memset(acReal6, '\0', 7*sizeof(char));

	memcpy(m_acRecName,acBuffLine,6);			CXYUtility::TrimChar(m_acRecName);
		
	memcpy( acSerial,	&(acBuffLine[6]),  5);	m_iSerial = atoi(acSerial);

	memcpy(m_acAtom,		&(acBuffLine[12]), 4);	CXYUtility::TrimChar(m_acAtom);
	memcpy(m_acAtomicSym,	&(acBuffLine[12]), 2);	CXYUtility::TrimChar(m_acAtomicSym);
	memcpy(m_acRemoteInd,	&(acBuffLine[14]), 1);	CXYUtility::TrimChar(m_acRemoteInd);
	memcpy(m_acBranchDes,	&(acBuffLine[15]), 1);	CXYUtility::TrimChar(m_acBranchDes);


	memcpy( m_acAltLoc,	&(acBuffLine[16]), 1);	CXYUtility::TrimChar(m_acAltLoc);
	memcpy( m_acResName,&(acBuffLine[17]), 3);	CXYUtility::TrimChar(m_acResName);
	memcpy( m_acChainID,&(acBuffLine[21]), 1);	CXYUtility::TrimChar(m_acChainID);
	
	memcpy( acResSeq,	&(acBuffLine[22]), 4);	m_iResSeq = atoi(acResSeq);
	memcpy( m_acICode,	&(acBuffLine[26]), 1);	CXYUtility::TrimChar(m_acICode);
	
	memcpy( acReal8,	&(acBuffLine[30]), 8);	m_fX = float(atof(acReal8));
	memcpy( acReal8,	&(acBuffLine[38]), 8);	m_fY = float(atof(acReal8));
	memcpy( acReal8,	&(acBuffLine[46]), 8);	m_fZ = float(atof(acReal8));
	m_point = CXYPoint3Df(m_fX, m_fY, m_fZ);
	
	memcpy( acReal6,	&(acBuffLine[54]), 6);	m_fOccupancy = float(atof(acReal6));
	memcpy( acReal6,	&(acBuffLine[60]), 6);	m_fTempFactor = float(atof(acReal6));

	memcpy( m_acSegID,	&(acBuffLine[72]), 4);	CXYUtility::TrimChar(m_acSegID);
	
	memcpy( m_acElement,&(acBuffLine[76]), 2);  CXYUtility::TrimChar(m_acElement);
	memcpy( m_acCharge,	&(acBuffLine[78]), 2);	CXYUtility::TrimChar(m_acCharge);
}
*/

//----------------------------------------------------------------------------
char * CXYPDBAtom::GetRecName()
{
	return m_acRecName;
}
void CXYPDBAtom::SetRecName(char *acRecName)
{
	if (char_traits<char>::length(acRecName) >0)
	{
		memcpy(m_acRecName,acRecName,7);
	}
}
//----------------------------------------------------------------------------
int	CXYPDBAtom::GetSerial()
{
	return m_iSerial;
}
void CXYPDBAtom::SetSerial(int iSerial)
{
	m_iSerial = iSerial;
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetName()
{
	return m_acAtom;
}
void CXYPDBAtom::SetName(char *acAtom)
{
	if (char_traits<char>::length(acAtom) >0)
	{
		memcpy(m_acAtom,acAtom,5);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetAtomicSym()
{
	return m_acAtomicSym;
}
void CXYPDBAtom::SetAtomicSym(char *acAtomicSym)
{
	if (char_traits<char>::length(acAtomicSym) >0)
	{
		memcpy(m_acAtomicSym,acAtomicSym,3);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetRemoteInd()
{
	return m_acRemoteInd;
}
void CXYPDBAtom::SetRemoteInd(char *acRemoteInd)
{
	if (char_traits<char>::length(acRemoteInd) >=0)
	{
		memcpy(m_acRemoteInd,acRemoteInd,2);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetBranchDes()
{
	return m_acBranchDes;
}
void CXYPDBAtom::SetBranchDes(char *acBranchDes)
{
	if (char_traits<char>::length(acBranchDes) >=0)
	{
		memcpy(m_acBranchDes,acBranchDes,2);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetAltLoc()
{
	return m_acAltLoc;
}
void CXYPDBAtom::SetAltLoc(char *acAltLoc)
{
	if (char_traits<char>::length(acAltLoc)>=0)
	{
		memcpy(m_acAltLoc,acAltLoc,2);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetResName()
{
	return m_acResName;
}
void CXYPDBAtom::SetResName(char *acResName)
{
	if (char_traits<char>::length(acResName) == 3)
	{
		memcpy(m_acResName,acResName,4);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetChainID()
{
	return m_acChainID;
}
void CXYPDBAtom::SetChainID(char *acChainID)
{
	if (char_traits<char>::length(acChainID) >= 0)
	{
		memcpy(m_acChainID,acChainID,2);
	}
}
//----------------------------------------------------------------------------
int CXYPDBAtom::GetResSeq()
{
	return m_iResSeq;
}
void CXYPDBAtom::SetResSeq(int iResSeq)
{
	m_iResSeq = iResSeq;
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetICode()
{
	return m_acICode;
}
void CXYPDBAtom::SetICode(char *acICode)
{
	if (char_traits<char>::length(acICode)>=0)
	{
		memcpy(m_acICode,acICode,2);
	}
}
//----------------------------------------------------------------------------
float CXYPDBAtom::GetX()
{
	return m_fX;
}
void CXYPDBAtom::SetX(float fX)
{
	m_fX = fX;
	m_point = CXYPoint3Df(m_fX,m_fY,m_fZ);
}
//----------------------------------------------------------------------------
float CXYPDBAtom::GetY()
{
	return m_fY;
}
void CXYPDBAtom::SetY(float fY)
{
	m_fY = fY;
	m_point = CXYPoint3Df(m_fX,m_fY,m_fZ);
}
//----------------------------------------------------------------------------
float CXYPDBAtom::GetZ()
{
	return m_fZ;
}
void CXYPDBAtom::SetZ(float fZ)
{
	m_fZ = fZ;
	m_point = CXYPoint3Df(m_fX,m_fY,m_fZ);
}
//----------------------------------------------------------------------------
float CXYPDBAtom::GetOccupancy()
{
	return m_fOccupancy;
}
void CXYPDBAtom::SetOccupancy(float fOccupancy)
{
	m_fOccupancy = fOccupancy;
}
//----------------------------------------------------------------------------
float CXYPDBAtom::GetTempFactor()
{
	return m_fTempFactor;
}
void CXYPDBAtom::SetTempFactor(float fTempFactor)
{
	m_fTempFactor = fTempFactor;
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetSegID()
{
	return m_acSegID;
}
void CXYPDBAtom::SetSegID(char *acSegID)
{
	if (char_traits<char>::length(acSegID)>=0)
	{
		memcpy(m_acSegID,acSegID,5);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetElement()
{
	return m_acElement;
}
void CXYPDBAtom::SetElement(char *acElement)
{
	if (char_traits<char>::length(acElement)>=0)
	{
		memcpy(m_acElement,acElement,3);
	}
}
//----------------------------------------------------------------------------
char *CXYPDBAtom::GetCharge()
{
	return m_acCharge;
}
void CXYPDBAtom::SetCharge(char *acCharge)
{
	if (char_traits<char>::length(acCharge)>=0)
	{
		memcpy(m_acCharge,acCharge,3);
	}
}
//----------------------------------------------------------------------------
CXYPoint3D<float> CXYPDBAtom::GetPoint() const
{
	return m_point;
}
void CXYPDBAtom::SetPoint(const CXYPoint3Df& rP)
{
	m_fX = rP.X();
	m_fY = rP.Y();
	m_fZ = rP.Z();
	m_point = rP;
}
//----------------------------------------------------------------------------
void CXYPDBAtom::ReadLine(const char *acBuffLine)
{
	char acSerial[6];
	char acResSeq[5];
	char acReal8[9];
	char acReal6[7];

	char acStr[9];

	memset(acSerial, '\0', 6*sizeof(char));
	memset(acResSeq, '\0', 5*sizeof(char));
	memset(acReal8,	'\0', 9*sizeof(char));
	memset(acReal6, '\0', 7*sizeof(char));

	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr,acBuffLine,6); sscanf(acStr,"%s",m_acRecName);
		
	memcpy( acSerial,	&(acBuffLine[6]),  5);	m_iSerial = atoi(acSerial);

	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[12]), 4); sscanf(acStr,"%s",m_acAtom);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[12]), 2); sscanf(acStr,"%s",m_acAtomicSym);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[14]), 1); sscanf(acStr,"%s",m_acRemoteInd);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[15]), 1); sscanf(acStr,"%s",m_acBranchDes);

	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[16]), 1); sscanf(acStr,"%s",m_acAltLoc); 
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[17]), 3); sscanf(acStr,"%s",m_acResName);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[21]), 1); sscanf(acStr,"%s",m_acChainID);

	memcpy( acResSeq,	&(acBuffLine[22]), 4);	m_iResSeq = atoi(acResSeq);

	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[26]), 1); sscanf(acStr,"%s",m_acICode);
	
	memcpy( acReal8,	&(acBuffLine[30]), 8);	m_fX = float(atof(acReal8));
	memcpy( acReal8,	&(acBuffLine[38]), 8);	m_fY = float(atof(acReal8));
	memcpy( acReal8,	&(acBuffLine[46]), 8);	m_fZ = float(atof(acReal8));
	m_point = CXYPoint3Df(m_fX, m_fY, m_fZ);
	
	memcpy( acReal6,	&(acBuffLine[54]), 6);	m_fOccupancy = float(atof(acReal6));
	memcpy( acReal6,	&(acBuffLine[60]), 6);	m_fTempFactor = float(atof(acReal6));

	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[72]), 4); sscanf(acStr,"%s",m_acSegID);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[76]), 2); sscanf(acStr,"%s",m_acElement);
	memset(acStr,'\0', 9*sizeof(char)); memcpy(acStr, &(acBuffLine[78]), 2); sscanf(acStr,"%s",m_acCharge);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
