#include "XYUtility.h"
#include <vector>
#include <algorithm>

//----------------------------------------------------------------------------
CXYUtility::CXYUtility(void)
{
}
//----------------------------------------------------------------------------
CXYUtility::~CXYUtility(void)
{
}
//----------------------------------------------------------------------------
/*
char *CXYUtility::TrimChar(char *pszSource)
{

   //字符串首指针
   char *pszHead = pszSource;
   //用于保存最后非空格类字符的位置的指针
   char *pszLast = &pszSource[strlen(pszSource)-1];
   //查找最后非空格指针的位置
   while (isSpaces(*pszLast))
      --pszLast;
   *(++pszLast) = '\0';
   //查找首个非空格类字符的位置，以便左移字符串
   for ( ; (*pszHead != '\0'); ++pszHead)
   {
      if (! isSpaces(*pszHead))
      {
         break;
      }
   }
   if ( pszSource != pszHead )
      strcpy(pszSource, pszHead);

   return pszSource;
}
*/
//----------------------------------------------------------------------------
float CXYUtility::Distance(CXYPoint3D<float> P1,CXYPoint3D<float> P2)
{
	//return CXYMath<float>::Sqrt(
	//	(P1.X()-P2.X())*(P1.X()-P2.X()) +
	//	(P1.Y()-P2.Y())*(P1.Y()-P2.Y()) +
	//	(P1.Z()-P2.Z())*(P1.Z()-P2.Z()) );
	return 1/CXYMath<float>::FastInvSqrt(
		(P1.X()-P2.X())*(P1.X()-P2.X()) +
		(P1.Y()-P2.Y())*(P1.Y()-P2.Y()) +
		(P1.Z()-P2.Z())*(P1.Z()-P2.Z()) );
}
//----------------------------------------------------------------------------
float CXYUtility::SquareLength(CXYPoint3D<float> P1,CXYPoint3D<float> P2)
{
	return ((P1.X()-P2.X())*(P1.X()-P2.X()) +
		(P1.Y()-P2.Y())*(P1.Y()-P2.Y()) +
		(P1.Z()-P2.Z())*(P1.Z()-P2.Z()) );
}
//----------------------------------------------------------------------------
char *CXYUtility::TrimChar(char *pszSource)
{
	char buf[100];
	sscanf(pszSource,"%s",buf);
	sprintf(pszSource,"%s",buf);
	return pszSource;
}

char* CXYUtility::BaseName(char *name, int with_suffix)
{
/*
	char *pdest;
	int result;
	int ch = '/';

	pdest = strrchr( fullFileName, ch );
	result = (int)(pdest - fullFileName + 1);
	if ( pdest != NULL )
		return (pdest+1);
	else
		return fullFileName;
*/
//	char *_name = (char *)(malloc (strlen (name) + 1));
//	if (_name == NULL) return NULL;
//	strcpy (_name,name);
	char *_name = strdup(name);

	char *_basename = strrchr(_name, '/');

	if (!_basename)
		_basename = _name;
	else
		_basename++;
	if (!with_suffix)
	{
		char *ext = strrchr(_basename, '.');

		if (ext)
			*ext = 0;
	}
	return _basename;
}

vector<int> CXYUtility::RandPerm( vector<int>& kV )
{
	int N = (int)(kV.size());
	vector<int> kVNew(N);
	int p[RANDOM_MAX];
  assert (N >= 0);
  for (int i = 0; i < N; ++i) {
    int j = rand() % (i + 1);
		p[i] = p[j];
		p[j] = i;
  }
	for (int i=0; i < N; ++i) {
		kVNew[i] = kV[p[i]];
	}
	return kVNew;
}

