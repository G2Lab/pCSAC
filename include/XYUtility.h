#ifndef XYUTILITY_H
#define	XYUTILITY_H
#pragma warning (disable: 4996)

#include <string.h>
#include <stdio.h>
#include <vector>
#include "XYPoint3D.h"

#define isSpaces(ch) (	ch == ' '  ||\
						ch == '\t' ||\
						ch == '\r' ||\
						ch == '\b' ||\
						ch == '\f' ||\
						ch == '\n')
const int RANDOM_MAX = 32767;
using namespace std;
class CXYUtility{
public:
	CXYUtility(void);
	~CXYUtility(void);
	//static char *TrimChar(char *pszSource);
	static float Distance(CXYPoint3D<float> P1,CXYPoint3D<float> P2);
	static float SquareLength(CXYPoint3D<float> P1,CXYPoint3D<float> P2);
	static char *trim_right( char *szSource );
	static char *trim_left( char *szSource );
	static char *TrimChar( char *szSource );
	//static char *BaseName(char *fullFileName);
	static char *BaseName(char *name, int with_suffix);
	static vector<int> RandPerm(vector<int>& kV);
private:
	
};
template<class T1, class T2 >
struct sort_pair_first_greater {
	bool operator()(const pair<T1,T2>&left, const pair<T1,T2>&right) {
		return left.first > right.first;
	}
};

template<class T1, class T2 >
struct sort_pair_first_less {
	bool operator()(const pair<T1,T2>&left, const pair<T1,T2>&right) {
		return left.first < right.first;
	}
};

template<typename T1,typename T2,typename T3>
class CXYTriple
{
public:
	CXYTriple(const T1 &t1,const T2 &t2,const T3 &t3):
	first(t1),
	second(t2),
	third(t3)
	{
	}
	T1 first;
	T2 second;
	T3 third;
};

#endif
