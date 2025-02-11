/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/* -----------------------------------------------------------------------------
 *  XYSO3Sequence.cpp
 *  CHROMATIN_SIS_COARSE
 *
 *  Created by Yun Xu on 10/11/11.
 *  Copyright 2011 __UIC__. All rights reserved.
 *
 *----------------------------------------------------------------------------- */
/* -----------------------------------------------------------------------------
 * Example usage:
 *	CXYSO3Sequence sO3sequence(72);
 *	sO3sequence.SetSO3Sequence();
 *	sO3sequence.printSO3Sequence();
 *----------------------------------------------------------------------------- */


#include "XYSO3Sequence.h"
//----------------------------------------------------------------------------

CXYSO3Sequence::CXYSO3Sequence(){
	m_iNumPoints = 72;
	SetSequenceBase();
	m_pMSpherepoints = new CXYMatrix<float> (m_iNumPoints,3);
}
//----------------------------------------------------------------------------
CXYSO3Sequence::CXYSO3Sequence(int iNumPoints){
	m_iNumPoints = iNumPoints;
	SetSequenceBase();
	m_pMSpherepoints = new CXYMatrix<float> (iNumPoints,3);
}
//----------------------------------------------------------------------------

CXYSO3Sequence::~CXYSO3Sequence(){
	delete m_pMSpherepoints;
}

//----------------------------------------------------------------------------
void CXYSO3Sequence::pix2ang_nest(long nside, long ipix, double *theta, double *phi){
  /*
	 c=======================================================================
	 subroutine pix2ang_nest(nside, ipix, theta, phi)
	 c=======================================================================
	 c     gives theta and phi corresponding to pixel ipix (NESTED) 
	 c     for a parameter nside
	 c=======================================================================
	 */
	
	int npix, npface, face_num;
	int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
	int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
	double z, fn, fact1, fact2;
	double piover2=0.5*M_PI;
	int ns_max=8192;
	
	static int pix2x[1024], pix2y[1024];
	//      common /pix2xy/ pix2x, pix2y
	
	int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
	//      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
	//      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
	jrll[0]=2;
	jrll[1]=2;
	jrll[2]=2;
	jrll[3]=2;
	jrll[4]=3;
	jrll[5]=3;
	jrll[6]=3;
	jrll[7]=3;
	jrll[8]=4;
	jrll[9]=4;
	jrll[10]=4;
	jrll[11]=4;
	jpll[0]=1;
	jpll[1]=3;
	jpll[2]=5;
	jpll[3]=7;
	jpll[4]=0;
	jpll[5]=2;
	jpll[6]=4;
	jpll[7]=6;
	jpll[8]=1;
	jpll[9]=3;
	jpll[10]=5;
	jpll[11]=7;
	
	
	if( nside<1 || nside>ns_max ) {
		fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
		exit(0);
	}
	npix = 12 * nside*nside;
	if( ipix<0 || ipix>npix-1 ) {
		fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
		exit(0);
	}
	
	/* initiates the array for the pixel number -> (x,y) mapping */
	if( pix2x[1023]<=0 ) mk_pix2xy(pix2x,pix2y);
	
	fn = 1.*nside;
	fact1 = 1./(3.*fn*fn);
	fact2 = 2./(3.*fn);
	nl4   = 4*nside;
	
	//c     finds the face, and the number in the face
	npface = nside*nside;
	
	face_num = ipix/npface;//  ! face number in {0,11}
	ipf = (int)fmod((float)ipix,(float)npface);//  ! pixel number in the face {0,npface-1}
	
	//c     finds the x,y on the face (starting from the lowest corner)
	//c     from the pixel number
	ip_low = (int)fmod((float)ipf,(float)1024);//       ! content of the last 10 bits
	ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
	ip_med = (int)fmod((float)ip_trunc,(float)1024);//  ! content of the next 10 bits
	ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits
	
	ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
	iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];
	
	//c     transforms this in (horizontal, vertical) coordinates
	jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
	jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}
	
	//c     computes the z coordinate on the sphere
	//      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
	jr =  jrll[face_num]*nside - jrt - 1;
	//      cout << "face_num=" << face_num << endl;
	//      cout << "jr = " << jr << endl;
	//      cout << "jrll(face_num)=" << jrll[face_num] << endl;
	//      cout << "----------------------------------------------------" << endl;
	nr = nside;//                  ! equatorial region (the most frequent)
	z  = (2*nside-jr)*fact2;
	kshift = (int)fmod(float(jr - nside),(float) 2);
	if( jr<nside ) { //then     ! north pole region
		nr = jr;
		z = 1. - nr*nr*fact1;
		kshift = 0;
	}
	else {
		if( jr>3*nside ) {// then ! south pole region
			nr = nl4 - jr;
			z = - 1. + nr*nr*fact1;
			kshift = 0;
		}
	}
	*theta = acos(z);
	
	//c     computes the phi coordinate on the sphere, in [0,2Pi]
	//      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
	jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
	if( jp>nl4 ) jp = jp - nl4;
	if( jp<1 )   jp = jp + nl4;
	
	*phi = (jp - (kshift+1)*0.5) * (piover2 / nr);
	
}
//----------------------------------------------------------------------------
void CXYSO3Sequence::mk_pix2xy(int *pix2x, int *pix2y){
  /* =======================================================================
   * subroutine mk_pix2xy
   * =======================================================================
   * constructs the array giving x and y in the face from pixel number
   * for the nested (quad-cube like) ordering of pixels
   *
   * the bits corresponding to x and y are interleaved in the pixel number
   * one breaks up the pixel number by even and odd bits
   * =======================================================================
   */
	
  int i, kpix, jpix, IX, IY, IP, ID;
  for (i = 0; i < 1023; i++) pix2x[i]=0;
  
  for( kpix=0;kpix<1024;kpix++ ) {
    jpix = kpix;
    IX = 0;
    IY = 0;
    IP = 1 ;//              ! bit position (in x and y)
    while( jpix!=0 ){// ! go through all the bits
      ID = (int)fmod((float)jpix,(float)2);//  ! bit value (in kpix), goes in ix
      jpix = jpix/2;
      IX = ID*IP+IX;
      
      ID = (int)fmod((float)jpix,(float)2);//  ! bit value (in kpix), goes in iy
      jpix = jpix/2;
      IY = ID*IP+IY;
      
      IP = 2*IP;//         ! next bit (in x and y)
    }
    
    pix2x[kpix] = IX;//     ! in 0,31
    pix2y[kpix] = IY;//     ! in 0,31
  }
  
  /* Later */
  return;
	
}
//----------------------------------------------------------------------------
vector<double> CXYSO3Sequence::find_point(int base_grid, long int point,long int level,long int healpix_point)
{
	int position=point%4;
	long int quo=0;
	double theta=0,phi=0;
	double vec[3];
	vector <double> Point;
	if(base_grid == 6 or base_grid == 7)
	{
		switch(position)//this switch statement translates between sequence of healpix and sequence for uniform points 
		{
			case 0:
				healpix_point+=3;
				break;
			case 1:
				healpix_point+=0;
				break;
			case 2: 
				healpix_point+=2;
				break;
			case 3:
				healpix_point+=1;
				break;
		}	
	}
	else if(base_grid == 3 or base_grid == 1 or base_grid == 9 or base_grid == 11)
	{
		switch(position)//this switch statement translates between sequence of healpix and sequence for uniform points 
		{
			case 0:
				healpix_point+=3;
				break;
			case 1:
				healpix_point+=0;
				break;
			case 2: 
				healpix_point+=1;
				break;
			case 3:
				healpix_point+=2;
				break;
		}	
	}
	else if(base_grid == 2 or base_grid == 0 or base_grid == 8 or base_grid == 10)
	{
		switch(position)//this switch statement translates between sequence of healpix and sequence for uniform points 
		{
			case 0:
				healpix_point+=0;
				break;
			case 1:
				healpix_point+=3;
				break;
			case 2: 
				healpix_point+=1;
				break;
			case 3:
				healpix_point+=2;
				break;
		}	
	}
	else if(base_grid == 4 or base_grid == 5)
	{
		switch(position)//this switch statement translates between sequence of healpix and sequence for uniform points 
		{
			case 0:
				healpix_point+=0;
				break;
			case 1:
				healpix_point+=3;
				break;
			case 2: 
				healpix_point+=2;
				break;
			case 3:
				healpix_point+=1;
				break;
		}	
	}
	
	quo=point/4;
	if(quo==0)
	{
		long int nside=(int) pow((float)2,(float)level);
		pix2ang_nest(nside,healpix_point,&theta,&phi);
		ang2vec(theta,phi,vec);
		Point.resize(0);
		Point.push_back(vec[0]);	
		Point.push_back(vec[1]);	
		Point.push_back(vec[2]);
		return Point;
	}
	else
	{
		return find_point(base_grid,quo-1,level+1,4*healpix_point);
	}
	
}



//----------------------------------------------------------------------------
void	CXYSO3Sequence:: SetSequenceBase(){
	static const int arr[12] = {
		6, 4, 1, 11, 9, 3, 5, 7, 10, 0, 2, 8
	};
	
	m_vSequence_base = vector<int> (arr, arr+sizeof(arr)/sizeof(int));
}
//----------------------------------------------------------------------------
vector<int >& CXYSO3Sequence::GetSequenceBase(){
	return m_vSequence_base;
}

//----------------------------------------------------------------------------
int CXYSO3Sequence::GetNumPoints(){
	return m_iNumPoints;
}
//----------------------------------------------------------------------------
void CXYSO3Sequence::ang2vec(double theta, double phi, double *vec) {
	
  double sz;
  double PI=M_PI;
	
  if( theta<0. || theta>PI) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
	
  sz = sin(theta);
	
  vec[0] = sz * cos(phi) ;
  vec[1] = sz * sin(phi) ;
  vec[2] = cos(theta)    ;
	
}

//----------------------------------------------------------------------------
int CXYSO3Sequence::SetSO3Sequence(){
	long int num_points=0;
	//	int num;
	vector < double > Points;
	vector <int> Sequence_base;
	ifstream input;
	ofstream output;
	int limit=0;
	double theta=0,phi=0;
	double vec[3];
	long int base_grid=0,cur_point=0;
	long int point_healpix=0;
	
	num_points = m_iNumPoints;
	Sequence_base = GetSequenceBase();
	
	
	//first twelve points are the base grid points;
//	output.open("data.pts");
	if(num_points<12)
		limit=num_points;
	else
		limit=12;
	for(int i=0;i<limit;i++)
	{
		pix2ang_nest(1,Sequence_base[i],&theta,&phi);
		//		cout << theta*180/M_PI << " " << phi*180/M_PI << endl; 
		ang2vec(theta,phi,vec);
//		output << vec[0] << "\t" << vec[1] << "\t" << vec[2] << endl;
		(*m_pMSpherepoints)[i][0] = vec[0];
		(*m_pMSpherepoints)[i][1] = vec[1];
		(*m_pMSpherepoints)[i][2] = vec[2];
		
	}
	
	for(int i=0;i<num_points-12;i++)
	{
		Points.resize(0);
		base_grid=i%12;
		cur_point=i/12;
		point_healpix=4*Sequence_base[base_grid];
		Points=find_point(Sequence_base[base_grid],cur_point,1,point_healpix); //current point value,level,current point in healpix
//		output << Points[0] << "\t" << Points[1] << "\t" << Points[2] << endl;
		(*m_pMSpherepoints)[i+12][0] = Points[0];
		(*m_pMSpherepoints)[i+12][1] = Points[1];
		(*m_pMSpherepoints)[i+12][2] = Points[2];
	}
//	output.close();
	return 0;
}


CXYMatrix<float>& CXYSO3Sequence::GetSO3Sequence(){
	return *m_pMSpherepoints;
}

void CXYSO3Sequence::printSO3Sequence(){
	CXYFile kFile;
	CXYMatrix<float>& kM = GetSO3Sequence();
	kFile.WriteMatrix((char*)"data.pts","w",kM);
}
