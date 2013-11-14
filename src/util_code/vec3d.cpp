//
// This file is released under the terms of the NASA Open Source Agreement (NOSA)
// version 1.3 as detailed in the LICENSE file which accompanies this software.
//

//******************************************************************************
//    
//   3D Vector Double Class
//   
//   J.R. Gloudemans - 7/7/93
//******************************************************************************

#include <math.h>
#include <float.h> //For DBL_EPSILON
#include <stdio.h>
#include "vec3d.h"
#include <cmath>

//****** Constructor:  vec3d x() ******//
 vec3d::vec3d()
{
     v[0] = v[1] = v[2] = 0.0;
}

//****** Constructor:  vec3d x(1.,0.,3.) ******//
 vec3d::vec3d(double xx, double yy, double zz)
{
    v[0] = xx;
    v[1] = yy;
    v[2] = zz;
}

//******* vec3d x = y ******//
 vec3d::vec3d(const vec3d& a)		
{
    v[0] = a.v[0];
    v[1] = a.v[1];
    v[2] = a.v[2];
}

//****** Equals:  x = y ******
 vec3d& vec3d::operator=(const vec3d& a)  
{
	if (this == &a)
	   return *this;

    v[0] = a.v[0];
    v[1] = a.v[1];
    v[2] = a.v[2];
    return *this;
}

//******* x = 35. ******//
 vec3d& vec3d::operator=(double a)
{
    v[0] = v[1] = v[2] = a;
    return *this;
}

//******* Set Point Values *******//
 vec3d& vec3d::set_xyz(double xx, double yy, double zz)
{
    v[0] = xx;
    v[1] = yy;
    v[2] = zz;
    return *this;
}

//******* Set Point Values *******//
 vec3d& vec3d::set_x(double xx)
{
    v[0] = xx;
    return *this;
}

//******* Set Point Values *******//
 vec3d& vec3d::set_y(double yy)
{
    v[1] = yy;
    return *this;
}

//******* Set Point Values *******//
 vec3d& vec3d::set_z(double zz)
{
    v[2] = zz;
    return *this;
}


//******* Transform *******//
 vec3d vec3d::transform(float mat[4][4]) const
{
    return( vec3d( (mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]),
                   (mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]),
                   (mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]) ) );

}

//******* Transform *******//
 vec3d vec3d::transform(double mat[4][4]) const
{
    return( vec3d( (mat[0][0]*v[0] + mat[1][0]*v[1] + mat[2][0]*v[2] + mat[3][0]),
                   (mat[0][1]*v[0] + mat[1][1]*v[1] + mat[2][1]*v[2] + mat[3][1]),
                   (mat[0][2]*v[0] + mat[1][2]*v[1] + mat[2][2]*v[2] + mat[3][2]) ) );

}

//************* x = a + b ******//
 vec3d operator+(const vec3d& a, const vec3d& b)
{
    vec3d  ret(a.v[0] + b.v[0],  a.v[1] + b.v[1],  a.v[2] + b.v[2]);
    return ret;
}

//************* x = a - b ******//
 vec3d operator-(const vec3d& a, const vec3d& b)
{
    vec3d  ret(a.v[0] - b.v[0],  a.v[1] - b.v[1],  a.v[2] - b.v[2]);
    return ret;
}

//******* x = a * b ******//
 vec3d operator*(const vec3d& a, double b)
{
    vec3d ret(a.v[0] * b, a.v[1] * b, a.v[2] * b);
    return ret;
}

//******* x = a * b ******//
 vec3d operator*(const vec3d& a, const vec3d& b)
{
    vec3d ret(a.v[0] * b.v[0], a.v[1] * b.v[1], a.v[2] * b.v[2]);
    return ret;
}

//******* x = a / b ******//
 vec3d operator/(const vec3d& a, double b)
{
    vec3d ret;
    //if (b != 0.0)
	if (!( b <= DBL_EPSILON && b >= 0.0))
      ret.set_xyz(a.v[0] / b, a.v[1] / b, a.v[2] / b);
    else
      ret.set_xyz(0.0, 0.0, 0.0);

    return ret;
}


//******* cout << a ******//
//ostream& operator<< (ostream& out, const vec3d& a)
//{
//    return ( out << "  " << a.v[0] << "  " <<
//	     a.v[1] << "  " << a.v[2] << "  " ) ;
//}

//******* distance between pnts ******//
double dist(const vec3d& a, const vec3d& b) 
{
    
    double xx = a.v[0] - b.v[0];   
    double yy = a.v[1] - b.v[1];   
    double zz = a.v[2] - b.v[2];
    return ( sqrt( xx*xx + yy*yy + zz*zz ));
}

//******* distance between pnts ******//
double dist_squared(const vec3d& a, const vec3d& b) 
{
    double xx = a.v[0] - b.v[0];   
    double yy = a.v[1] - b.v[1];   
    double zz = a.v[2] - b.v[2];
    return ( xx*xx + yy*yy + zz*zz );
}

//******* Magnitude:   x = a.mag() ******//
 double vec3d::mag() const
{
    return( sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) );
} 

//****** Normalize:   a.normalize()  ******//
 void vec3d::normalize()
{

    double length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (length <= 0.0)
      {
        v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
      }
    else
      {
        v[0] /= length; v[1] /= length; v[2] /= length;
      }
}


int vec3d::major_comp() const
{
	int i = 0;
	double c = abs( v[i] );

	if( abs(v[1]) > c )
	{
		i = 1;
		c = abs( v[i] );
	}

	if( abs(v[2]) > c )
	{
		i = 2;
		c = abs( v[i] );
	}
	return i;
}

int vec3d::minor_comp() const
{
	int i = 0;
	double c = abs( v[i] );

	if( abs(v[1]) < c )
	{
		i = 1;
		c = abs( v[i] );
	}

	if( abs(v[2]) < c )
	{
		i = 2;
		c = abs( v[i] );
	}
	return i;
}

//******* Dot Product:  x = a.dot(b) ******//
 double dot(const vec3d& a, const vec3d& b)
{
    return (a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2]) ;
}

//******* Cross Product:  a = cross(a, b) ******//
 vec3d cross(const vec3d& a, const vec3d& b)
{
    return vec3d(a.v[1]*b.v[2] - a.v[2]*b.v[1], 
                 a.v[2]*b.v[0] - a.v[0]*b.v[2], 
                 a.v[0]*b.v[1] - a.v[1]*b.v[0]);
}


//******* Angle Between Vectors ******//
double angle(const vec3d& a, const vec3d& b)
{
    double angle = dot(a, b)/(a.mag()*b.mag());

    if ( angle < -1.0 )			angle = -1.0;
    else if ( angle > 1.0 )		angle = 1.0;

    return(acos(angle));

}
		
double signed_angle(const vec3d& a, const vec3d& b, const vec3d& ref )
{
	double ang = angle( a, b );
	vec3d c = cross( a, b );

	double d = dot( c, ref );
	if ( d < 0 )
		ang = -ang;

	return ang;
}



//v1.Normalize();
//v2.Normalize();
//Vec3 c=v1.CrossProduct(v2);
//angle=std::atan2(c.Magnitude(),v1.Dot(v2));
//angle=c.Dot(Vec3(1.f,0.f,0.f)) < 0.f ? -angle : angle;
//std::cout<<"angle:"<<angle<<std::endl;
//if (angle>M_PI_2+0.0001f||angle<-M_PI_2-0.0001f)
//	return true;
//
//
//	double cosine = dot(a, b)/(a.mag()*b.mag());
//	if (cosine > 1) cosine = 1;
//	else if (cosine < -1) cosine = -1;
// 
//  if ((data[0]*v[1] - data[1]*v[0]) < 0)
//    return -acos(cosine);
//  else
//    return acos(cosine);
//}

//******* Cosine of Angle Between Vectors ******//
double cos_angle(const vec3d& a, const vec3d& b)
{
   double angle = dot(a, b)/(a.mag()*b.mag());

   if (angle < -1.0 )
	   return -1.0;
   else if ( angle > 1.0)
	   return 1.0;

   return angle;
}

//******* Radius of Circle Passing Thru 3 Points ******//
double radius_of_circle(vec3d& p1, vec3d& p2, vec3d& p3) 
{
   vec3d a = p1 - p3;
   vec3d b = p2 - p3;

   double denom = 2.0*cross(a, b).mag();
 
   //if (denom == 0.0)
   if (denom <= DBL_EPSILON)
     return(1.0e06);

   else
     {
       return((a.mag()*b.mag()*(a - b).mag())/denom);
     }
}

//******* Center And Radius of Circle Passing Thru 3 Points ******//
void center_of_circle(vec3d& p1, vec3d& p2, vec3d& p3, vec3d& center) 
{
   vec3d a = p1 - p3;
   vec3d b = p2 - p3;

   double temp = cross(a, b).mag();
   double denom = 2.0*temp*temp;
 
   //if (denom == 0.0)
   if( denom <= DBL_EPSILON)
     {
       center = p1 + vec3d(1.0e06, 1.0e06, 1.0e06);
     }

   else
     {

       double a_mag  = a.mag();
       double a_sqr  = a_mag*a_mag;
       double b_mag  = b.mag();
       double b_sqr  = b_mag*b_mag;
       double ab_dot = dot(a,b);
       center = a*(b_sqr*(a_sqr - ab_dot)) + b*(a_sqr*(b_sqr - ab_dot));
       center = center * (1.0/denom) + p3;
 
     }
}


//******* Dist Between Point And Plane ******//
double dist_pnt_2_plane(vec3d& org, vec3d& norm, vec3d& pnt) 
{
  //===== NORM SHOULD BE NORMALIZED ====//
  double d = dot((pnt-org), norm);

  return(fabs(d)); 

}

//******* Dist Between Point And Line ******//
double dist_pnt_2_line(vec3d& line_pt1, vec3d& line_pt2, vec3d& pnt) 
{
  vec3d A_B = pnt - line_pt1;
  vec3d C_B = line_pt2 - line_pt1;

  double denom = C_B.mag();

  if (denom < 0.0 )
    return( A_B.mag() );

  return( cross(A_B, C_B).mag()/denom );
}

//******* Distance Between Point And Line Segment ******//
/*
double dist_pnt_2_line_seg(vec3d& line_pt1, vec3d& line_pt2, vec3d& pnt) 
{
  vec3d p_ln1 = pnt - line_pt1;
  vec3d ln2_ln1 = line_pt2 - line_pt1;

  if ( cos_angle( p_ln1, ln2_ln1 ) < 0.0 )
    return ( dist(pnt, line_pt1) );

  vec3d p_ln2 = pnt - line_pt2;
  vec3d ln1_ln2 = line_pt1 - line_pt2;
  
  if ( cos_angle( p_ln2, ln1_ln2 ) < 0.0 )
    return ( dist(pnt, line_pt2) );

  double denom = ln2_ln1.mag();

  if (denom < 0.0 )
    return( dist(pnt, line_pt1 ));

  return( cross(p_ln1, ln2_ln1).mag()/denom );
}
*/

//******* Project Pnt Onto Line Seg ******//
vec3d proj_pnt_on_line_seg(const vec3d& line_pt1, const vec3d& line_pt2, const vec3d& pnt) 
{
  vec3d p_ln1 = pnt - line_pt1;

  if ( p_ln1.mag() <= 0.0000000001 )
	  return line_pt1;
 
  vec3d ln2_ln1 = line_pt2 - line_pt1;

  if ( cos_angle( p_ln1, ln2_ln1 ) <= 0.0 )
    return ( line_pt1 );

  vec3d p_ln2 = pnt - line_pt2;

  if ( p_ln2.mag() <= 0.0000000001 )
	  return line_pt2;

  vec3d ln1_ln2 = line_pt1 - line_pt2;
  
  if ( cos_angle( p_ln2, ln1_ln2 ) <= 0.0 )
    return ( line_pt2 );

  double denom = ln2_ln1.mag();

  if (denom <= 0.0 )
    return ( line_pt1 );

  double numer =  cos_angle( p_ln1, ln2_ln1 )*p_ln1.mag();

  return(line_pt1 + ln2_ln1*(numer/denom));

}

//******* Project Pnt Onto Ray ******//
vec3d proj_pnt_on_ray(const vec3d& line_pt1, const vec3d& line_pt2, const vec3d& pnt) 
{
  vec3d p_ln1 = pnt - line_pt1;
  vec3d ln2_ln1 = line_pt2 - line_pt1;

  vec3d p_ln2 = pnt - line_pt2;
  vec3d ln1_ln2 = line_pt1 - line_pt2;
  
  double denom = ln2_ln1.mag();

  if (denom <= 0.000000000012 )
    return ( line_pt1 );

  double numer =  cos_angle( p_ln1, ln2_ln1 )*p_ln1.mag();

  return(line_pt1 + ln2_ln1*(numer/denom));

}

//******* Project Pnt Onto Line ******//
vec3d proj_pnt_on_line(const vec3d& line_pt1, const vec3d& line_pt2, const vec3d& pnt) 
{
  vec3d p_ln1 = pnt - line_pt1;
  vec3d ln2_ln1 = line_pt2 - line_pt1;

  vec3d p_ln2 = pnt - line_pt2;
  vec3d ln1_ln2 = line_pt1 - line_pt2;
  
  double denom = ln2_ln1.mag();

  if ( fabs(denom) <= 0.000000000012 )
    return ( line_pt1 );

  double p_ln1_mag = p_ln1.mag();
  if ( fabs(p_ln1_mag) <= 0.000000000012 )
    return ( line_pt1 );


  double numer =  cos_angle( p_ln1, ln2_ln1 )*p_ln1_mag;

  return(line_pt1 + ln2_ln1*(numer/denom));

}


//******* Project Line To Plane******//
//======= NOT TESTED !!!!!!!!!!!!!!! ===//
vec3d proj_pnt_to_plane(vec3d& org, vec3d& plane_ln1, vec3d& plane_ln2, vec3d& pnt) 
{
  vec3d normal = cross(plane_ln1, plane_ln2);

  vec3d proj_pnt = proj_pnt_on_ray(org, org+normal, pnt);

  vec3d proj_vec = org - proj_pnt;

  return(pnt + proj_vec); 

}

//******* Find The Point On Line AB nearest to Line CD******//
//======= NOT TESTED !!!!!!!!!!!!!!! ===//
int ray_ray_intersect(vec3d& A, vec3d& B, vec3d& C, vec3d& D, vec3d& int_pnt1, vec3d& int_pnt2)
{
  vec3d line1 = B - A;
  vec3d line2 = C - D;

  vec3d normal = cross(line1, line2);

  if ( normal.mag() <= 0.0 )
    {
      //===== Paralle Lines =====//
      return(0);
    }
  else
    {
      double t = 0.0;
      if ( plane_ray_intersect(A, line1, normal, C, line2, t) )
        {
          int_pnt2 = C + (line2*t);
        }
      else
        {
          // What UP?
          cout << " RAY RAY INTERSECT - WHAT UP 1? " << endl;
        }
      if ( plane_ray_intersect(C, line2, normal, A, line1, t) )
        {
          int_pnt1 = A + (line1*t);
        }
      else
        {
          // What UP?
          cout << " RAY RAY INTERSECT - WHAT UP 2? " << endl;
        }
    }
  return(1);
}

//******* Triangle - Line Segment Intersection ******//
// ==== Triangle - Line Segment Intersection ====//
// ==== A - Base Point on Triangle
// ==== B - Vector for one   Side of Tri
// ==== C - Vector for other Side of Tri
// ==== D - Base Point for Line Seg
// ==== E - Vector for Line Seg
// ==============================================//
int tri_seg_intersect(vec3d& A, vec3d& B, vec3d& C, vec3d& D, vec3d& E, 
        double& u, double& w, double& t)
{
   double zero = -1.0e-08;
   double one  = 1.0 - zero;

   vec3d cs = cross(B, C);
   double denom = dot(cs, E);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   t = (dot(cs, A) -  dot(cs, D))/denom;

   if ( (t < zero) || (t > one) ) return(0);

   cs = cross(C, E);
   denom = dot(cs, B);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   u = (dot(cs, D) - dot(cs, A))/denom;
       
   if ( (u < zero) || (u > one) ) return(0);

   cs = cross(B, E);
   denom = dot(cs, C);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   w = (dot(cs, D) - dot(cs, A))/denom;

   if ( (w < zero) || (w > one) ) return(0);

   if ( (w+u) > one) return(0);

   return(1);
}

//******* Triangle - Ray Intersection ******//
// ==== Triangle - Line Segment Intersection ====//
// ==== A - Base Point on Triangle
// ==== B - Vector for one   Side of Tri
// ==== C - Vector for other Side of Tri
// ==== D - Base Point for Ray
// ==== E - Vector for Ray
// ==============================================//
int tri_ray_intersect(vec3d& A, vec3d& B, vec3d& C, vec3d& D, vec3d& E, 
        double& u, double& w, double& t)
{
   double zero = -1.0e-08;
   double one  = 1.0 - zero;

   vec3d cs = cross(B, C);
   double denom = dot(cs, E);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   t = (dot(cs, A) -  dot(cs, D))/denom;

   cs = cross(C, E);
   denom = dot(cs, B);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   u = (dot(cs, D) - dot(cs, A))/denom;
       
   if ( (u < zero) || (u > one) ) return(0);

   cs = cross(B, E);
   denom = dot(cs, C);

   //if ( fabs(denom) == 0.0 ) return(0);
   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   w = (dot(cs, D) - dot(cs, A))/denom;

   if ( (w < zero) || (w > one) ) return(0);

   if ( (w+u) > one) return(0);

   return(1);
}


//******* Plane - Ray Intersection ******//
// ==== Plane - Line Segment Intersection ====//
// ==== A - Base Point on Plane
// ==== B - Vector for one   Side of Plane
// ==== C - Vector for other Side of Plane
// ==== D - Base Point for Ray
// ==== E - Vector for Ray
// ==============================================//
int plane_ray_intersect(vec3d& A, vec3d& B, vec3d& C, vec3d& D, vec3d& E, double& t)
{
   vec3d cs = cross(B, C);
   double denom = dot(cs, E);

   if ( fabs(denom) <= DBL_EPSILON ) return(0);

   t = (dot(cs, A) -  dot(cs, D))/denom;

   return(1);
}
//******* Plane - Ray Intersection ******//
// ==== Plane - Line Segment Intersection ====//
// ==== orig - Orig of Plane
// ==== norm - Normal of Plane
// ==== D - Base Point for Ray
// ==== E - Vector for Ray
// ==============================================//
int plane_ray_intersect(vec3d& orig, vec3d& norm, vec3d& D, vec3d& E, double& t)
{
	double denom = dot(norm, E);

	if ( fabs(denom) <= DBL_EPSILON ) return(0);

	t = (dot(norm, orig) -  dot(norm, D))/denom;

	return(1);
}

//******* Signed Volume Of Tetrahedron Defined By ******//
//*******       Three Vecs From Common Pnt        ******//
double tetra_volume(vec3d& A, vec3d& B, vec3d& C)
{
  double determ = A.v[0]*B.v[1]*C.v[2] + B.v[0]*C.v[1]*A.v[2]
                + C.v[0]*A.v[1]*B.v[2] - C.v[0]*B.v[1]*A.v[2]
                - B.v[0]*A.v[1]*C.v[2] - A.v[0]*C.v[1]*B.v[2];

  return( determ/6.0 );
}


//******* Signed Area Of Tri Defined By           ******//
//*******       Three Vecs From Common Pnt        ******//
double area_squared(vec3d& A, vec3d& B, vec3d& C)
{

  double mBA = (B - A).mag();
  double mCA = (C - A).mag();
  double mCB = (C - B).mag();

  double s = 0.5*(mBA + mCA + mCB);

  double area = s*(s-mBA)*(s-mCA)*(s-mCB);

  return( area );

}


//******* Signed Area Of Tri Defined By           ******//
//*******       Three Vecs From Common Pnt        ******//
double area(vec3d& A, vec3d& B, vec3d& C)
{

  double mBA = (B - A).mag();
  double mCA = (C - A).mag();
  double mCB = (C - B).mag();

  double s = 0.5*(mBA + mCA + mCB);

  double area = s*(s-mBA)*(s-mCA)*(s-mCB);

  if ( area <= 0.0 ) 
    return(0.0);
 
//  else if ( area <= 0.0 )
//    return( -sqrt(-area) );

  return( sqrt(area) );

}

// dist3D_Segment_to_Segment based on code by Dan Sunday
// http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
//
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
//
// dist3D_Segment_to_Segment():
double dist3D_Segment_to_Segment( vec3d& S1P0, vec3d& S1P1, vec3d& S2P0, vec3d& S2P1)
{
	double SMALL_NUM = 0.0000001;

    vec3d   u = S1P1 - S1P0;
    vec3d   v = S2P1 - S2P0;
    vec3d   w = S1P0 - S2P0;
    double    a = dot(u,u);        // always >= 0
    double    b = dot(u,v);
    double    c = dot(v,v);        // always >= 0
    double    d = dot(u,w);
    double    e = dot(v,w);
    double    D = a*c - b*b;       // always >= 0
    double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    vec3d   dP = w + (u * sc) - (v* tc );  // = S1(sc) - S2(tc)

    return dP.mag();   // return the closest distance
}

// dist3D_Segment_to_Segment based on code by Dan Sunday
// http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
//
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
//
// dist3D_Segment_to_Segment():
double dist3D_Segment_to_Segment( vec3d& S1P0, vec3d& S1P1, vec3d& S2P0, vec3d& S2P1, 
								 double* Lt, vec3d* Ln, double* St, vec3d* Sn )
{
	double SMALL_NUM = 0.0000001;

    vec3d   u = S1P1 - S1P0;
    vec3d   v = S2P1 - S2P0;
    vec3d   w = S1P0 - S2P0;
    double    a = dot(u,u);        // always >= 0
    double    b = dot(u,v);
    double    c = dot(v,v);        // always >= 0
    double    d = dot(u,w);
    double    e = dot(v,w);
    double    D = a*c - b*b;       // always >= 0
    double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
    double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    vec3d   dP = w + (u * sc) - (v* tc );  // = S1(sc) - S2(tc)

	*Ln = S1P0 + u*sc;
	*Sn = S2P0 + v*tc;

	*Lt = sc;
	*St = tc;

    return dP.mag();   // return the closest distance
}

//==== Find Nearest Points On 2 Line Segs - return dist ====//
double nearSegSeg(const vec3d& L0,const vec3d& L1,const vec3d& S0,const vec3d& S1, 
				  double* Lt, vec3d* Ln, double* St, vec3d* Sn )
{
	vec3d u = L1 - L0;
	vec3d v = S1 - S0;
	vec3d w = L0 - S0;

	double a = dot(u,u);
	double b = dot(u,v);
	double c = dot(v,v);
	double d = dot(u,w);
	double e = dot(v,w);

	double D = a*c - b*b;
	double sc = 0.0;
	double sN = 0.0;
	double sD = D; 
	double tc = 0.0;
	double tN = 0.0;
	double tD = D; 

    // compute the line parameters of the two closest points
    if (D < 0.0000001) { // the lines are almost parallel
        sN = 0.0;
        tN = e;
        tD = c;
		sD = c;
    }
    else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = sN / sD;
    tc = tN / tD;

	*Ln = L0 + u*sc;
	*Sn = S0 + v*tc;

	*Lt = sc;
	*St = tc;

	return dist( *Ln, *Sn );

 
}


double pointLineDistSquared( vec3d& X0, vec3d& X1, vec3d& X2, double* t )
{
	vec3d X10 = X1 - X0;
	vec3d X21 = X2 - X1;

	double denom = dist_squared( X2, X1 );

	if ( denom < 0.000000001 )
		*t = 0.0;
	else
		*t = -dot(X10, X21)/dist_squared( X2, X1 );

	vec3d Xon = X1 + X21 * (*t);

	return dist_squared( Xon, X0 );
}
	
double pointSegDistSquared( vec3d& p, vec3d& sp0, vec3d& sp1, double* t )
{
	double dSqr = pointLineDistSquared( p, sp0, sp1, t );

	if ( *t < 0 )
	{
		*t = 0;
		vec3d vec = p - sp0;
		dSqr = dist_squared( p, sp0 );
	}
	else if ( *t > 1 )
	{
		*t = 1;
		vec3d vec = p - sp1;
		dSqr = dist_squared( p, sp1 );
	}

	return dSqr;

}

vec2d MapToPlane( vec3d & P, vec3d & B, vec3d & e0, vec3d & e1 )
{
	vec2d result;
	vec3d BmP = B - P;
	vec3d zero;
	vec3d me1 = zero - e1;
	double a = dot(e0, e0);
	double b = dot(e0, e1);
	double c = dot(e1, e1);
	double d = dot(e0, BmP);
	double e = dot(e1, BmP);
	double f = dot(BmP, BmP);

	double s = 0;
	double t = 0;
	double denom = a*c - b*b;
	if ( denom )
	{
		s = (b*e - c*d)/denom;
		t = (b*d - a*e)/denom;
	}

	result.set_xy( s, t );
	
	return result;
}

vec3d MapFromPlane( vec2d & uw, vec3d & B, vec3d & e0, vec3d & e1 )
{
	vec3d result = B + e0*uw.x() + e1*uw.y();
	return result;
}

 int plane_half_space( vec3d & planeOrig, vec3d & planeNorm, vec3d & pnt )
 {
	double od = dot( planeNorm, planeOrig );
	double pd = dot( planeNorm, pnt );

	if ( pd > od )
		return 1;
	else
		return -1;
 }

// pint12 = p1 + s*(p2-p1)
// pint34 = p3 + t*(p4-p3)
bool line_line_intersect( vec3d & p1, vec3d & p2, vec3d & p3, vec3d & p4, double* s, double* t )
{
	vec3d p13 = p1 - p3;
	vec3d p43 = p4 - p3;
	double d4343 = dot( p43, p43 );
	if ( d4343 < DBL_EPSILON ) return false;

	vec3d p21 = p2 - p1;
	double d2121 = dot( p21, p21 );
	if ( d2121 < DBL_EPSILON ) return false;

	double d1343 = dot( p13, p43 );
	double d4321 = dot( p43, p21 );
	double d1321 = dot( p13, p21 );

	double denom = d2121*d4343 - d4321*d4321;
	if ( fabs(denom) < DBL_EPSILON ) return false;

	double numer = d1343 * d4321 - d1321 * d4343;

	*s = numer/denom;
	*t = (d1343 + d4321 * (*s)) / d4343;

	return true;
}


////	Rotate a point p by angle theta around an arbitrary axis r
////	Return the rotated point.
////	Positive angles are anticlockwise looking down the axis
////	towards the origin.
////	Assume right hand coordinate system.
vec3d RotateArbAxis(vec3d & p, double theta, vec3d & r)		// Radians
{
   vec3d q( 0, 0, 0 );
   double costheta,sintheta;

   r.normalize();
   costheta = cos(theta);
   sintheta = sin(theta);

   q[0] += (costheta + (1 - costheta) * r[0] * r[0]) * p[0];
   q[0] += ((1 - costheta) * r[0] * r[1] - r[2] * sintheta) * p[1];
   q[0] += ((1 - costheta) * r[0] * r[2] + r[1] * sintheta) * p[2];

   q[1] += ((1 - costheta) * r[0] * r[1] + r[2] * sintheta) * p[0];
   q[1] += (costheta + (1 - costheta) * r[1] * r[1]) * p[1];
   q[1] += ((1 - costheta) * r[1] * r[2] - r[0] * sintheta) * p[2];

   q[2] += ((1 - costheta) * r[0] * r[2] - r[1] * sintheta) * p[0];
   q[2] += ((1 - costheta) * r[1] * r[2] + r[0] * sintheta) * p[1];
   q[2] += (costheta + (1 - costheta) * r[2] * r[2]) * p[2];

   return(q);
}

//==== Find The Area of a Concave Polygon ===//
double poly_area( vector< vec3d > & pnt_vec, vec3d& center )
{
	if ( pnt_vec.size() < 3 )
		return 0.0;

	double total_area = 0.0;
	for ( int i = 0 ; i < (int)(pnt_vec.size()-1) ; i++ )
	{
		total_area += area( center, pnt_vec[i], pnt_vec[i+1] );
	}

	if ( dist( pnt_vec[0], pnt_vec.back() ) > 0.0000001 )
	{
		total_area += area( center, pnt_vec.back(), pnt_vec[0] );
	}

	return total_area;
}
