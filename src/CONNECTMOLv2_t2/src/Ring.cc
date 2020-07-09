// 
// File:   Ring.cc
// Author: vidhu
//
// Created on August 28, 2006, 11:23 AM
//

#include "Ring.h"


//
// Constructor
///
bool Ring::isPlanar()const{
    if(mVect.size()<3)
        return false;
    if(mVect.size()==3)
        return true;
    
    bool two_pi;
    Point_3d pt1 = mPdbRef[mVect[0]].getXYZ();
    Point_3d pt2 = mPdbRef[mVect[1]].getXYZ();
    Point_3d pt3 = mPdbRef[mVect[2]].getXYZ();
    Point_3d pt4 = mPdbRef[mVect[3]].getXYZ();
//    mPdbRef[mVect[0]].print();
    double angle1;
    angle1=(180/PI)*torsion_angle(pt1,pt2,pt3,pt4,two_pi);
//    cout << "Angle1 : "<< angle1 ;
        
    if(mVect.size()>4){
        Point_3d pt5 = mPdbRef[mVect[4]].getXYZ();
        double angle2=(180/PI)*torsion_angle(pt2,pt3,pt4,pt5,two_pi);
      
  //      cout <<"  Angle2 : "<< angle2 <<endl;
        return ( (angle1 < 6) && (angle2 < 6) );
    }
    else
        return (angle1 < 5);  
    return true;
}


double Ring::torsion_angle(const Point_3d &p1,const Point_3d &p2,
                    const Point_3d &p3,const Point_3d &p4, bool &two_pi)const
{
	Point_3d e=(p3-p2).cross_with(p1-p2);
	Point_3d f=(p4-p3).cross_with(p2-p3);
	two_pi=false;
	double temp=divide(e.dot_with(f), e.mod()*f.mod());
	double phi=my_acos(temp);

/*	if(e.dot_with(p4-p3) < 0.0){
		two_pi=true;
	phi=2*PI-phi;
	}
*/
	return phi;
}

inline double Ring::my_acos(double t)const{

	return( (abs(t)>=1.0) ? ( (t<0.0)? PI : 0.0) : acos(t) ) ;
}

bool Ring::isInRing(unsigned atomNo)const{
    for(unsigned i=0;i<mVect.size();i++)
       if(atomNo==mVect[i]){
           return true;          
       }
    return false;
}

unsigned Ring::getAtomNo(int i)const{
    return mVect[i];
}

int Ring::countAtom(string symbol)const{
    int nAtom=0;
    for(int i=0;i<size();i++)
    {
        if(mPdbRef[mVect[i]].getSymbol()==symbol)
            nAtom++;
    }
    return nAtom;   
}

//
// Destructor
//
Ring::~Ring()
{
    for(unsigned i=0;i<mVect.size();i++)
        mPdbRef[mVect[i]].unsetRingPointer();
  
}

