// 
// File:   Ring.h
// Author: vidhu
//
// Created on August 28, 2006, 11:23 AM
//

#ifndef _Ring_H
#define	_Ring_H
#include "Pdb.hpp"
#include "vector"
#include "constants.hpp"
#include "Graph.hpp"
#include "PdbRecord.hpp"
#include "Point_3d.hpp" 

using namespace std;


class Ring {
        
public:
	Ring(const vector<unsigned> vect,Pdb& rPdb)
            :mVect(vect),mPdbRef(rPdb){
            for(unsigned i=0;i<mVect.size();i++)
               mPdbRef[mVect[i]].setRingPointer(this);                 
            }
        bool isPlanar()const;
        int size()const{return mVect.size();}
        bool isInRing(unsigned atomNo)const;
        unsigned getAtomNo(int i)const;
        int countAtom(string symbol)const;
	~Ring();
                 
protected:

private:
        vector<unsigned> mVect;
        Pdb& mPdbRef;
        double torsion_angle(const Point_3d &p1,const Point_3d &p2,\
                    const Point_3d &p3,const Point_3d &p4, bool &two_pi)const;
        inline double my_acos(double t)const;
        inline double divide(double Nr, double Dr)const{
            return((abs(Dr)<VERYSMALL) ?( (Dr<0.0) ? \
            -Nr/VERYSMALL : Nr/VERYSMALL ) : Nr/Dr );
        }
   
};

#endif	/* _Ring_H */

