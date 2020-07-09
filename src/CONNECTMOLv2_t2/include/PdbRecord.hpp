//// dt 14-03-05////////
#ifndef _INCLUDED_PDBRECORD_HPP
#define _INCLUDED_PDBRECORD_HPP

#include <string>
#include <iostream>
#include <ctype.h>
#include "Point_3d.hpp"
#include "constants.hpp"
#include <cstdlib>
#include <cstdio>
//#include "Ring.h"

using namespace std;
class Ring;
class PdbRecord{
    
    friend class Ring;

private:
	enum ATOM_SYMBOL{ H,C,N,O,S,P,Cl,Br,I };
	bool isValidAtomFormat(char* record);
	bool mainChain_;
        short valency;
	
        // a standard Pdb record
        string recordType;      ////ATOM or HETATM			1-6
        unsigned serial;        ///atom No or serial no			7-11
        string name;            ///atom name				13-16
        char altLoc;            ///altername location indicator		17
        string resName;         ///amino or residue name		18-20
        char chainId;           ///chain identifier			22
        unsigned resSeq;        ///Residue sequence number		23-26
        char iCode;             ///Code for insertion of residues	27
	Point_3d coordinates;            ///cartesian coordinates x,y,z		31-38,39-46,47-54
	Ring* mRingPointer;	
        string atomType;
        double charge;
	double Ep; ///epsilon
	double R;  ///Radii of the atom
	int hy; ///hydrophobic flag //@surojit
        void setValency();
        void setRingPointer(Ring* ring){ mRingPointer=ring;}
        void unsetRingPointer(){mRingPointer =0;}
public:
	//////////member function/////
        PdbRecord(){ mRingPointer=0;}
	bool read(char* );
	void write(char* );
	void print()const;

        unsigned getSerial() const{return serial; }
        string getName() const;
        char getAltLoc() const{return altLoc; }
        string getResName() const{ return resName; }
        char getChainId() const{ return chainId; }
        unsigned getResSeq() const{return resSeq; }
        char getICode() const{return iCode; }
        const Point_3d& getXYZ()const{return coordinates; }
        string getAtomType() const{return atomType;}
        double getCharge() const{return charge;}
        string getSymbol() const{ 
		
		char tmpChar[3];
                int j=0;
                for(unsigned i=0;i<name.length();i++)
                    if(!isdigit(name[i])){
                        tmpChar[j]=name[i];
                        j++;
                     } 
		if(tmpChar[0]=='H')
			tmpChar[1]='\0';
                tmpChar[j]='\0';
                //cout << "******"<<tmpChar<<"**********" <<endl;
                return(string(tmpChar));
	}
        short getValency()const{
            return valency;
        }
	int operator==(const PdbRecord&  p) const{ return ( getSerial()==p.getSerial());}
	double getRadius()const {

		char atomSymbol;
		atomSymbol=isdigit(name[0])? 'H' : name[0];

		switch(atomSymbol){
			case 'H':
				return 1.0;
			case 'C': 
				return 1.7;
			case 'N':
				return 1.5;
			case 'O':
				return 1.4;
			case 'S':
				return 1.8;
			default:
				cout<<"getAtomRadius, shouldn't have come to this pt"<<endl;
				exit(-1);
		}
	}
	double getAtomWt() const{

		char atomSymbol;
		atomSymbol=isdigit(name[0])? 'H' : name[0];

		switch(atomSymbol){
			case 'H':
				return  1.0;
			case 'C': 
				return 12.0;
			case 'N':
				return 14.0;
			case 'O':
				return 16.0;
			case 'S':
				return 32.0;
			default:
				cout<<"getAtom, shouldn't have come to this pt"<<endl;
				exit(-1);
		}
	}
	
	bool isMainChain() const{ return mainChain_; }
       	void setXYZ(const Point_3d& crPdb){
		coordinates=crPdb;
	}
       	void setXYZ(double x, double y, double z){
		coordinates.setX(x); coordinates.setY(y); coordinates.setZ(z);
	}
        bool isInRing()const{ return mRingPointer;}
        Ring* getRingPointer()const{ return (mRingPointer );}
        
        ~PdbRecord(){ mRingPointer=0;}
	
	string getAtomSymbol() const;
	float getRadii()const;
};
#endif
