#ifndef _INCLUDED_PDB_HPP
#define _INCLUDED_PDB_HPP

#include <vector>
#include <string>
#include <cassert>

#include "constants.hpp"
#include "PdbRecord.hpp"
#include "Graph.hpp"
class PdbRecord;
using namespace std;
class Pdb{
	private:
		vector<PdbRecord> pdbVector;
		unsigned sz;	///size -> no of atoms
	public:
		double  dihedralEnergy,iDihedralEnergy, bondEnergy, angleEnergy,hydroEnergy, \
			vanderEnergy,electroEnergy,vander1_4Energy,electro1_4Energy; ///@surojit hydroEnergy
		
		Pdb(){};
		Pdb(char *filename){read(filename);}
		Pdb(string& filename){read(filename);}
		void read(char* );
		void read(string&);
		void write(char* );
		void write(string&);
		unsigned size() const { return sz; } //size of Pdb ie total no of atoms
		unsigned residueSize()const { return (*this)[size()].getResSeq(); } //no of residue
		PdbRecord& operator[](unsigned i) { 
			//assert(i>0);
	       		return pdbVector[i-1];			
		}
		const PdbRecord& operator[](unsigned i)const {
			//assert(i>0);
	       		return pdbVector[i-1];			
		}
};

#endif
