#include "PdbRecord.hpp"

string PdbRecord::getAtomSymbol() const{
    int len=name.length();
    if(len ==1 && !isdigit(name[0]) && isupper(name[0])){
	return name;
    }
    char c[4];
    string s;
    int j=0;
	
    for(int i=0; i<len-1;++i)
    {			
	 if(!isdigit((name[i])))
	 {	 
		 c[j++] = name[i];
		 char nxt=name[i+1];
		 //addition by shailesh
		 if(isupper(nxt))
		 {
		    if(nxt=='L'|| nxt=='R')
			c[j++]=nxt;
		    //if()
		 }
		 //end
		 if(islower(nxt))
		     c[j++] = nxt;
		 break;
	 }	 
    }
    c[j]='\0';
    s=c;
    return (s);
}

float PdbRecord::getRadii()const{
	string atmName=getAtomSymbol();
	if(atmName == "H")
		return 0.30;
	if(atmName == "D")
		 return 0.30;
	if(atmName == "C")
		 return 0.85;
	if(atmName == "N")
		 return 0.78;
	if(atmName == "O")
		 return 0.85;
	if(atmName =="F")
		return 0.53;
	if(atmName == "Na")
		return 0.97;
	if(atmName == "Mg")
		return 0.66;
	if(atmName == "P")
		return 1.08;
	if(atmName == "S")
		return 1.04;
	if(atmName == "Cl" || atmName == "CL")
		return 0.92;
	if(atmName == "K")
		return 1.33;
	if(atmName == "Ca")
		return 0.99;
	if(atmName == "Fe")
		return 0.74;
	if(atmName == "Mn")
		return 0.80;
	if(atmName == "Cu")
		return 0.96;
	if(atmName == "Zn")
		return 0.74;
	if(atmName == "Br" || atmName == "BR")
		return 1.09;
	if(atmName == "I")
		return 1.50;
	if(atmName == "Si")
		return 1.11;
	if(atmName == "AL")
		return 0.64;
	if(atmName == "Li")
		return 0.85;
	if(atmName == "Be")
		return 0.56;
	if(atmName == "B")
		return 0.39;
	if(atmName == "Sc")
		return 0.84;
	if(atmName == "Ti")
		return 0.77;
	if(atmName == "V")
		return 0.74;
	if(atmName == "Cr")
		return 0.72;
	if(atmName == "Co")
		return 0.79;	
	if(atmName == "Ni")
		return 0.79;
	if(atmName == "Ga")
		return 0.72;
	if(atmName == "Ge")
		return 0.82;
	if(atmName == "As")
		return 1.13;
	if(atmName == "Se")
		return 1.13;
	if(atmName == "Sr")
		return 1.12;
	if(atmName == "Ba")
		return 1.41;
	if(atmName == "Ru")
		return 0.78;
	if(atmName == "Rh")
		return 0.76;
	if(atmName == "Pd")
		return 0.95;
	if(atmName == "Ag")
		return 1.02;
	if(atmName == "Cd")
		return 1.03;
	if(atmName == "Pt")
		return 0.89;
	if(atmName == "Au")
		return 1.42;
	if(atmName == "Hg")
		return 1.26;
	if(atmName == "Tl")
		return 1.55;
	if(atmName == "Pb")
		return 1.26;
	cout << "Atom symbol not found for " << atmName <<endl;
	exit(1);
}

string PdbRecord::getName() const{
	if(isdigit(name[0]))
	{
		char ch;
		const unsigned sz = name.size(); 
		string atmString(sz,' ');
		ch=name[0];
		for(unsigned i=1;i<name.size();i++)
		{
			atmString[i-1]=name[i];
		}
		atmString[name.size()-1]=ch;	
			return atmString;
	}
        return name;
}

/*
   
1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C
*/

	
bool PdbRecord::read(char *argData)// reads data from a string in Pdb format
{  
        char recordTypeI[8],readLine[8],nameI[8],resNameI[8];
	char atmType[4];
	double I,J,K;
                                                                                                                           
        //sscanf(argData,"%s%u%s%s%u%lf%lf%lf",recordTypeI,&serial,nameI,resNameI,&resSeq,&I,&J,&K);
        sscanf(argData,"%6s",readLine);
	readLine[7]='\0';
	//coordinates.setXYZ(I,J,K);
        recordType=readLine;
//	cout << recordType <<endl;	
        //name = nameI;
        //resName = resNameI;
	
//	if( recordType=="ATOM" && isValidAtomFormat(argData) ){
if( recordType=="ATOM" || recordType=="HETATM"){
		sscanf(argData,"%6s%5u%*c%4s%c%3s%*c%c%4u%c%*3c%8lf%8lf%8lf%*3c%s%lf%lf%lf",\
			recordTypeI,&serial,nameI,&altLoc,resNameI,&chainId,&resSeq,&iCode,&I,&J,&K,atmType,&R,&Ep,&charge);
		
		coordinates.setXYZ(I,J,K);
		recordType=recordTypeI;
		name = nameI;
		resName = resNameI;
		atomType = atmType;
			
		//set mainchain flag
		string atmName = this->getName();
		if(atmName == "CA" || atmName == "N" || atmName == "O" || atmName == "C" || atmName == "CB" ){
			mainChain_ = true;
		}
		else 
			mainChain_ = false;
		setValency();	
		return true;
	}
        
	return false;

}
/*
void PdbRecord::write(char *argData)// writes data to a string in Pdb format
{ 
        const char *recordTypeI,*nameI,*resNameI,*atomTypeI;
                                                                                                                           
        recordTypeI=recordType.c_str();
        nameI = name.c_str();
        resNameI = resName.c_str();
        atomTypeI = atomType.c_str();
////do the correcting of Hydrogen as used by pankaj/////
//	ckout << nameI << "  " << getAtomName().size() << "     " << getResidueName().size() << " " << getResidueName() <<endl; 
	if(getAtomName().size()==4)
		sprintf(argData,"%-6s%5d %-4s %3s %4d%12.3lf%8.3lf%8.3lf\n",recordTypeI,serial,nameI,resNameI,resSeq,coordinates.x,coordinates.y,coordinates.z);
	else
		sprintf(argData,"%-6s%5d  %-3s %3s %4d%12.3lf%8.3lf%8.3lf\n",recordTypeI,serial,nameI,resNameI,resSeq,coordinates.x,coordinates.y,coordinates.z);
                                                                                                                           
}
*/
void PdbRecord::write(char *argData)// writes data to a string in Pdb format
{ 
        const char *recordTypeI,*nameI,*resNameI,*atomTypeI;
                                                                                                                           
        recordTypeI=recordType.c_str();
        nameI = name.c_str();
        resNameI = resName.c_str();
        atomTypeI = atomType.c_str();
	////do the correcting of Hydrogen as used by pankaj/////
	//	ckout << nameI << "  " << getAtomName().size() << "     " << getResidueName().size() << " " << getResidueName() <<endl; 
	if(getName().size()==4)
		sprintf(argData,"%-6s%5u %-4s%c%3s %c%4u%c   %8.3lf%8.3lf%8.3lf	0 %-2s%8.4lf%8.4lf%8.4lf\n",recordTypeI,serial,nameI,altLoc,resNameI,chainId,resSeq,iCode,coordinates.getX(),coordinates.getY(),coordinates.getZ(),atomTypeI,R,Ep,charge);
	else
		sprintf(argData,"%-6s%5u  %-3s%c%3s %c%4u%c   %8.3lf%8.3lf%8.3lf 0 %-2s%8.4lf%8.4lf%8.4lf\n",recordTypeI,serial,nameI,altLoc,resNameI,chainId,resSeq,iCode,coordinates.getX(),coordinates.getY(),coordinates.getZ(),atomTypeI,R,Ep,charge);
                                                                                                                           
}
                                                                                                                          
void PdbRecord::print()const
{
        const char *recordTypeI,*nameI,*resNameI;
                                                                                                                           
        recordTypeI=recordType.c_str();
        nameI = name.c_str();
        resNameI = resName.c_str();
                                                                                                                           
       printf("%-6s%5d %-4s %3s  %4d    %8.3lf%8.3lf%8.3lf\n",recordTypeI,serial,nameI,resNameI,resSeq,coordinates.getX(),coordinates.getY(),coordinates.getZ());
       cout << mainChain_ << endl;
       //ckprintf("%5d %-4s %3s  %4d  %f  %s \n",serial,nameI,resNameI,resSeq,charge,atomType.c_str());

}

bool PdbRecord::isValidAtomFormat(char* record){

//1         2         3         4         5         6         7         8
//12345678901234567890123456789012345678901234567890123456789012345678901234567890
//ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N

	int i;
	// check col 7-11 serial no. should be integer
	for(i=6;i<11;i++){

		if( !isdigit(record[i]) && !isspace(record[i])  ){
			cout << "ATOM format invalid ! [column 7-11]" << endl;
			return false;
		}
	}
	
	// check col 12 - should be blank
	if( !isspace(record[11]) ){
		cout << "ATOM format invalid ! [column 12 should be blank]" << endl;
		return false;
	}

	//check col 13-16 -  atom name either alphanumeric or blank
	for(i=12;i<16;i++){
		if(!isalnum(record[i]) && !isspace(record[i])){
			cout << "ATOM format invalid ! [column 13-16]"<< endl;
			return false;
		}
	}

	//check col 17 - atlLoc should be alphabet
	if( !isalpha(record[16]) && !isspace(record[16]) ){
		cout << "ATOM format invalid ! [column 17]" << endl;
		return false;
	}
	
	//check col 18-20 - residue name 
	for(i=17;i<20;i++){
		if( !isalnum(record[i]) && !isspace(record[i])){
			cout << "ATOM format invalid ! [column 18-20]" << endl;
			return false;
		}
	}
	
	//check col 21 - should be blank
	if( !isspace(record[20]) ){
		cout << "ATOM format invalid ! [column 21 should be blank]" << endl;
		return false;
	}

	//check col 22 - chain ID integer
	if( !isdigit(record[21]) && !isspace(record[i])){
		cout << "ATOM format invalid ! [column 22]" << endl;
		return false;
	}
	
	//check col 23-26 - residue sequence 
	for(i=22;i<26;i++){
		if( !isdigit(record[i]) && !isspace(record[i])){
			cout << "ATOM format invalid ! [column 23-26]" << endl;
			return false;
		}
	}
	
	//check col 27 - iCode
	if( !isalpha(record[26]) && !isspace(record[26])  ){
		cout << "ATOM format invalid ! [column 27]" << endl;
		return false;
	}

	//check col 28,29,30 - blank
	for(i=27;i<30;i++){
		if( !isspace(record[i]) ){
			cout << "ATOM format invalid ! [column 28,29,30 must be whitespaces]" << endl;
			return false;
		}
	}		

	//check col 31 - 38 x coordinate float
	for(i=30;i<38;i++){
		if( !isdigit(record[i]) && record[i]!='.'&& record[i]!='-' && record[i]!='+'  && !isspace(record[i]) ){
			cout << "ATOM format invalid ! [column 31-38]" << endl;
			return false;
		}
	}		
	
	//check col 39 - 46 y coordinate float
	for(i=38;i<46;i++){
		if( !isdigit(record[i]) && record[i]!='.' && record[i]!='-' && record[i]!='+'  && !isspace(record[i]) ){
			cout << "ATOM format invalid ! [column 39-46]" << endl;
			return false;
		}
	}		
	
	//check col 47 - 54 z coordinate float
	for(i=46;i<54;i++){
		if( !isdigit(record[i]) && record[i]!='.' && record[i]!='-' && record[i]!='+'  && !isspace(record[i]) ){
			cout << "ATOM format invalid ! [column 47-54]" << endl;
			return false;
		}
	}		
	
	return true;
		
		
}

void PdbRecord::setValency(){
    //cout << getAtomSymbol()<<endl;
    if(getAtomSymbol()=="H")
        valency=1;
    else if (getAtomSymbol()=="C")
        valency = 4;
    else if(getAtomSymbol()=="N")
        valency = 3;
    else if(getAtomSymbol()=="O")
        valency = 2;
    else if(getAtomSymbol()=="Cl" || getAtomSymbol()=="CL")
        valency = 1;
    else if(getAtomSymbol()=="Br" || getAtomSymbol()=="BR")
        valency = 1;
    else if(getAtomSymbol()=="I")
        valency = 1;
    else if(getAtomSymbol()=="F")
        valency = 1;
    else if(getAtomSymbol()=="S")
        valency = 6;
    else if(getAtomSymbol()=="P")
        valency = 5;
    else if(getAtomSymbol()=="B")
	valency = 3;
}
