
#include <iostream>
#include "Graph.hpp"
#include "Pdb.hpp"
#include "Ring.h"
#include <cstdio>
using namespace std;

/*
    Author	    :   Vidhu S. Pandey
    Description	    :   Main module for generating connectivity and bond order in
			small molecules
    Date	    :	10/12/2006
    Revised	    :	09/06/2007
    
    Version 2.0
    Antechamber bond distance is incorporated (no need for bonddist.txt)
    
    
    Copyright. All rights reserved.
*/

void appendConnectInfo(Graph<unsigned>& graph,char *fileName);
bool isAromatic(Graph<unsigned>& graph,Pdb& pdb,Ring* ring);
bool isResonatingRing(Graph<unsigned>& graph,Pdb& pdb,Ring* ring);
int getTotalBondOrder(Graph<unsigned>& graph,unsigned atomNo);
int getNSingleBond(Graph<unsigned>& graph,unsigned atomNo);
int getNDoubleBond(Graph<unsigned>& graph,unsigned atomNo);
int getNTripleBond(Graph<unsigned>& graph,unsigned atomNo);
int hasResonatingOxygenBond(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo);
bool isAmideNitrogen(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo);
int isAmideCarbon(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo);
bool isValidBondLength(float dist, float ri, float rj);
void printUsage(){
    cout <<"\nUSAGE : ./executable   arg1(input file)   arg2(output file) [OPTION](0 or 1)"<<endl;
    cout << "OPTION :\n"
         << "0(default) -  prints normal bond orders without resonating bonds\n"
	 << "1          -  prints resonating bonds as 4 \n"<<endl;

}

int main(int argc,char** argv){

    if(argc==2){
	cout<< "\nPlease specify output file name"<<endl; 
	printUsage();
	exit(0);
    }
    else if(argc==1){
	printUsage();
	exit(0);
    }
    int makeResonating=0;
    if(argc==4){
	if(atoi(argv[3])==1)
		makeResonating = true;
    }
    
    
    FILE *fa;
    int errorFlag=0;
//    cout << "Before reading "<<endl;
    Pdb ligand(argv[1]);
//    cout << "After reading "<<endl;
    Graph<unsigned> graph;
    vector<vector<unsigned> > rings;
    vector<Ring*> ringVect;
    
    fa = fopen("output.log","a");
    fprintf(fa,"\n***** %s \n",argv[1]);

    //adding atom serial no.s in graph
    for(unsigned i=1;i<=ligand.size();i++){
        graph.addVertex(ligand[i].getSerial());
    }

    //connecting atoms based on bond distances
    for(unsigned i=1;i<ligand.size();++i)
    {
	Point_3d pt=ligand[i].getXYZ();
	float ri=ligand[i].getRadii();
	for(unsigned j=i+1;j<=ligand.size();++j)
	{
	    double dist=pt.distance(ligand[j].getXYZ());
	    float rj=ligand[j].getRadii();
	    if(isValidBondLength(dist,ri,rj)){
		 graph.addEdge(i,j);
	    }
	} 
    }	 

    //setting the degree of each node in graph
    graph.setDegree();
    //graph.print();
    
    
    //print warning if any atom is not  connected
    for(unsigned i=1;i<ligand.size();i++){
       if(graph.getDegree(i)==0){
          fprintf(fa,"ERROR ! Atom %d not connected\n",i);
          errorFlag=1;
       }
    }


    //finding rings
    unsigned st=1;
    rings=graph.findCycle(st);
    //cout << rings.size()<<endl;
    for(unsigned i=0;i<rings.size();i++){
        if(rings[i].size()<9){		
            cout <<rings[i].size() <<  " MEMBERED RING" << endl;
            for( unsigned j=0;j<rings[i].size();j++)
                cout << rings[i][j] << " " ;
                Ring* r = new Ring(rings[i],ligand); 
                ringVect.push_back(r);
            //    cout <<"Planarity "<< r->isPlanar()<<"  "<<endl;
                cout << endl << endl;
        }
    }
   
//cout <<"**************"<<endl;	    
    rings.clear();
 
    //singe bonds & single degree atoms
    for(unsigned i=1;i<=ligand.size();i++){
        vector<unsigned> tmpVect;
        int degree=graph.getDegree(i);
        int valency=ligand[i].getValency();
        graph.getAdjacencyList(i,tmpVect);
        
	if(graph.getDegree(i)==1){
           if(ligand[i].getAtomSymbol()=="O"){ 
                //cout << ligand[tmpVect[0]].getAtomSymbol()<<endl;
               //carboxylic RCOO  charged or RNOO	
               if(ligand[tmpVect[0]].getAtomSymbol()=="C" ||
                  ligand[tmpVect[0]].getAtomSymbol()=="N"){
                   
                    vector<unsigned> tmpVect2;
                    graph.getAdjacencyList(tmpVect[0],tmpVect2);
                    int flag=0;
                    //check for 1 existing C=O or N=O bond
                    for(unsigned j=0;j<tmpVect2.size()&& flag==0;j++){
                       if(ligand[tmpVect2[j]].getAtomSymbol()=="O"&&
                          graph.getWeight(tmpVect[0],tmpVect2[j])==2)
                          flag=1;
                    }
                    if(flag==1){
                      graph.setWeight(i,tmpVect[0],1);
                    }
                    else{ 	      
                      graph.setWeight(i,tmpVect[0],2);  
                    }
               }
	       //sulphonyl oxygen SO3-
               else if(ligand[tmpVect[0]].getAtomSymbol()=="S"){
                    vector<unsigned> tmpVect2;
                    graph.getAdjacencyList(tmpVect[0],tmpVect2);
                    int flag=0;
                    //check for 2 S=O bond
                    for(unsigned j=0;j<tmpVect2.size();j++){
                       if(ligand[tmpVect2[j]].getAtomSymbol()=="O"&&
                          graph.getWeight(tmpVect[0],tmpVect2[j])==2)
                          flag++;
                       //cout <<"***"<< ligand[tmpVect2[j]].getAtomSymbol()<<endl; 
                       //cout <<"***"<<graph.getWeight(tmpVect[0],tmpVect2[j])<<endl;
                    }
                    if(flag>=2)
                      graph.setWeight(i,tmpVect[0],1);  
                    else 	      
                      graph.setWeight(i,tmpVect[0],2);  
               }
	    }	       
            else if(ligand[i].getAtomSymbol()=="S") 
               graph.setWeight(i,tmpVect[0],2); //S=           
            else
               graph.setWeight(i,tmpVect[0],valency);// O=,N||| 
        }
	//2 degree bond    
	// -N-S
	else if(degree==2 && (ligand[i].getAtomSymbol()=="N")
                &&(!ligand[i].isInRing()) ){
                    
		if(ligand[tmpVect[0]].getAtomSymbol()=="S"||
		   ligand[tmpVect[1]].getAtomSymbol()=="S"){
			graph.setWeight(i,tmpVect[0],1);
			graph.setWeight(i,tmpVect[1],1);
		}
	}
	
        //for -S-
        else if( degree==2 && ligand[i].getAtomSymbol()=="S"){
            for(unsigned j=0;j<tmpVect.size();j++){
                if((graph.getWeight(i,tmpVect[j])==0))                
                    graph.setWeight(i,tmpVect[j],1);
            }
        }
	
	//degree == valency
        //for C,N,-O-
        else if( (degree==valency)&& degree>1)
        {
           for(unsigned j=0;j<tmpVect.size();j++){
                if((graph.getWeight(i,tmpVect[j])==0))                
                    graph.setWeight(i,tmpVect[j],1);
           }               
        }
        
    }
    
    //6 membered ring aromatic
    for(unsigned i=0;i<ringVect.size();i++){
        if(ringVect[i]->size()==6){
            if(isAromatic(graph,ligand,ringVect[i])){    
                Ring *tmpRing=ringVect[i];
		int flag=-1;
		for(unsigned j=0;j<6 && flag==-1;j++){
                 //   cout<<j <<" "<<getTotalBondOrder(graph,tmpRing->getAtomNo(j))<<endl;
                    if(getTotalBondOrder(graph,tmpRing->getAtomNo(j))>1){
                            flag=j;
                    }
		}
		if(flag!=-1){
                    unsigned atomNo1;
                    unsigned atomNo2;
                    unsigned atomNo3;
                    if(flag==0)
                            atomNo1=tmpRing->getAtomNo(5);
                    else
                            atomNo1=tmpRing->getAtomNo(flag-1);

                    if(flag==5)
                            atomNo3=tmpRing->getAtomNo(0);
                    else
                            atomNo3=tmpRing->getAtomNo(flag+1);

                    atomNo2=tmpRing->getAtomNo(flag);

                    if(graph.getWeight(atomNo1,atomNo2)==0)		
                            graph.setWeight(atomNo1,atomNo2,1);		
                    if(graph.getWeight(atomNo2,atomNo3)==0)		
                            graph.setWeight(atomNo2,atomNo3,1);
                }
                else{
                    graph.setWeight(tmpRing->getAtomNo(0),tmpRing->getAtomNo(1),2);
                    graph.setWeight(tmpRing->getAtomNo(1),tmpRing->getAtomNo(2),1);
                    graph.setWeight(tmpRing->getAtomNo(2),tmpRing->getAtomNo(3),2);
                    graph.setWeight(tmpRing->getAtomNo(3),tmpRing->getAtomNo(4),1);
                    graph.setWeight(tmpRing->getAtomNo(4),tmpRing->getAtomNo(5),2);
                    graph.setWeight(tmpRing->getAtomNo(5),tmpRing->getAtomNo(0),1);
                }
            }
        }
    }
    
   
    //5 membered ring
   /* for(int i=0;i<ringVect.size();i++){
        Ring *tmpRing=ringVect[i];
        unsigned atomNo;
        if(tmpRing->size()==5){
            for(int j=0;j<tmpRing->size();j++){
                
                atomNo = tmpRing->getAtomNo(j);
                vector<unsigned> tmpVect;
                int degree=graph.getDegree(atomNo);        
                int valency=ligand[atomNo].getValency();
                graph.getAdjacencyList(atomNo,tmpVect);
                int nBond=0;//total bonds made
                int wt=0;
                int nVal=0; //total valency satisfied
                for(int k=0;k<tmpVect.size();k++){
                 wt=graph.getWeight(atomNo,tmpVect[k]);
                    if(wt!=0)
                        nBond++;
                    nVal += wt;  
                }
                
                //Carbon C
                if(ligand[atomNo].getAtomSymbol()=="C"){
                    //if degree is 3 with 1 double bond and 1 single
                    if(degree==3 && nVal==3){
                        for(int k=0;k<tmpVect.size();k++){
                            if(tmpRing->isInRing(tmpVect[k]) &&
                                (graph.getWeight(atomNo,tmpVect[k])==0))
                               graph.setWeight(atomNo,tmpVect[k],1);
                        }
                    }
                    
                    //if one double bond found 
                    if(degree==3 && nVal==2){
                        for(int k=0;k<tmpVect.size();k++){
                            if(tmpRing->isInRing(tmpVect[k]) &&
                                (graph.getWeight(atomNo,tmpVect[k])==0))
                               graph.setWeight(atomNo,tmpVect[k],1);
                        }
                    }
                }
                
                //Nitrogen N
                else if(ligand[atomNo].getAtomSymbol()=="N"){
                    if(degree==2 && (nVal==1||nVal==2)){
                        for(int k=0;k<tmpVect.size();k++){
                            if(tmpRing->isInRing(tmpVect[k]) &&
                                (graph.getWeight(atomNo,tmpVect[k])==0))
                               graph.setWeight(atomNo,tmpVect[k],valency-nVal);
                        }
                    }      
                    
                } 
               
            }
        }        
    }
 //   graph.print();
  */
    for(unsigned i=1;i<=ligand.size();i++){
        vector<unsigned> tmpVect;
        int valency=ligand[i].getValency();
	cout<<ligand[i].getAtomSymbol()<<"\n";
        int nVal=0;
        int diffVal;
        int degree=graph.getDegree(i);
        int wt; 
        int nBond=0;
        graph.getAdjacencyList(i,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            wt=graph.getWeight(i,tmpVect[j]);
            if(wt!=0)
                nBond++;
            nVal += wt;            
        }       
     
	//if only one bond remains
        diffVal=valency-nVal;
        if(degree-nBond ==1){
            for(unsigned j=0;j<tmpVect.size();j++){
     
                if(graph.getWeight(i,tmpVect[j])==0){
                    graph.setWeight(i,tmpVect[j],diffVal);
                }
            }
        }
	
        else if((degree-nBond)==diffVal){
            for(unsigned j=0;j<tmpVect.size();j++){     
                if(graph.getWeight(i,tmpVect[j])==0){
                    graph.setWeight(i,tmpVect[j],1);
                }
            } 
        }
    }
    
    for(unsigned i=1;i<=ligand.size();i++){
        vector<unsigned> tmpVect;
        int valency=ligand[i].getValency();
        int nVal=0;
        int diffVal;
        int degree=graph.getDegree(i);
        int wt; 
        int nBond=0;
        graph.getAdjacencyList(i,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            wt=graph.getWeight(i,tmpVect[j]);
            if(wt!=0)
                nBond++;
            nVal += wt;            
        }       
     
        diffVal=valency-nVal;
        if(degree-nBond ==1){
            for(unsigned j=0;j<tmpVect.size();j++){
     
                if(graph.getWeight(i,tmpVect[j])==0){
                    graph.setWeight(i,tmpVect[j],diffVal);
                }
            }
        }
    }

    for(unsigned i=1;i<=ligand.size();i++){
        vector<unsigned> tmpVect;
        int valency=ligand[i].getValency();
        int nVal=0;
        int diffVal;
        int degree=graph.getDegree(i);
        int wt; 
        int nBond=0;
        graph.getAdjacencyList(i,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            wt=graph.getWeight(i,tmpVect[j]);
            if(wt!=0)
                nBond++;
            nVal += wt;            
        }       
     
        diffVal=valency-nVal;
        if(degree-nBond ==1){
            for(unsigned j=0;j<tmpVect.size();j++){
     
                if(graph.getWeight(i,tmpVect[j])==0){
                    graph.setWeight(i,tmpVect[j],diffVal);
                }
            }
        }
    }
   
   //print warnings if valency not satisfied
   for(unsigned i=1;i<=ligand.size();i++){
        string symbol=ligand[i].getAtomSymbol();
        int valency=ligand[i].getValency();
        int nBO=getTotalBondOrder(graph,i);
        /*if(symbol=="S"){
            if(nBO!=2 && nBO!=6)
                fprintf(fa,"\nWARNING! ATOM %d %s VALENCY NOT SATISFIED EXPECTED 2 OR 6   SHOWN %d",i,symbol.c_str(),nBO);
	}
        else if(nBO != valency){
                fprintf(fa,"\nWARNING! ATOM %d %s VALENCY NOT SATISFIED EXPECTED %d       SHOWN %d",i,symbol.c_str(),valency,nBO);
        }*/    
        if( (nBO != valency && symbol=="C") || nBO==0){
		if(symbol=="C" && nBO == 3){
		    vector<unsigned> adj;
		    graph.getAdjacencyList(i,adj);
		    for(unsigned j=0;j<adj.size();j++)
			if(ligand[adj[j]].getAtomSymbol()=="N"){
			    graph.setWeight(i,adj[j],2);
			    break;
			}
		}
		else{
		    
		    fprintf(fa,"\nWARNING! ATOM %d %s VALENCY NOT SATISFIED EXPECTED %d       SHOWN %d",i,symbol.c_str(),valency,nBO);
		    errorFlag=1;
		}
	}
   }
  
   //assigning 4 to resonating bonds in rings if makeResonating is set true
   if(makeResonating){
	  for(unsigned i=0;i<ringVect.size();i++){
	      
	    Ring *tmpRing=ringVect[i];
	    if(ringVect[i]->isPlanar()){
		if(isResonatingRing(graph,ligand,tmpRing)){
		    if(ringVect[i]->size()==6){
			graph.setWeight(tmpRing->getAtomNo(0),tmpRing->getAtomNo(1),4);
			graph.setWeight(tmpRing->getAtomNo(1),tmpRing->getAtomNo(2),4);
			graph.setWeight(tmpRing->getAtomNo(2),tmpRing->getAtomNo(3),4);
			graph.setWeight(tmpRing->getAtomNo(3),tmpRing->getAtomNo(4),4);
			graph.setWeight(tmpRing->getAtomNo(4),tmpRing->getAtomNo(5),4);
			graph.setWeight(tmpRing->getAtomNo(5),tmpRing->getAtomNo(0),4);
		    }
		    else if(ringVect[i]->size()==5 ){
			graph.setWeight(tmpRing->getAtomNo(0),tmpRing->getAtomNo(1),4);
			graph.setWeight(tmpRing->getAtomNo(1),tmpRing->getAtomNo(2),4);
			graph.setWeight(tmpRing->getAtomNo(2),tmpRing->getAtomNo(3),4);
			graph.setWeight(tmpRing->getAtomNo(3),tmpRing->getAtomNo(4),4);
			graph.setWeight(tmpRing->getAtomNo(4),tmpRing->getAtomNo(0),4);

		    }
		}
	    }
	}
       
       //assigning 4 to aliphatic bonds 
       /*for(unsigned i=1;i<=ligand.size();i++){
	    if(!ligand[i].isInRing()){
		 if(isAmideCarbon(graph,ligand,i)>1){
		    vector<unsigned> tmpVect;
		    graph.getAdjacencyList(i,tmpVect);
		    for(unsigned j=0;j<tmpVect.size();j++){
		       if(ligand[tmpVect[j]].getAtomSymbol()=="N")
			    graph.setWeight(i,tmpVect[j],4);	   
		    }
		 }	
		 if(hasResonatingOxygenBond(graph,ligand,i)>1){
		    vector<unsigned> tmpVect;
		    graph.getAdjacencyList(i,tmpVect);
	    //		cout << i<<endl;
		    for(unsigned j=0;j<tmpVect.size();j++){
		      //cout << ligand[tmpVect[j]].getAtomSymbol() <<endl;
		      if(ligand[tmpVect[j]].getAtomSymbol()=="O"){
			    graph.setWeight(i,tmpVect[j],4);
		      }
		    }	   
		 }
	    }
       }
       */
   } 
   graph.print();
   if(!errorFlag){   
      ligand.write(argv[2]);
      appendConnectInfo(graph,argv[2]);
   }
   if(errorFlag){  
      char fname[50];
      strcpy(fname,argv[2]);
      //strcat(fname,".err");     //Change by Goutam so that no error file will be generated
      ligand.write(fname);
      appendConnectInfo(graph,fname);
      FILE *ef;
      ef = fopen("error.log","a");
      fprintf(ef,"%s\n",argv[1]);
      fclose(ef);
   }
   for(unsigned i=0;i<ringVect.size();i++)
	   delete ringVect[i];

   fprintf(fa,"\nCONNECTION DONE. \n");
   fclose(fa);

   //cout << "Connection Done. Output written in output.pdb."<< endl;
   return 0;

}

void appendConnectInfo(Graph<unsigned>& graph,char *fileName){

    FILE *fp;
    fp = fopen(fileName,"a+");
    
    for (unsigned i=0;i<graph.size();i++){
        vector<unsigned> tmpVect; 
        
        graph.getAdjacencyList(i+1,tmpVect);
    
        char line[80]="";
        char strNum[10]="";
        strcat(line,"CONNECT  ");     //change CONECT to CONNECT by Goutam
        sprintf(strNum,"%u   ",i+1);
        strcat(line,strNum);
        for(unsigned j=0;j<tmpVect.size();j++){   
       	    char  strNum1[10]="";		
            sprintf(strNum1,"%u  ",tmpVect[j]);
            strcat(line,strNum1);
        }
        fprintf(fp,"%s\n",line);    
    }
    fprintf(fp,"END\n");
    
    for (unsigned i=0;i<graph.size();i++){
        vector<unsigned> tmpVect;
                
	graph.getAdjacencyList(i+1,tmpVect);
        
        for(unsigned j=0;j<tmpVect.size();j++){   
	    if(i+1<tmpVect[j]){
		char line[80]="";
		char  strNum1[80]="";
		sprintf(strNum1,"%03u%03u%03u   0  0  0",i+1,tmpVect[j],
			graph.getWeight(i+1,tmpVect[j]));
		strcat(line,strNum1);
		fprintf(fp,"%s\n",line);  
	    }	
        }
    }
    
	fprintf(fp,"EOF");
    fclose(fp);
}

bool isAromatic(Graph<unsigned>& graph,Pdb& pdb,Ring* ring){
  
    for(int i=0;i<ring->size();i++)
    {
	unsigned atomNo = ring->getAtomNo(i);
        int degree=graph.getDegree(atomNo);
        string symbol=pdb[atomNo].getAtomSymbol();
	//cout<< symbol<<"\n";
        if((symbol=="C" && degree!=3  )){
//            cout <<atomNo<<" "<< symbol << "  " <<degree <<
//                 "  "<<getTotalBondOrder(graph,atomNo)<<endl;
            return false;
	}
        else if((symbol=="N" && degree!=2 ) ){
	    return false;
	}
    }
    return true;
}

//check for ring C attached to O with double bond
bool isResonatingRing(Graph<unsigned>& graph,Pdb& pdb,Ring* ring){
  
    if(ring->isPlanar()){
       for(int i=0;i<ring->size();i++){
	  unsigned atomNo = ring->getAtomNo(i);
       	  vector<unsigned> tmpVect;
          graph.getAdjacencyList(atomNo,tmpVect);
	  string symbol=pdb[atomNo].getAtomSymbol();
	  int degree = graph.getDegree(atomNo);
	  if((symbol=="C")){
	     if(degree!=3)
	   	return false;	     
             for(unsigned j=0;j<tmpVect.size();j++)
                if(pdb[tmpVect[j]].getAtomSymbol()=="O"&&
		   graph.getWeight(atomNo,tmpVect[j])==2)
			return false;	
	  }		
       }
    }
    return true;	    
}

/*
bool isResonating6(Graph<unsigned>& graph,Pdb& pdb,Ring* ring){
    for(int i=0;i<ring->size();i++)
    {
	unsigned atomNo = ring->getAtomNo(i);
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        int degree=graph.getDegree(atomNo);
        string symbol=pdb[atomNo].getAtomSymbol();
        if((symbol=="C" && degree!=3  )){
//            cout <<atomNo<<" "<< symbol << "  " <<degree <<
//                 "  "<<getTotalBondOrder(graph,atomNo)<<endl;
            return false;
	}
        else if((symbol=="N" && degree!=2 ) ){
	    return false;
	}
    }
    return true;
}
*/

//returns no. of valency satisfied
int getTotalBondOrder(Graph<unsigned>& graph,unsigned atomNo){
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            nVal += graph.getWeight(atomNo,tmpVect[j]);
        }       
    	return nVal; 
}

//returns no of single bonds formed
int getNSingleBond(Graph<unsigned>& graph,unsigned atomNo){
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            if(graph.getWeight(atomNo,tmpVect[j])==1);
		    nVal++;
        }       
    	return nVal; 
}

//returns no of double bonds formed
int getNDoubleBond(Graph<unsigned>& graph,unsigned atomNo){
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            if(graph.getWeight(atomNo,tmpVect[j])==2);
		    nVal++;
        }       
    	return nVal; 
}

//returns no of triple bonds formed
int getNTripleBond(Graph<unsigned>& graph,unsigned atomNo){
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            if(graph.getWeight(atomNo,tmpVect[j])==3);
		    nVal++;
        }       
    	return nVal; 
}

//returns the no of Oxygens bound to atomNo if a charge O- oxygen is found else 
//return 0 (for COO-ve & N02)
int hasResonatingOxygenBond(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo){
        vector<unsigned> tmpVect;
        int nVal=0,flag=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            if(ligand[tmpVect[j]].getAtomSymbol()=="O"){
		nVal++;
		//check for -ve charged O
                if(graph.getDegree(tmpVect[j])==1 &&
                   graph.getWeight(atomNo,tmpVect[j])==1)
                   flag=1;
            }
        }       
        if(flag==1)
            return nVal;
    	return 0;  
}

//return true if the atom no is N of amide group NH2 
bool isAmideNitrogen(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo){
        if(ligand[atomNo].getAtomSymbol()!="N")
            return false;
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++){
            if(ligand[tmpVect[j]].getAtomSymbol()=="H"){
		nVal++;
            }
        }       
        if (nVal==2)
            return true;
        return false;    
}

//return no of NH2 attached to the atom no C 
int isAmideCarbon(Graph<unsigned>& graph,Pdb& ligand,unsigned atomNo){
        if(ligand[atomNo].getAtomSymbol()!="C")
            return 0;
        vector<unsigned> tmpVect;
        int nVal=0;
        graph.getAdjacencyList(atomNo,tmpVect);
        for(unsigned j=0;j<tmpVect.size();j++)
            if(isAmideNitrogen(graph,ligand,tmpVect[j]))
		nVal++;
        return nVal;    
}

bool isValidBondLength(float dist, float ri, float rj){
    float dlimit=0.0;
    float radius=ri+rj;
    if(dist<=1.5)
	dlimit=radius+0.15;
    else if((dist<=1.9 && dist>1.5))
	dlimit=radius+0.11;
    else if((dist<=2.05 && dist>1.9))
	dlimit=radius+0.09;
    else if(dist>2.05)
	dlimit=radius+0.08;

    if(dist>(radius/2)){
	if(dist<dlimit)
	{
	    return true;
	}
    }
    return false;
}    
