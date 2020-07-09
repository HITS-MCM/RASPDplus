#include "Pdb.hpp"
#include "cstring"

void Pdb::read(char* filename)
{
	PdbRecord pdbNode; ///one unit
	unsigned count=0;
      	char name[10];
	char store[500];
        FILE *fp1;
        fp1=fopen(filename,"r");
	//unsigned atomNo,residueNo,
        while(fgets(store,500,fp1)!=NULL)
        {
		sscanf(store,"%6s",name);
                if(strcmp(name,"TER")==0)
			break;
			
		if(pdbNode.read(store)){
			//graph.addVertex(pdbNode);
			pdbVector.push_back(pdbNode);
        	        count++; ///total atoms
		}
        }
        fclose(fp1);
	sz=count;
	pdbVector.resize(sz);
}
void Pdb::read(string& filename)
{
	const char* flName;
	flName = filename.c_str();
	PdbRecord pdbNode; ///one unit
	unsigned count=0;
      	char name[10];
	char store[120];
        FILE *fp1;
        fp1=fopen(flName,"r");
	//unsigned atomNo,residueNo,
        while(fgets(store,120,fp1)!=NULL)
        {
		sscanf(store,"%6s",name);
		if(strcmp(name,"TER")==0)
			break;
			
		if(pdbNode.read(store)){
			//graph.addVertex(pdbNode);
			pdbVector.push_back(pdbNode);
        	        count++; ///total atoms
		}
        }
        fclose(fp1);
	sz=count;
	pdbVector.resize(sz);
}
/*
void Pdb::read(char* filename)
{
	PdbRecord pdbNode; ///one unit
	unsigned count=0;
        char name[10],store[120];
        FILE *fp1;
        fp1=fopen(filename,"r");
        while(fgets(store,120,fp1)!=NULL)
        {
		sscanf(store,"%s",name);
                if(strcmp(name,"ATOM")==0)
                {
			pdbNode.read(store);
			pdbVector.push_back(pdbNode);
			//graph.addVertex(pdbNode);
                        count++; ///total atoms
                }
        }
        fclose(fp1);
//	graph.buildIndex();
	sz=count;
	pdbVector.resize(sz);
}
*/
void Pdb::write(char* filename)
{
	FILE *fpw;
        char store[120];
        fpw=fopen(filename,"w");
	assert(fpw!=NULL);
        for(unsigned i=1;i<=size();i++) ///i=atomno
	{
        	(*this)[i].write(store); ///PdbRecord write function
//                graph[i].print(); ///PdbRecord print function
                fputs(store,fpw);
        }
	fclose(fpw);
}
void Pdb::write(string& filename)
{
	FILE *fpw;
        char store[120];
	const char* flName;
	flName = filename.c_str();
        fpw=fopen(flName,"w");
	assert(fpw!=NULL);
        for(unsigned i=1;i<=size();i++) ///i=atomno
	{
        	(*this)[i].write(store); ///PdbRecord write function
//                graph[i].print(); ///PdbRecord print function
                fputs(store,fpw);
        }
	fclose(fpw);
}
