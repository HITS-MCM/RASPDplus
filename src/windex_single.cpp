/*
#>  Copyright (c) 2013 2014 2015 2016 2017 2018 2019 2020
#>  Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
#>
#>  Supercomputing Facility for Bioinformatics and Computational Biology, IIT Delhi (http://www.scfbio-iitd.res.in/)
#>  Indian Institute of Technology
#>  Hauz Khas, New Delhi - 110016, India
#>
#>  Please send your contact address to get information on updates and
#>  new features to "mcmsoft@h-its.org". Questions will be
#>  answered as soon as possible.
#>  References:
#>  A rapid identification of hit molecules for target proteins via physico-chemical descriptors.
#>  (2013) Phys. Chem. Chem. Phys., 15, 9107-9116.
#>  Authors: Goutam Mukherjee and B. Jayaram
#>  Version 1.0 (April 2013)

#>  RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features
#>  Authors: Stefan Holderbach, Lukas Adam, B. Jayaram, Rebecca C. Wade, Goutam Mukherjee
#>  Version 1.0 (June 2020)
#>  Manuscript in preparation, 2020

*To calculate the Wiener index of a molecule from a PDB file.
WORKING:
      1.Searches for lines starting with "ATOM" or "HETATM" to calculate the number of atom in molecule
      2.Creates a graph with number of nodes = number of atoms.
      3.Searches for lines starting with "CONECT" and then adds edges between relevant atoms using functions processline()  and addEdge().
      4.Stores Bond Order between atoms in array bond[][];
      5.Runs Dijkstra's algorithm to establish shortest path between atoms and stores this information in array windex[][]
      6.Calculates and prints the Wiener index

USAGE : Usage: ./windex_list.exe <PDB file>
WARNINGS: Program takes a specific format of PDB file, example file appended at the end
DEPENDENCIES:PDB file
Author: Abhinav Singh      Created: 06-07-2014     Modified: 10-07-2014
Status: Completed
**/


#include<algorithm>
#include<unistd.h>
#include<fstream>
#include<stdio.h>
//#include <conio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <iostream>
#include <string>
#include<float.h>
#include<string.h>
#define SBLENGTH 2147483647

using namespace std;

class Atmnolist
{
    public:
    long int Z;
    const char *atm;
};
	Atmnolist atmlist[19];


typedef struct AdjListNode
{
    long int dest;
    struct AdjListNode* next;
}AdjListNode;


typedef struct AdjList
{
     AdjListNode *head;
}AdjList;


typedef struct Graph
{
    long int total,V;
    AdjList* arrays;
}Graph;

AdjListNode* newAdjListNode(long int dest)
{
    AdjListNode* newNode =( AdjListNode*) calloc(1,sizeof( AdjListNode));
    if (newNode == NULL)
    {
        fprintf(stderr, "mem allocation error\n");
        exit(EXIT_FAILURE);
    }
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

Graph* createGraph(long int v)
{

    Graph* graph = (Graph*) calloc(1,sizeof(Graph));
    if (graph == NULL) {
  fprintf(stderr, "mem allocation error\n");
  exit(EXIT_FAILURE);
}
    graph->V = v;


    graph->arrays = (AdjList*) calloc(v, sizeof(AdjList));
    if (graph->arrays == NULL) {
  fprintf(stderr, "mem allocation error\n");
  exit(EXIT_FAILURE);
}


    long int i;
    for (i = 0; i < v; ++i)
        graph->arrays[i].head = NULL;
    return graph;
}


void addEdge(Graph* graph, long int src, long int dest)
{
    AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->arrays[src].head;
    graph->arrays[src].head = newNode;

}


long int minDistance(Graph *graph,long double *dist, bool *sptSet)
{

   long double min = SBLENGTH;
   long int min_index,v;
   for (v = 0; v < graph->V; v++)
     {
         if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
     }

   return min_index;
}

long double* calcweight(long double *weight,long int u,Graph *graph,long int *atmno,long int **bond)
{
    long int m;

    AdjListNode *p;
    p=(graph->arrays[u]).head;
    double temp;
    while(p)
    {
        m=p->dest;

        if(bond[u][m]!=4)
        {
            weight[m]=(((double)(6*6)/(double)(atmno[u]*atmno[m]))/bond[u][m]);
        }

        else
        {
            weight[m]=(((double)(6*6)/(double)(atmno[u]*atmno[m]))/1.5);
        }
        p=p->next;

    }


    return weight;
}

long double printSolution(long double *dist,Graph *graph,long int src,long double **windex)
{
    long int i,j;
   for (i = 0; i <graph->V; i++)
      {
          windex[src][i]=dist[i];
      }


}


void dijkstra(Graph *graph, long int src, long double **windex, long int *atmno,long int **bond)
{
    long int i=0,counter=0,j=0,u=0;
    long double *weight=NULL;
    weight= (long double*) calloc ((graph->V),sizeof(long double));
    std::fill_n(weight,graph->V,0);

    long double *dist = (long double*) calloc((graph->V),(sizeof(long double)));     // The output arrays.  dist[i] will hold the shortest
                                                                            // distance from src to i
    std::fill_n(dist,graph->V,0);

    bool *sptSet =(bool*) calloc((graph->V),(sizeof(bool))); // sptSet[i] will true if vertex i is included in shortest
                                                             // path tree or shortest distance from src to i is finalized

    std::fill_n(sptSet,graph->V,0);
                                                             // Initialize all distances as INFINITE and stpSet[] as false
    for (i = 0; i < graph->V; i++)
        {
            dist[i] = LDBL_MAX, sptSet[i] = false;
            weight[i]=0;
        }

                                                            // Distance of source vertex from itself is always 0
    dist[src] = 0;

                                                            // Find shortest path for all vertices
    for (counter = 0; counter < (graph->V)-1; counter++)
    {

                                                            // Pick the minimum distance vertex from the set of vertices not

                                                            // yet processed. u is always equal to src in first iteration.
        u = minDistance(graph,dist, sptSet);

        weight = calcweight(weight,u,graph,atmno,bond);

                                                            // Mark the picked vertex as processed
        sptSet[u] = true;


                                                            // Update dist value of the adjacent vertices of the picked vertex.
        for (j = 0; j < graph->V; j++)

        {                                                   // Update dist[j] only if is not in sptSet, there is an edge from
                                                            // u to j, and total weight of path from src to  j through u is
                                                            // smalle0r than current value of dist[j]

            if (!sptSet[j] && weight[j] && dist[u] != LDBL_MAX&& dist[u]+weight[j] < dist[j])
            {
                dist[j] = dist[u] + weight[j];

            }

         }
     }


     printSolution(dist,graph,src,windex);
}

long double** make2darray(long int );
long int** make2iarray(long int );
void processline (string,string , Graph *,long int []);

long int getatmno (char *);


int main(int argc, char* argv[])
{

        atmlist[0].atm="H";atmlist[0].Z=1;
        atmlist[1].atm="B";atmlist[1].Z=5;
        atmlist[2].atm="C";atmlist[2].Z=6;
        atmlist[3].atm="N";atmlist[3].Z=7;
        atmlist[4].atm="O";atmlist[4].Z=8;
        atmlist[5].atm="F";atmlist[5].Z=9;
        atmlist[6].atm="Na";atmlist[6].Z=11;
        atmlist[7].atm="Mg";atmlist[7].Z=12;
        atmlist[8].atm="P";atmlist[8].Z=15;
        atmlist[9].atm="S";atmlist[9].Z=16;
        atmlist[10].atm="Cl";atmlist[10].Z=17;
        atmlist[11].atm="K";atmlist[11].Z=19;
        atmlist[12].atm="Ca";atmlist[12].Z=20;
        atmlist[13].atm="Br";atmlist[13].Z=35;
        atmlist[14].atm="I";atmlist[14].Z=53;
        atmlist[15].atm="BR";atmlist[15].Z=35;
        atmlist[16].atm="NA";atmlist[16].Z=11;
        atmlist[17].atm="CL";atmlist[17].Z=17;
        atmlist[18].atm="CA";atmlist[18].Z=20;


        Graph *graph=NULL;
        char type[10]={'a','a','a','a','a','a','a','a','a','a'},atm[30]={'a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a'},tempchar1[4]={'a','a','a','a'},tempchar2[4]={'a','a','a','a'},tempchar3[4]={'a','a','a','a'},wastechar[200];
        long int *atmno=NULL,*Hpresent=NULL,var[8],breakflag=0,foundend=0,foundcnt=0,bo=0,src=0,dest=0,total=0,V=0,j=0,i=0,z=1,mno=0,xyz=0,waste;
        long int tempint1=0,tempint2=0;
        long double wi=0;
        fstream myfile,eol,pdblist;
        const char *c,*d;

        std::fill_n(wastechar,200,'a');
        string line,str="CONNECT",filename;
        std::fill_n(var,8, 32000);

                myfile.open(argv[1],ios::in);
                if(myfile.is_open())
                {
                    while (getline (myfile,line))
                    {
                        long int foundatm = line.find("ATOM");
                        long int foundhetatm = line.find("HETATM");

                        if(foundatm==0||foundhetatm==0)
                        {
                            V++;

                        }

                        if (foundatm!=0&&foundhetatm!=0)
                        {
                            break;
                        }
                    }

                myfile.close();
                }

                else
                {
                    cout<<"file not found \n";
                    return 0;
                }


                long int **bond = make2iarray(V);
                for(i=0;i<V;i++)
                {
                    std::fill_n(bond[i],V,1);
                }
                long double **windex = make2darray(V);
                atmno = (long int*) calloc(V,sizeof(long int));

                std::fill_n(atmno,V,0);


                Hpresent = new long int[V*sizeof(long int)];
                std::fill_n(Hpresent,V,0);


                myfile.open( argv[1],ios::in);
                if (myfile.is_open())
                    {
                        while (getline (myfile,line))
                        {

                            long int foundatm = line.find("ATOM");
                            long int foundhetatm = line.find("HETATM");

                            if(foundatm==0||foundhetatm==0)
                            {
                                breakflag=1;


                                c = line.c_str();
                                sscanf (c,"%s %ld %*s",type,&src);

                                char chars[] = "()-0123456789.";

                                for (unsigned int i = 0; i < strlen(chars); ++i)
                                {
                                    line.erase (std::remove(line.begin(), line.end(), chars[i]), line.end());
                                }

                                sscanf (c,"%s %s %*s",type,atm);

                                if(strcmp(atm,"H")!=0)
                                {


                                    atmno[src-1] = getatmno(atm);

                                }

                                else
                                {
                                    atmno[src-1] = getatmno(atm);
                                    Hpresent[src-1]= 1;
                                }


                            }

                            else if (foundatm!=0&&foundhetatm!=0)
                            {
                                break;
                            }

                        }




                        graph = createGraph(V);

                        processline (line,str,graph,var);

                        while (getline (myfile,line))
                        {
                            foundcnt = line.find(str);
                            long int foundend =line.find("END");
                            if (foundcnt!=std::string::npos)
                            {

                            usleep(1);
                            processline (line,str,graph,var);
                            }


                            else
                            {
                                int rubbish;

                                std::fill_n(tempchar3,4,'a');

                                //ASSUMING BOND ORDER IS RIGHT AFTER 'CONECT'
                                const char *c = line.c_str();
                                sscanf(c,"%3c %3c %3c %*c",tempchar1,tempchar2,tempchar3);
                                src = atol(tempchar1);
                                dest=atol(tempchar2);
                                bo= atol(tempchar3);

                                rubbish=bond[dest-1][src-1];
                                bond[src-1][dest-1]=bo;
                                bond[dest-1][src-1]=bo;

                                while (getline (myfile,line))
                                {
                                    std::fill_n(tempchar3,4,'a');

                                    const char *c = line.c_str();

                                    sscanf(c,"%3c %3c %3c %*c",tempchar1,tempchar2,tempchar3);
                                    src = atol(tempchar1);
                                    dest=atol(tempchar2);
                                    bo= atol(tempchar3);
                                    rubbish=bond[dest-1][src-1];
                                    bond[src-1][dest-1]=bo;
                                    bond[dest-1][src-1]=bo;
                                }

                            }
                        }

                    }


                myfile.close();


                for (i=0;i<V;i++)
                {

                    dijkstra(graph,i,windex,atmno,bond);

                }
                for(i=0;i<V;i++)
                {

                    if(Hpresent[i]!=1)
                    {
                        windex[i][i]=(float)1-((float)6/(float)atmno[i]);


                        for(j=i;j<V;j++)
                        {
                            if(Hpresent[j]!=1)
                            {
                                wi=windex[i][j]+wi;
                            }
                        }
                    }
                    else
                    {
                        windex[i][i]=0;
                    }

                }
                cout<<filename<<"   "<<wi<<"\n";

return 0;
}



//method to add edges between required nodes after checking if line has word 'CONECT' in it.
void processline (string line,string str, Graph *graph,long int var[])
{

    long int src=0,i=0;
    char type[10];
    std::fill_n(type,10,'a');
    long int foundcnt = line.find(str);
    if (foundcnt!=std::string::npos)
            {
                std::fill_n(var,8, 32000);
                const char *c = line.c_str();
                sscanf (c,"%s %ld  %ld  %ld  %ld  %ld  %ld  %ld  %ld  %ld",type,&src,&var[0],&var[1],&var[2],&var[3],&var[4],&var[5],&var[6],&var[7]);
                i=0;
                while(var[i]!=32000)
                    {
                        addEdge(graph, src-1,var[i]-1);
                        i++;
                    }
            }

}


//method to obtain atomic number from symbol read from file

long int getatmno (char *c)
{
    int i=0;
    for (i=0;i<19;i++)
    {
        if(strcmp(c,atmlist[i].atm)==0)
        {
            return (atmlist[i].Z);

        }

    }

    cout<<"some atom not in list check source code";
    exit(-1);

}

long double** make2darray(long int V)
{
    long double **array2;
    long int i;
        array2 = (long double**) calloc (V,sizeof(long double*));
    if (array2== NULL)
        {
            fprintf(stderr, "mem allocation error\n");
            exit(EXIT_FAILURE);
        }
    for(i=0;i<V;i++)
    {

        array2[i] = (long double*) calloc (V,sizeof(long double));
        if (array2[i]== NULL)
        {
            fprintf(stderr, "mem allocation error\n");
            exit(EXIT_FAILURE);
        }
        std::fill_n(array2[i],V,0);
    }
    return array2;

}

long int** make2iarray(long int V)
{
    long int **array2,i;
        array2 = (long int**) calloc (V,sizeof(long int*));
    if (array2== NULL)
        {
            fprintf(stderr, "mem allocation error\n");
            exit(EXIT_FAILURE);
        }
    for(i=0;i<V;i++)
    {

        array2[i] = (long int*) calloc (V,sizeof(long int));
        if (array2[i]== NULL)
        {
            fprintf(stderr, "mem allocation error\n");
            exit(EXIT_FAILURE);
        }
        std::fill_n(array2[i],V,0);
    }
    return array2;

}




/**
A sample PDB file  0HK.pdb  its Wiener index is 1286.99
ATOM      1  C   DRG     1      -0.305   0.285  -0.555 0 .00  0.0000  0.0000  0.0000
ATOM      2  C   DRG     1      -2.394   1.171   0.400 0 .00  0.0000  0.0000  0.0000
ATOM      3  C   DRG     1      -2.007  -1.259   0.327 0 .00  0.0000  0.0000  0.0000
ATOM      4  C   DRG     1      -2.345   2.349   1.023 0 .00  0.0000  0.0000  0.0000
ATOM      5  C   DRG     1      -3.414   3.190   0.765 0 .00  0.0000  0.0000  0.0000
ATOM      6  C   DRG     1      -4.331   2.694  -0.068 0 .00  0.0000  0.0000  0.0000
ATOM      7  C   DRG     1      -1.593  -2.398   0.884 0 .00  0.0000  0.0000  0.0000
ATOM      8  C   DRG     1      -2.354  -3.510   0.571 0 .00  0.0000  0.0000  0.0000
ATOM      9  C   DRG     1      -3.386  -3.275  -0.240 0 .00  0.0000  0.0000  0.0000
ATOM     10  C   DRG     1       5.128   0.447  -1.493 0 .00  0.0000  0.0000  0.0000
ATOM     11  C   DRG     1       3.741  -1.076  -0.183 0 .00  0.0000  0.0000  0.0000
ATOM     12  C   DRG     1       2.905  -1.117  -1.473 0 .00  0.0000  0.0000  0.0000
ATOM     13  C   DRG     1       1.891   0.030  -1.444 0 .00  0.0000  0.0000  0.0000
ATOM     14  C   DRG     1       2.618   1.356  -1.201 0 .00  0.0000  0.0000  0.0000
ATOM     15  C   DRG     1       3.451   1.226   0.084 0 .00  0.0000  0.0000  0.0000
ATOM     16  C   DRG     1       2.574   0.586   1.185 0 .00  0.0000  0.0000  0.0000
ATOM     17  C   DRG     1       2.765  -0.935   1.008 0 .00  0.0000  0.0000  0.0000
ATOM     18  C   DRG     1       5.505   0.198   0.880 0 .00  0.0000  0.0000  0.0000
ATOM     19  C   DRG     1      -1.353   0.086   0.511 0 .00  0.0000  0.0000  0.0000
ATOM     20  N   DRG     1       4.497   0.205  -0.189 0 .00  0.0000  0.0000  0.0000
ATOM     21  O   DRG     1      -0.580   0.884  -1.567 0 .00  0.0000  0.0000  0.0000
ATOM     22  O   DRG     1       3.219  -0.230   2.165 0 .00  0.0000  0.0000  0.0000
ATOM     23  O   DRG     1       0.934  -0.198  -0.376 0 .00  0.0000  0.0000  0.0000
ATOM     24  O   DRG     1      -0.738   0.146   1.799 0 .00  0.0000  0.0000  0.0000
ATOM     25  S   DRG     1      -3.854   1.079  -0.576 0 .00  0.0000  0.0000  0.0000
ATOM     26  S   DRG     1      -3.438  -1.565  -0.649 0 .00  0.0000  0.0000  0.0000
ATOM     27  H   DRG     1      -1.534   2.627   1.679 0 .00  0.0000  0.0000  0.0000
ATOM     28  H   DRG     1      -3.500   4.173   1.204 0 .00  0.0000  0.0000  0.0000
ATOM     29  H   DRG     1      -5.229   3.205  -0.382 0 .00  0.0000  0.0000  0.0000
ATOM     30  H   DRG     1      -0.729  -2.448   1.531 0 .00  0.0000  0.0000  0.0000
ATOM     31  H   DRG     1      -2.129  -4.495   0.955 0 .00  0.0000  0.0000  0.0000
ATOM     32  H   DRG     1      -4.085  -4.020  -0.590 0 .00  0.0000  0.0000  0.0000
ATOM     33  H   DRG     1      -0.297   0.985   1.987 0 .00  0.0000  0.0000  0.0000
ATOM     34  H   DRG     1       5.597   1.431  -1.492 0 .00  0.0000  0.0000  0.0000
ATOM     35  H   DRG     1       5.885  -0.315  -1.678 0 .00  0.0000  0.0000  0.0000
ATOM     36  H   DRG     1       4.372   0.406  -2.276 0 .00  0.0000  0.0000  0.0000
ATOM     37  H   DRG     1       4.388  -1.947  -0.085 0 .00  0.0000  0.0000  0.0000
ATOM     38  H   DRG     1       3.561  -1.003  -2.336 0 .00  0.0000  0.0000  0.0000
ATOM     39  H   DRG     1       2.379  -2.070  -1.537 0 .00  0.0000  0.0000  0.0000
ATOM     40  H   DRG     1       1.366   0.075  -2.398 0 .00  0.0000  0.0000  0.0000
ATOM     41  H   DRG     1       3.273   1.594  -2.039 0 .00  0.0000  0.0000  0.0000
ATOM     42  H   DRG     1       1.878   2.148  -1.078 0 .00  0.0000  0.0000  0.0000
ATOM     43  H   DRG     1       3.878   2.179   0.395 0 .00  0.0000  0.0000  0.0000
ATOM     44  H   DRG     1       1.614   1.052   1.410 0 .00  0.0000  0.0000  0.0000
ATOM     45  H   DRG     1       1.942  -1.645   1.091 0 .00  0.0000  0.0000  0.0000
ATOM     46  H   DRG     1       5.919   1.199   0.996 0 .00  0.0000  0.0000  0.0000
ATOM     47  H   DRG     1       5.040  -0.113   1.815 0 .00  0.0000  0.0000  0.0000
ATOM     48  H   DRG     1       6.303  -0.498   0.622 0 .00  0.0000  0.0000  0.0000
ATOM     49  H   DRG     1       6.303  -0.498   0.622 0 .00  0.0000  0.0000  0.0000
CONECT  1   19  21  23
CONECT  2   4  19  25
CONECT  3   7  19  26
CONECT  4   2  5  27
CONECT  5   4  6  28
CONECT  6   5  25  29
CONECT  7   3  8  30
CONECT  8   7  9  31
CONECT  9   8  26  32
CONECT  10   20  34  35  36
CONECT  11   12  17  20  37
CONECT  12   11  13  38  39
CONECT  13   12  14  23  40
CONECT  14   13  15  41  42
CONECT  15   14  16  20  43
CONECT  16   15  17  22  44
CONECT  17   11  16  22  45
CONECT  18   20  46  47  48  49
CONECT  19   1  2  3  24
CONECT  20   10  11  15  18
CONECT  21   1
CONECT  22   16  17
CONECT  23   1  13
CONECT  24   19  33
CONECT  25   2  6
CONECT  26   3  9
CONECT  27   4
CONECT  28   5
CONECT  29   6
CONECT  30   7
CONECT  31   8
CONECT  32   9
CONECT  33   24
CONECT  34   10
CONECT  35   10
CONECT  36   10
CONECT  37   11
CONECT  38   12
CONECT  39   12
CONECT  40   13
CONECT  41   14
CONECT  42   14
CONECT  43   15
CONECT  44   16
CONECT  45   17
CONECT  46   18
CONECT  47   18
CONECT  48   18
CONECT  49   18
001019001   0  0  0
001021002   0  0  0
001023001   0  0  0
002004004   0  0  0
002019001   0  0  0
002025004   0  0  0
003007004   0  0  0
003019001   0  0  0
003026004   0  0  0
004005004   0  0  0
004027001   0  0  0
005006004   0  0  0
005028001   0  0  0
006025004   0  0  0
006029001   0  0  0
007008004   0  0  0
007030001   0  0  0
008009004   0  0  0
008031001   0  0  0
009026004   0  0  0
009032001   0  0  0
010020001   0  0  0
010034001   0  0  0
010035001   0  0  0
010036001   0  0  0
011012001   0  0  0
011017001   0  0  0
011020001   0  0  0
011037001   0  0  0
012013001   0  0  0
012038001   0  0  0
012039001   0  0  0
013014001   0  0  0
013023001   0  0  0
013040001   0  0  0
014015001   0  0  0
014041001   0  0  0
014042001   0  0  0
015016001   0  0  0
015020001   0  0  0
015043001   0  0  0
016017001   0  0  0
016022001   0  0  0
016044001   0  0  0
017022001   0  0  0
017045001   0  0  0
018020001   0  0  0
018046001   0  0  0
018047001   0  0  0
018048001   0  0  0
018049001   0  0  0
019024001   0  0  0
024033001   0  0  0
**/
