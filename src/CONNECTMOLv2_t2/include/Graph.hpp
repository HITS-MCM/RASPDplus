#ifndef __GRAPH_HPP
#define __GRAPH_HPP 
#include <vector>
#include <stack>
#include <iostream>
#include <queue>

using namespace std;

#include "DArray"
/*
    Author	: Vidhu S. Pandey
    Description	: A template class implementation of graph
    Date	: 1/11/2005
    
    Copyright. All rights reserved.
*/

template <class T>
class edge;  //forward declaration

template <class T>
class vertex{
public:
	T info;
	int coloured;		// 1 for coloured 0 for uncoloured
	int distance;		// distance from a particular node
	int degree; 		// no of edges incident on the vertex
	vertex<T> *next;
	edge<T> *edges;
};

template <class T>
class edge{
public:
	int weight;            //single bond double bond etc
	vertex<T> *end;
	edge<T> *next;
};

template <class T>
class Graph{
private:
	unsigned nNodes_;////no of nodes
	vertex<T> *startVertex_,*endVertex_,*pointer_;
	vertex<T> *makeVertex();
	edge<T> *makeEdge();
	void resetColour()const; 			// uncolours all the vertices
	vertex<T> *findVertex(const T&)const;
	void addEdge(vertex<T>*,vertex<T>*,int=0);
	void removeEdge(vertex<T>*,vertex<T>*); 
	void findSubsetRing(vector<vector<T> >&)const;
	bool isSubsetRing(const vector<T>&,const vector<T>&)const;
public:
	
	Graph();
	Graph(const Graph<T> &); 
	Graph<T>& operator=(const Graph<T>&);
	void addVertex(T);
	void addEdge(const T&,const T&,int=0);
	void removeEdge(const T&,const T&); 
	void initializeDistance(int);
	int shortestDistance(int,int) const ;
	int shortestDistance2(int,int) const ;
	vector<int> shortestPath(int,int) const; 
	vector<vector<T> > findCycle(T&) const ;
	void setDegree();
        int getDegree(const T&) const;
        int getAdjacencyList(const T&,vector<T>&) const;
        void setWeight(const T&,const T&,int);
        int getWeight(const T&,const T&)const;
	vector<int> BFS(int) const;
	vector<int> DFS(int) const ;
	void print() const;
	void printVertex() const;           		//print()
	void printEdge(vertex<T>*) const;   		// print()
	unsigned size() const { return nNodes_; } 
	void clear();
        
	~Graph();
};

template <class T>
Graph<T>::Graph()
{
	startVertex_ = endVertex_ = pointer_ = NULL;
	nNodes_=0;
//	cout << "GRAPH CONSTRUCTOR CALLED" << endl;
}



template <class T>
Graph<T>::Graph(const Graph<T>& gr) 
{
	startVertex_ = endVertex_ = pointer_ = NULL;
	nNodes_=0;
	vertex<T>* tmpVer;
	edge<T>* tmpEdge;
	for(tmpVer=gr.startVertex_;tmpVer!=NULL;tmpVer=tmpVer->next)  
		addVertex(tmpVer->info);	
	
	for(tmpVer=gr.startVertex_;tmpVer!=NULL;tmpVer=tmpVer->next)  
		for(tmpEdge=tmpVer->edges;tmpEdge!=NULL;tmpEdge=tmpEdge->next)
			addEdge(tmpVer->info,tmpEdge->end->info);
	
}


template <class T>
Graph<T>& Graph<T>::operator=(const Graph<T>& gr){
	
	if(&gr == this) return *this;
	
	clear();
	vertex<T>* tmpVer;
	edge<T>* tmpEdge;
	for(tmpVer=gr.startVertex_;tmpVer!=NULL;tmpVer=tmpVer->next)  
		addVertex(tmpVer->info);	
	
	for(tmpVer=gr.startVertex_;tmpVer!=NULL;tmpVer=tmpVer->next)  
		for(tmpEdge=tmpVer->edges;tmpEdge!=NULL;tmpEdge=tmpEdge->next)
			addEdge(tmpVer->info,tmpEdge->end->info);
	
	return *this;
}

template <class T>
vertex<T>* Graph<T>::makeVertex()
{
	vertex<T>* tmp;
	tmp = new vertex<T>;
	if(!tmp) cout << "INSUFFICIENT MEMORY - Graph::makeVertex()" << endl;
	tmp->next = NULL;
	tmp->edges = NULL;
	return tmp;
}

template <class T>
edge<T>* Graph<T>::makeEdge()  
{
	edge<T>* tmp;
	tmp = new edge<T>;
	if(!tmp) cout << "INSUFFICIENT MEMORY - Graph::makeEdge()" << endl;
	tmp->next = NULL;
	tmp->end =NULL;
	return tmp;
}

template <class T>
void Graph<T>::resetColour()const{ // uncolours all the vertices

	vertex<T>* tmp;

	for(tmp=startVertex_;tmp!=NULL;tmp=tmp->next)
		tmp->coloured = 0;
}

template <class T>
void Graph<T>::initializeDistance(int argn){

	vertex<T>* tmp;

	for(tmp=startVertex_;tmp!=NULL;tmp=tmp->next)
		tmp->distance =argn;

}

template <class T>
void Graph<T>::setDegree(){

	vertex<T>* tmpVert;
	edge<T>* tmpEdge;
	int countEdge;

	for(tmpVert=startVertex_;tmpVert != NULL;tmpVert=tmpVert->next){

		countEdge=0;
		for(tmpEdge = tmpVert->edges;tmpEdge != NULL;tmpEdge = tmpEdge->next) 
			countEdge++;

		tmpVert->degree = countEdge;
	}
}

template <class T>
int Graph<T>::getDegree(const T& item)const{
    vertex<T>* tmpVert;
    tmpVert=findVertex(item);
    if(tmpVert)
        return tmpVert->degree;
    return 0;
}

template <class T>
int Graph<T>::getAdjacencyList(const T& item,vector<T>& vect)const{
    vertex<T>* tmpVert;
    edge<T>* tmp;
    vect.clear();
    tmpVert=findVertex(item);
    //cout << item <<"  ";
    if(tmpVert){
	for(tmp=tmpVert->edges;tmp!=NULL;tmp=tmp->next){
            vect.push_back(tmp->end->info);
      //      cout<<"  "<< tmp->end->info <<"    ";
        }
        //cout <<endl;
        return 1;
    }
    return 0;
}



template <class T>
void Graph<T>::addVertex(T argItem)
{
	vertex<T>*  temp;
	temp=makeVertex();
	temp->info=argItem;
	if(startVertex_==NULL)
	{
		startVertex_=endVertex_=temp;
	}
	else
	{	
		endVertex_->next=temp;
		endVertex_=temp;
	}
	nNodes_++;
}

template <class T>
void Graph<T>::addEdge(vertex<T>* argVertex1,vertex<T>* argVertex2,int wt) 
{
	if((argVertex1==NULL)||(argVertex2==NULL))
		return;
	edge<T>* temp,*endEdge;

	temp = makeEdge();	
	temp->weight = wt;
	temp->end=argVertex2;

	if(argVertex1->edges==NULL)
		argVertex1->edges=temp;
	
	else{	
		for(endEdge=argVertex1->edges;endEdge->next!=NULL;endEdge=endEdge->next);
			endEdge->next=temp;	
	}
}

template <class T>
void Graph<T>::addEdge(const T& info1,const T& info2,int wt)  
{
	vertex<T> *ver1, *ver2;
	
	ver1=findVertex(info1);
	ver2=findVertex(info2);
	addEdge(ver1,ver2,wt);
        addEdge(ver2,ver1,wt);

}

template <class T>
void Graph<T>::setWeight(const T& item1,const T& item2,int wt)  
{
	vertex<T> *ver1, *ver2;
	edge<T> *tmpEdge=0;
        
	ver1=findVertex(item1);
	ver2=findVertex(item2);
	for(tmpEdge=ver1->edges;tmpEdge!=0;tmpEdge=tmpEdge->next)
            if(tmpEdge->end->info == item2){
                tmpEdge->weight=wt;
                break;
            }
                
        for(tmpEdge=ver2->edges;tmpEdge!=0;tmpEdge=tmpEdge->next)
            if(tmpEdge->end->info == item1){
                tmpEdge->weight=wt;
                break;
            }
            
}

template <class T>
int Graph<T>::getWeight(const T& item1,const T& item2)const{
    vertex<T> *ver1;
    edge<T> *tmpEdge=0;

    ver1=findVertex(item1);
    for(tmpEdge=ver1->edges;tmpEdge!=0;tmpEdge=tmpEdge->next)
        if(tmpEdge->end->info == item2){
            return tmpEdge->weight;
        }
    return 0;
}

template <class T>
void Graph<T>::removeEdge(vertex<T>* ver1,vertex<T>* ver2) 
{
	edge<T> *tmpEdge1,*tmpEdge2;
	
	if(ver1->edges->end==ver2)
	{
		tmpEdge1=ver1->edges;
		ver1->edges=tmpEdge1->next;
		delete tmpEdge1;
		return;
	}

	else
	{
		for(tmpEdge2=ver1->edges;tmpEdge2->next!=NULL;tmpEdge2=tmpEdge2->next)
		{
			if(tmpEdge2->next->end==ver2)
			{
				tmpEdge1=tmpEdge2->next;
				tmpEdge2->next=tmpEdge1->next;
				delete tmpEdge1;
				return;
			}
		}
	}
}


template <class T>
vertex<T>* Graph<T>::findVertex(const T& item)const
{
	vertex<T>* tmp;
	for(tmp=startVertex_;tmp!=NULL;tmp=tmp->next)
	{	
		if(tmp->info==item)
			return tmp;
	}
	return NULL;
}
 
template <class T>
vector<int>  Graph<T>::shortestPath(int argAtomNo1,int argAtomNo2) const{
	vertex<T> *tmpVert1, *tmpVert2;
	edge<T>   *tmpEdge;
	std::queue <vertex<T>* > que;	//que -> holds the vertex element using the properties enqueue(to insert),deque 
	vector<int> originList,queList,pathList;	//list->traversal linklist,originList->maintain 
	stack<int> stk;//to reverse
	unsigned tmp,flag;//flag to break out of  the outer loop
	
	if(argAtomNo1==argAtomNo2){
		pathList.push_back(argAtomNo1);
		return pathList;
	}

	resetColour();//uncolours the graph
	initializeDistance(-1);//sets distances of all nodes to -1

	tmpVert1 = findVertex(argAtomNo1);

	que.push(tmpVert1); ////insert first time
	tmpVert1->coloured=1;
	tmpVert1->distance=0;
	queList.push_back(tmpVert1->info.getAtomSerial());
	originList.push_back(0);


//	cout << "BREADTH FIRST SEARCH" << endl;

	while (!que.empty()){
		
		tmpVert1=que.front();
		que.pop();

	//	list.add(tmpVert1->info.getAtomSerial());///not used in this program as yet 

		flag=0;
		for(tmpEdge=tmpVert1->edges;tmpEdge != NULL;tmpEdge=tmpEdge->next){
		
			tmpVert2 = tmpEdge->end;
		
			if(tmpVert2->coloured == 0){ 

				tmpVert2->distance = tmpVert1->distance+1;

				if((tmpVert2->info.getAtomSerial() == (unsigned)argAtomNo2)||(tmpVert2->distance==4)){
				
	//				list.add(tmpVert2->info.getAtomSerial());
					queList.push_back(tmpVert2->info.getAtomSerial());
					originList.push_back(tmpVert1->info.getAtomSerial());
					flag=1;
			//		queList.print();
			//		originList.print();
					break;
				}


				que.push(tmpVert2);
				tmpVert2->coloured=1;
				queList.push_back(tmpVert2->info.getAtomSerial());
				originList.push_back(tmpVert1->info.getAtomSerial());
			}
		}
		if(flag==1)
			break;
	}
		
	//queList.print();
	//originList.print();

	tmp = queList[queList.size()-1];

	while(tmp!=(unsigned)argAtomNo1){

		stk.push(tmp);

		for(unsigned i=0;i<queList.size();i++)
			if((unsigned)queList[i] == tmp)
				tmp = originList[i];
	}
	stk.push(tmp);

	while(!stk.empty()){
		pathList.push_back(stk.top());
		stk.pop();	
	}

	
	return pathList;

}

template <class T>
int  Graph<T>::shortestDistance(int argAtomNo1,int argAtomNo2) const{

	vertex<T>* tmpVert1,*tmpVert2;
	edge<T>* tmpEdge;
	std::queue< vertex<T>* > que; 
	int tmp,i,flag;//flag to break out of  the outer loop
	
	if(argAtomNo1==argAtomNo2){
		return 0;
	}

	resetColour();//uncolours the graph
	initializeDistance(-1);//sets distances of all nodes to -1

	tmpVert1 = findVertex(argAtomNo1);

	que.push(tmpVert1); ////insert first time
	tmpVert1->coloured=1;
	tmpVert1->distance=0;


	while (!que.empty()){
		
		tmpVert1=que.front();
		que.pop();

		flag=0;
		for(tmpEdge=tmpVert1->edges;tmpEdge != NULL;tmpEdge=tmpEdge->next){
		
			tmpVert2 = tmpEdge->end;
		
			if(tmpVert2->coloured == 0){ 

				tmpVert2->distance = tmpVert1->distance+1;

				if((tmpVert2->info.getAtomSerial() == argAtomNo2)){ 
				
					flag=1;
					break;
				}

				que.push(tmpVert2);
				tmpVert2->coloured=1;
			}
		}
		if(flag==1)
			break;
	}
		
	return tmpVert2->distance;

}

template <class T>
vector<int> Graph<T>::BFS(int argAtomNo) const{

	vertex<T>* tmp1,*tmp3;
	edge<T>* tmp2;
	std::queue< vertex<T>*> que;
	vector<int> list;
//	vector<vertex<T>*> priorityList;

	resetColour(); //uncolours the graph

	tmp1 = findVertex(argAtomNo);

	que.push(tmp1);
	tmp1->coloured=1;

//	cout << "BREADTH FIRST SEARCH" << endl;

	while (!que.empty()){
		
		tmp1=que.front();
		que.pop();

		list.push_back(tmp1->info.getAtomSerial());

	//	cout << tmp1->info.getAtomSerial() << "   ";

		for(tmp2=tmp1->edges;tmp2 != NULL;tmp2=tmp2->next){
			tmp3 = tmp2->end;
			if(tmp3->coloured== 0){ 
				que.push(tmp3);
				tmp3->coloured=1;
			}
		}
	}
//	cout <<  endl << endl;
	return list;
}


template <class T>
vector<vector<T> > Graph<T>::findCycle(T& argAtomNo) const{

	vertex<T>* tmp1,*tmp3;
	edge<T>* tmp2;
	stack<vertex<T>*> stk;
	vector<T> ringList;
	vector<vector<T> > rings;
	DArray<T> fatherList(nNodes_+1);
	T lp1,lp2,father;
	resetColour();
	const T STARTATOM = startVertex_->info; 

	tmp1 = findVertex(argAtomNo);

//cout << tmp1->info.getResidueName() << " " <<tmp1->info.getResidueSequence() << "   " << n << endl;
	//colour 0 - unvisited, 1 - discovered, 2 - visited

	stk.push(tmp1);
	tmp1->coloured = 1; // discovered
	fatherList[0] = -1;

	while(!stk.empty()){

		tmp1 = stk.top();
		stk.pop();
		tmp1->coloured = 2; //visited

		//list.push_back(tmp1->info.getAtomSerial());

		for(tmp2 = tmp1->edges;tmp2 != NULL;tmp2=tmp2->next){
			tmp3 = tmp2->end;
			if(tmp3->coloured == 0){
				stk.push(tmp3);
				tmp3->coloured = 1; // discovered
				fatherList[tmp3->info-STARTATOM] = tmp1->info;
			//	cout << tmp3->info.getAtomSerial() << "  " << fatherList[tmp3->info.getAtomSerial()] << endl;
			}
			else if(tmp3->coloured==1  ){
	//			cout << "\nCYCLE" << "   ";
				lp1 =  tmp1->info;
				lp2 =  tmp3->info;
			//	cout << lp2 << "   " << lp1 <<endl;
				
				father = lp1;
		//		cout<< fatherList[lp2-STARTATOM]<< "  "  << lp2 << "  ";
				
				ringList.clear();
				ringList.push_back(fatherList[lp2-STARTATOM]);	
				ringList.push_back(lp2);	
				
				while(father!=fatherList[lp2-STARTATOM]){
			//		cout << father << "  ";
					ringList.push_back(father);	
					father = fatherList[father-STARTATOM];
				}
				rings.push_back(ringList);
		//		cout << endl;
			}
		}
	}
	findSubsetRing(rings);
	return rings;
}

template <class T>
void Graph<T>::findSubsetRing(vector<vector<T> >& rings)const{

    vector<T> rng1,rng2;
    vector<vector<T> > tmpVect=rings;
    vector<T> newRing;
    for(unsigned i=0;i<tmpVect.size();i++)
        for(unsigned j=i+1;j<tmpVect.size();j++){
            if((tmpVect[i][0]==tmpVect[j][0])&&(tmpVect[i][1]==tmpVect[j][1])){
                rng1 = tmpVect[i];
                rng2 = tmpVect[j];
		if(isSubsetRing(rng1,rng2)){
			if(rng1.size()>rng2.size()){
			    unsigned k;
			    for(k=1;k<rng1.size()&&(rng1[k]!=rng2[2]);k++)
				newRing.push_back(rng1[k]);
				newRing.push_back(rng1[k]);
			    rings.erase(rings.begin()+(i),rings.begin()+(i+1));
			    rings.push_back(newRing);
			}
			else if(rng2.size()>rng1.size()){
			    unsigned k;
			    for(k=1;k<rng2.size()&&(rng2[k]!=rng1[2]);k++)
				newRing.push_back(rng2[k]);
				newRing.push_back(rng2[k]);
			    rings.erase(rings.begin()+(j),rings.begin()+(j+1));
			    rings.push_back(newRing);
			}
		}
            }
        }
   //for(unsigned i=0;i<newRing.size();i++)
   //	 cout << newRing[i] << "  " ;
}

template <class T>
bool Graph<T>::isSubsetRing(const vector<T>& ring1,const vector<T>& ring2)const{

	// if ring 2 is subset of ring 1
	if(ring1.size()==ring2.size())
		return false;
	if(ring1.size()>ring2.size()){
		for(unsigned i=0;i<ring2.size();i++){
			int flag=0;
			for(unsigned j=0;j<ring1.size();j++){
				if(ring2[i]==ring1[j]){
					flag=1;
					break;
				}
			}
			if(flag==0)
				return false;
		}
	}
	// if ring 1 is subset of ring 2
	else{
		for(unsigned i=0;i<ring1.size();i++){
			int flag=0;
			for(unsigned j=0;j<ring2.size();j++){
				if(ring1[i]==ring2[j]){
					flag=1;
					break;
				}
			}
			if(flag==0)
				return false;
		}
	}
	return true;
}


template <class T>
vector<int> Graph<T>::DFS(int argAtomNo) const{

	vertex<T>* tmp1,*tmp3;
	edge<T>* tmp2;
	stack<vertex<T>*> stk;
	vector<int> list;

	resetColour();

	tmp1 = findVertex(argAtomNo);

	stk.push(tmp1);
	tmp1->coloured = 1;

//	cout << "DEPTH FIRST SEARCH" << endl;

	while(!stk.empty()){

		tmp1 = stk.top();
		stk.pop();

		list.push_back(tmp1->info.getAtomSerial());

//		cout << tmp1->info.getAtomSerial() << "   ";

		for(tmp2 = tmp1->edges;tmp2 != NULL;tmp2=tmp2->next){
			tmp3 = tmp2->end;
			if(tmp3->coloured == 0){
				stk.push(tmp3);
				tmp3->coloured = 1;
			}
		}
	}
//	cout <<  endl << endl;
	return list;
}


template <class T>
void Graph<T>::printVertex() const // prints all the vertices
{
	vertex<T>* tmp;
	for(tmp=startVertex_;tmp!=NULL;tmp=tmp->next)
	{
		cout<<tmp->info<<endl;
	}
}
template <class T>
void Graph<T>::printEdge(vertex<T>* argVertex) const
{
	edge<T>* tmp;
	for(tmp=argVertex->edges;tmp!=NULL;tmp=tmp->next)
//		cout<<"  "<< tmp->end->info.getAtomSerial() << "(" << tmp->weight << ") "  ;
		cout<<"  "<< tmp->end->info << "("<< tmp->weight <<")    ";
}

template <class T>
void Graph<T>::print() const // prints the whole graph (vertex followed by its adjoining vertices )
{
	vertex<T>* tmp;
//	cout<<"First "<<startVertex_->info.<<"            "<<startVertex_->edges->end->info.getAtomSerial()<<endl;
	for(tmp=startVertex_;tmp!=NULL;tmp=tmp->next)
	{
//		cout<<tmp->info.getAtomSerial() << "{" << tmp->degree << "}" <<"    ";
		cout<<tmp->info << "["<<tmp->degree<<"]   " ;
		printEdge(tmp);
		cout<<endl;
	}
}


template <class T>
Graph<T>::~Graph(){
	clear();
//	cout << "GRAPH DESTRUCTOR CALLED"<<endl;
}

template <class T>
void Graph<T>::clear(){

	vertex<T> *tmpVert,*storeVert;
	edge<T> *tmpEdge,*storeEdge;

	tmpVert = startVertex_;
	while(tmpVert){
		tmpEdge=tmpVert->edges;
		
		while(tmpEdge){
			storeEdge = tmpEdge->next;
			delete tmpEdge;
			tmpEdge = storeEdge;
		}
		
		storeVert = tmpVert->next;
		delete tmpVert;
		tmpVert = storeVert;
	}
	startVertex_ = endVertex_ = pointer_ = NULL;
	nNodes_=0;
}


template <class T>
int  Graph<T>::shortestDistance2(int argAtomNo1,int argAtomNo2) const{

	vertex<T>* tmpVert1,*tmpVert2;
	edge<T>* tmpEdge;
	std::queue< vertex<T>* > que; 
	int tmp,i,flag;//flag to break out of  the outer loop
	
	if(argAtomNo1==argAtomNo2){
		return 0;
	}

	resetColour();//uncolours the graph
	initializeDistance(-1);//sets distances of all nodes to -1

	tmpVert1 = findVertex(argAtomNo1);

	que.push(tmpVert1); ////insert first time
	tmpVert1->coloured=1;
	tmpVert1->distance=0;


	while (!que.empty()){
		
		tmpVert1=que.front();
		que.pop();

		flag=0;
		for(tmpEdge=tmpVert1->edges;tmpEdge != NULL;tmpEdge=tmpEdge->next){
		
			tmpVert2 = tmpEdge->end;
		
			if(tmpVert2->coloured == 0){ 

				tmpVert2->distance = tmpVert1->distance+1;

				if((tmpVert2->info.getAtomSerial() == argAtomNo2)){ 
				
					flag=1;
					break;
				}

				que.push(tmpVert2);
				tmpVert2->coloured=1;
			}
		}
		if(flag==1)
			break;
	}
		
	return tmpVert2->distance;

}

#endif	
