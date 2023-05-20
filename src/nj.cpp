// $Id: nj.cpp 9948 2011-10-23 15:53:03Z cohenofi $

// version 1.00
// last modified 3 Nov 2002


#include "nj.h"
#include "errorMsg.h"
#include "logFile.h"
#include "treeUtil.h"
#include "threads.h"
#include <cassert>
#include <algorithm>
#include <map>
using namespace std;


//------------------------------------------
// general outline:
// we follow Swofford's book, "Molecular Systematics" pg489.
// currentNodes is the vector of the nodes that are "in process".
// in the beggining, these are all the leaves. Once, 2 leaves are separeted, 
// they are excluded from currentNodes, and their father is added to currentNodes.
// we (almost) finish the algorithm when currentNodes's size is 3. (i.e., we know the topology).
// thus when we start from an evolutionary tree, all we do, is to construct a star (start) tree
//------------------------------------------




//------------------------------------------
// constructor and start
//------------------------------------------
tree NJalg::computeTree(VVdouble distances,const vector<string>& names, const tree * const constriantTree /*= NULL*/){
	assert(distances.size() == names.size());
	tree resTree = startingTree(names);
	if (distances.size()<3) return resTree;
	vector<tree::nodeP> currentNodes;
	resTree.getAllLeaves(currentNodes,resTree.getRoot());
	if (constriantTree) {
		njConstraint njc(resTree, *constriantTree);
		while (currentNodes.size() >= 3) NJiterate(resTree,currentNodes,distances, njc);
	} else {
		while (currentNodes.size() >= 3) NJiterate(resTree,currentNodes,distances);
	}
	resTree.create_names_to_internal_nodes();
	resTree.makeSureAllBranchesArePositive();
	LOGDO(5,resTree.output(myLog::LogFile()));
	return resTree;
}

tree NJalg::startingTree(const vector<string>& names) {
	return starTree(names);
}

tree NJalg::startingTree(const tree& inTree) {
	tree et;
	et.createRootNode();
	vector<tree::nodeP> allLeaves;
	inTree.getAllLeaves(allLeaves,inTree.getRoot());
	
	vector<string> names(allLeaves.size());
	for (size_t k = 0 ; k < allLeaves.size(); ++k)
	  names[k]=allLeaves[k]->name();
	  
	return startingTree(names);
}

void NJalg::updateBranchDistance(const VVdouble& distanceTable,
								 const Vdouble& rValues,
								 tree::nodeP nodeNew,
								 tree::nodeP nodeI,
								 tree::nodeP nodeJ,
								 int Iplace,
								 int Jplace) {
	MDOUBLE dis= (Iplace<Jplace) ? distanceTable[Iplace][Jplace] : distanceTable[Jplace][Iplace];
	MDOUBLE DisI_new = dis/2.0;
	MDOUBLE tmp = rValues[Iplace] - rValues[Jplace];
	tmp/= (  2.0*(distanceTable.size()-2)  );
	DisI_new = DisI_new+ tmp;
	MDOUBLE DisJ_new = dis - DisI_new;
	if (DisI_new<tree::SHORT_LENGTH_VALUE) DisI_new=tree::SHORT_LENGTH_VALUE; // no negative..
	if (DisJ_new<tree::SHORT_LENGTH_VALUE) DisJ_new=tree::SHORT_LENGTH_VALUE; // no negative..
	nodeI->setDisToFather(DisI_new);
	nodeJ->setDisToFather(DisJ_new);
}

void NJalg::NJiterate(tree& et,
					  vector<tree::nodeP>& currentNodes,
							 VVdouble& distanceTable) {
	Vdouble	rVector = calc_r_values(currentNodes,distanceTable);//CHECK2
	
	if (currentNodes.size() == 3) {
		update3taxaLevel(distanceTable,rVector,currentNodes);
		currentNodes.clear();
		return;
	}
	
	int minRaw,minCol;
	calc_M_matrix(currentNodes,distanceTable,rVector,minRaw,minCol);//CHECK3
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(et,nodeI,nodeJ);
	//CHECK4
	updateBranchDistance(distanceTable,rVector,theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6
		et.create_names_to_internal_nodes();
	UpdateDistanceTableAndCurrentNodes(currentNodes,distanceTable,nodeI,nodeJ,theNewNode,minRaw,minCol);
}

void NJalg::NJiterate(tree& et,
		      vector<tree::nodeP>& currentNodes,
		      VVdouble& distanceTable,
		      njConstraint& njc) {
	Vdouble	rMatrix = calc_r_values(currentNodes,distanceTable);//CHECK2
	
	if (currentNodes.size() == 3) {
		update3taxaLevel(distanceTable,rMatrix,currentNodes);
		currentNodes.clear();
		return;
	}
	
	int minRaw,minCol;
	calc_M_matrix(currentNodes,distanceTable,rMatrix,minRaw,minCol, njc);//CHECK3
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(et,nodeI,nodeJ);
	njc.join(nodeI, nodeJ, theNewNode);
	//CHECK4
	updateBranchDistance(distanceTable,rMatrix,theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6
		et.create_names_to_internal_nodes();
	UpdateDistanceTableAndCurrentNodes(currentNodes,distanceTable,nodeI,nodeJ,theNewNode,minRaw,minCol);
	LOGDO(15,et.output(myLog::LogFile(),tree::ANCESTORID));

}



Vdouble NJalg::calc_r_values(vector<tree::nodeP>& currentNodes,
							 const VVdouble& distanceTable) {
	Vdouble r_values(currentNodes.size(),0.0);
	for (size_t i=0; i <r_values.size();++i) {
		for (size_t j =0; j < r_values.size();++j) {
			MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
			r_values[i] += dis;
		}
	}
	return r_values;
}

void NJalg::calc_M_matrix(vector<tree::nodeP>& currentNodes,
						  const VVdouble& distanceTable,
						  const Vdouble & r_values,
						  int& minRaw,int& minCol){
	MDOUBLE min = VERYBIG;
	for (size_t i=0; i < currentNodes.size();++i){
		for (size_t j =i+1; j < currentNodes.size();++j) {
			MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
			MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
			if (tmp<min) {minRaw = i;minCol=j;min=tmp;}
			
		}
	}
}

void NJalg::calc_M_matrix(vector<tree::nodeP>& currentNodes,
			  const VVdouble& distanceTable,
			  const Vdouble & r_values,
			  int& minRaw,int& minCol, 
			  const njConstraint& njc){
	MDOUBLE min = VERYBIG;
	MDOUBLE min_noc =  VERYBIG;
	int minRaw_noc=-1,minCol_noc=-1;
	for (size_t i=0; i < currentNodes.size();++i){
	  for (size_t j =i+1; j < currentNodes.size();++j) {
	    if (njc.isCompatible(currentNodes[i],currentNodes[j])) {
	      MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
	      MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
	      if (tmp<min) {minRaw = i;minCol=j;min=tmp;}
	    }
	    LOGDO(10,{
		    MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
		    MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
		    if (tmp<min_noc) {minRaw_noc = i;minCol_noc=j;min_noc=tmp;}
		  });
	      
	      }
	}
	LOGDO(10, {if (min_noc != min) {myLog::LogFile() 
		    << "NJ-constratin changes outcome " <<
		    currentNodes[minRaw_noc]->name()<<","<<currentNodes[minCol_noc]->name() <<"-> " <<
		    currentNodes[minRaw]    ->name()<<","<<currentNodes[minCol]    ->name()<< 
		    "  ("<<min-min_noc<<")"<<endl; 
		  njc.isCompatible(currentNodes[minRaw_noc], currentNodes[minCol_noc], true);
		  myLog::LogFile() << njc <<endl;
		}
	      });
}

tree::nodeP NJalg::SeparateNodes(tree& et, tree::nodeP node1,
											 tree::nodeP node2) {
	if (node1->father() != node2->father()) 
	 errorMsg::reportError(" error in function NJalg::SeparateNodes - nodes don't have the same father");

	tree::nodeP fatherNode = node1->father();

	tree::nodeP theNewNode = et.createNode(fatherNode,et.getNodesNum());
	node1->setFather(theNewNode);
	theNewNode->setSon(node1);
	node2->setFather(theNewNode);
	theNewNode->setSon(node2);

	// remove from son list of father node.
	fatherNode->removeSon(node1); 

	fatherNode->removeSon(node2); 
	return theNewNode;
}

void NJalg::update3taxaLevel(VVdouble& distanceTable,Vdouble & r_values,
							 vector<tree::nodeP>& currentNodes) {
	// update the distance of the 3 taxa that are left in the end, to the root.
	
	MDOUBLE dis0root = distanceTable[0][1]/2+0.5*(r_values[0]-r_values[1]);
	MDOUBLE dis1root = distanceTable[0][1]/2+0.5*(r_values[1]-r_values[0]);
	MDOUBLE dis2root = distanceTable[0][2]/2+0.5*(r_values[2]-r_values[0]);
	if (dis0root<tree::SHORT_LENGTH_VALUE) dis0root=tree::SHORT_LENGTH_VALUE; // no negative..
	if (dis1root<tree::SHORT_LENGTH_VALUE) dis1root=tree::SHORT_LENGTH_VALUE; // no negative..
	if (dis2root<tree::SHORT_LENGTH_VALUE) dis2root=tree::SHORT_LENGTH_VALUE; // no negative..
	currentNodes[0]->setDisToFather(dis0root);
	currentNodes[1]->setDisToFather(dis1root);
	currentNodes[2]->setDisToFather(dis2root);
}

void NJalg::UpdateDistanceTableAndCurrentNodes(vector<tree::nodeP>& currentNodes,
											   VVdouble& distanceTable,
											   tree::nodeP nodeI,
											   tree::nodeP nodeJ,
											   tree::nodeP theNewNode,
											   int Iplace,
											   int Jplace) {
	//	Iplace is the place of i in the "old" currentNodes vector
	size_t i,j;
	//	updating currentNodes
	vector<tree::nodeP> newCurrentNode= currentNodes;

	vector<tree::nodeP>::iterator vec_iter1=remove(
		newCurrentNode.begin(),newCurrentNode.end(),nodeI );
	newCurrentNode.erase(vec_iter1,newCurrentNode.end());

	vector<tree::nodeP>::iterator vec_iter2=remove(
	newCurrentNode.begin(),newCurrentNode.end(),nodeJ );
	newCurrentNode.erase(vec_iter2,newCurrentNode.end());
	
	newCurrentNode.push_back(theNewNode);

	map<tree::nodeP,int> nodeIntMap1;
	for (size_t z=0; z<currentNodes.size();++z) {
		nodeIntMap1.insert(map<tree::nodeP,int>::value_type(currentNodes[z],z));
	}

	VVdouble newDisTable;
	newDisTable.resize(newCurrentNode.size());
	for (size_t z1=0;z1<newDisTable.size();++z1) newDisTable[z1].resize(newCurrentNode.size(),0.0);

// updatine the table
	for (i=0; i < newCurrentNode.size(); i++) {
		for (j=i+1; j < newCurrentNode.size() ; j++) {
			if ((i!=newCurrentNode.size()-1) && (j!=newCurrentNode.size()-1)) {// both old nodes
				int oldI = nodeIntMap1[newCurrentNode[i]];
				int oldJ = nodeIntMap1[newCurrentNode[j]];
				MDOUBLE dis= (oldI<oldJ) ? distanceTable[oldI][oldJ] : distanceTable[oldJ][oldI];
				newDisTable[i][j] = dis;
			} //else if (i==newCurrentNode.size()-1) { // i is new
			//	newDisTable[i][j] = (dis(Iplace,NewOldPlaces[j])+dis(Jplace,NewOldPlaces[j])-dis(Iplace,Jplace))/2.0;
			//}
			else if (j==newCurrentNode.size()-1) { // j is new
				int oldI = Iplace;
				int oldJ = Jplace;
				int oldK = nodeIntMap1[newCurrentNode[i]];
				MDOUBLE disIK= (oldI<oldK) ? distanceTable[oldI][oldK] : distanceTable[oldK][oldI];
				MDOUBLE disIJ= (oldI<oldJ) ? distanceTable[oldI][oldJ] : distanceTable[oldJ][oldI];
				MDOUBLE disJK= (oldJ<oldK) ? distanceTable[oldJ][oldK] : distanceTable[oldK][oldJ];
				newDisTable[i][j] = 0.5*(disIK+disJK-disIJ); //EQ. 43 SWOFFORD PAGE 489.
			}
		}
	}

	currentNodes=newCurrentNode;
	distanceTable=newDisTable;
}

/*
* This function does what UpdateDistanceTableAndCurrentNodes does but more efficiently
* without allocating new memory. It is also written in a more elegant way that is understandable:
* we want to replace two nodes with a new one using an equation so we will iterate over
* each row in the table that holds the distances between the nodes.
* In each row we will add the new node distance using an equation, The equation changes depending if 
* the row index is smaller than the nodes that should be replaced between them or larger than them.
* after we added the new node distance to the row we will erase the old nodes distances.
*/
void NJalg::UpdateDistanceTableAndCurrentNodesInPlace(vector<tree::nodeP>& currentNodes,
											   VVdouble& distanceTable,
											   tree::nodeP nodeI,
											   tree::nodeP nodeJ,
											   tree::nodeP theNewNode,
											   int Iplace,
											   int Jplace) {
	//	Iplace and Jplace are the indexes of i and j in the input currentNodes vector
	int k;
	int low_index, high_index;
	int distance_table_length = currentNodes.size();

	//for several reasons it is much easier to understand this algorithm if you know which index is larger Iplace
	//or Jplace so we'll put them in new variables in the correct order
	if (Iplace < Jplace) {
		low_index = Iplace;
		high_index = Jplace;
	} else {
		high_index = Iplace;
		low_index = Jplace;
	}

	//erase the nodes from the node vector and add the new one in their stead
	currentNodes.erase(currentNodes.begin()+high_index);
	currentNodes.erase(currentNodes.begin()+low_index);
	currentNodes.push_back(theNewNode);

	MDOUBLE disIJ = distanceTable[low_index][high_index];
	for (k=0; k < low_index; k++){          /*disIK                 disIJ   disJK*/
		distanceTable[k].push_back(0.5*(distanceTable[k][low_index]-disIJ+distanceTable[k][high_index]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}
	for (++k; k < high_index; k++) {
		distanceTable[k].push_back(0.5*(distanceTable[low_index][k]-disIJ+distanceTable[k][high_index]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}
	for (++k; k < distance_table_length; k++){
		distanceTable[k].push_back(0.5*(distanceTable[low_index][k]-disIJ+distanceTable[high_index][k]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}
	distanceTable.erase(distanceTable.begin() + high_index);
	distanceTable.erase(distanceTable.begin() + low_index);	
}

struct node_vector_data {
	vector<tree::nodeP> &currentNodes;
	tree::nodeP aNewNode;
	int high_index;
	int low_index;

	node_vector_data(vector<tree::nodeP> &currentNodes, tree::nodeP aNewNode,int high_index,int low_index) :
		currentNodes(currentNodes), aNewNode(aNewNode), high_index(high_index), low_index(low_index){}
};

struct node_vector_data {
	vector<tree::nodeP> &currentNodes;
	tree::nodeP aNewNode;
	int high_index;
	int low_index;

	node_vector_data(vector<tree::nodeP> &currentNodes, tree::nodeP aNewNode,int high_index,int low_index) :
		currentNodes(currentNodes), aNewNode(aNewNode), high_index(high_index), low_index(low_index){}
};

void *replace_nodes_with_new(void *args)
{
	struct node_vector_data *nodes = (struct node_vector_data *)args;
	vector<tree::nodeP>& currentNodes = nodes->currentNodes;
	currentNodes.erase(currentNodes.begin()+nodes->high_index);
	currentNodes.erase(currentNodes.begin()+nodes->low_index);
	currentNodes.push_back(nodes->aNewNode);

	return NULL;
}

struct distance_table_data {
	VVdouble *distanceTable;
	int start_index;
	int last_index;
	int high_index_to_remove;
	int low_index_to_remove;
    void update_distance_table_fields(VVdouble *distanceTable, int start_index, int last_index, int high_index_to_remove, int low_index_to_remove)
	{
          distanceTable = distanceTable;
          start_index = start_index;
          last_index = last_index;
          high_index_to_remove = high_index_to_remove;
          low_index_to_remove = low_index_to_remove;
	}
	distance_table_data() = default;
	~distance_table_data() = default;
};


#define MIN(x, y) ((x) < (y) ? (x) : (y))

void *replace_nodes_with_new_in_distance_table(void *args)
{
	struct distance_table_data *distance_table_data = (struct distance_table_data *)args;
	VVdouble& distanceTable = *(distance_table_data->distanceTable);
	int distance_table_length = distanceTable[0].size();
	int start_index = distance_table_data->start_index;
	int last_index = distance_table_data->last_index;
	MDOUBLE disIJ = distanceTable[distance_table_data->low_index_to_remove][distance_table_data->high_index_to_remove];

	int low_index = MIN(last_index, distance_table_data->low_index_to_remove);
	int high_index = MIN(last_index, distance_table_data->high_index_to_remove);
	last_index = MIN(last_index, distance_table_length);//we can't have the last index be bigger then the table size
	int k;

	for (k=start_index; k < low_index; k++){          /*disIK                 disIJ   disJK*/
		distanceTable[k].push_back(0.5*(distanceTable[k][low_index]-disIJ+distanceTable[k][high_index]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}
	for (; k < high_index; k++) {
		distanceTable[k].push_back(0.5*(distanceTable[low_index][k]-disIJ+distanceTable[k][high_index]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}
	for (; k < last_index; k++){
		distanceTable[k].push_back(0.5*(distanceTable[low_index][k]-disIJ+distanceTable[high_index][k]));//EQ. 43 SWOFFORD PAGE 489.
		distanceTable[k].erase(distanceTable[k].begin() + high_index);
		distanceTable[k].erase(distanceTable[k].begin() + low_index);
	}

	return NULL;
}

/*
 * A multicore lockless algorithm to improve the UpdateDistanceTableAndCurrentNodesInPlace
 * function by dividing the work between different cores. should scale linearly by
 * the number of cores.
 */
void NJalg::UpdateDistanceTableAndCurrentNodesMulticore(vector<tree::nodeP>& currentNodes,
											   VVdouble& distanceTable,
											   tree::nodeP nodeI,
											   tree::nodeP nodeJ,
											   tree::nodeP theNewNode,
											   int Iplace,
											   int Jplace) {
	int low_index_to_remove, high_index_to_remove;
	int distance_table_length = currentNodes.size();
	int nodes_per_core;

	//for several reasons it is much easier to understand this algorithm if you know which index is larger Iplace
	//or Jplace so we'll put them in a snew variables in the correct order
	if (Iplace < Jplace) {
		low_index_to_remove = Iplace;
		high_index_to_remove = Jplace;
	} else {
		high_index_to_remove = Iplace;
		low_index_to_remove = Jplace;
	}

	//multicore implementation dividing the work between workers:
	struct distance_table_data distance_table_data[MAX_THREADS];//data each worker will need for his work
	size_t number_of_cores = NUM_OF_WORKERS-1;//there is one core which is saved for replacing the nodes in currentNodes vector
											 //while all the others will work in adjusting the distance table
	size_t number_of_nodes_per_core = 0;

	number_of_nodes_per_core = distance_table_length/number_of_cores;
	//In case we have more cores than nodes we will reduce the number of cores to the number of nodes
	//as two cores can't work on 1 node effectively with lockless algorithm
	if (0  == number_of_nodes_per_core){
		number_of_nodes_per_core = 1;
		number_of_cores = distance_table_length;
	}

	number_of_cores++;//should keep one core for the reminder nodes: In C integer math 5/2=2
		                  // while 2*2=4 so if we use only 2 cores we will need another core for the fifth node
	//initialize the structs for each core and set it to work
	int start_index = 0;
	int last_index = number_of_nodes_per_core;
	int core_num;
	for (core_num=0; core_num<number_of_cores; ++core_num) {
		/*we iniate the struct with it's Ctor for each worker:*/
		distance_table_data[core_num].update_distance_table_fields(&distanceTable, start_index, last_index, high_index_to_remove, low_index_to_remove);
		last_index += number_of_nodes_per_core;//we let the workers themselves check if it is out of bounds
		start_index = last_index;
		schedule_work(replace_nodes_with_new_in_distance_table, &distance_table_data[core_num], core_num);
	}


	struct node_vector_data node_vector_data(currentNodes, theNewNode, high_index_to_remove, low_index_to_remove);
	schedule_work(replace_nodes_with_new, &node_vector_data, core_num);
	check_that_work_is_done(0, core_num);

	distanceTable.erase(distanceTable.begin() + high_index_to_remove);
	distanceTable.erase(distanceTable.begin() + low_index_to_remove);
}
/*
NJalg::NJalg(){
	_myET = NULL;
}



//-----------------------------
// The algorithm
//-----------------------------

void NJalg::GetDisTable(const sequenceContainer& sd,const vector<MDOUBLE>  * weights) {
	
	VVresize(_startingDistanceTable,distanceTable.size(),distanceTable.size());// for printing stuff later.
	VVresize(LTable,distanceTable.size(),distanceTable.size());// for printing stuff later.

	int i,j;
	_nodeNames.resize(currentNodes.size());
	for ( i=0; i < currentNodes.size(); i++) {
		_nodeNames[i] =(currentNodes[i]->name()); 
		for ( j=i+1; j < currentNodes.size(); j++) {
			MDOUBLE tempDis = -2000.0;
			MDOUBLE resLikelihood;
			int seqnodeI_ID = sd.getId(currentNodes[i]->name());
			int seqnodeJ_ID = sd.getId(currentNodes[j]->name());
			const sequence& snodeI = *sd.getSeqPtr(seqnodeI_ID,true);
			const sequence& snodeJ = *sd.getSeqPtr(seqnodeJ_ID,true);
			tempDis = _cd->giveDistance(snodeI,snodeJ,weights,&resLikelihood);
			distanceTable[i][j] = tempDis;
			LTable[i][j] = resLikelihood;
		}
	}
	if (myLog::LogLevel()>4) {
		for (i=0; i < currentNodes.size(); i++) {
			for (j=i+1; j < currentNodes.size(); j++) {
				LOG(100,<<"nj distance ["<<i<<"]["<<j<<"] ="<<distanceTable[i][j]<<endl);
			}
		}
	}
	//if (myLog::LogLevel()>4) {
	//	for (i=0; i < currentNodes.size(); i++) {
	//		for (j=i+1; j < currentNodes.size(); j++) {
	//			LOG(4,<<"nj likelihood for distance["<<i<<"]["<<j<<"] ="<<LTable[i][j]<<endl);
	//		}
	//	}
	//}
	// for printing stuff later.
	for (int tmp1=0; tmp1<distanceTable.size();++tmp1)
	for (int tmp2=0; tmp2<distanceTable.size();++tmp2) 
	_startingDistanceTable[tmp1][tmp2] = distanceTable[tmp1][tmp2];
}






void NJalg::NJiterate() {
	getMmatrixFromDistanceTable();
	int minRaw,minCol;
	findMinM(minRaw,minCol);
	
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(nodeI,nodeJ);

	//CHECK4

	updateBranchDistance(theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6

	UpdateDistanceTableAndCurrentNodes(nodeI,nodeJ,theNewNode,minRaw,minCol);
}

		
		












//CHECK1
//cout<<"\n-----------------------------------------------"<<endl;
//for (int h=0; h < currentNodes.size(); h++) cout<<currentNodes[h]->name()<<" = "<<h<<endl;

//CHECK2
//	for (int i =0; i < r_values.size();++i) cout<<"r["<<i<<"] = "<<r_values[i]<<endl;

//CHECK3
//	for (i =0; i < currentNodes.size();++i) 
//		for (int j =i+1; j <currentNodes.size();++j) 
//			cout<<"M["<<i<<"]["<<j<<"] = "<<Mmatrix[i][j]<<endl;

//CHECK4
//	string htuname = "HTU";
//	char k = 'a'+currentNodes.size();
//	htuname+=k;
//	theNewNode->SetName(htuname);		
	
//CHECK5
//_myET->getRoot()->SetName("RootOfStar");		
	
//CHECK6
//	et.output(cout,et.getRoot(),tree::ANCESTOR);

	



*/
