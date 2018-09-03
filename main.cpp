/// \file main.cpp
/// \brief This file works as the front-end interface for the user. Takes the matrix as input and outputs the corresponding H-Matrix.
/// \author Sahil Gupta (ae13b046@smail.iitm.ac.in)

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <queue>
#include <stack>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "tree.h"
#include "graph_cluster.h"
#include "block_cluster.h"
#include "h_mat.h"

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMat;
typedef MatrixXd Mat;

/// \brief This function inputs a matrix in '.csv' format.
///
/// \param Reference to the sparse matrix which holds the matix A
/// \return void
///
///
void input_matrix(SpMat&);
/// \brief This function executes the graph coarsening process based on modified HEM algorithm.
///
/// \param 'v' A vector of pointers to graph_cluster objects, which store the information about graphs at each step of the coarsening process.
/// \param 'n_cols' Number of columns (/rows) in the matrix.
/// \return void
///
//!< NOTE: The HEM algorithm is only applicable to symmetric systems.
void generate_graphs(std::vector<graph_cluster*>&, int); // function for graph coarsening process
/// \brief This function reorders the original input matrix ('A') as per the index set, using a permutation matrix.
///
/// \param 's1' the original matrix ('A')
/// \param 'idx_set' the index set as computed from the index tree.
/// \return void
///
///
void reorder_matrix(SpMat&, std::vector<unsigned int>&);
/// \brief This function creates the graphs again based on the reordered matrix. The process is not computationally intensive because priority groups need not be found again. This process is important because graphs will be needed while creating block cluster tree.
///
/// \param 'graphs' vector containing graphs from previous coarsening process.
/// \param 'bt'
/// \return
///
///
void reorder_graphs(std::vector<graph_cluster*>&, tree&);

//void generate_block_cluster_tree(bct_node*, int, tree&, tree&, std::vector<graph_cluster*>&);

int main()
{
	// create dummy matrix for testing
	//Mat m(6,6);
//	m<<0,0,0,0,
//        0,0,0,0,
//        0,0,0,0,
//        0,0,0,0;
//	 m<<1,4,0,0,0,7,
//	 	4,1,0,0,0.5,2,
//	 	0,0,1,3,0,1,
//	 	0,0,3,1,0,0,
//	 	0,0.5,0,0,1,3,
//	 	7,2,1,0,3,1;
//	 m<< 1,3,0,1,
//	 	3,1,1,0,
//	 	0,1,1,0,
//	 	1,0,0,1;
//    Mat m(8,8);
//	m<< 1,2,0,0,0.1,0.05,0.5,0.25,
//		2,1,1,1,0.1,0.05,0.5,0.25,
//		0,1,1,2,0.2,0.3,0,0,
//		0,1,2,1,0.2,0.3,0,1,
//		0.1,0.01,0.3,0.2,1,2,1,0,
//		0.1,0.01,0.1,0.1,2,1,0,0,
//		0.2,0.01,0,0,1,0,1,2,
//		0.2,0.02,0,1,0,0,2,1;
	//cout<<"Dense Matrix:\n"<<m<<endl;
	cout<<"Converting to sparse matrix!"<<endl;
//	SpMat s1 = m.sparseView();
	// convert dense matrix to sparse format
	SpMat s1(804,804);
	input_matrix(s1);
    s1.cwiseAbs();
	graph_cluster g1(&s1);

	//coarsening process starts here
	cout<<"Graph coarsening started. Step 1"<<endl;
	std::vector<graph_cluster*> graphs;
	graphs.push_back(&g1);
	cout<<"Graph coarsening process. Step 2"<<endl;
	generate_graphs(graphs,s1.cols());
	std::cout<<"Graph coarsening completed. Step 3"<<std::endl;
	// coarsening process completes here

	// create tree from indices of clusters
	std::vector<unsigned int> dum_v;
	dum_v.push_back(0);
	tree bt(dum_v);
	bt.graphs_to_tree(graphs);
	cout<<"-----------------------------------------------------"<<endl;
	cout<<"Binary Tree corresponding to coarsened graphs completed."<<endl;
	//cout<<bt;

	// index mapping
	cout<<"Reordering process started."<<endl;
	std::vector<unsigned int> idx_set; // INDEX SET: this set corresponds to leaves of the above tree from left to right
	bt.map_index(graphs, idx_set);
	cout<<"-----------------------------------------------------"<<endl;
	cout<<"Reordering of Binary tree completed."<<endl;
	//cout<<bt;
	cout<<"-----------------------------------------------------"<<endl;

	cout<<"-----------------------------------------------------"<<endl;
	// permute the matrix as per the index set
	reorder_matrix(s1,idx_set);
	cout<<"Reordering of matrix completed."<<endl;
	//cout<<MatrixXd(s1)<<endl;
	cout<<"-----------------------------------------------------"<<endl;

	cout<<"-----------------------------------------------------"<<endl;
	// generate graphs from reordered matrix
	reorder_graphs(graphs, bt);
	cout<<"Reordering of graphs completed. "<<endl;
//	for(std::vector<graph_cluster* >::iterator itr=graphs.begin();itr!=graphs.end();++itr)
//	{
//		cout<<*(*itr);
//	}
//	cout<<endl;
	cout<<"-----------------------------------------------------"<<endl;

	//update binary index so that cluster and binary index tree are stored in same object
	cout<<"-----------------------------------------------------"<<endl;
	cout<<"Index Tree created."<<endl;
	// creating cluster tree
	bt.update_bt_idx();
	//bt.index_tree(); // prints the index tree to console
	cout<<"-----------------------------------------------------"<<endl;


	cout<<"-----------------------------------------------------"<<endl;
	cout<<"Cluster Tree created."<<endl;
	// creating cluster tree
	bt.cluster_tree(s1.cols());
	//cout<<bt;
	cout<<"-----------------------------------------------------"<<endl;

	cout<<"-----------------------------------------------------"<<endl;
	cout<<"Block Cluster Tree created. "<<endl;
	bctree bct;
	int leaf_size = 80;
	bct.block_cluster(bt, graphs, leaf_size);
	//cout<<bct<<endl;
	bct.output();
	cout<<"-----------------------------------------------------"<<endl;

	cout<<"-----------------------------------------------------"<<endl;
	cout<<"H-Matrix successfully created. "<<endl;
   	hmat hMatrix(bct, &s1, 1);
   	cout<<"-----------------------------------------------------"<<endl;
}

void input_matrix(SpMat& sm)
{
    ifstream ip("matrix.csv");

    if(!ip.is_open()) cout<<"Error: file open"<<endl;

    int row, col;
    double val;
    string line;
    int line_n=1;
    while(ip.good())
    {
        getline(ip,line);
        //cout<<line_n<<endl;
        istringstream ss(line);
        string substr;
        int i=0;
        while(getline(ss,substr,'\t'))
        {
            if (i==0)
                row= stoi(substr)-1;
            else if(i==1)
                col = stoi(substr)-1;
            else if(i==2)
                val = stod(substr);
            i+=1;

        }
        sm.insert(row,col) = val;
        line_n+=1;
    }

    std::cout<<"Input success: Matrix dimensions "<<sm.rows()<<","<<sm.cols()<<std::endl;
    ip.close();
}

void generate_graphs(std::vector<graph_cluster*>& v, int n_cols)
{
    //! PRE-PROCESSING: The input vector contains only the original matrix and no priority groups exists. So, the Index Set is split into two parts at the middle based on the number of rows (/cols). These groups serve as priority groups for the 1st iteration. The Priority Match algorithm is executed using these groups as input.
	std::vector<unsigned int> gp1,gp2;
	if(n_cols%2==0)
	{
		for(unsigned int i=0;i<n_cols/2;i++)
		{
			gp1.push_back(i);
		}
		for(unsigned int i=n_cols/2;i<n_cols;i++)
		{
			gp2.push_back(i);
		}
	}
	else
	{
		for(unsigned int i=0;i<1+n_cols/2;i++)
		{
			gp1.push_back(i);
		}
		for(unsigned int i=1+n_cols/2;i<n_cols;i++)
		{
			gp2.push_back(i);
		}
	}
	v.back()->priority_match(gp1,gp2);
	//cout<<*v.back();

	std::vector<SpMat*> matrices;

    //! Process:
    //! We keep iterating until the number of clusters left in the graph_cluster object is '1'.
    //! Memory is allocated using the new operator to store the coarsen graph at next step.
    //! convert_to_coarser_graph is used to compute the coarsened graph.
    //! Memory is allocated using the new operator for a new graph_cluster object.
    //! This object is added to the vector passed as an input to the generate_graphs function.
	int n=0;
	while(n!=2)
	{
		n = v.back()->get_n_clusters(); // get number of clusters from the last graph
		//cout<<"Graph Clustering Process Step 2. Number of nodes: "<<n<<endl;
		SpMat* s = new SpMat(n,n); // initialize a new sparse matrix to store the next coarse graph
		v.back()->convert_to_coarser_graph(*s); // computes the next coarse graph: in form of matrix
		graph_cluster* g = new graph_cluster; // initialize a new graph cluster object
		g->set_matrix(s); // add matrix to the newly initialized object
//		if(n<20)
//            cout<<*s<<endl;
		gp1 = v.back()->get_priority_group1();
		gp2 = v.back()->get_priority_group2();
		g->priority_match(gp1,gp2);
		v.push_back(g);
		//cout<<*g;
	}
}

void reorder_matrix(SpMat& s1, std::vector<unsigned int>& idx_set)
{
	// generate permutation matrix
	SpMat pMat(s1.cols(),s1.cols());
	pMat.reserve(s1.cols());
	int dum_ctr=0;
	for(std::vector<unsigned int>::iterator itr=idx_set.begin();itr!=idx_set.end();++itr)
	{
		pMat.insert(dum_ctr,*itr) = 1;
		dum_ctr+=1;
	}

	s1 = pMat*s1*pMat.transpose();
}

void reorder_graphs(std::vector<graph_cluster*>& graphs, tree& bt)
{
	// iterate over the original graphs to calculate reordered graphs
	// FIVE properties to be updated: 1. Matrix(/graph); 2. Clusters; 3. n_clusters;

	// reverse BFS: use queue and stack
	queue <node*> bt_nodes;
	bt_nodes.push(bt.get_root());
	stack <node*> bt_nodes_reverse;

	node* current_node = new node;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();

        //cout<<"DB1"<<endl;
		if(current_node!=NULL)
		{
		    //cout<<"DB2"<<endl;
			if(current_node->left==NULL && current_node->right==NULL)
			{
			    //cout<<"inside if statement"<<endl;
			    break;
			}
            //cout<<"DB3"<<endl;
			bt_nodes_reverse.push(current_node);
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
		}
	}
	//cout<<"successful exit"<<endl;
	// the nodes of cluster tree have been collected in reverse order (one level before leaves to root) in the stack
	std::vector<std::vector<unsigned int> > dum_clusters;
	int current_level = graphs.size()-1;
	//cout<<"current_level: "<<current_level<<endl;
	//cout<<bt_nodes_reverse.size()<<endl;
	for(std::vector<graph_cluster*>::iterator itr= graphs.begin();itr!=std::prev(graphs.end());)
	{
		// collect clusters for this level
		dum_clusters.clear();
		current_node = bt_nodes_reverse.top();
		bt_nodes_reverse.pop();
		while(current_node->level == current_level)
		{
			//cout<<"current_level: "<<current_level<<endl;
			std::vector<unsigned int> dum_v;
			unsigned int dum_el =0;
			if(current_node->left!=NULL)
                dum_v.push_back(current_node->left->data.at(dum_el));
			if(current_node->right!=NULL)
				dum_v.push_back(current_node->right->data.at(dum_el));
			dum_clusters.push_back(dum_v);
			current_node = bt_nodes_reverse.top();
			bt_nodes_reverse.pop();
		}
		bt_nodes_reverse.push(current_node);
		std::reverse(dum_clusters.begin(),dum_clusters.end());
		// cout<<"dum_clusters-size: "<<dum_clusters.size()<<endl;
		// cout<<"//////"<<endl;
		// for(std::vector<std::vector<int> >::iterator itr1=dum_clusters.begin();itr1!=dum_clusters.end();++itr1)
		// {
		// 	cout<<" <- ";
		// 	for(std::vector<int>::iterator itr2=(*itr1).begin();itr2!=(*itr1).end();++itr2)
		// 	{
		// 		cout<<" "<<*itr2<<" ";
		// 	}
		// 	cout<<" -> ";
		// }
		// cout<<"///////"<<endl;

		(*itr)->set_clusters(dum_clusters);
		// coarsening of graph
		SpMat* dum_mat = new SpMat(dum_clusters.size(),dum_clusters.size());
		(*itr)->convert_to_coarser_graph(*dum_mat, dum_clusters);
		//cout<<"error"<<endl;
		std::advance(itr,1);
		(*itr)->set_matrix(dum_mat);
		current_level-=1;
	}
	current_node=NULL;
	delete current_node;
}
