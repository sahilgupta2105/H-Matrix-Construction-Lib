/// \file block_cluster.cpp
/// \brief Class for storage and manipulation of block cluster trees (quad-tree).

#include "block_cluster.h"
#include <queue>
#include <iostream>
#include <fstream>

bctree::bctree()
{
	bct_node* ptr = new bct_node;
	ptr->cluster1 = NULL;
	ptr->cluster2 = NULL;
	ptr->left = NULL;
	ptr->left_left = NULL;
	ptr->right = NULL;
	ptr->right_right = NULL;
	ptr->type=3;
	root = ptr;
}

void bctree::block_cluster(tree& bt, std::vector<graph_cluster*>& graphs, int leaf_size)
{
    int verbose=0;
	// initialize a queue for traversal of block cluster tree as it is created
	root->cluster1 = bt.get_root();
	root->cluster2 = bt.get_root();

	std::queue <bct_node*> bct_nodes;
	bct_nodes.push(root);
	bct_node* current_node = new bct_node;
	while(!bct_nodes.empty())
	{
		if (verbose)
            std::cout<<"----->>---------->>----------"<<std::endl;
		current_node = bct_nodes.front();
		bct_nodes.pop();
		// clus1 and clus2 are the two nodes of cluster tree for current node in bct
		node* clus1 = current_node->cluster1;
		node* clus2 = current_node->cluster2;


		// check the size clusters: leaf-size condition
		int size1 = clus1->data.size();
		int size2 = clus2->data.size();
		//std::cout<<"DB: size1: "<<size1<<" size2: "<<size2<<std::endl;
		if(size1 <= leaf_size || size2 <= leaf_size)
		{
			// this is a dense node: skip it
			if (verbose)
                std::cout<<"Dense Block!"<<std::endl;
			//std::cout<<"bt_idx1: "<<clus1->bt_idx<<" bt_idx2: "<<clus2->bt_idx<<std::endl;
			current_node->type = 2;
		}
		else
		{
			// Admissibility Condition
			// Check if connected in graph or not
			if(clus1->level!=clus2->level)
			{
				std::cout<<"Block_Cluster: How can the two nodes be from different levels?"<<std::endl;
				break;
			}

			int connection=0;

			if(clus1->bt_idx == clus2->bt_idx)
				connection = 1;
			else
			{
			    unsigned int dum_el = graphs.size() - clus1->level;
				graph_cluster* current_graph = graphs.at(dum_el);
				connection = current_graph->edge_weight(clus1->bt_idx, clus2->bt_idx);
			}
			if (verbose)
                std::cout<<"connection: "<<connection<<std::endl;
			if(connection)
			{
				// Inadmissible Block

				// retrieve children of both the nodes
                if (verbose)
                    std::cout<<"Inadm. block!"<<std::endl;
				current_node->type=3;
				if (verbose)
                    std::cout<<"bt_idx1: "<<clus1->bt_idx<<" bt_idx2: "<<clus2->bt_idx<<std::endl;

				std::vector<node*> clus1_child;
				clus1_child.push_back(clus1->left);
				clus1_child.push_back(clus1->right);
				std::vector<node*> clus2_child;
				clus2_child.push_back(clus2->left);
				clus2_child.push_back(clus2->right);
				//std::cout<<"size(bct_nodes): "<<bct_nodes.size()<<std::endl;
				// cartesian product of these children
				// as cluster tree is a binary tree, every node must have atmost 2 children
				if (verbose)
                    std::cout<<"child_size1: "<<clus1_child.size()<<" child_size2: "<<clus2_child.size()<<std::endl;
				for(std::vector<node*>::iterator c2 = clus2_child.begin(); c2!= clus2_child.end(); ++c2)
				{
					for(std::vector<node*>::iterator c1 = clus1_child.begin(); c1!=clus1_child.end(); ++c1)
					{
						//std::cout<<"cp"<<std::endl;
						//std::cout<<"	size(bct_nodes): "<<bct_nodes.size()<<std::endl;
						bct_node* dum_ptr = new bct_node;
						dum_ptr->cluster1 = *c1;
						dum_ptr->cluster2 = *c2;
						dum_ptr->left_left = NULL;
						dum_ptr->left = NULL;
						dum_ptr->right = NULL;
						dum_ptr->right_right = NULL;
						dum_ptr->type=3;
						if(*c1==NULL || *c2==NULL)
							continue;
						if(current_node->left_left == NULL)
						{
                            //std::cout<<"1"<<std::endl;
							current_node->left_left = dum_ptr;
							bct_nodes.push(current_node->left_left);
						}
						else if(current_node->left == NULL)
						{
                            //std::cout<<"2"<<std::endl;
							current_node->left = dum_ptr;
							bct_nodes.push(current_node->left);
						}
						else if(current_node->right == NULL)
						{
                            //std::cout<<"3"<<std::endl;
							current_node->right = dum_ptr;
							bct_nodes.push(current_node->right);
						}
						else if(current_node->right_right == NULL)
						{
                            //std::cout<<"4"<<std::endl;
							current_node->right_right = dum_ptr;
							bct_nodes.push(current_node->right_right);
						}
					}
				}
				if (verbose)
                    std::cout<<"Cartesian product completed!"<<std::endl;
				//std::cout<<"size(bct_nodes): "<<bct_nodes.size()<<std::endl;
			}
			else
			{
				// Admissible Block
				if (verbose)
                    std::cout<<"Adm. block!"<<std::endl;
				//std::cout<<"bt_idx1: "<<clus1->bt_idx<<" bt_idx2: "<<clus2->bt_idx<<std::endl;
				current_node->type=1;
			}

		}
	}
	current_node=NULL;
	delete current_node;
}

// output bct to a file
void bctree::output(){

    std::queue <bct_node*> bct_nodes;
	bct_nodes.push(root);
	bct_node* current_node = new bct_node;
	std::ofstream myfile;
	myfile.open("block.txt");

	while(!bct_nodes.empty()){
        current_node = bct_nodes.front();
		bct_nodes.pop();
        if(current_node!=NULL){
            int node_type = current_node->type;
            if(node_type==1){
                // rk block
                myfile<<node_type<<",";
                std::vector<unsigned int> clus1, clus2;
                clus1 = current_node->cluster1->data;
                clus2 = current_node->cluster2->data;
                for(std::vector<unsigned int>::iterator itr = clus1.begin();itr!=clus1.end();++itr){
                    myfile<<*itr<<",";
                }
                myfile<<-1077<<",";
                for(std::vector<unsigned int>::iterator itr = clus2.begin();itr!=clus2.end();++itr){
                    myfile<<*itr<<",";
                }
                myfile<<"\n";
            }else if (node_type==2){
                // full block

                myfile<<node_type<<",";
                std::vector<unsigned int> clus1, clus2;
                clus1 = current_node->cluster1->data;
                clus2 = current_node->cluster2->data;
                for(std::vector<unsigned int>::iterator itr = clus1.begin();itr!=clus1.end();++itr){
                    myfile<<*itr<<",";
                }
                myfile<<-1077<<",";
                for(std::vector<unsigned int>::iterator itr = clus2.begin();itr!=clus2.end();++itr){
                    myfile<<*itr<<",";
                }
                myfile<<"\n";

            }else{
                // internal node
                bct_nodes.push(current_node->left);
                bct_nodes.push(current_node->left_left);
                bct_nodes.push(current_node->right);
                bct_nodes.push(current_node->right_right);
            }
        }
	}
	myfile.close();
	current_node = NULL;
	delete current_node;

}

// overload print function
std::ostream& operator<<(std::ostream& os, bctree& bt)
{
	os<<"-----------------------------------------------------"<<"\n";
	std::queue <bct_node*> bt_nodes;
	bt_nodes.push(bt.get_root());
	os<<"BFS of block cluster tree: \n";
	bct_node* current_node = new bct_node;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		if(current_node!=NULL)
		{
			std::vector<unsigned int> current_data = current_node->cluster1->data;
            std::vector<unsigned int> current_data1 = current_node->cluster2->data;
			//os<<" <- level: "<< current_node->level<<" - ";
			os<<" C1: <- ";
			for(std::vector<unsigned int>::iterator itr = current_data.begin(); itr!= current_data.end(); ++ itr)
			{
				os<<" "<<*itr<<" ";
			}
			os<<"-> ";
            os<<" C2: <- ";
			for(std::vector<unsigned int>::iterator itr = current_data1.begin(); itr!= current_data1.end(); ++ itr)
			{
				os<<" "<<*itr<<" ";
			}
			os<<"-> ";
			os<<"<- Type= ";
			os<<current_node->type;
			os<<"->";
			os<<" | ";
			bt_nodes.push(current_node->left_left);
            bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
            bt_nodes.push(current_node->right_right);
		}
		else
		{
			os<<" <- "<<" NULL "<<" -> ";
		}
	}
	os<<"\n";
	current_node=NULL;
	delete current_node;
	return os;
}

bct_node* bctree::get_root(void)
{
    return root;
}
