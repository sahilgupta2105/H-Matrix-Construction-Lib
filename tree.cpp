/// \file tree.cpp
/// \brief Class for storage and manipulation of binary trees during H-Matrix build process.

#include "tree.h"
#include <queue>
#include <iostream>
#include <stack>

//constructor

tree::tree()
{
	node* ptr = new node;
	root=ptr;
}

tree::tree(std::vector<unsigned int>& x)
{
	//custom constructor with 'x' as the root of binary tree
	node* ptr = new node;
	ptr->data=x;
	ptr->left=NULL;
	ptr->right=NULL;
	ptr->level=0;
	ptr->bt_idx = 0;
	root = ptr;
}

// returns pointer to root of tree
node* tree::get_root()
{
	return root;
}

// graphs to binary tree
void tree::graphs_to_tree(std::vector<graph_cluster*>& graphs)
{
	// this algorithm creates a binary tree from the set of coarsened graphs
	// root of the binary tree is already initiated
	// use BFS to fill the tree
	std::cout<<"graph to cluster loaded!"<<std::endl;
	std::queue<node*> bt_nodes;
	bt_nodes.push(root);
	node* current_node = new node;

	int current_level;
	int n_graphs = graphs.size()-1;
	unsigned int n_cluster; // cluster number
	graph_cluster* current_graph;

	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		//std::cout<<"current_node loaded!"<<std::endl;
		// check for NULL POINTERS
		if(current_node!=NULL)
		{
			//std::cout<<"loading current_level!"<<std::endl;
			current_level = current_node->level;
			//std::cout<<"n_graphs: "<<n_graphs<<std::endl;
			//std::cout<<"current_level: "<<current_level<<std::endl;
			if(current_level == n_graphs + 1)
            {
                //std::cout<<"Last level debug: "<<current_node->data.at(0)<<std::endl;
                continue;
            }
            unsigned int dum_el = n_graphs - current_level;
            unsigned int dum_el1;
			current_graph = graphs.at(dum_el); // children of this current node are found in this graph
			dum_el=0;
			n_cluster = (current_node->data).at(dum_el);
			//std::cout<<"n_clusters: "<<n_cluster<<std::endl;
			std::vector<unsigned int> current_cluster = current_graph->get_cluster(n_cluster);
			if(current_cluster.size()==2)
			{
				node* dum_node1 = new node;
				node* dum_node2 = new node;
				std::vector<unsigned int> dum_v;
				dum_el1 = current_cluster.at(dum_el);
				dum_v.push_back(dum_el1);
				//std::cout<<"data1: "<<dum_el<<"; data1 (dum_v): "<<dum_el1<<std::endl;
				dum_node1->data = dum_v;
				dum_node1->left = NULL;
				dum_node1->right = NULL;
				dum_node1->level = current_level + 1;
				dum_v.clear();
				dum_el=1;
				dum_el1 = current_cluster.at(dum_el);
				dum_v.push_back(dum_el1);
				//std::cout<<"data2: "<<dum_el<<"; data2 (dum_v): "<<dum_el1<<std::endl;
				dum_node2->data = dum_v;
				dum_node2->left = NULL;
				dum_node2->right = NULL;
				dum_node2->level = current_level + 1;
				current_node->left = dum_node1;
				current_node->right = dum_node2;
				bt_nodes.push(dum_node1);
				bt_nodes.push(dum_node2);
			}
			else
			{
				node* dum_node1 = new node;
				std::vector<unsigned int> dum_v;
				dum_el=0;
				dum_el1=current_cluster.at(dum_el);
				dum_v.push_back(dum_el1);
                //std::cout<<"data(single_node): "<<dum_el<<"; data(single_node) (dum_v): "<<dum_el1<<std::endl;
				dum_node1->data = dum_v;
				dum_node1->left = NULL;
				dum_node1->right = NULL;
				dum_node1->level = current_level + 1;
				current_node->left = dum_node1;
				current_node->right = NULL;
				bt_nodes.push(dum_node1);
			}

		}
	}
	current_node = NULL;
	delete current_node;
}

// maps leaf nodes and updates index set
void tree::map_index(std::vector<graph_cluster*>& graphs, std::vector<unsigned int>& v)
{
	// this function creates an implicit mapping of the indices as the location (index) of elements are similar to mapping
	// BFS of tree
	std::queue <node*> bt_nodes;
	bt_nodes.push(root);
	node* current_node = new node;
    //std::cout<<"DB: index mapping"<<std::endl;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
        //std::cout<<"DB1"<<std::endl;
		if(current_node!=NULL)
		{
            //std::cout<<"DB2: "<<current_node->left<<std::endl;
		    //std::cout<<"DB3: "<<current_node->right<<std::endl;
			if(current_node->left==NULL && current_node->right==NULL)
			{
				std::vector<unsigned int> current_data = current_node->data;
				unsigned int dum_el=0;
				//std::cout<<"DB4: "<<current_data.at(dum_el)<<std::endl;
				v.push_back(current_data.at(dum_el));
			}
			else
			{
				bt_nodes.push(current_node->left);
				bt_nodes.push(current_node->right);
			}
		}
	}

	current_node=NULL;
	delete current_node;

	std::cout<<"-----------------------------------------------------"<<std::endl;
	std::cout<<"Index mapping set: ";

	for(std::vector<unsigned int>::iterator itr = v.begin();itr!=v.end();++itr)
	{
        std::cout<<" << ";
		std::cout<<" "<<*itr<<" ";
        std::cout<<" >> ";
	}

	std::cout<<std::endl;

	unsigned int level_count=0; // used to maintain count of level nodes; will be useful for filtering NULL nodes
	int current_level=0;

	std::vector<unsigned int> dum_v;

	bt_nodes.push(root);

	current_node = NULL;

	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		if(current_node!=NULL)
		{
			if(current_node->level != current_level)
			{
				current_level+=1;
				level_count=0;
			}
			dum_v.clear();
			dum_v.push_back(level_count);
			current_node->data = dum_v;
			level_count+=1;
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
		}
		//std::cout<<"current_level: "<<current_level<<" level_count: "<<level_count<<std::endl;

	}
	current_node=NULL;
	delete current_node;
}

// generates cluster tree using binary tree
void tree::cluster_tree(int sz_mat)
{
	std::queue <node*> bt_nodes;
	std::stack <node*> bt_nodes_reverse;
	bt_nodes.push(root);

	node* current_node = new node;
    //std::cout<<"stack collection started."<<std::endl;
	// collect all node pointers in a stack
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes_reverse.push(current_node);
		bt_nodes.pop();

		if(current_node != NULL)
		{
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
		}
	}
	//std::cout<<"stack collection completed."<<std::endl;
	// build cluster tree using stack
	//std::cout<<"cluster tree build started."<<std::endl;
	while(!bt_nodes_reverse.empty())
	{
		current_node = bt_nodes_reverse.top();
		bt_nodes_reverse.pop();

		if(current_node!=NULL)
		{
			if(current_node->left==NULL && current_node->right==NULL)
			{
				// leaf node
				// dont do anything
				//std::cout<<"leaf node."<<std::endl;
				continue;
			}
			else
			{
			    //std::cout<<"std node."<<std::endl;
				 // collect children in a vector
				std::vector<unsigned int> child1, child2;
				if(current_node->left!=NULL)
					child1 = current_node->left->data;
				if(current_node->right!=NULL)
					child2 = current_node->right->data;
				if(child1.empty())
				{
					current_node->data = child2;
				}
				else if(child2.empty())
				{
					current_node->data = child1;
				}
				else
				{
					child1.insert(child1.end(),child2.begin(),child2.end());
					current_node->data = child1;
				}

			}
		}
	}
	std::cout<<"cluster tree build completed."<<std::endl;
	current_node=NULL;
	delete current_node;
}

//updates index corresponding to binary tree
void tree::update_bt_idx(void)
{
	std::queue <node*> bt_nodes;
	bt_nodes.push(root);
	node* current_node = new node;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		if(current_node!=NULL)
		{
			current_node->bt_idx = current_node->data.at(0);
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
		}
	}
	current_node = NULL;
	delete current_node;
}

void tree::index_tree(void)
{
	std::cout<<"-----------------------------------------------------"<<"\n";
	std::queue <node*> bt_nodes;
	bt_nodes.push(root);
	std::cout<<"BFS of binary tree: \n";
	node* current_node = new node;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		if(current_node!=NULL)
		{
			std::cout<<" <- ";
			std::cout<<" "<<current_node->bt_idx<<" ";
			std::cout<<"-> ";
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
		}
		else
		{
			std::cout<<" <- "<<" NULL "<<" -> ";
		}
	}
	std::cout<<"\n";
	current_node=NULL;
	delete current_node;
}

// overload print function
std::ostream& operator<<(std::ostream& os, tree& bt)
{
	os<<"-----------------------------------------------------"<<"\n";
	std::queue <node*> bt_nodes;
	bt_nodes.push(bt.get_root());
	os<<"BFS of binary tree: \n";
	node* current_node = new node;
	while(!bt_nodes.empty())
	{
		current_node = bt_nodes.front();
		bt_nodes.pop();
		if(current_node!=NULL)
		{
		    //std::cout<<"dbb"<<std::endl;
		    if(current_node->data.empty())
                os<<"data is NULL"<<std::endl;
			std::vector<unsigned int> current_data;
            current_data= current_node->data;
			//std::cout<<"dbb1"<<std::endl;
			os<<" <- level: "<< current_node->level<<" - ";
			os<<" <- ";
			for(std::vector<unsigned int>::iterator itr = current_data.begin(); itr!= current_data.end(); ++ itr)
			{
				os<<" "<<*itr<<" ";
			}
			os<<"-> ";
			bt_nodes.push(current_node->left);
			bt_nodes.push(current_node->right);
			current_data.clear();
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
