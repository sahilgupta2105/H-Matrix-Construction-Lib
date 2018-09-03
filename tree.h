//class for general tree: depending on application can be used as binary or quad-tree
//! This class can be used for a binary tree.
#ifndef TREE_H
#define TREE_H

#include <vector>
#include <iostream>
#include "graph_cluster.h"

/// "node" represents a node in the tree and it consists of the following attributes:
/// 1. data: vector to hold data (in this case subsets of the index set).
/// 2. left: pointer to the left child node in the tree.
/// 3. right: pointer to the right child node in the tree.
/// 4. level: level number of the node.
/// 5. bt_idx: index number of the node from index tree; used to store cluster and index tree in the same object.
struct node
{
	std::vector<unsigned int> data;
	node* left;
	node* right;
	int level;
	int bt_idx;
};

/// Class to store and manipulate index tree and cluster tree.
class tree
{
	node* root;
public:
    /// Default constructor for the tree.
	tree();
	/// Custom constructor with the input vector as the root.
	tree(std::vector<unsigned int>&);
	/// Helper function; returns a pointer to the root of the tree.
	node* get_root(void);
    /// This method uses the vector of graphs from the coarsening process to create binary tree, based on DFS traversal.
	void graphs_to_tree(std::vector<graph_cluster*>&);
	/// This method maps the leaf nodes from left to right in the input index set. Mapping set is used to permute the original matrix so that clustered rows and cols are together, which is important when constructing the block cluster tree.
	/// It also modifies the index tree according to the index set to create an ascending tree at each level.
	void map_index(std::vector<graph_cluster*>&, std::vector<unsigned int>&);
	/// Generates the cluster tree from bottom to top; stores in the same object.
	void cluster_tree(int);
	/// Updates the 'bt_idx' attribute at every node to accommodate both index and cluster tree in one object.
	void update_bt_idx(void);
	/// Prints the index tree on the console.
	void index_tree(void);
	friend std::ostream& operator<<(std::ostream& os, tree& gc);
};

#endif
