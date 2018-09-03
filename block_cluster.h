// struct for holding block cluster tree: as the tree may be sparse, vectors are not space efficient
//! This class can be used for a block cluster tree.
#ifndef BLOCK_H
#define BLOCK_H

#include "tree.h"
#include "graph_cluster.h"
#include <vector>

/// "bct_node" represents a node in the block cluster tree and it consists of the following attributes:
/// 1. cluster1: vector to hold data from the 1st set used for cartesian product.
/// 2. cluster2: vector to hold data from the 2nd set used for cartesian product.
/// 3. left_left: a pointer to the left most child.
/// 4. left: a pointer to the left child.
/// 5. right: a pointer to the right child.
/// 6. right_right: a pointer to the right most child.
/// 7. type: rk or full or to be split.
struct bct_node
{
	node* cluster1;
	node* cluster2;
	bct_node* left_left;
	bct_node* left;
	bct_node* right;
	bct_node* right_right;
	int type;
};

class bctree
{
private:
	bct_node* root;
public:
	bctree();
	/// Creates the block cluster tree using cluster tree, graphs and number of cols as input.
	void block_cluster(tree&, std::vector<graph_cluster*>&, int);
	/// Helper function; returns a pointer to the root of the block cluster tree.
	bct_node* get_root(void);
	void output();
	friend std::ostream& operator<<(std::ostream& os, bctree& gc);
};

#endif
