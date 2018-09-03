// class for H-Matrix
// 3 structs: a) R-k Matrix b) Full Matrix c) SuperMatrix
//! This class can be used for a block cluster tree.
#ifndef HMAT_H
#define HMAT_H

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <vector>
#include "block_cluster.h"

/// Three structs for handling the blocks during the partition process. The structs are described below:
/// rkmat: used for handling R-K Matrix blocks.
/// fullmat: used for handling leaves of the block cluster tree which are stored as full matrices.
/// supermat: used for handling the blocks in process of partitioning.

///rkmat:
/// k: expected rank of rk block (same as input rank 'r').
/// kt: rank of the rk block
/// a: stores the column vectors.
/// b: stores the row vectors.
/// Note: both 'a' and 'b' are stored as row vectors (due to some error in Eigen library). Care to be taken when performing operations on rk blocks.
struct rkmat
{
	int k;
	int kt;
	std::vector<Eigen::VectorXd > a;
	std::vector<Eigen::VectorXd > b;
};

///fullmat:
/// Eigen sparse matrix for holding dense blocks
struct fullmat
{
	Eigen::SparseMatrix<double>* m;
};

struct supermat
{
	int type; // 1 == rk- matrix; 2 == full matrix; 3 == supermatrix (internal node)
	int rows,cols; // rows and cols of this supermatrix
	// depending on the type, other two pointers are set to NULL
	rkmat* r;
	fullmat* f;
	std::vector<supermat*> s; // contains children of a particular node
};

/// Helper for custom_sort. Tracks the original index of the element after sorting the array.
struct new_el
{
    double val;
    int idx;
};

/// The class contains a pointer to the root of the block cluster tree and consequently, the H-Matrix is constructed by recursively traveling down the block cluster tree.
class hmat
{
private:
	supermat* root;
public:
	hmat();
	/// Custom constructor which uses block cluster tree and matrix to build the H-Matrix.
	hmat(bctree&, Eigen::SparseMatrix<double>*, int);
	/// Helper function for constructing H-Matrix. We traverse the block cluster tree and mark each node as R-K, Full or Super matrix.
	supermat* create_hmat(bctree&, Eigen::SparseMatrix<double>*, int);
	void CA_partial_pivot(Eigen::MatrixXd&, rkmat*, int);
};

/// Helper function for Cross-Approximation partial pivoting algorithm.
/// Finds the next row maximizer index based on the Partition set.
int find_index(Eigen::VectorXd&, std::vector<int>&);
bool custom_sort(struct new_el, struct new_el);

#endif
