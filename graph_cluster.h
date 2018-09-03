
#ifndef GRAPHCLUSTER_H
#define GRAPHCLUSTER_H

#include <vector>
#include <Eigen/SparseCore>

/// "priority_groups" holds two vectors as 'group1' and 'group2', where group1 is the higher priority (i.e. the unmatched nodes) and group2 is the lower priority group.
/// This ensures a more balanced cluster tree as explained in Fang Yang journal paper.
struct priority_groups
{
	std::vector<unsigned int> group1;
	std::vector<unsigned int> group2;
};

/// "graph_cluster" class is used for handling the graphs built during the coarsening process.
/// This class has methods which handle the coarsening process based on Modified Heavy Edge Matching (HEM) algorithm.
/// An object of 'graph_cluster' class holds the following attributes:
/// 1. A pointer to the row compressed matrix (object of SparseMatrix, Eigen) of the input matrix 'A'.
/// 2. The clusters in the current graph_cluster object.
/// 3. Total number of clusters in the object.
/// 4. Number of singleton clusters present.
/// 5. A struct object holding priority groups to be used in the next iteration.
class graph_cluster
{
private:
	Eigen::SparseMatrix<double>* mat_ptr;
	std::vector<std::vector<unsigned int> > clusters;
	int n_clusters;
	int n_single_clusters;
	priority_groups groups;
public:
    /// Default constructor for 'graph_cluster' class; creates an empty object.
	graph_cluster(void);
	/// Custom constructor for 'graph_cluster' class.
	graph_cluster(Eigen::SparseMatrix<double>* dum_ptr);
	/// Method for executing the Priority Match algorithm (source: Fang Yang journal). This method uses 'match' and 'create_priority_groups' methods of a graph_cluster object.
	void priority_match(std::vector<unsigned int>, std::vector<unsigned int>);
	/// Method finds the maximum edge adjacent to the input index ('s'). If no match is found (e.g. zero weight), then the function returns negative index.
	int match(unsigned int s, std::vector<unsigned int> idx_set, unsigned int*);
	/// Function overloading for directly printing the graph_cluster object to the console.
	friend std::ostream& operator<<(std::ostream& os, graph_cluster& gc);
	/// Method to convert the current graph to a coarser graph using the clusters obtained using HEM algorithm.
	void convert_to_coarser_graph(Eigen::SparseMatrix<double>&);
	void convert_to_coarser_graph(Eigen::SparseMatrix<double>&,std::vector<std::vector<int unsigned> >);
	/// Helper function to assign matrix to the graph_cluster object.
	void set_matrix(Eigen::SparseMatrix<double>*);
	int get_n_clusters(void);
	std::vector<unsigned int> get_cluster(unsigned int);
	/// Method to divide clusters into priority groups, with singleton sets with higher priority.
	void create_priority_groups(void);
	std::vector<unsigned int> get_priority_group1(void);
	std::vector<unsigned int> get_priority_group2(void);
	void set_clusters(std::vector<std::vector<unsigned int> >&);
	double edge_weight(int, int);
};

// function to remove element from vectors during priority matching
/// Helper function to remove element from the input vector 'v' located at 'x1' and 'x2'.
std::vector<unsigned int> remove_el(std::vector<unsigned int> v, int x1, int x2=-100);

#endif
