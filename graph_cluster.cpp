/// \file graph_cluster.cpp
/// \brief Class for storage and manipulation of graphs during H-Matrix build process.
#include <iostream>
#include <cmath>
#include "graph_cluster.h"

using namespace std;
using namespace Eigen;

typedef vector<unsigned int> Vt;

// removes element from vector
Vt remove_el(Vt v, int x1,int x2)
{
	//cout<<"Currently inside remove_el: "<<x1<<" "<<x2<<endl;
	Vt dum_v;
	dum_v = v;
	dum_v.erase(remove(dum_v.begin(), dum_v.end(), x1), dum_v.end());
	if (x2>=0)
		dum_v.erase(remove(dum_v.begin(), dum_v.end(), x2), dum_v.end());
	//cout<<"Exiting remove_el!"<<endl;
	return dum_v;
}


graph_cluster::graph_cluster(void)
{
	mat_ptr=NULL;
	n_clusters=0;
	n_single_clusters=0;
}


graph_cluster::graph_cluster(SparseMatrix<double>* dum_ptr)
{
	mat_ptr = dum_ptr;
	n_clusters = 0;
	n_single_clusters=0;
}

// priority match algorithm for coarsening process
void graph_cluster::priority_match(Vt gp1,Vt gp2)
{
	//cout<<"Priority matching loaded!"<<endl;
	Vt index_set1,index_set2;

	index_set1 = gp1;
	index_set2 = gp2;

	n_clusters=0;

	while(!index_set1.empty())
	{
		unsigned int s = *(index_set1.begin()); // pick a s in V1
		// union of V1 and V2
		Vt union_index_set;
		union_index_set = index_set1;
		union_index_set.insert(union_index_set.end(),index_set2.begin(),index_set2.end());

		unsigned int t;
		// find a cluster corresponding to s
		//cout<<"Removing task started! "<<s<<" "<<t<<endl;
		if (!this->match(s,union_index_set,&t))
		{
			// means no match found for this s
			Vt dum_set1;
			dum_set1.push_back(s);
			clusters.push_back(dum_set1);
			index_set1 = remove_el(index_set1,s);
		}
		else
		{
			// push both indicies to the clusters vector
			Vt dum_set1;
			dum_set1.push_back(s);
			dum_set1.push_back(t);
			clusters.push_back(dum_set1);
			index_set1 = remove_el(index_set1,s,t);
			index_set2 = remove_el(index_set2,t);
		}

		n_clusters+=1;
	}
	// next process the index_set2
	while(!index_set2.empty())
	{
		unsigned int s = *(index_set2.begin()); // pick a s in V1

		unsigned int t;
		// find a cluster corresponding to s

		if (!this->match(s,index_set2,&t))
		{
			// means no match found for this s
			Vt dum_set1;
			dum_set1.push_back(s);
			clusters.push_back(dum_set1);
			index_set2 = remove_el(index_set2,s);
		}
		else
		{
			// push both indicies to the clusters vector
			Vt dum_set1;
			dum_set1.push_back(s);
			dum_set1.push_back(t);
			clusters.push_back(dum_set1);
			index_set2 = remove_el(index_set2,s,t);
		}

		n_clusters+=1;
	}
	if(n_clusters==1)
		std::sort(clusters.begin()->begin(),clusters.begin()->end());

	this->create_priority_groups();
}

// HEAVY EDGE MATCHING ALGORITHM

int graph_cluster::match(unsigned int s, Vt idx_set, unsigned int* t)
{
	// if match is not found this function returns negative values
	unsigned int col = s;
	double dum_max=0.0;
	unsigned int dum_idx;
	int dum=1;
	for(SparseMatrix<double>::InnerIterator it(*mat_ptr,col);it;++it)
	{
		if(find(idx_set.begin(),idx_set.end(),it.index())!=idx_set.end())
		{
			// checks if the current index search in graph is in the input index set

			if(it.index()==s)
			{
				// do nothing
			}
			else
			{
				if(abs(it.value())>=dum_max)
				{
					dum_idx = it.index();
					dum_max = abs(it.value());
					//std::cout<<"dum_max: "<<dum_max<<std::endl;
					//std::cout<<"dum_idx: "<<dum_idx<<std::endl;
				}
			}
		}
		else
		{
			// do nothing because the current index is not in the input index set, which means that pairing is not possible
		}
	}

    if(dum_max==0.0)
        dum=0;

    *t = dum_idx;

	return dum;

}

// creates priority groups to be used in the next graph's priority matching

void graph_cluster::create_priority_groups(void)
{
	// compute sizes of both priority groups
	int first_cluster_size;
	if (n_clusters%2==0)
		first_cluster_size = n_clusters/2;
	else
		first_cluster_size = n_clusters/2 + 1;
	int second_cluster_size = n_clusters/2;

	// compute cluster sizes
	std::vector<int> cluster_sizes;
	for (std::vector<Vt>::iterator itr=clusters.begin();itr!=clusters.end();++itr)
	{
		cluster_sizes.push_back((*itr).size());
	}

	// find cluster with size=1
	int idx=-1;
	for (std::vector<int>::iterator itr=cluster_sizes.begin();itr!=cluster_sizes.end();++itr)
	{
		if (*itr==1)
		{
			idx = itr - cluster_sizes.begin();
			break;
		}
	}

	if (idx<0)
	{
		// no single element cluster
		// process first priority group index set
		for (unsigned int i=0;i<clusters.size();i++)
		{
			if(groups.group1.size()<first_cluster_size)
				groups.group1.push_back(i);
			else
				groups.group2.push_back(i);
		}
	}
	else
	{
		// first add single element cluster index to group1
		groups.group1.push_back(idx);
		for (unsigned int i=0;i<clusters.size();i++)
		{
			if(i!=idx)
			{
				if(groups.group1.size()<first_cluster_size)
					groups.group1.push_back(i);
				else
					groups.group2.push_back(i);
			}
		}
	}
}

// creates new graph after clustering is completed

void graph_cluster::convert_to_coarser_graph(Eigen::SparseMatrix<double>& dum_mat)
{
	// initialize a sparse matrix which holds the new coarser graph
	SparseMatrix<double> new_graph(n_clusters,n_clusters);

	if(!clusters.empty())
	{
		//iterate over each cluster to fill the sparse matrix
		for(unsigned int i=0;i<n_clusters;i++)
		{
			//cout<<"i: "<<i<<endl;
			Vt dum_set_col = clusters.at(i);
			for(unsigned int j=0;j<n_clusters;j++)
			{
				//cout<<"  it= "<<j<<endl;
				if(j==i)
				{
					// diagonal element
					new_graph.coeffRef(j,i)=1;
				}
				else
				{
					Vt dum_set_row = clusters.at(j);
					double dum_coeff =0;
					for(Vt::iterator it_col = dum_set_col.begin();it_col!=dum_set_col.end();++it_col)
					{
						for(Vt::iterator it_row = dum_set_row.begin();it_row!=dum_set_row.end();++it_row)
						{
							dum_coeff += mat_ptr->coeff(*it_row,*it_col);
						}
					}
					//cout<<"DB: dum_coeff= "<<dum_coeff<<endl;
					new_graph.coeffRef(j,i)=dum_coeff;
				}
			}
		}
		//cout<<"Operation success"<<endl;
		//cout<<Eigen::MatrixXd(new_graph)<<endl;
	}
	else
	{
		cout<<"Error in convert_to_coarser_graph function! Empty clusters!"<<endl;
		return;

	}
	dum_mat = new_graph;
}

//function overloading for convert to coarser graph
void graph_cluster::convert_to_coarser_graph(Eigen::SparseMatrix<double>& dum_mat,std::vector<std::vector<unsigned int> > clusters_v)
{
	// initialize a sparse matrix which holds the new coarser graph
	SparseMatrix<double> new_graph(clusters_v.size(),clusters_v.size());

	if(!clusters_v.empty())
	{
		//iterate over each cluster to fill the sparse matrix
		for(unsigned int i=0;i<clusters_v.size();i++)
		{
			//cout<<"i: "<<i<<endl;
			Vt dum_set_col = clusters_v.at(i);
			for(unsigned int j=0;j<clusters_v.size();j++)
			{
				//cout<<"  j: "<<j<<endl;
				if(j==i)
				{
					// diagonal element
					new_graph.coeffRef(j,i)=1;
				}
				else
				{
					Vt dum_set_row = clusters_v.at(j);
					double dum_coeff =0;
					for(Vt::iterator it_col = dum_set_col.begin();it_col!=dum_set_col.end();++it_col)
					{
						for(Vt::iterator it_row = dum_set_row.begin();it_row!=dum_set_row.end();++it_row)
						{
							dum_coeff += mat_ptr->coeff(*it_row,*it_col);
						}
					}
					//cout<<"DB: dum_coeff= "<<dum_coeff<<endl;
					new_graph.coeffRef(j,i)=dum_coeff;
				}
			}
		}
		//cout<<"Operation success"<<endl;
		//cout<<Eigen::MatrixXd(new_graph)<<endl;
	}
	else
	{
		cout<<"Error in convert_to_coarser_graph function! Empty clusters!"<<endl;
		return;

	}
	dum_mat = new_graph;
}

// sets matrix pointer for a graph cluster object
void graph_cluster::set_matrix(Eigen::SparseMatrix<double>* dum_ptr)
{
	//cout<<"DB: pointer updated successfully"<<endl;
	//cout<<"DB1------>before->\n"<<Eigen::MatrixXd(*dum_ptr)<<endl;
	mat_ptr = dum_ptr;
	//cout<<"DB1------>after->\n"<<Eigen::MatrixXd(*mat_ptr)<<endl;
}

// gets number of cluster of graph_cluster object

int graph_cluster::get_n_clusters(void)
{
	return n_clusters;
}

// return cluster number: x
Vt graph_cluster::get_cluster(unsigned int x)
{
	return clusters.at(x);
}
// returns priority groups

Vt graph_cluster::get_priority_group1(void)
{
	return groups.group1;
}

Vt graph_cluster::get_priority_group2(void)
{
	return groups.group2;
}

// set clusters of graph_cluster object
void graph_cluster::set_clusters(std::vector<std::vector<unsigned int> >& x)
{
	clusters = x;
}

//edge weight of two given verticies
double graph_cluster::edge_weight(int i1, int i2)
{
	return mat_ptr->coeff(i1,i2);
}

// print overloading

ostream& operator<<(ostream& os, graph_cluster& gc)
{
	os<<"-----------------------------------------------------"<<"\n";
	os<<"Matrix_Graph:\n";
	os<<Eigen::MatrixXd(*(gc.mat_ptr))<<"\n";
	os<<"Clusters "<<"(n = "<<gc.n_clusters<<"): ";
	if(!gc.clusters.empty())
	{
		for(std::vector<Vt>::iterator itr=gc.clusters.begin();itr!=gc.clusters.end();++itr)
		{
			Vt dum_set1 = *itr;
			os<<"{";
			for(Vt::iterator itr1 = dum_set1.begin();itr1!=dum_set1.end();++itr1)
			{
				os<<*itr1<<",";
			}
			os<<"\b";
			os<<"},";
		}
	}
	os<<"\b"<<" ";
	os<<"\n";
	os<<"Priority Groups \n";
	os<<"Group 1  ";
	if(!gc.clusters.empty())
	{
		for(Vt::iterator itr2=gc.groups.group1.begin();itr2!=gc.groups.group1.end();++itr2)
		{
			os<<*itr2<<" ";
		}
		os<<"\n";
		os<<"Group 2  ";
		for(Vt::iterator itr2=gc.groups.group2.begin();itr2!=gc.groups.group2.end();++itr2)
		{
			os<<*itr2<<" ";
		}
		os<<"\n";
	}
	os<<"-----------------------------------------------------"<<"\n";
	return os;
}
