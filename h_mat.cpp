/// \file h_mat.cpp
/// \brief Class for storage and manipulation of Hierarchical Matrices.

#include "h_mat.h"
#include <queue>
#include <cstdlib>
#include <algorithm>

//deafult constructor
hmat::hmat()
{
	supermat* s = new supermat;
	root =s;
}

hmat::hmat(bctree& bct, Eigen::SparseMatrix<double>* mat, int r=10)
{
	root = create_hmat(bct,mat,r);
}

supermat* hmat::create_hmat(bctree& bct, Eigen::SparseMatrix<double>* mat, int r)
{
	// a queue is needed for traversal of block cluster tree
	std::queue<bct_node*> bct_nodes;
	bct_nodes.push(bct.get_root());
	bct_node* current_node = new bct_node;

	// this is the supermat pointer which will be returned to the calling function
	supermat* dum_root = new supermat;

	dum_root->type = bct.get_root()->type;
	dum_root->rows = (bct.get_root()->cluster1->data).size();
	dum_root->cols = (bct.get_root()->cluster2->data).size();

	// a queue is needed for traversal of set of H-Matrices
	std::queue<supermat*> hmat_nodes;
	hmat_nodes.push(dum_root);
	supermat* current_block = new supermat;
	// this stores the current matrix block information
	//Eigen::SparseMatrix<double>* current_mat = mat;
	while(!hmat_nodes.empty())
	{
		current_node = bct_nodes.front();
		bct_nodes.pop();

		current_block = hmat_nodes.front();
		hmat_nodes.pop();

        //std::cout<<"DB: "<<current_node->type<<","<<current_node->cluster1->data.at(0)<<std::endl;

		if(current_node->type==1)
		{
			// this is R-k Leaf
			unsigned int dum_el =0;
			int start_row = current_node->cluster1->data.at(dum_el);
			int start_col = current_node->cluster2->data.at(dum_el);
			int n_rows = current_node->cluster1->data.size();
			int n_cols = current_node->cluster2->data.size();
			current_block->f = NULL;
			current_block->s.clear();
			current_block->rows = n_rows;
			current_block->cols = n_cols;
			current_block->type = 1;
			Eigen::MatrixXd dum_mat = Eigen::MatrixXd(mat->block(start_row,start_col,n_rows,n_cols));
			rkmat* dum_rk = new rkmat;
			current_block->r= dum_rk;
			CA_partial_pivot(dum_mat, current_block->r, r);
			/////////////////////////////////////////////////////////////////////////////////////////////////
			//debug rk-block
			//std::cout<<"start_row, start_col, n_rows, n_cols: "<<start_row<<", "<<start_col<<", "<<n_rows<<", "<<n_cols<<std::endl;
			//std::cout<<"block being approximated:"<<std::endl;
			//std::cout<<dum_mat<<std::endl;
			//std::cout<<"rk approximation: "<<std::endl;
			//std::cout<<"a-vec"<<std::endl;
			//for(std::vector<Eigen::VectorXd >::iterator itr=current_block->r->a.begin();itr!=current_block->r->a.end();++itr)
            //{
           //     std::cout<<*itr<<std::endl;
            //}
           // std::cout<<"b-vec"<<std::endl;
			//for(std::vector<Eigen::VectorXd >::iterator itr=current_block->r->b.begin();itr!=current_block->r->b.end();++itr)
           // {
           //     std::cout<<*itr<<std::endl;
           // }
            // debug ends
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
		else if(current_node->type==2)
		{
			// this is a dense node
			// for this case we need to partition the matrix and store the block in full matrix pointer 'f'
			unsigned int dum_el =0;
			int start_row = current_node->cluster1->data.at(dum_el);
			int start_col = current_node->cluster2->data.at(dum_el);
			int n_rows = current_node->cluster1->data.size();
			int n_cols = current_node->cluster2->data.size();
			//std::cout<<"start_row, start_col, n_rows, n_cols: "<<start_row<<","<<start_col<<","<<n_rows<<","<<n_cols<<std::endl;
			fullmat* dum_f = new fullmat;
			Eigen::SparseMatrix<double>* dum_mat = new Eigen::SparseMatrix<double>;
			*dum_mat = mat->block(start_row,start_col,n_rows,n_cols);
            dum_f->m = dum_mat;
			current_block->f= dum_f;
			current_block->s.clear();
			current_block->r=NULL;
			current_block->type=2;
			current_block->rows = n_rows;
			current_block->cols = n_cols;
		}
		else if(current_node->type==3)
		{
			// this is supermatrix node
			std::vector<bct_node*> child_current_node;
            child_current_node.clear();

			if(current_node->left_left!=NULL)
            {
				child_current_node.push_back(current_node->left_left);
                //std::cout<<"DB: (type) "<<current_node->left_left->type<<std::endl;
            }

			if(current_node->left!=NULL)
            {
				child_current_node.push_back(current_node->left);
                //std::cout<<"DB: (type) "<<current_node->left->type<<std::endl;
            }

			if(current_node->right!=NULL)
            {
				child_current_node.push_back(current_node->right);
				//std::cout<<"DB: (type) "<<current_node->right->type<<std::endl;
            }

			if(current_node->right_right!=NULL)
            {
                child_current_node.push_back(current_node->right_right);
                // debug

                // BUG/Issue: the current node right_right node type value keeps changing randomly
                // No possible error detected as of yet.
                // To make the code work, run the program multiple times; works on a random iteration.

                //std::cout<<"Debug output begins..."<<std::endl;
                std::vector<unsigned int> current_data = current_node->cluster1->data;
                std::vector<unsigned int> current_data1 = current_node->cluster2->data;
                //std::cout<<" C1: <- ";
                //for(std::vector<unsigned int>::iterator itr = current_data.begin(); itr!= current_data.end(); ++ itr)
                //{
               //     std::cout<<" "<<*itr<<" ";
               // }
               // std::cout<<"-> ";
               // std::cout<<" C2: <- ";
               // for(std::vector<unsigned int>::iterator itr = current_data1.begin(); itr!= current_data1.end(); ++ itr)
               // {
                //    std::cout<<" "<<*itr<<" ";
               // }
               // std::cout<<"-> ";
               // std::cout<<"DB: (type) "<<current_node->right_right->type<<std::endl;
               // std::cout<<"Debug output ends..."<<std::endl;
            }


			current_block->r=NULL;
			current_block->f=NULL;
			current_block->type=3;
			// create supermatrix nodes corresponding to the non-NULL children of current node of bct
			for(std::vector<bct_node*>::iterator itr=child_current_node.begin(); itr!= child_current_node.end(); ++itr)
			{
				supermat* new_sp = new supermat;
				new_sp->type = (*itr)->type;
				new_sp->rows = (*itr)->cluster1->data.size();
				new_sp->cols = (*itr)->cluster2->data.size();
				current_block->s.push_back(new_sp);
				bct_nodes.push(*itr);
				hmat_nodes.push(new_sp);
			}
		}
		else
		{
			std::cout<<"Error in create_hmat: unknown type of matrix block!"<<std::endl;
			break;
		}
	}
    delete current_block,current_node;

	return dum_root;
}

__attribute__((force_align_arg_pointer)) void hmat::CA_partial_pivot(Eigen::MatrixXd& dum_mat, rkmat* rk, int r)
{
	// Cross Approximation with partial pivoting
	// input: required rank
	// ouput: 'rk' matrix is filled
	bool db=0;

	if (db)
        std::cout<<"Inside cross approx."<<std::endl;
    rk->k=r;
    //Eigen::MatrixXd dum_mat = Eigen::MatrixXd(dum_mat1);
	unsigned int current_i=0;
	unsigned int current_j=0;
	std::vector<int> collected_indicies;

	int mu = 1;
    Eigen::MatrixXd::Index max_index;
    if (db)
        std::cout<<"DB1"<<std::endl;
	while(mu<=r)
    {
        dum_mat.row(current_i).maxCoeff(&max_index);
        current_j = int(max_index);
        if (db)
            std::cout<<"DB2"<<std::endl;
        double rk_sum=0.0;

        for(unsigned int i=1;i<mu;i++)
        {
            rk_sum = rk_sum + (rk->a.at(i-1)(current_i))*(rk->b.at(i-1)(current_j));
        }
        if (db)
            std::cout<<"DB3"<<std::endl;
        double delta = dum_mat(current_i,current_j) - rk_sum;
        if (db)
            std::cout<<"DB: special"<<std::endl;
        Eigen::VectorXd a_vec = Eigen::VectorXd::Zero(dum_mat.rows(),1);
        if (db)
            std::cout<<"DB: special1"<<std::endl;
        Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(dum_mat.cols(),1);
        if (db)
            std::cout<<"DB4"<<std::endl;
        if (db)
            std::cout<<"delta: "<<delta<<std::endl;
        if (delta==0)
        {
            if (db)
                std::cout<<"inside delta==0"<<std::endl;
            if(collected_indicies.size()==dum_mat.rows()-1)
                break;
            else
                mu=mu-1;
        }
        else
        {
            if (db)
                std::cout<<"DB6"<<std::endl;
            for(unsigned int i=1;i<mu;i++)
            {
                a_vec = a_vec + (rk->a.at(i-1))*(rk->b.at(i-1)(current_j));
                b_vec = b_vec + (rk->a.at(i-1)(current_i))*(rk->b.at(i-1));
            }
            if (db)
                std::cout<<"DB7"<<std::endl;
            a_vec = dum_mat.block(0,current_j,dum_mat.rows(),1) - a_vec;
            b_vec = (dum_mat.block(current_i,0,1,dum_mat.cols()).transpose() - b_vec)/delta;
            if (db)
                std::cout<<"DB8"<<std::endl;
            rk->a.push_back(a_vec);
            rk->b.push_back(b_vec);
        }
        if (db)
            std::cout<<"DB5"<<std::endl;
        collected_indicies.push_back(current_i);
        current_i = find_index(a_vec,collected_indicies);

        if(current_i<0)
            break;

        mu = mu+1;
    }

    rk->kt = mu;
}


int find_index(Eigen::VectorXd& vec, std::vector<int>& collected_idx)
{
    vec = vec.cwiseAbs();

    struct new_el arr[vec.size()];

    for(int i=0;i<vec.size();i++)
    {
        arr[i].val = vec(i);
        arr[i].idx = i;
    }

    std::sort(arr,arr+vec.size(),custom_sort);

    int j=0;
    int next_idx=-100;

    while(1)
    {
        if(std::find(collected_idx.begin(),collected_idx.end(),arr[j].idx) == collected_idx.end())
        {
            next_idx = arr[j].idx;
            break;
        }
        else
            j=j+1;
    }
    if (next_idx<0)
        std::cout<<"Error in find_index algorithm."<<std::endl;
    return next_idx;

}

bool custom_sort(struct new_el e1, struct new_el e2)
{
    return e1.val > e2.val;
}
