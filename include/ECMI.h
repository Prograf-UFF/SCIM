#pragma once

#include <fstream>
#include <vector>
#include <set>

#pragma GCC system_header
#include <Eigen/Core>

#include <unordered_map>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "leitorbasenumerica.h"

#include "Cluster.h"
#include "FpTree.h"
#include "CFITree.h"

// CONSTANT PARAMETERS 
const double LIMITE_INFERIOR_DELTA = 1e-8; // using for eliminate irrelevant dimension
const double MAXIMUM_ANGLE_DEGREE = 90.0; // using to find worts formation of clusters give a item
typedef double real_type;


typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>	dynamic_size_matrix;
typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1>					column_vector_type;
typedef Eigen::Matrix<real_type, 1, Eigen::Dynamic>					row_vector_type;

class ECMI
{
public:
    ECMI(const std::string& _strFile, const std::string& _strFileOut);
	~ECMI();

	void runAlgorithm(real_type _correlation_precision);

	int getQtdItemsetGerados() const { return qtdItemset; }

	dynamic_size_matrix getXprojected() { return x_projected; }

	dynamic_size_matrix getDistanceMatrix() { return distance_matrix; }

	dynamic_size_matrix getAngleMatrix() { return angle_par_item; }

	real_type getRadiusDistanceCluster(size_t id) const { return radius_distance_cluster[id]; }
	
	real_type getRadiusWortsDistanceCluster(size_t id) const { return worts_case_cluster[id]; }

	std::vector<size_t> getSortByDistance(size_t idItem) const ;

	real_type getDistanceParItem(size_t idItemA, size_t idItemB) const { return distance_matrix(idItemA, idItemB); }

	real_type getAngleParItem(size_t idItemA, size_t idItemB) const { return angle_par_item(idItemA, idItemB); }
	
private:
	
    void calculateDualScaling(BASE_NUM& _dado);

	void calculateDistanceMatrixChiSquare();

	void calculateAngularMatrix();

	void generateOverlapCluster();

    void generateFpTree(BASE_NUM& _dado);

    void fpGrowth(std::shared_ptr<FPTree> _tree, std::set<int>& _prefix, const std::map<int, int>& _mapSupport);

    void allCombinationsOfPrefixPath(std::shared_ptr<FPTree> _tree, std::set<int>& _prefix, const std::map<int, int>& _mapSupport);

    void saveItemset(std::vector<int> _itemset, int _countOcurrence);

    void sortedOriginalOrder(std::vector<int>& _itemset) const;

private:

    std::shared_ptr<FPTree> m_pTree;

    std::shared_ptr<CFITree> m_pCFITree;

    std::ofstream m_fileOutput;

     std::map<int,int> mapIdSupport;

    int nTransactions;

    real_type correlation_precision;

	// dualscaling atributes
	dynamic_size_matrix x_normed, x_projected;
	dynamic_size_matrix y_normed, y_projected;
	row_vector_type rho, delta, fc;
	column_vector_type fr;
	
	// distance metric (chi_square) , distance for each par of item
	dynamic_size_matrix distance_matrix;
	
	// radius distance of each cluster ( each item generates a cluster)
	row_vector_type radius_distance_cluster;

	// angle of pair of item
	dynamic_size_matrix angle_par_item;
	row_vector_type worts_case_cluster;

	// item cluster generate
    Clusters vector_cluster;

    long int qtdItemset; // quantity  itemset

     map<int, int> 		m_cabecalhoBaseIdReal;

     BASE_NUM base;

     int countRecursao;
};

