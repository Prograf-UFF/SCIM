/*
Copyright 2018 Altobelli de Brito Mantuan.
Distributed under the GNU General Public License, Version 3.
Last update : August 25th, 2018.
*/

#include "ECMI.h"
#include "dual_scaling.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>
#include <fstream>

typedef ds::multiple_choice_data<dynamic_size_matrix>	mcd;

string printVectorToString(std::vector<size_t>& r) {
    string out;
    for (size_t i = 0; i < r.size(); i++)
        out += to_string( r[i] + 1) + ",";

    if(!out.empty())
        out.pop_back();

    return out;
}

void print(const std::vector<size_t>& r) {
    for (size_t i = 0; i < r.size(); i++)
        std::cout << "\t" << r[i];
    std::cout << std::endl;
}

void print(const std::vector<int>& r) {
    for (size_t i = 0; i < r.size(); i++)
        std::cout << "\t" << r[i] + 1;
    std::cout << std::endl;
}

void print(const std::set<int>& r) {
    for (auto i = r.begin(); i != r.end(); ++i)
        std::cout << "\t" << *i + 1;
    std::cout << std::endl;
}


ECMI::ECMI(const std::string& _strFile, const std::string& _strFileOut)
   : m_pTree(new FPTree()), m_pCFITree(new CFITree()), m_fileOutput(_strFileOut), qtdItemset(0){

    countRecursao = 0;

    auto start = std::chrono::high_resolution_clock::now();

    base = BASE_NUM();
    if(!LeitorBaseNumerico::obterDadoArquivo(_strFile, base)){
        std::cout << "Erro read file!" << std::endl;
        return;
    }

    auto finish = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = finish - start;

    nTransactions = base.getSizeTransation();

    calculateDualScaling(base);

    calculateDistanceMatrixChiSquare();

	calculateAngularMatrix();


}

ECMI::~ECMI(){
    m_fileOutput.close();
}

void ECMI::calculateDualScaling(BASE_NUM& _dado){

    ds::dual_scaling(mcd(_dado.getMatrix()), x_normed, y_normed, x_projected, y_projected, rho, delta, fc, fr, LIMITE_INFERIOR_DELTA);
}

void ECMI::calculateDistanceMatrixChiSquare() {

	const int qtdItens = x_projected.rows();
	const int qtdDimension = x_projected.cols();

	distance_matrix.resize(qtdItens, qtdItens);
	distance_matrix = Eigen::MatrixXd::Zero(qtdItens, qtdItens);

	for (int i = 0; i < qtdItens; i++) {
		for (int j = 1 + i; j < qtdItens; j++) {
			double quadradoDistancia = 0.0;
			for (int aux = 0; aux < qtdDimension; aux++) {
				double delta = (x_projected(i, aux) / sqrt(fc(i)/fr.size())) - (x_projected(j, aux) / sqrt(fc(j)/fr.size()));
				quadradoDistancia += (delta * delta)*rho(aux);
			}
            distance_matrix(i, j) = quadradoDistancia;
			distance_matrix(j, i) = distance_matrix(i, j);
		}
	}

	radius_distance_cluster.resize(qtdItens);

    for (int i = 0; i < qtdItens; i++) {
		double quadradoDistancia = 0.0;
		for (int aux = 0; aux < qtdDimension; aux++) {
			double delta = x_projected(i, aux) / sqrt(fc(i) / fr.size());
			quadradoDistancia = quadradoDistancia + (delta * delta)*rho(aux);
		}
		
        radius_distance_cluster(i) = quadradoDistancia;
	}

}

std::vector<size_t> sort_indexes(const row_vector_type &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

   // sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

template<typename Lambda>
std::vector<size_t> findItemSortIndex(const std::vector<size_t>& _sortIndex, Lambda _func ) {
	std::vector<size_t> index;
	for (auto i : _sortIndex)
		if (_func(i))
			index.push_back(i);
		else
			return index;
	
	return index;
}

double acosd(real_type x){

	if (x < -1.0) x = -1.0;
	else if (x > 1.0) x = 1.0;

	return acos(x) / M_PI * 180;
}

void ECMI::calculateAngularMatrix() {
	const int qtdItens = x_projected.rows();
	
	// calculate the angle pair of item
	angle_par_item.resize(qtdItens, qtdItens);
	angle_par_item = Eigen::MatrixXd::Zero(qtdItens, qtdItens);

	for (int i = 0; i < qtdItens; i++) {
		const auto& a = x_projected.row(i);
		for (int j = 1 + i; j < qtdItens; j++) {
			const auto& b = x_projected.row(j);
			
            angle_par_item(i, j) = acosd((a.cwiseProduct(b)).sum() / (a.norm()*b.norm()));
			angle_par_item(j, i) = angle_par_item(i, j);
		}
	}

    // calculate de maximun distance of each cluster, using radius <= 90º
	worts_case_cluster.resize(qtdItens);

	for (int i = 0; i < qtdItens; i++) {
		const row_vector_type& col_angle = angle_par_item.row(i);
		
		const std::vector<size_t>& I = sort_indexes(col_angle);				// matlab -> [Y I] = sort(col_angle)

        const std::vector<size_t> IndexitensWorstCaseCLuster = findItemSortIndex(I,
			[&col_angle](size_t i1) {return col_angle[i1] <= MAXIMUM_ANGLE_DEGREE; });  // matlab -> [NI] = I(Y <= 90)
		
        real_type valueDistance = 0;
		for (auto j : IndexitensWorstCaseCLuster) {							// matlab -> distance_pair_item = M(i,IndexitensWorstCaseCLuster)
			if (valueDistance < distance_matrix(i, j))
				valueDistance = distance_matrix(i, j);						// matlab -> max(distance_pair_item)
		}
		worts_case_cluster(i) = valueDistance;
	}
}

void ECMI::generateOverlapCluster(){

    if(!vector_cluster.empty())
        vector_cluster.clear();

    for (int i = 0; i < x_projected.rows(); i++) {

        real_type distanceAutomatic = radius_distance_cluster(i);
        real_type error_distance_range = worts_case_cluster(i) - distanceAutomatic;
        real_type cut_distance = distanceAutomatic + (error_distance_range > 0 ? (error_distance_range*(1 - correlation_precision)) : 0.0);

        const row_vector_type& distancesPairItem = distance_matrix.row(i);
        const std::vector<size_t>& I = sort_indexes(distancesPairItem);										// matlab -> [Y I] = sort(distance_matrix)

        std::vector<size_t> IndexItemAutomaticCluster = findItemSortIndex(I,[&distancesPairItem, &distanceAutomatic](size_t i1) {return distancesPairItem[i1] <= distanceAutomatic; });
        const std::vector<size_t> IndexItemsCluste = findItemSortIndex(I,[&distancesPairItem, &cut_distance](size_t i1) {return distancesPairItem[i1] <= cut_distance; });   // matlab -> [IndexItemCluster] = I(distancesPairItem < cut_distance)

        //remove reference item cluster
        IndexItemAutomaticCluster.erase(IndexItemAutomaticCluster.begin());

        std::set<int> ItemsCluster(IndexItemsCluste.begin(),IndexItemsCluste.end());
        std::set<int> AutomaticCluster(IndexItemAutomaticCluster.begin(),IndexItemAutomaticCluster.end());

        vector_cluster.push_back(Cluster{ItemsCluster,AutomaticCluster});
    }
}

void ECMI::generateFpTree(BASE_NUM& _dado){

   for(signed i = 0; i < fc.cols(); i++ )
       mapIdSupport[i] = fc(i);

   auto listTransation = _dado.getTransation();

   auto cabecalho = _dado.getCabecalho();

   m_cabecalhoBaseIdReal = _dado.getCabecalhoIdReal();

   for (auto it = listTransation.begin(); it != listTransation.end(); ++it){
        auto idTransation = *it;

        std::vector<int> transation;
        for(auto id: idTransation)
              transation.push_back(cabecalho[id]);

        std::vector<int> newtransation = vector_cluster.getRelevantItens(transation);

        std::sort(newtransation.begin(), newtransation.end(), [this](const int& lhs, const int& rhs) {
            if(mapIdSupport[lhs] > mapIdSupport[rhs]) return true ;
            if(mapIdSupport[lhs] < mapIdSupport[rhs]) return false ;

            if(lhs > rhs) return true;
            if(lhs < rhs) return false;

            return false;
        });

       m_pTree->addTransaction(newtransation);
   }

    m_pTree->createHeaderList(mapIdSupport);
}

void ECMI::runAlgorithm(real_type _correlation_precision){

    correlation_precision = _correlation_precision;

    qtdItemset = 0;

    generateOverlapCluster();

    generateFpTree(base);

    if (!m_pTree->empty()){
        std::set<int> prefixItemset;
        fpGrowth(m_pTree, prefixItemset, mapIdSupport);
    }

}

void ECMI::fpGrowth(std::shared_ptr<FPTree> _tree, std::set<int>& _prefix, const std::map<int, int>& _mapSupport) {

    if (_tree->contains_single_path()) {
        allCombinationsOfPrefixPath(_tree, _prefix, _mapSupport);
    }
    else {
        for (int i = (int)_tree->getHeader().size() -1; i >=0 ; i--) {

            int prefItem = _tree->getHeader().at(i);
            int frequence = _mapSupport.find(prefItem)->second;

            std::set<int> newprefix = _prefix;
            newprefix.insert(prefItem);

            // pruning using clusters
            if (!vector_cluster.existItemset(newprefix))
                continue;

            std::vector<std::vector<std::shared_ptr<const FPNode>>> prefixPaths;
            std::shared_ptr<const FPNode> path = _tree->getFirstNodeTID(prefItem);
            std::map<int, int > mapSupportBeta;

            while (path != nullptr) {
                if (!path->isRootNode()) {
                    std::vector<std::shared_ptr<const FPNode>> prefixPath;
                    prefixPath.push_back(path);

                    int pathCount = path->getFrequency();
                    std::shared_ptr<const FPNode> parent = path->getParent();

                    while (!parent->isRootNode()) {
                        prefixPath.push_back(parent);
                        mapSupportBeta[parent->getItemID()] += pathCount;
                        parent = parent->getParent();
                    }
                    prefixPaths.push_back(prefixPath);
                }
                path = path->getNodeLink();
            }

            std::vector<int> itemsetPrefix;
            itemsetPrefix.insert(std::end(itemsetPrefix), std::begin(newprefix), std::end(newprefix));
            sortedOriginalOrder(itemsetPrefix);

             //if(m_pCFITree->passSubsetChecking(itemsetPrefix,frequence)){

                    std::shared_ptr<FPTree> pTreeBeta = std::make_shared<FPTree>(prefixPaths);

                    if (!pTreeBeta->empty()){
                          pTreeBeta->createHeaderList(mapIdSupport);
                          countRecursao+=1;
                          fpGrowth(pTreeBeta, newprefix, mapSupportBeta);
                    }

                    if(m_pCFITree->passSubsetChecking(itemsetPrefix,frequence) && vector_cluster.obeyFormationRules(set<int>(itemsetPrefix.begin(),itemsetPrefix.end())) )
                        saveItemset(itemsetPrefix,frequence);
           // }
        }
    }
}

struct cmpStruct {
  bool operator() (std::set<int> const & lhs, std::set<int> const & rhs) const
  {
    if( lhs.size() > rhs.size() ) return true;
    if( lhs.size() < rhs.size() ) return false;

    if( lhs > rhs ) return true;
    if( lhs < rhs ) return false;

    return false;
  }
};

void ECMI::allCombinationsOfPrefixPath(std::shared_ptr<FPTree> _tree, std::set<int>& _prefix, const std::map<int, int>& _mapSupport) {

    std::vector<int> idItens;
    std::shared_ptr<const FPNode> node = _tree->getRoot();

    if(node == nullptr) return;

    node = node->getFirstChildren();
    while (node != nullptr) {
        idItens.push_back(node->getItemID());
        node = node->getFirstChildren();
    }
    unsigned int auxRange = 0;
    for (size_t r = 0; r < idItens.size(); r++) {

        if(r == idItens.size() -1){
             std::set<std::set<int>, cmpStruct> allclosedItemsets;

             std::set<int> closeditemset;
             closeditemset.insert(idItens.begin(),idItens.end());
             closeditemset.insert(std::begin(_prefix), std::end(_prefix));

             for(auto it : closeditemset){
                 if(!vector_cluster[it].existItemset(_prefix))
                     continue;

                     std::set<int> closedSubPath = vector_cluster[it].getPossibleItemset(std::set<int>(idItens.begin(),idItens.end()), _prefix);
                     if(!closedSubPath.empty()){
                         bool bAdd = false;
                         for(auto item: closedSubPath){
                             unsigned int rpos = std::find(idItens.begin(), idItens.end(), item) - idItens.begin();
                             if(auxRange <= rpos && rpos < idItens.size()){
                                bAdd  = true;
                                break;
                             }
                         }
                        if(bAdd)
                            allclosedItemsets.insert(closedSubPath);
                     }
             }

             for(auto itAux: allclosedItemsets){

                 std::vector<int> itemset(itAux.begin(),itAux.end());
                 itemset.insert(std::end(itemset), std::begin(_prefix), std::end(_prefix));

                 sortedOriginalOrder(itemset);

                 int pos = 0;
                 for(auto item: itAux){
                     int rpos = std::find(idItens.begin(), idItens.end(), item) - idItens.begin();
                     if(pos < rpos)
                         pos = rpos;
                 }

                 if(m_pCFITree->passSubsetChecking(itemset,_mapSupport.at(idItens[pos])))
                     saveItemset(itemset,_mapSupport.at(idItens[r]));
             }

        }
        else
            if(_mapSupport.at(idItens[r]) != _mapSupport.at(idItens[r+1])){

                std::set<std::set<int>, cmpStruct> allclosedItemsets;

                std::set<int> closeditemset;
                closeditemset.insert(idItens.begin(),idItens.begin()+(r+1));
                closeditemset.insert(std::begin(_prefix), std::end(_prefix));

                //print(closeditemset);

                for(auto it : closeditemset){
                    if(!vector_cluster[it].existItemset(_prefix))
                        continue;

                    std::set<int> closedSubPath = vector_cluster[it].getPossibleItemset(std::set<int>(idItens.begin(),idItens.begin()+(r+1)), _prefix);
                    if(!closedSubPath.empty()){
                        bool bAdd = false;
                        for(auto item: closedSubPath){
                            unsigned int rpos = std::find(idItens.begin(), idItens.end(), item) - idItens.begin();
                            if(auxRange >= rpos && rpos < idItens.size()){
                               bAdd  = true;
                               break;
                            }
                        }
                       if(bAdd)
                             allclosedItemsets.insert(closedSubPath);
                    }
                }

                for(auto itAux: allclosedItemsets){

                    std::vector<int> itemset(itAux.begin(),itAux.end());
                    itemset.insert(std::end(itemset), std::begin(_prefix), std::end(_prefix));

                    sortedOriginalOrder(itemset);

                    if(m_pCFITree->passSubsetChecking(itemset,_mapSupport.at(idItens[r])))
                        saveItemset(itemset,_mapSupport.at(idItens[r]));
                 }

                auxRange = r+1;
            }
       }
}

void ECMI::saveItemset(std::vector<int> _itemset, int _countOcurrence) {

    m_pCFITree->addCFI(_itemset, _countOcurrence);

    float fSupport = (float)_countOcurrence / (float)nTransactions;

    if (_itemset.size() <= 1) return;

    int nMaxOcurrence = 0;
    for (auto item : _itemset)
            mapIdSupport[item] > nMaxOcurrence ? nMaxOcurrence = mapIdSupport[item] : nMaxOcurrence;

    float fHconfidence = fSupport / ((float)nMaxOcurrence / (float)nTransactions);


    std::sort(_itemset.begin(),_itemset.end());

    std::string strItemset;
    for (auto item : _itemset)
        strItemset += std::to_string(m_cabecalhoBaseIdReal[item]) + ",";

    strItemset.pop_back();

    m_fileOutput << std::to_string(fSupport) + " " + std::to_string(fHconfidence) + " " + std::to_string(_itemset.size()) + " " + strItemset << std::endl;

    qtdItemset++;

}

void ECMI::sortedOriginalOrder(std::vector<int>& _itemset) const
{
    std::sort(_itemset.begin(), _itemset.end(), [this](const int& lhs, const int& rhs) {
        if(mapIdSupport.at(lhs) > mapIdSupport.at(rhs)) return true ;
        if(mapIdSupport.at(lhs) < mapIdSupport.at(rhs)) return false ;

        if(lhs > rhs) return true;
        if(lhs < rhs) return false;

        return false;
    });
}

std::vector<size_t> ECMI::getSortByDistance(size_t idItem) const {
	const row_vector_type& distancesPairItem = distance_matrix.row(idItem);
	return sort_indexes(distancesPairItem);
}
