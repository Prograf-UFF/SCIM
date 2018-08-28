/*
Copyright 2018 Altobelli de Brito Mantuan.
Distributed under the GNU General Public License, Version 3.
Last update : August 25th, 2018.
*/

#define __USE_MINGW_ANSI_STDIO 0
#include "Cluster.h"
#include <algorithm>
#include <iostream>
#include <fstream>

Cluster::Cluster(const std::set<int>& _itemsCluster, const std::set<int>& automatic_region)
{
    items = _itemsCluster;
    itemsAutomaticRegion = automatic_region;
}

Cluster::~Cluster()
{
}

bool Cluster::existItemset(const std::set<int>& _itemset) const {

    unsigned int numElement = 0;
    for (int i : _itemset)
		if (items.find(i) != items.end()) ++numElement;

	return numElement == _itemset.size();
}

bool Cluster::obeyFormationRules(const std::set<int>& _itemset) const {

	bool hasItemAutomaticRegion = false;

	if (itemsAutomaticRegion.size() == 0) {
		hasItemAutomaticRegion = true;
	}
	else {
        for (int i : _itemset)
			if (itemsAutomaticRegion.find(i) != itemsAutomaticRegion.end() ) {
				hasItemAutomaticRegion = true;
				break;
			}
	}

	return hasItemAutomaticRegion && existItemset(_itemset);
}

std::set<int> Cluster::getPossibleItemset(const std::set<int>& _cloasePath, const std::set<int>& prefix){

    std::set<int> v_intersection;

    for (int i : _cloasePath)
        if (items.find(i) != items.end() )
            v_intersection.insert(i);

    std::set<int> itemset = v_intersection;
    itemset.insert(prefix.begin(),prefix.end());

    if(obeyFormationRules(itemset))
        return v_intersection;

    return std::set<int>();
}

std::set<int> Cluster::getIntersection(const std::set<int>& _itemset) const{

    std::vector<int> v(_itemset.size());

    auto it = std::set_intersection(_itemset.begin(), _itemset.end(), items.begin(), items.end(), v.begin());

    v.resize(it-v.begin());

    std::set<int> subTransation(v.begin(),v.end());

    if(subTransation.size() > 1 && obeyFormationRules(subTransation))
        return subTransation;

    return  std::set<int>();

}

bool Clusters::existItemset(const std::set<int>& _itemset) const{

	for (const Cluster cluster : *this) {
        if (cluster.existItemset(_itemset)) return true;
	}
	return false;
}

bool Clusters::obeyFormationRules(const std::set<int>& _itemset) const{

    for (int item: _itemset)
        if (this->at(item).obeyFormationRules(_itemset))
			return true;

	return false;
}

std::vector<int> Clusters::getRelevantItens(const std::vector<int>& _transation) const{

    std::set<int> setTransation(_transation.begin(),_transation.end());
    std::set<int> newTransation;
    for (int item: _transation) {
       auto subTransation = this->at(item).getIntersection(setTransation);
       newTransation.insert(subTransation.begin(), subTransation.end());
    }

    return std::vector<int>(newTransation.begin(),newTransation.end());
}

void Clusters::saveFile(const std::string & _file) {
    std::ofstream file;
    file.open(_file);
        std::string lineCluster;
        for( int i = 0; i < this->size(); i++){
            lineCluster = std::to_string(i);
            lineCluster += "-";
            for (int item: this->at(i).getItemsAutomaticRegion()){
                lineCluster += std::to_string(item) + ",";
            }
            if(this->at(i).getItemsAutomaticRegion().size() != 0)
                lineCluster.pop_back();

            lineCluster += "-";
            for (int item: this->at(i).getItemsPrecisionCorrelation(i)){
                lineCluster += std::to_string(item) + ",";
            }
            lineCluster.pop_back();
        file << lineCluster << std::endl;
        }
    file.close();
}
