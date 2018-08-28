#pragma once

#include <set>
#include <vector>
#include <string>

/*
                    {automaticRegion} {Correlation_precision}
    Itens cluster = |a|b|c|d|e|f|g|h  |i|j|
*/

class Cluster
{
public:
    Cluster(const std::set<int>& _items, const std::set<int>& automatic_region);
	~Cluster();

    bool existItemset(const std::set<int>& _itemset) const;
    bool obeyFormationRules(const std::set<int>& _itemset) const;
    std::set<int> getPossibleItemset(const std::set<int>& _cloasePath, const std::set<int>& prefix);
    std::set<int> getIntersection(const std::set<int>& _itemset) const;

    const std::set<int> getItemsPrecisionCorrelation(int _ref) const {
        std::set<int> precisionCorr;

        for (int i : items)
            if( (i != _ref) && (itemsAutomaticRegion.find(i) == itemsAutomaticRegion.end()) )
                precisionCorr.insert(i);

        return precisionCorr;
    }

    const std::set<int> getItemsAutomaticRegion() const {return itemsAutomaticRegion; }

private:
    std::set<int> items;
    std::set<int> itemsAutomaticRegion;
};

class Clusters : public std::vector<Cluster> {

public:
    bool existItemset(const std::set<int>& _itemset) const;
    bool obeyFormationRules(const std::set<int>& _itemset) const ;
    std::vector<int> getRelevantItens(const std::vector<int>& _transation) const;
    void saveFile(const std::string & _file) ;
};
