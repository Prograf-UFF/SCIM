#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <utility>


class CFINode {
    friend class CFITree;

public:
    CFINode(int _id);
    ~CFINode();

    std::shared_ptr<const CFINode> getChildId(int _id) const;
    bool isRootNode() const { return itemID == -1; }
    std::shared_ptr<const CFINode> getParent() const { return parent; }
    std::shared_ptr<const CFINode> getNodeLink() const { return node_link; }
    int getFrequency() const { return counter; }
    int getItemID() const { return itemID; }
    std::shared_ptr<const CFINode> getFirstChildren() const { return children.size() == 1 ? children.front() : nullptr; }
    std::string toString(std::string indent) const;

private:
    std::shared_ptr<CFINode> getChildId( int _id);

private:
    int itemID;

    int counter; // frequency counter  (a.k.a. support)

    unsigned int level;  // at which level in the CFI tree this node appears

    std::shared_ptr<CFINode> parent; // the parent node of that node or null if it is the root

    std::shared_ptr<CFINode> node_link; // link to next node with same item id

    std::list<std::shared_ptr<CFINode>> children;

};

class CFITree {
public:
    CFITree();
    ~CFITree();

    void addCFI(const std::vector<int>& _itemset, int counter);

    bool passSubsetChecking(const std::vector<int>& headWithP, int headWithPSupport);

    std::string toString() const;

private:

    void fixNodeLinks(std::shared_ptr<CFINode> _newNode);

    void removelink(std::shared_ptr<CFINode>& _nodechild);

    bool issASubsetOfPrefixPath(const std::vector<int>& headWithP, std::shared_ptr<CFINode> node);


private:
    std::shared_ptr<CFINode> root;

    std::map<int, std::shared_ptr<CFINode>> mapItemNode; // List of pairs (item, frequency) of the header table

    std::map<int, std::shared_ptr<CFINode>> mapItemLastNode; // Map that indicates the last node for each item using the node links. key: item   value: an fp tree node

    std::shared_ptr<CFINode> lastAddedItemsetNode; // last added itemset

};

