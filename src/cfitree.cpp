/*
Copyright 2018 Altobelli de Brito Mantuan.
Distributed under the GNU General Public License, Version 3.
Last update : August 25th, 2018.
*/

#include "CFITree.h"

CFINode::CFINode(int _id) :
    itemID( _id ), counter( 1 ), level(0) ,parent( nullptr ) ,node_link( nullptr )
{
}

CFINode::~CFINode() {
    children.clear();
    node_link = nullptr;
    parent = nullptr;
}

std::shared_ptr<const CFINode> CFINode::getChildId(int _id) const {
    for (std::shared_ptr<CFINode> child : children)
        if (child->itemID == _id)
            return child;

    return nullptr;
}

std::shared_ptr<CFINode> CFINode::getChildId(int _id) {
    for (std::shared_ptr<CFINode> child : children)
        if (child->itemID == _id)
            return child;

    return nullptr;
}

std::string CFINode::toString(std::string indent) const {
    std::string output;
    output += " " + std::to_string(itemID);
    output += " (count=" + std::to_string(counter) + "level=" + std::to_string(level) + ")\n";
    std::string newIndent = indent + "  ";
    for (std::shared_ptr<const CFINode> child : children)
        output += newIndent + child->toString(newIndent);

    return output;
}

CFITree::CFITree():
 root(std::make_shared<CFINode>((int)-1)), lastAddedItemsetNode(nullptr) {
}


CFITree::~CFITree(){
    lastAddedItemsetNode = nullptr;
    mapItemNode.clear();
    mapItemLastNode.clear();
    removelink(root);
}


void CFITree::removelink(std::shared_ptr<CFINode>& _nodechild){
    auto it = _nodechild->children.begin();
        while (it != _nodechild->children.end()) {
            removelink(*it);
            it++;
        }
        _nodechild = nullptr;
}


void CFITree::addCFI(const std::vector<int> &_itemset, int counter)
{
    std::shared_ptr<CFINode> currentNode = root;
        for (size_t i = 0; i < _itemset.size(); i++) {
            std::shared_ptr<CFINode> child = currentNode->getChildId(_itemset[i]);
            if (child) {
                if(child->counter < counter)
                    child->counter = counter;

                currentNode = child;
            }
            else {

                std::shared_ptr<CFINode> newNode = std::make_shared<CFINode>(_itemset[i]);
                newNode->parent = currentNode;
                newNode->level = i+1;
                newNode->counter = counter;
                currentNode->children.push_back(newNode);

                currentNode = newNode;
                fixNodeLinks(newNode);
            }
        }
      lastAddedItemsetNode = currentNode;
}

void CFITree::fixNodeLinks(std::shared_ptr<CFINode> _newNode) {
    auto ItLastNode = mapItemLastNode.find(_newNode->itemID);
    if (ItLastNode != mapItemLastNode.end()) {
        auto node = ItLastNode->second;
        node->node_link = _newNode;
    }

    mapItemLastNode[_newNode->itemID] = _newNode;

    auto ItNode = mapItemNode.find(_newNode->itemID);
    if (ItNode == mapItemNode.end()) {
        mapItemNode[_newNode->itemID] = _newNode;
    }

}

bool CFITree::passSubsetChecking(const std::vector<int>& headWithP, int headWithPSupport){

    if(lastAddedItemsetNode != nullptr && lastAddedItemsetNode->counter == headWithPSupport){
        if(issASubsetOfPrefixPath(headWithP,lastAddedItemsetNode))
            return false;
    }

    int firstItem = headWithP[headWithP.size()-1];

    auto Itnode = mapItemNode.find(firstItem);
    if(Itnode == mapItemNode.end())
        return true;

    auto node = Itnode->second;

    do {
        // for a node, we will check if "headwithP" is a subset of the path ending at node
        // if it is a subset, then "headWithP" is in the CFI-tree, we return false
       if(node->counter == headWithPSupport && issASubsetOfPrefixPath(headWithP, node))
            return false;

        // go to the next itemset to test
        node = node->node_link;
       }while(node != nullptr);

    // the itemset is not in the CFI-TREE.  Itemset passed the test!
    return true;

}

bool CFITree::issASubsetOfPrefixPath(const std::vector<int>& headWithP, std::shared_ptr<CFINode> node){

        // optimization proposed in the fpmax* paper: if there is less than itemset node in that branch, we don't need to check it
        if(node->level >= headWithP.size()) {
            // check if "itemset" is contained in the prefix path ending at "node"
            // We will start comparing from the parent of "node" in the prefix path since
            // the last item of itemset is "node".
            std::shared_ptr<CFINode> nodeToCheck = node;
            int positionInItemset = headWithP.size()-1;
            int itemToLookFor = headWithP[positionInItemset];
            // for each item in itemset
            do {
                if(nodeToCheck->itemID == itemToLookFor) {
                    positionInItemset--;
                    // we found the itemset completely, so the subset check test is failed
                    if(positionInItemset < 0)
                        return true;

                    itemToLookFor = headWithP[positionInItemset];
                }
                nodeToCheck = nodeToCheck->parent;
            }while(nodeToCheck != nullptr);
        }
        return false;
}

std::string CFITree::toString() const {
    std::string temp = "FCI\n";

    temp += root->toString("");

    return temp;
}

