/*
Copyright 2018 Altobelli de Brito Mantuan.
Distributed under the GNU General Public License, Version 3.
Last update : August 25th, 2018.
*/

#define __USE_MINGW_ANSI_STDIO 0

#include "FpTree.h"
#include <algorithm>
#include <cassert>
#include <utility>
#include <iterator>
#include <string>

FPNode::FPNode(int _id) :
	itemID( _id ), frequency( 1 ),  parent( nullptr ) ,node_link( nullptr )
{
}

FPNode::~FPNode() {
	children.clear();
	node_link = nullptr;
	parent = nullptr;
}

std::shared_ptr<const FPNode> FPNode::getChildId(int _id) const {
	for (std::shared_ptr<FPNode> child : children)
		if (child->itemID == _id)
			return child;

	return nullptr;
}

std::shared_ptr<FPNode> FPNode::getChildId(int _id) {
	for (std::shared_ptr<FPNode> child : children)
		if (child->itemID == _id)
			return child;

	return nullptr;
}

std::string FPNode::toString(std::string indent) const {
	std::string output;
    output += " " + std::to_string(itemID + 1);
	output += " (count=" + std::to_string(frequency) + ")\n";
	std::string newIndent = indent + "  ";
	for (std::shared_ptr<const FPNode> child : children)
		output += newIndent + child->toString(newIndent);

	return output;
}

FPTree::FPTree() :
    root( std::make_shared<FPNode>((int)-1) ) {

}

FPTree::FPTree(std::vector<std::vector<std::shared_ptr<const FPNode>>>& _paths) :
    root(std::make_shared<FPNode>((int)-1)) {

	for (auto path : _paths)
        addPrefixPath(path);

}

void FPTree::removelink(std::shared_ptr<FPNode>& _nodechild) {
	auto it = _nodechild->children.begin();
	while (it != _nodechild->children.end()) {
		removelink(*it);
		it++;
	}
	_nodechild = nullptr;
}

FPTree::~FPTree() {
	mapItemNode.clear();
	mapItemLastNode.clear();
	removelink(root);
}

void FPTree::fixNodeLinks(std::shared_ptr<FPNode> _newNode) {
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

void FPTree::addTransaction(std::vector<int> _trans) {

	std::shared_ptr<FPNode> currentNode = root;
    for (int item : _trans) {
		std::shared_ptr<FPNode> child = currentNode->getChildId(item);
		if (child) {
			child->frequency++;
			currentNode = child;
		}
		else {
			std::shared_ptr<FPNode> newNode = std::make_shared<FPNode>(item);
			newNode->parent = currentNode;
			currentNode->children.push_back(newNode);

			currentNode = newNode;
			fixNodeLinks(newNode);
		}
	}
}

void FPTree::addPrefixPath(std::vector< std::shared_ptr<const FPNode>> _path) {
    int pathFrequence = _path[0]->frequency;
	std::shared_ptr<FPNode> currentNode = root;

	for (size_t i = _path.size() - 1; i >= 1; i--) {
		std::shared_ptr<const FPNode> pathItem = _path[i];

		std::shared_ptr<FPNode> child = currentNode->getChildId(pathItem->itemID);
		if (child) {
			child->frequency += pathFrequence;
			currentNode = child;
		}
		else {
			std::shared_ptr<FPNode> newNode = std::make_shared<FPNode>(pathItem->itemID);
			newNode->parent = currentNode;
			newNode->frequency = pathFrequence;
			currentNode->children.push_back(newNode);

			currentNode = newNode;
			fixNodeLinks(newNode);
		}
	}
}

void FPTree::createHeaderList( std::map<int, int> _mapSupport) {

    std::transform(mapItemNode.begin(), mapItemNode.end(), back_inserter(header), [](std::pair<int, std::shared_ptr<FPNode>> p) { return p.first; });

    std::sort(header.begin(), header.end(), [&_mapSupport](const int& lhs, const int& rhs)
	{
        if(_mapSupport[lhs] > _mapSupport[rhs]) return true ;
        if(_mapSupport[lhs] < _mapSupport[rhs]) return false ;

        if(lhs > rhs) return true;
        if(lhs < rhs) return false;

        return false;

	});    
}

bool FPTree::empty() const {
    assert( root );
    return root->children.size() == 0;
}

bool FPTree::contains_single_path(const std::shared_ptr<FPNode>& fpnode) const{
	assert(fpnode);
	if (fpnode->children.size() == 0) { return true; }
	if (fpnode->children.size() > 1) { return false; }
	return contains_single_path(fpnode->children.front());
}

bool FPTree::contains_single_path() const{
	if (empty()) return true;
	return contains_single_path(root);
}

std::shared_ptr<const FPNode> FPTree::getFirstNodeTID(int _id) const {
	auto it = mapItemNode.find(_id);
	if (it != mapItemNode.end())
		return it->second;

	return nullptr;
}

std::shared_ptr<const FPNode> FPTree::getLastNodeTID(int _id) const {
	auto it = mapItemLastNode.find(_id);
	if (it != mapItemLastNode.end())
		return it->second;

	return nullptr;
}

std::string FPTree::toString() const {
	std::string temp = "F\n";
    temp += "HeaderList: [ " ;

    for (int item : header)
        temp += std::to_string(item + 1) + " ,";

	temp.pop_back();
	temp += " ] \n";

	temp += root->toString("");

	return temp;
}
