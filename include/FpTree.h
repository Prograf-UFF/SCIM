#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <utility>


using Item = std::string;
using Transaction = std::vector<Item>;
using TransformedPrefixPath = std::pair<std::vector<Item>, int>;
using Pattern = std::pair<std::set<Item>, int>;

class FPTree;

class FPNode {
	friend class FPTree;

public:
    FPNode(int _id);
	~FPNode();

    std::shared_ptr<const FPNode> getChildId(int _id) const;
    bool isRootNode() const { return itemID == -1; }
	std::shared_ptr<const FPNode> getParent() const { return parent; }
	std::shared_ptr<const FPNode> getNodeLink() const { return node_link; }
    int getFrequency() const { return frequency; }
    int getItemID() const { return itemID; }
	std::shared_ptr<const FPNode> getFirstChildren() const { return children.size() == 1 ? children.front() : nullptr; }
	std::string toString(std::string indent) const;

private:
    std::shared_ptr<FPNode> getChildId( int _id);

private:
    int itemID;

    int frequency;

	std::shared_ptr<FPNode> parent;

	std::shared_ptr<FPNode> node_link; // link to next node with same item id

	std::list<std::shared_ptr<FPNode>> children;

};

class FPTree {
public:
    FPTree();
	~FPTree();

    FPTree(std::vector<std::vector<std::shared_ptr<const FPNode>>>& _paths);

    void addTransaction(std::vector<int> _trans);

    bool empty() const;

	bool contains_single_path() const ;

    std::vector<int> getHeader() const { return header; }

    std::shared_ptr<const FPNode> getFirstNodeTID(int _id) const;

    std::shared_ptr<const FPNode> getLastNodeTID(int _id) const;

	std::shared_ptr<const FPNode> getRoot() const { return root; }

    void createHeaderList(std::map<int, int> _mapSupport);

	std::string toString() const;

private:
	void fixNodeLinks(std::shared_ptr<FPNode> _newNode);

	bool contains_single_path(const std::shared_ptr<FPNode>& fpnode) const;

	void addPrefixPath(std::vector< std::shared_ptr<const FPNode>> _path);

	void removelink(std::shared_ptr<FPNode>& _nodechild);


private:
	std::shared_ptr<FPNode> root;

    std::vector<int> header;

    std::map<int, std::shared_ptr<FPNode>> mapItemNode;

    std::map<int, std::shared_ptr<FPNode>> mapItemLastNode;

};

