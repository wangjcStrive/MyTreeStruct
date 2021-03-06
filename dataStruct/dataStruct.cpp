// dataStruct.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

//#include "pch.h"
#include <iostream>
#include "tree.h"

using namespace std;

int main()
{
	Tree ins;
	TreeNode* root;
	TreeNode* anotherRoot;
	TreeNode* newTree;
	//vector<int> nodeArray = { 0, 1, 2, 3, 4, 5 };
	vector<int> nodeArray;
	for (int i = 1; i < 9; i++)
	{
		nodeArray.push_back(i);
	}
	//nodeArray[3] = 0;
	root = ins.createTreeRecursion(nodeArray);


	vector<int> anotherNodeArray;
	for (int i = 6; i < 15; i++)
	{
		anotherNodeArray.push_back(i);
	}
	anotherNodeArray[3] = 0;
	anotherNodeArray[5] = 0;
	anotherRoot = ins.createTreeRecursion(anotherNodeArray);

	cout << endl;
	cout << ins.maxDepth(root) << endl;
	cout << ins.maxDepth(anotherRoot) << endl;
	cout << ins.isValidBST(root) << endl;
	cout << ins.isValidBST(anotherRoot) << endl;
	////newTree = ins.mergeTrees(root, anotherRoot);
	//newTree = ins.mergeTreesIterative(root, anotherRoot);
	//ins.iterativePreOrderTraverse(newTree);
	

	//ins.iterativePreOrderTraverse(root);
	//cout << endl;
	//ins.iterativeMidOrderTraverse(root);
	//cout << endl;
	//ins.iterativePostTraverse(root);
	//cout << endl;
	//cout << "LevelTraverse: ";
	//ins.LevelTraverse(root);
	//cout << endl;

	//cout << ins.getCountUniqueBinarySearchTree(3) << endl;
}