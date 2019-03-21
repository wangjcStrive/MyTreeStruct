#include <iostream>
#include <stack>
#include <vector>
#include <queue>
#include <climits>
//#include <algorithm>
using namespace std;


typedef struct BinaryTreeNode
{
	int val;
	struct BinaryTreeNode* left;
	struct BinaryTreeNode* right;
	BinaryTreeNode() :val(0), left(nullptr), right(nullptr) {}
	BinaryTreeNode(int val) :val(val), left(nullptr), right(nullptr) { }
}TreeNode;


class Tree
{
public:
	/*
	 ** 1. traverse遍历
	 */
	void preOrderTraverse(TreeNode* root)
	{
		if (root == nullptr)
		{
			return;
		}

		cout << root->val << " ";
		preOrderTraverse(root->left);
		preOrderTraverse(root->right);
	}



	void midOrderTraverse(TreeNode* root)
	{
		if (root == nullptr)
			cout << "midOrder, empty tree" << endl;

		midOrderTraverse(root->left);
		root->val;
		midOrderTraverse(root->right);
	}

	void postOrderTraverse(TreeNode* root)
	{
		if (root == nullptr)
			cout << "midOrder, empty tree" << endl;

		midOrderTraverse(root->left);
		midOrderTraverse(root->right);
		root->val;
	}


	//非递归遍历
	//先序遍历：先输出节点值 -> 入栈 -> 遍历左子树 -> 退栈时遍历栈顶节点的右子树
	void iterativePreOrderTraverse(TreeNode* root)
	{
		stack<TreeNode*> s;
		TreeNode* tempRoot = root;

		if (root == nullptr)
			return;


		while (tempRoot || !s.empty())
		{
			if (tempRoot)
			{
				cout << tempRoot->val << " ";
				s.push(tempRoot);
				tempRoot = tempRoot->left;

			}
			else
			{
				tempRoot = s.top();
				s.pop();
				tempRoot = tempRoot->right;
			}
		}
	}

	//中序遍历
	//先进栈 -> 遍历左子树 -> 出栈时输出节点 -> 遍历右子树
	void iterativeMidOrderTraverse(TreeNode* root)
	{
		stack<TreeNode*> s;
		TreeNode* tempRoot = root;
		if (root == nullptr)
			return;

		while (tempRoot || !s.empty())
		{
			if (tempRoot)
			{
				s.push(tempRoot);
				tempRoot = tempRoot->left;
			}
			else
			{
				tempRoot = s.top();
				s.pop();
				cout << tempRoot->val << " ";
				tempRoot = tempRoot->right;
			}
		}
	}
	//后序遍历
	//用一个变量标识节点的右子树状态。true为已经遍历，false表示没遍历
	void iterativePostTraverse(TreeNode* root)
	{
		stack< pair<TreeNode*, bool> > s;
		TreeNode* tempRoot = nullptr;
		tempRoot = root;
		while (tempRoot || !s.empty())
		{
			if (tempRoot)
			{
				s.push(make_pair(tempRoot, false));
				tempRoot = tempRoot->left;
			}
			else
			{
				if (s.top().second == false)
				{
					s.top().second = true;
					tempRoot = s.top().first->right;
				}
				else
				{
					cout << s.top().first->val << " ";
					s.pop();
				}
			}
		}
	}

	//分层遍历 广度遍历
	//如果队列不为空：每次弹出一个节点 -> 如果左右不为空则压入队列.
	void LevelTraverse(TreeNode* root)
	{
		if (root == nullptr)
			return;

		queue<TreeNode*> q;
		q.push(root);
		while (!q.empty())
		{
			TreeNode* tempNode = q.front();
			q.pop();
			cout << tempNode->val << " ";

			if (tempNode->left != nullptr)
				q.push(tempNode->left);
			if (tempNode->right != nullptr)
				q.push(tempNode->right);
		}
	}

	vector<vector<int>> levelOrder(TreeNode* root)
	{
		vector<vector<int>> result;
		if (root == nullptr)
			return result;
		queue<TreeNode*> q;
		q.push(root);
		while (!q.empty())
		{
			vector<int> oneLevel;
			int size = q.size();
			for (int i = 0; i < size; i++)
			{
				TreeNode* node = q.front();
				q.pop();
				oneLevel.push_back(node->val);
				if (node->left != nullptr)
					q.push(node->left);
				if (node->right != nullptr)
					q.push(node->right);
			}
			result.push_back(oneLevel);
		}
		return result;
	}



	//get第K层节点个数
	int getNodeQuantityInKLevel(TreeNode* root, int k)
	{
		if (root == nullptr || k < 0)
			return 0;
		if (k = 0)
			return 1;
		int numLeft = getNodeQuantityInKLevel(root->left, k - 1);
		int numRight = getNodeQuantityInKLevel(root->right, k - 1);
		return (numLeft + numRight);
	}


	//LeetCode 96. Unique Binary Search Trees
	//递归
	//以下方法不正确，当n为3,for循环的i为2时计算不正确
	int getCountUniqueBinarySearchTree(int n)
	{
		return getCountUniqueBinarySearchTreeA(1, n);
	}
	int getCountUniqueBinarySearchTreeA(int start, int end)
	{
		int count = 0;
		int nLeft = 0;
		int nRight = 0;
		if (start <= end)
			return 0;
		if (start == end - 1)
			return 1;
		for (int i = start; i <= end; i++)
		{
			nLeft = getCountUniqueBinarySearchTreeA(start, i);
			nRight = getCountUniqueBinarySearchTreeA(i + 1, end);
			count += nLeft + nRight;
		}
		return count;
	}
	//卡特兰数公式
	int getCountUniqueBinarySearchTree_web(int n)
	{
		if (n <= 0)
			return 0;
		std::vector<int> v(n + 1);
		v[0] = v[1] = 1;
		for (int i = 2; i <= n; i++)
			for (int j = 0; j < i; j++)
				v[i] += v[j] * v[i - j - 1];

		return v[n];
	}

	//leetcode 100: same tree
	//需要判断节点内的值
	bool isSameTree(TreeNode* p, TreeNode* q)
	{
		if (p == nullptr && q == nullptr)
			return true;
		else if (p == nullptr || q == nullptr)
			return false;
		else if (p->val != q->val)
			return false;
		bool leftResult = isSameTree(p->left, q->left);
		bool rightResult = isSameTree(p->right, q->right);
		return (leftResult && rightResult);
	}

	//leetcode 101: Symmetric Tree
	//[1,2,2,3,4,4,3]
	//递归
	bool isSymmetric(TreeNode* root)
	{
		if (root == nullptr)
			return false;
		return isSymmetricSubTree(root->left, root->right);
	}
	bool isSymmetricSubTree(TreeNode* left, TreeNode* right)
	{
		if (left == nullptr && right == nullptr)
			return true;
		else if (left == nullptr || right == nullptr)
			return false;
		else if (left->val != right->val)
			return false;

		bool leftResult = isSymmetricSubTree(left->left, right->right);
		bool rightResult = isSymmetricSubTree(left->right, right->left);

		return (leftResult && rightResult);
	}
	//迭代
	bool isSymmetricIterator(TreeNode* root)
	{
		if (!root)
			return true;
		queue<TreeNode*> leftTree, rightTree;
		leftTree.push(root->left);
		rightTree.push(root->right);

		while (!leftTree.empty() && !rightTree.empty())
		{
			TreeNode* node1 = leftTree.front();
			TreeNode* node2 = rightTree.front();
			leftTree.pop();
			rightTree.pop();
			if ((node1 && !node2) || (!node1 && !node2))
				return false;
			if (node1)
			{
				if (node1->val != node2->val)
					return false;
				leftTree.push(node1->left);
				leftTree.push(node1->right);
				rightTree.push(node2->right);
				rightTree.push(node2->left);
			}
		}
	}
	//98. Validate Binary Search Tree
	//需要考虑隔代节点的大小
	//solution 1.很好理解的方案。但是没有accept，感觉像是Lcode的例子有问题似的....
	bool isValidBST(TreeNode* root)
	{
		if (root == NULL)
			return true;

		return dfs(root, INT_MIN, INT_MAX);
	}
	bool dfs(TreeNode* root, int low, int up)
	{
		if (root == NULL)
			return true;

		if (root->val >= up || root->val <= low)
			return false;
		return dfs(root->left, low, root->val) && dfs(root->right, root->val, up);
	}

	//solution 2. from LCode. 递归
	bool isValidBST_recursive(TreeNode* root)
	{

	}

	//104. Maximum Depth of Binary Tree
	int maxDepth(TreeNode* root)
	{
		int depth = 0;
		if (root == NULL)
			return depth;
		return 1 + max( maxDepth(root->left), maxDepth(root->right) );
	}



	//105. Construct Binary Tree from Preorder and Inorder Traversal
	//recursive
	//前序遍历第一个节点是树的根节点，在中序遍历中根据该几点区分出左子树和右子树
	TreeNode* buildTree(vector<int>& preorder, vector<int>& inorder)
	{
		reBuildTree(preorder, 0, preorder.size(), inorder, 0, inorder.size());
	}
	//
	TreeNode* reBuildTree(vector<int>& preorder, int i, int j, vector<int>& inorder, int ii, int jj)
	{
		if (i >= j || ii >= j)
			return nullptr;

		int mid = preorder[i];
		auto f = find(inorder.begin() + ii, inorder.begin() + jj, mid);

		int dis = f - inorder.begin() - ii;

		TreeNode* root = new TreeNode(mid);
		root->left = reBuildTree(preorder, i + 1, i + 1 + dis, inorder, ii, ii + dis);
		root->right = reBuildTree(preorder, i + 1 + dis, j, inorder, ii + dis + 1, jj);
		return root;
	}


	~Tree()
	{
		//TODO. delete all node and set m_pRoot to nullptr
	}

	TreeNode* createTreeRecursion(const vector<int> inputArray)
	{
		return createTree(inputArray, 0);
	}

	// LeetCode 617. Merge Two Binary Trees
	TreeNode* mergeTrees(TreeNode* t1, TreeNode* t2)
	{
		TreeNode* tempNewNode = nullptr;

		if (t1 || t2)
		{
			int leftValue = t1 == nullptr ? 0 : t1->val;
			int rightValue = t2 == nullptr ? 0 : t2->val;
			tempNewNode = new TreeNode(leftValue + rightValue);
			tempNewNode->left = mergeTrees(t1 == nullptr ? nullptr : t1->left, t2 == nullptr ? nullptr : t2->left);
			tempNewNode->right = mergeTrees(t1 == nullptr ? nullptr : t1->right, t2 == nullptr ? nullptr : t2->right);
		}
		return tempNewNode;
	}
	//Iterative Mode. not correct， 如果t1树上节点为空，没有处理好.
	TreeNode* mergeTreesIterative(TreeNode* t1, TreeNode* t2)
	{
		if (t1 == nullptr && t2 == nullptr)
			return nullptr;

		stack< pair<TreeNode*, TreeNode*> > s;
		s.push(make_pair(t1, t2));
		while (!s.empty())
		{
			//pair<TreeNode*, TreeNode*> tempPair = s.top();
			auto tempPair = s.top();
			s.pop();
			auto p1 = tempPair.first;
			auto p2 = tempPair.second;

			if (p1 != nullptr)
				p1->val += p2 ? p2->val : 0;
			else
			{
				if (p2 != nullptr)
					p1 = new TreeNode(p2->val);
			}

			auto leftSubNode1 = p1 ? p1->left : nullptr;
			auto leftSubNode2 = p2 ? p2->left : nullptr;
			auto rightSubNode1 = p1 ? p1->right : nullptr;
			auto rightSubNode2 = p2 ? p2->right : nullptr;

			if(leftSubNode1 || leftSubNode2)
				s.push(make_pair(leftSubNode1, leftSubNode2));
			if (rightSubNode1 || rightSubNode2)
				s.push(make_pair(rightSubNode1, rightSubNode2));
		}
		return t1;
	}

private:
	/*
	 * 0
	 * 1 2
	 * 3 4 5 6
	 * left:2i+1
	 * right: 2i+2
	 */
	TreeNode* createTree(const vector<int> root, size_t index)
	{
		size_t leftIndex, rightIndex;
		try
		{
			if (index >= root.size())
				return nullptr;
			TreeNode* rootNode = new TreeNode(root[index]);

			leftIndex = index * 2 + 1;
			if (leftIndex < root.size() && root[leftIndex] != 0)
			{
				rootNode->left = createTree(root, leftIndex);
			}
			else
				rootNode->left = nullptr;

			rightIndex = index * 2 + 2;
			if (rightIndex < root.size() && root[rightIndex] != 0)
				rootNode->right = createTree(root, rightIndex);
			else
				rootNode->right = nullptr;

			return rootNode;
		}
		catch (const std::exception& e)
		{
			cout << e.what() << endl;
		}
	}
};
