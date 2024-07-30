#include <iostream>
#include <fstream>
using namespace std;

enum Boolean {FALSE, TRUE};
template <class Type>
class Element {
public:
    Type key;
};
template <class Type>
class TBST;
template <class Type>
class TBSTNode {
friend class TBST<Type>;
public:
private:
   Element<Type> data;
   TBSTNode *LeftChild, *RightChild;
   void display(int i);
   void treeprint(int i);
};

template <class Type>
class TBST {
public:
   TBST(TBSTNode<Type> *init = 0) {root = init;};

//   TBST& operator=(const TBST&);

   Boolean Insert(const Element<Type>& x);
//   Boolean Delete(const Element<Type>&);
   TBSTNode<Type>* Search(TBSTNode<Type>*, const Element<Type>&);
   TBSTNode<Type>* Search(const Element<Type>&);
   TBSTNode<Type>* IterSearch(const Element<Type>&);

   void treeprint() { cout << "\n"; root->treeprint(1); }

   void display() {cout << "\n";
		   if (root) root->display(1);
		   else cout << "0\n";};
   Element<Type>* Split(Type i, TBST& B, TBST& C, Element<Type> &x);
private:
   TBSTNode<Type> *root;
};

template <class Type>
void TBSTNode<Type>::display(int i)
{
   cout << "Position " << i << ": data.key " << data.key << "\n";
   if (LeftChild) LeftChild->display(2*i);
   if (RightChild) RightChild->display(2*i + 1);
};

template <class Type>
Element<Type>* TBST<Type>::Split(Type i, TBST<Type>& B, TBST<Type>& C,
Element<Type> &x)
// Split the binary search tree with respect to key @i@
{
    if (!root) {B.root = C.root = 0; return 0;} // empty tree
    TBSTNode<Type> *Y = new TBSTNode<Type>; TBSTNode<Type> *L = Y;
    TBSTNode<Type> *Z = new TBSTNode<Type>; TBSTNode<Type> *R = Z;
    TBSTNode<Type> *t = root;
    while (t)
	if (i == t->data.key) {  // split at @t@
	    L->RightChild = t->LeftChild;
	    R->LeftChild = t->RightChild;
	    x = t->data;
	    B.root = Y->RightChild; delete Y;
	    C.root = Z->LeftChild; delete Z;
	    return &x;
	}
	else if (i < t->data.key) {
	    R->LeftChild = t;
	    R = t; t = t->LeftChild;
	}
	else {
	    L->RightChild = t;
	    L = t; t = t->RightChild;
	}
    L->RightChild = R->LeftChild = 0;
    B.root = Y->RightChild; delete Y;
    C.root = Z->LeftChild; delete Z;
    return 0;
}

template <class Type>
void TBSTNode<Type>::treeprint(int l)
{
   if (this) {
      (this->RightChild)->treeprint(l+1);
      for (int i = 1; i <= l; i++) cout << "   ";
      cout << this->data.key << "\n";
      (this->LeftChild)->treeprint(l+1);
   };
};

template <class Type>
TBSTNode<Type>* TBST<Type>::Search(TBSTNode<Type>* b, const Element<Type> &x)
{
   if (!b) return 0;
   if (x.key == b->data.key) return b;
   if  (x.key < b->data.key) return Search(b->LeftChild,x);
   return Search(b->RightChild,x);
};

template <class Type>
TBSTNode<Type>* TBST<Type>::Search(const Element<Type>& x)
{
   return Search(root, x);
}

template <class Type>
TBSTNode<Type>* TBST<Type>::IterSearch(const Element<Type>& x)
{
   for (TBSTNode<Type> *t = root;  t; )
   {
      if (x.key == t->data.key) return t;
      if (x.key < t->data.key) t = t->LeftChild;
      else t = t->RightChild;
   }
   return 0;
}

template <class Type>
Boolean TBST<Type>::Insert(const Element<Type>& x)
{
   TBSTNode<Type> *p = root;  TBSTNode<Type> *q = 0;
   while (p) {
      q = p;
      if (x.key == p->data.key) return FALSE; //x.key is already in t
      if (x.key < p->data.key) p = p->LeftChild;
      else p = p->RightChild;
   };

   p = new TBSTNode<Type>;
   p->LeftChild = p->RightChild = 0; p->data = x;
   if (!root) root = p;
   else if (x.key < q->data.key) q->LeftChild = p;
	else q->RightChild = p;
   return TRUE;
}
template <class Type>
class BST
{
	TBST<Type> internal_bst;
public:
	void Insert(int data)
	{
		//return;
		Element<Type> element;
		element.key = data;
		this->internal_bst.Insert(element);
	}
	int IsExist(Type data)
	{
		//return 0;
		Element<Type> element;
		element.key = data;
		if(this->internal_bst.Search(element)) return 1;
		return 0;
	}
};