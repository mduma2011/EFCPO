#ifndef _HEURISTIC_H
#define _HEURISTIC_H

#include "Graph.h"
#include "Tour.h"
#include <assert.h>
#include "Time_Random.h"
#include <iostream>
#include <fstream>
#include "mex.h"
using namespace std;

class Heuristics
{
	struct Node;
	struct Candidate;
	struct Segment;
	struct SSegment;
	struct SwapRecord;
	struct Node 
	{
		int Id = -1;
		int Subproblem;
		Node *Pred, *Suc;
		Node *BestPred, *BestSuc;
		Candidate *CandidateSet;		
		double X, Y, Z;
		double Xc, Yc, Zc;
		char Axis;
		//int BestPi;
		Node *Nearest;
		Node *Next;
		int LastV;
		int V;
		int Degree;
		Node *Tail;
		int Cost;
		Node *OldPred, *OldSuc;
		Node *FixedTo1, *FixedTo2;
		int Level;
		int mark;
		int Loc;
		int Rank;
		char OldPredExcluded, OldSucExcluded;  
		Segment *Parent;
		short IsActivated;
		int PredCost, SucCost;
		Node *Dad;
		Node *InputSuc;
	};
	struct Candidate
	{
		Node *To;
		int Cost = 0;
		int Alpha = 0;
	};
	struct Segment {
		char Reversed; 
		Node *First, *Last; 
		Segment *Pred, *Suc;
		int Rank; 
		int Size; 
		SSegment *Parent;   
	};
	struct SSegment {
		char Reversed;  
		Segment *First, *Last; 
		SSegment *Pred, *Suc; 
		int Rank;   
		int Size;   
	};
	struct SwapRecord
	{
		Node *t1, *t2, *t3, *t4;
	};

	class Heap
	{
		int Size;
		Node **Nodes;
		int HeapCount;    /* Its current number of elements */
		int HeapCapacity; /* Its capacity */
	public:		
		Heap(int Size)
		{
			this->Size = Size;
			Nodes = (Node **) malloc((Size + 1) * sizeof(Node *));
			HeapCapacity = Size;
			HeapCount = 0;
		}
		~Heap()
		{
			delete Nodes;			
		}
		void HeapSiftUp(Node * N)
		{
			int Loc = N->Loc, Parent = Loc / 2;

			while (Parent && N->Rank < Nodes[Parent]->Rank) {
				Nodes[Loc] = Nodes[Parent];
				Nodes[Loc]->Loc = Loc;
				Loc = Parent;
				Parent /= 2;
			}
			Nodes[Loc] = N;
			N->Loc = Loc;
		}
		void HeapSiftDown(Node * N)
		{
			int Loc = N->Loc, Child;

			while (Loc <= HeapCount / 2) {
				Child = 2 * Loc;
				if (Child < HeapCount && Nodes[Child + 1]->Rank < Nodes[Child]->Rank)
					Child++;
				if (N->Rank <= Nodes[Child]->Rank)
					break;
				Nodes[Loc] = Nodes[Child];
				Nodes[Loc]->Loc = Loc;
				Loc = Child;
			}
			Nodes[Loc] = N;
			N->Loc = Loc;
		}
		Node *HeapDeleteMin()
		{
			Node *Remove;

			if (!HeapCount)
				return 0;
			Remove = Nodes[1];
			Nodes[1] = Nodes[HeapCount--];
			Nodes[1]->Loc = 1;
			HeapSiftDown(Nodes[1]);
			Remove->Loc = 0;
			return Remove;
		}
		void HeapInsert(Node * N)
		{
			HeapLazyInsert(N);
			HeapSiftUp(N);
		}
		void HeapDelete(Node * N)
		{
			int Loc = N->Loc;
			if (!Loc)
				return;
			Nodes[Loc] = Nodes[HeapCount--];
			Nodes[Loc]->Loc = Loc;
			if (Nodes[Loc]->Rank > N->Rank)
				HeapSiftDown(Nodes[Loc]);
			else
				HeapSiftUp(Nodes[Loc]);
			N->Loc = 0;
		}
		void HeapLazyInsert(Node * N)
		{   
			if(HeapCount < HeapCapacity)
            {    
			  Nodes[++HeapCount] = N;
			  N->Loc = HeapCount;
            }
		}
		void Heapify()
		{
			int Loc;
			for (Loc = HeapCount / 2; Loc >= 1; Loc--)
				HeapSiftDown(Nodes[Loc]);
		}
	};

	Segment *FirstSegment;
	Segment *SegmentSet;
	SSegment *FirstSSegment;
	Node** NodeSet;
	Node* FirstNode;
	int Dimension;
	int MaxSwaps;
	Graph* graph;
	SwapRecord *SwapStack;
	Node **KDTree;
	Candidate *CandidateSet;
	double *XMin, *XMax, *YMin, *YMax, *ZMin, *ZMax;
	int cutoff;
	int Swaps;
	int cand_num;
	int GroupSize;
	int Groups;
	int SGroupSize;
	int SGroups;
	short Reversed;
	short OldReversed;
	enum CoordTypes { TWOD_COORDS, THREED_COORDS, NO_COORDS };
#define Link(a, b) { (a)->Suc = (b); (b)->Pred = (a); }
#define Coord(N, axis) (axis == 0 ? (N)->X : axis == 1 ? (N)->Y : (N)->Z)
#define Free(s) { free(s); s = 0; }
#define Fixed(a, b) ((a)->FixedTo1 == (b) || (a)->FixedTo2 == (b))
#define FixedOrCommon(a, b) Fixed(a, b)
#define Swap1(a1,a2,a3)\
	FLIP(a1,a2,a3,0)
#define Swap2(a1,a2,a3, b1,b2,b3)\
	(Swap1(a1,a2,a3), Swap1(b1,b2,b3))
#define Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3)\
	(Swap2(a1,a2,a3, b1,b2,b3), Swap1(c1,c2,c3))
#define Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3)\
	(Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3), Swap1(d1,d2,d3))
#define Swap5(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3, e1,e2,e3)\
	(Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3), Swap1(e1,e2,e3))

#define SPLIT_CUTOFF 0.75
#define InBestTour(a, b) ((a)->BestSuc == (b) || (b)->BestSuc == (a))
#define InInputTour(a, b) ((a)->InputSuc == (b) || (b)->InputSuc == (a))
#define InInitialTour(a, b)\
	((a)->InitialSuc == (b) || (b)->InitialSuc == (a))
#define Near(a, b)\
	((a)->BestSuc ? InBestTour(a, b) : (a)->Dad == (b) || (b)->Dad == (a))
#define InInputTour(a, b) ((a)->InputSuc == (b) || (b)->InputSuc == (a))
#define Pred(a) (Reversed == (a)->Parent->Reversed ? (a)->Pred : (a)->Suc)
#define Suc(a) (Reversed == (a)->Parent->Reversed ? (a)->Suc : (a)->Pred)
#define BETWEEN(a, b, c) Between_SL(a, b, c)
#define FLIP(a, b, c, d) Flip_SL(a, b, c)
#define Link(a, b) { (a)->Suc = (b); (b)->Pred = (a); }
#define Follow(b, a)\
	{ Link((b)->Pred, (b)->Suc); Link(b, b); Link(b, (a)->Suc); Link(a, b); }
#define Precede(a, b)\
	{ Link((a)->Pred, (a)->Suc); Link(a, a); Link((b)->Pred, a); Link(a, b); }
#define SLink(a, b) { (a)->Suc = (b); (b)->Pred = (a); }

	typedef int (Heuristics::*ContainsFunction) (Node * T, int Q, Node * N);
	typedef int (Heuristics::*BoxOverlapsFunction) (Node * T, int Q, Node * N);
	int Candidates, Radius;
	ContainsFunction Contains;
	BoxOverlapsFunction BoxOverlaps;
	int Level;
	Node *FirstActive, *LastActive; 

	void Exclude(Node * ta, Node * tb)
	{
		if (ta == tb->Pred || ta == tb->Suc)
			return;
		if (ta == tb->OldPred)
			tb->OldPredExcluded = 1;
		else if (ta == tb->OldSuc)
			tb->OldSucExcluded = 1;
		if (tb == ta->OldPred)
			ta->OldPredExcluded = 1;
		else if (tb == ta->OldSuc)
			ta->OldSucExcluded = 1;
	}
	int SegmentSize(Node * ta, Node * tb)
	{
		Segment *Pa, *Pb;
		int n, nLeft, nMid, nRight;

		Pa = ta->Parent;
		Pb = tb->Parent;
		if (Pa == Pb) {
			n = Reversed == Pa->Reversed ? tb->Rank - ta->Rank :
				ta->Rank - tb->Rank;
		return (n < 0 ? n + Dimension : n) + 1;
		}
		nLeft =
			Reversed ==
			Pa->Reversed ? Pa->Last->Rank - ta->Rank : ta->Rank -
			Pa->First->Rank;
		if (nLeft < 0)
			nLeft += Pa->Size;
		nMid = !Reversed ? Pb->Rank - Pa->Rank : Pa->Rank - Pb->Rank;
		if (nMid < 0)
			nMid += Groups;
		nMid = nMid == 2 ? (!Reversed ? Pa->Suc : Pa->Pred)->Size
			: (nMid - 1) * GroupSize;
		nRight =
			Reversed ==
			Pb->Reversed ? tb->Rank -
			Pb->First->Rank : Pb->Last->Rank - tb->Rank;
		if (nRight < 0)
			nRight += Pb->Size;
		return nLeft + nMid + nRight + 2;
	}
	int Excludable(Node * ta, Node * tb)
	{
		if (ta == tb->OldPred)
			return !tb->OldPredExcluded;
		if (ta == tb->OldSuc)
			return !tb->OldSucExcluded;
		return 0;
	}
	int Between_SL(const Node * ta, const Node * tb, const Node * tc)
	{
		const Segment *Pa, *Pb, *Pc;

		if (tb == ta || tb == tc)
			return 1;
		if (ta == tc)
			return 0;
		Pa = ta->Parent;
		Pb = tb->Parent;
		Pc = tc->Parent;
		if (Pa == Pc) {
			if (Pb == Pa)
				return (Reversed == Pa->Reversed) ==
				(ta->Rank < tc->Rank ?
				tb->Rank > ta->Rank && tb->Rank < tc->Rank :
			tb->Rank > ta->Rank || tb->Rank < tc->Rank);
			return (Reversed == Pa->Reversed) == (ta->Rank > tc->Rank);
		}
		if (Pb == Pc)
			return (Reversed == Pb->Reversed) == (tb->Rank < tc->Rank);
		if (Pa == Pb)
			return (Reversed == Pa->Reversed) == (ta->Rank < tb->Rank);
		return Reversed !=
			(Pa->Rank < Pc->Rank ?
			Pb->Rank > Pa->Rank && Pb->Rank < Pc->Rank :
		Pb->Rank > Pa->Rank || Pb->Rank < Pc->Rank);
	}
	void SplitSegment(Node * t1, Node * t2)
	{
		Segment *P = t1->Parent, *Q;
		Node *t, *u;
		int i, Count, Temp;

		if (t2->Rank < t1->Rank) {
			t = t1;
			t1 = t2;
			t2 = t;
		}
		Count = t1->Rank - P->First->Rank + 1;
		if (2 * Count < P->Size) {
			/* The left part of P is merged with its neighbouring segment, Q */
			Q = P->Reversed ? P->Suc : P->Pred;
			t = P->First->Pred;
			i = t->Rank;
			if (t == Q->Last) {
				if (t == Q->First && t->Suc != P->First) {
					u = t->Suc;
					t->Suc = t->Pred;
					t->Pred = u;
					Q->Reversed ^= 1;
					Temp = t->SucCost;
					t->SucCost = t->PredCost;
					t->PredCost = Temp;
				}
				for (t = P->First; t != t2; t = t->Suc) {
					t->Parent = Q;
					t->Rank = ++i;
				}
				Q->Last = t1;
			} else {
				for (t = P->First; t != t2; t = u) {
					t->Parent = Q;
					t->Rank = --i;
					u = t->Suc;
					t->Suc = t->Pred;
					t->Pred = u;
					Temp = t->SucCost;
					t->SucCost = t->PredCost;
					t->PredCost = Temp;
				}
				Q->First = t1;
			}
			P->First = t2;
		} else {
			/* The right part of P is merged with its neighbouring segment, Q */
			Q = P->Reversed ? P->Pred : P->Suc;
			t = P->Last->Suc;
			i = t->Rank;
			if (t == Q->First) {
				if (t == Q->Last && t->Pred != P->Last) {
					u = t->Suc;
					t->Suc = t->Pred;
					t->Pred = u;
					Q->Reversed ^= 1;
					Temp = t->SucCost;
					t->SucCost = t->PredCost;
					t->PredCost = Temp;
				}
				for (t = P->Last; t != t1; t = t->Pred) {
					t->Parent = Q;
					t->Rank = --i;
				}
				Q->First = t2;
			} else {
				for (t = P->Last; t != t1; t = u) {
					t->Parent = Q;
					t->Rank = ++i;
					u = t->Pred;
					t->Pred = t->Suc;
					t->Suc = u;
					Temp = t->SucCost;
					t->SucCost = t->PredCost;
					t->PredCost = Temp;
				}
				Q->Last = t2;
			}
			Count = P->Size - Count;
			P->Last = t1;
		}
		P->Size -= Count;
		Q->Size += Count;
	}
	void Flip(Node * t1, Node * t2, Node * t3)
	{
		Node *s1, *s2, *t4;
		int R, Temp, Ct2t3, Ct4t1;

		assert(t1->Pred == t2 || t1->Suc == t2);
		if (t3 == t2->Pred || t3 == t2->Suc)
			return;
		t4 = t1->Suc == t2 ? t3->Pred : t3->Suc;
		if (t1->Suc != t2) {
			s1 = t1;
			t1 = t2;
			t2 = s1;
			s1 = t3;
			t3 = t4;
			t4 = s1;
		}
		/* Find the segment with the smallest number of nodes */
		if ((R = t2->Rank - t3->Rank) < 0)
			R += Dimension;
		if (2 * R > Dimension) {
			s1 = t3;
			t3 = t2;
			t2 = s1;
			s1 = t4;
			t4 = t1;
			t1 = s1;
		}
		Ct2t3 = graph->D(t2->Id, t3->Id);
		Ct4t1 = graph->D(t4->Id, t1->Id);
		/* Swap segment (t3 --> t1) */
		R = t1->Rank;
		t1->Suc = 0;
		s2 = t3;
		while ((s1 = s2)) {
			s2 = s1->Suc;
			s1->Suc = s1->Pred;
			s1->Pred = s2;
			s1->Rank = R--;
			Temp = s1->SucCost;
			s1->SucCost = s1->PredCost;
			s1->PredCost = Temp;
		}
		(t3->Suc = t2)->Pred = t3;
		(t4->Suc = t1)->Pred = t4;
		t3->SucCost = t2->PredCost = Ct2t3;
		t1->PredCost = t4->SucCost = Ct4t1;
		SwapStack[Swaps].t1 = t1;
		SwapStack[Swaps].t2 = t2;
		SwapStack[Swaps].t3 = t3;
		SwapStack[Swaps].t4 = t4;
		Swaps++;
	}
	void Flip_SL(Node * t1, Node * t2, Node * t3)
	{
		Node *t4, *a, *b, *cc, *d;
		Segment *P1, *P2, *P3, *P4, *Q1, *Q2;
		Node *s1, *s2;
		int i, Temp;
		assert(t1->Pred == t2 || t1->Suc == t2);
		if (t3 == t2->Pred || t3 == t2->Suc)
			return;
		if (Groups == 1) {
			Flip(t1, t2, t3);
			return;
		}
		t4 = t2 == Suc(t1) ? Pred(t3) : Suc(t3);
		P1 = t1->Parent;
		P2 = t2->Parent;
		P3 = t3->Parent;
		P4 = t4->Parent;
		/* Split segments if needed */
		if (P1 != P3 && P2 != P4) {
			if (P1 == P2) {
				SplitSegment(t1, t2);
				P1 = t1->Parent;
				P2 = t2->Parent;
			}
			if (P3 == P4 && P1 != P3 && P2 != P4) {
				SplitSegment(t3, t4);
				P3 = t3->Parent;
				P4 = t4->Parent;
			}
		} else
			if ((P1 == P3
				&& abs(t3->Rank - t1->Rank) > SPLIT_CUTOFF * GroupSize)
				|| (P2 == P4
				&& abs(t4->Rank - t2->Rank) > SPLIT_CUTOFF * GroupSize)) {
					if (P1 == P2) {
						SplitSegment(t1, t2);
						P1 = t1->Parent;
						P2 = t2->Parent;
						P3 = t3->Parent;
						P4 = t4->Parent;
					}
					if (P3 == P4) {
						SplitSegment(t3, t4);
						P1 = t1->Parent;
						P2 = t2->Parent;
						P3 = t3->Parent;
						P4 = t4->Parent;
					}
			}
			/* Check if it is possible to flip locally within a segment */
			b = 0;
			if (P1 == P3) {
				/* Either the t1 --> t3 path or the t2 --> t4 path lies 
				within one segment */
				if (t1->Rank < t3->Rank) {
					if (P1 == P2 && P1 == P4 && t2->Rank > t1->Rank) {
						a = t1;
						b = t2;
						cc = t3;
						d = t4;
					} else {
						a = t2;
						b = t1;
						cc = t4;
						d = t3;
					}
				} else {
					if (P1 == P2 && P1 == P4 && t2->Rank < t1->Rank) {
						a = t3;
						b = t4;
						cc = t1;
						d = t2;
					} else {
						a = t4;
						b = t3;
						cc = t2;
						d = t1;
					}
				}
			} else if (P2 == P4) {
				/* The t2 --> t4 path lies within one segment */
				if (t4->Rank < t2->Rank) {
					a = t3;
					b = t4;
					cc = t1;
					d = t2;
				} else {
					a = t1;
					b = t2;
					cc = t3;
					d = t4;
				}
			}
			if (b) {
				int Cbc = graph->D(b->Id, cc->Id), Cda = graph->D(d->Id, a->Id);
				/* Flip locally (b --> d) within a segment */
				i = d->Rank;
				d->Suc = 0;
				s2 = b;
				while ((s1 = s2)) {
					s2 = s1->Suc;
					s1->Suc = s1->Pred;
					s1->Pred = s2;
					s1->Rank = i--;
					Temp = s1->SucCost;
					s1->SucCost = s1->PredCost;
					s1->PredCost = Temp;
				}
				d->Pred = a;
				b->Suc = cc;
				d->PredCost = Cda;
				b->SucCost = Cbc;
				if (a->Suc == b) {
					a->Suc = d;
					a->SucCost = d->PredCost;
				} else {
					a->Pred = d;
					a->PredCost = d->PredCost;
				}
				if (cc->Pred == d) {
					cc->Pred = b;
					cc->PredCost = b->SucCost;
				} else {
					cc->Suc = b;
					cc->SucCost = b->SucCost;
				}
				if (b->Parent->First == b)
					b->Parent->First = d;
				else if (d->Parent->First == d)
					d->Parent->First = b;
				if (b->Parent->Last == b)
					b->Parent->Last = d;
				else if (d->Parent->Last == d)
					d->Parent->Last = b;
			} else {
				int Ct2t3, Ct4t1;
				/* Reverse a sequence of segments */
				if (P1->Suc != P2) {
					a = t1;
					t1 = t2;
					t2 = a;
					a = t3;
					t3 = t4;
					t4 = a;
					Q1 = P1;
					P1 = P2;
					P2 = Q1;
					Q1 = P3;
					P3 = P4;
					P4 = Q1;
				}
				/* Find the sequence with the smallest number of segments */
				if ((i = P2->Rank - P3->Rank) < 0)
					i += Groups;
				if (2 * i > Groups) {
					a = t3;
					t3 = t2;
					t2 = a;
					a = t1;
					t1 = t4;
					t4 = a;
					Q1 = P3;
					P3 = P2;
					P2 = Q1;
					Q1 = P1;
					P1 = P4;
					P4 = Q1;
				}
				Ct2t3 = graph->D(t2->Id, t3->Id);
				Ct4t1 = graph->D(t4->Id, t1->Id);
				/* Reverse the sequence of segments (P3 --> P1). 
				Mirrors the corresponding code in the Flip function */
				i = P1->Rank;
				P1->Suc = 0;
				Q2 = P3;
				while ((Q1 = Q2)) {
					Q2 = Q1->Suc;
					Q1->Suc = Q1->Pred;
					Q1->Pred = Q2;
					Q1->Rank = i--;
					Q1->Reversed ^= 1;
				}
				P3->Suc = P2;
				P2->Pred = P3;
				P1->Pred = P4;
				P4->Suc = P1;
				if (t3->Suc == t4) {
					t3->Suc = t2;
					t3->SucCost = Ct2t3;
				} else {
					t3->Pred = t2;
					t3->PredCost = Ct2t3;
				}
				if (t2->Suc == t1) {
					t2->Suc = t3;
					t2->SucCost = Ct2t3;
				} else {
					t2->Pred = t3;
					t2->PredCost = Ct2t3;
				}
				if (t1->Pred == t2) {
					t1->Pred = t4;
					t1->PredCost = Ct4t1;
				} else {
					t1->Suc = t4;
					t1->SucCost = Ct4t1;
				}
				if (t4->Pred == t3) {
					t4->Pred = t1;
					t4->PredCost = Ct4t1;
				} else {
					t4->Suc = t1;
					t4->SucCost = Ct4t1;
				}
			}
			SwapStack[Swaps].t1 = t1;
			SwapStack[Swaps].t2 = t2;
			SwapStack[Swaps].t3 = t3;
			SwapStack[Swaps].t4 = t4;
			Swaps++;
			//Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
			//    (Rand[t3->Id] * Rand[t4->Id]) ^
			//    (Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
	}
	void AllocateStructures()
	{    
		int i, K = 5;
		SwapStack = new SwapRecord[MaxSwaps + 6 * K];
	}
	void Activate(Node * N)
	{
		if (N->Next != 0)
			return;
		if (FirstActive == 0)
			FirstActive = LastActive = N;
		else
			LastActive = LastActive->Next = N;
		LastActive->Next = FirstActive;
	}
	Node* RemoveFirstActive()
	{
		Node *N = FirstActive;
		if (FirstActive == LastActive)
			FirstActive = LastActive = 0;
		else 
			LastActive->Next = FirstActive = FirstActive->Next;
		if (N)
			N->Next = 0;
		return N;
	}
	int Contains2D(Node * T, int Q, Node * N)
	{
		switch (Q) {
		case 1:
			return T->X >= N->X && T->Y >= N->Y;
		case 2:
			return T->X <= N->X && T->Y >= N->Y;
		case 3:
			return T->X <= N->X && T->Y <= N->Y;
		case 4:
			return T->X >= N->X && T->Y <= N->Y;
		default:
			return 1;
		}
	}
	int Contains3D(Node * T, int Q, Node * N)
	{
		switch (Q) {
		case 1:
			return T->X >= N->X && T->Y >= N->Y && T->Z >= N->Z;
		case 2:
			return T->X <= N->X && T->Y >= N->Y && T->Z >= N->Z;
		case 3:
			return T->X <= N->X && T->Y <= N->Y && T->Z >= N->Z;
		case 4:
			return T->X >= N->X && T->Y <= N->Y && T->Z >= N->Z;
		case 5:
			return T->X >= N->X && T->Y >= N->Y && T->Z <= N->Z;
		case 6:
			return T->X <= N->X && T->Y >= N->Y && T->Z <= N->Z;
		case 7:
			return T->X <= N->X && T->Y <= N->Y && T->Z <= N->Z;
		case 8:
			return T->X >= N->X && T->Y < N->Y && T->Z <= N->Z;
		default:
			return 1;
		}
	}
	int BoxOverlaps2D(Node * T, int Q, Node * N)
	{
		switch (Q) {
		case 1:
			return XMax[T->Id] >= N->X && YMax[T->Id] >= N->Y;
		case 2:
			return XMin[T->Id] <= N->X && YMax[T->Id] >= N->Y;
		case 3:
			return XMin[T->Id] <= N->X && YMin[T->Id] <= N->Y;
		case 4:
			return XMax[T->Id] >= N->X && YMin[T->Id] <= N->Y;
		default:
			return 1;
		}
	}
	int BoxOverlaps3D(Node * T, int Q, Node * N)
	{
		switch (Q) {
		case 1:
			return XMax[T->Id] >= N->X && YMax[T->Id] >= N->Y &&
				ZMax[T->Id] >= N->Z;
		case 2:
			return XMin[T->Id] <= N->X && YMax[T->Id] >= N->Y &&
				ZMax[T->Id] >= N->Z;
		case 3:
			return XMin[T->Id] <= N->X && YMin[T->Id] <= N->Y &&
				ZMax[T->Id] >= N->Z;
		case 4:
			return XMax[T->Id] >= N->X && YMin[T->Id] <= N->Y &&
				ZMax[T->Id] >= N->Z;
		case 5:
			return XMax[T->Id] >= N->X && YMax[T->Id] >= N->Y &&
				ZMin[T->Id] <= N->Z;
		case 6:
			return XMin[T->Id] <= N->X && YMax[T->Id] >= N->Y &&
				ZMin[T->Id] <= N->Z;
		case 7:
			return XMin[T->Id] <= N->X && YMin[T->Id] <= N->Y &&
				ZMin[T->Id] <= N->Z;
		case 8:
			return XMax[T->Id] >= N->X && YMin[T->Id] <= N->Y &&
				ZMin[T->Id] <= N->Z;
		default:
			return 1;
		}
	}
	Node** BuildKDTree(int Cutoff)
	{
		int i;
		Node *N;

		cutoff = Cutoff >= 1 ? Cutoff : 1;
		KDTree = new Node*[Dimension];
		for (i = 0, N = FirstNode; i < Dimension; i++, N = N->Suc)
			KDTree[i] = N;
		BuildSubKDTree(0, Dimension - 1);
		return KDTree;
	}
	void BuildSubKDTree(int start, int end)
	{
		if (end - start + 1 > cutoff) {
			int mid = (start + end) / 2;
			char axis = FindMaxSpread(start, end);
			Partition(start, end, mid, axis);
			KDTree[mid]->Axis = axis;
			BuildSubKDTree(start, mid - 1);
			BuildSubKDTree(mid + 1, end);
		}
	}
	void Partition(int start, int end, int k, int axis)
	{
		while (start < end) {
			int i = start, j = end - 1, mid = (start + end) / 2;
			double pivot;
			if (Coord(KDTree[mid], axis) < Coord(KDTree[start], axis))
				Swap(start, mid);
			if (Coord(KDTree[end], axis) < Coord(KDTree[start], axis))
				Swap(start, end);
			if (Coord(KDTree[end], axis) < Coord(KDTree[mid], axis))
				Swap(mid, end);
			if (end - start <= 2)
				return;
			Swap(mid, j);
			pivot = Coord(KDTree[j], axis);
			while (1) {
				while (Coord(KDTree[++i], axis) < pivot);
				while (pivot < Coord(KDTree[--j], axis));
				if (i >= j)
					break;
				Swap(i, j);
			}
			Swap(i, end - 1);
			if (i >= k)
				end = i - 1;
			if (i <= k)
				start = i + 1;
		}
	}
	void Swap(int i, int j)
	{
		Node *T = KDTree[i];
		KDTree[i] = KDTree[j];
		KDTree[j] = T;
	}
	char FindMaxSpread(int start, int end)
	{
		int i, axis;
		Node *N;
		double Min[3], Max[3];

		N = KDTree[start];
		Min[0] = Max[0] = N->X;
		Min[1] = Max[1] = N->Y;
		Min[2] = Max[2] = N->Z;
		for (i = start + 1; i <= end; i++) {
			for (axis = 1; axis >= 0; axis--) {
				N = KDTree[i];
				if (Coord(N, axis) < Min[axis])
					Min[axis] = Coord(N, axis);
				else if (Coord(N, axis) > Max[axis])
					Max[axis] = Coord(N, axis);
			}
		}
		if (Max[0] - Min[0] > Max[1] - Min[1])
			return 0;
		return 1;
	}
	int AddCandidate(Node * From, Node * To, int Cost, int Alpha)
	{
		int Count;
		Candidate *NFrom;

		if (From->Subproblem != FirstNode->Subproblem)
			return 0;
		if (From->CandidateSet == 0)
			From->CandidateSet = (Candidate *) calloc(this->Number_Of_Candidates, sizeof(Candidate));
		if (From == To || To->Subproblem != FirstNode->Subproblem)
			return 0;
		Count = 0;
		for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To; NFrom++)
			Count++;
		if (NFrom->To)
			return 0;
		NFrom->Cost = Cost;
		NFrom->Alpha = Alpha;
		NFrom->To = To;
		From->CandidateSet =(Candidate *) realloc(From->CandidateSet,(Count + 2) * sizeof(Candidate));
		From->CandidateSet[Count + 1].To = 0;
		return 1;
	}	
	void NQN(Node * N, int Q, int start, int end, int K)
	{
		double diff;
		int mid = (start + end) / 2, d, i;
		Node *T = KDTree[mid], P;
		int axis = T->Axis;

		if (start <= end && T != N && (this->*Contains)(T, Q, N) &&
			!InCandidateSet(N, T) &&        
			(d = this->graph->D(N->X, N->Y, T->X, T->Y)) <= Radius) {
				i = Candidates;
				while (--i >= 0 && d < CandidateSet[i].Cost)
					CandidateSet[i + 1] = CandidateSet[i];
				CandidateSet[i + 1].To = T;
				CandidateSet[i + 1].Cost = d;
				if (Candidates < K)
					Candidates++;
				if (Candidates == K)
					Radius = CandidateSet[Candidates - 1].Cost;
		}
		if (start < end && (this->*BoxOverlaps)(T, Q, N)) {
			P.X = axis == 0 ? T->X : N->X;
			P.Y = axis == 1 ? T->Y : N->Y;
			P.Z = axis == 2 ? T->Z : N->Z;

			diff = Coord(T, axis) - Coord(N, axis);
			if (diff >= 0) {
				if (Overlaps(Q, diff, 0, axis))
					NQN(N, Q, start, mid - 1, K);
				if (Overlaps(Q, diff, 1, axis) &&                
					this->graph->D(N->X, N->Y, (&P)->X, (&P)->Y)  <= Radius)
					NQN(N, Q, mid + 1, end, K);
			} else {
				if (Overlaps(Q, diff, 1, axis))
					NQN(N, Q, mid + 1, end, K);
				if (Overlaps(Q, diff, 0, axis) &&                
					this->graph->D(N->X, N->Y, (&P)->X, (&P)->Y)  <= Radius)
					NQN(N, Q, start, mid - 1, K);
			}
		}
	}
	void NearestQuadrantNeighbors(Node * N, int Q, int K)
	{
		Candidates = 0;
		Radius = INT_MAX;
		NQN(N, Q, 0, Dimension - 1, K);
	}
	int Overlaps(int Q, double diff, int High, int axis)
	{
		switch (Q) {
		case 1:
			return High || diff >= 0;
		case 2:
			return axis == 0 ? !High || diff <= 0 : High || diff >= 0;
		case 3:
			return axis <= 1 ? !High || diff <= 0 : High || diff >= 0;
		case 4:
			return axis == 1 ? !High || diff <= 0 : High || diff >= 0;
		case 5:
			return axis <= 1 ? High || diff >= 0 : !High || diff <= 0;
		case 6:
			return axis == 1 ? High || diff >= 0 : !High || diff <= 0;
		case 7:
			return !High || diff <= 0;
		case 8:
			return axis == 0 ? High || diff >= 0 : !High || diff <= 0;
		default:
			return 1;
		}
	}
	int IsCandidate(Node * ta, Node * tb)
	{
		Candidate *Nta;

		for (Nta = ta->CandidateSet; Nta && Nta->To; Nta++)
		{
			if (Nta->To == tb)
				return 1;
		}
		return 0;
	}
	int InCandidateSet(Node * N, Node * T)
	{
		int i;
		for (i = 0; i < Candidates; i++)
			if (CandidateSet[i].To == T)
				return 1;
		return IsCandidate(N, T);
	}
	void FreeCandidateSets()
	{
		Node *N = FirstNode;
		if (!N)
			return;
		do {
			Free(N->CandidateSet);
		}
		while ((N = N->Suc) != FirstNode);
	}
	void ComputeBounds(int start, int end)
	{
		if (start <= end) 
        {
			int mid = (start + end) / 2, i;
 			Node *T = KDTree[mid];
           
 			XMin[T->Id] = YMin[T->Id] = DBL_MAX;
 			XMax[T->Id] = YMax[T->Id] = -DBL_MAX;
 			for (i = start; i <= end; i++) {
 				Node *N = KDTree[i];
 				if (N == T)
 					continue;
 				if (N->X < XMin[T->Id])
 					XMin[T->Id] = N->X;
 				if (N->X > XMax[T->Id])
 					XMax[T->Id] = N->X;
 				if (N->Y < YMin[T->Id])
 					YMin[T->Id] = N->Y;
 				if (N->Y > YMax[T->Id])
 					YMax[T->Id] = N->Y;            
 			}
			ComputeBounds(start, mid - 1);
			ComputeBounds(mid + 1, end);
		}
	}
	void CreateNearestNeighborCandidateSet(int K)
	{
		Node *From, *To;
		int i;


		KDTree = BuildKDTree(1);
		XMin =	(double *) malloc((1 + Dimension) * sizeof(double));
		XMax =  (double *) malloc((1 + Dimension) * sizeof(double));
		YMin =	(double *) malloc((1 + Dimension) * sizeof(double));
		YMax =	(double *) malloc((1 + Dimension) * sizeof(double));
 		ComputeBounds(0, Dimension - 1);
		Contains = &Heuristics::Contains2D;
		BoxOverlaps =&Heuristics::BoxOverlaps2D;
		CandidateSet =	(Candidate *) malloc((K + 1) * sizeof(Candidate));

		From = FirstNode;
		do {
			NearestQuadrantNeighbors(From, 0, K);
			for (i = 0; i < Candidates; i++) {
				To = CandidateSet[i].To;
				AddCandidate(From, To, this->graph->D(From->X,From->Y, To->X,To->Y), 1);
			}
		} while ((From = From->Suc) != FirstNode);

		free(CandidateSet);
		free(KDTree);
		free(XMin);
		free(XMax);
		free(YMin);
		free(YMax);
	}
	void CreateNodes()
	{
		NodeSet = (Node**)malloc(sizeof(Node*) * Dimension);
		for (int i = 0 ; i < Dimension ; i++ )
		{
			NodeSet[i] = new Node();
			NodeSet[i]->Id = i;
			NodeSet[i]->X = this->graph->X(i);
			NodeSet[i]->Y = this->graph->Y(i);
		}
		Node *Prev = 0, *N = 0;

		if (Dimension <= 0)
			printf("DIMENSION is not positive (or not specified)");

		for (int i = 0; i < Dimension; i++, Prev = N) {
			N = NodeSet[i];
			if (i == 0)
				FirstNode = N;
			else
				Link(Prev, N);
			N->Id = i; 
		}
		Link(N, FirstNode);
	}
	void Make2OptMove(Node * t1, Node * t2, Node * t3, Node * t4)
	{
		Swap1(t1, t2, t3);
	}
	void Make3OptMove(Node * t1, Node * t2, Node * t3, Node * t4, Node * t5, Node * t6, int Case)
	{
		switch (Case) {
		case 1:
		case 2:
			Swap2(t1, t2, t3, t6, t5, t4);
			return;
		case 5:
			Swap3(t1, t2, t4, t6, t5, t4, t6, t2, t3);
			return;
		case 6:
			Swap2(t3, t4, t5, t1, t2, t3);
			return;
		default:
			printf("Make3OptMove: Internal error");
		}
	}
	void Make4OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
		Node * t5, Node * t6, Node * t7, Node * t8, int Case)
	{
		if (Suc(t1) != t2)
			Reversed ^= 1;
		switch (Case) {
		case 1:
		case 2:
			Swap3(t1, t2, t3, t6, t5, t4, t7, t8, t1);
			return;
		case 3:
		case 4:
			Swap3(t1, t2, t3, t8, t7, t6, t5, t8, t1);
			return;
		case 5:
			if (!BETWEEN(t2, t7, t3))
				Swap3(t5, t6, t7, t2, t1, t4, t1, t4, t5);
			else if (BETWEEN(t2, t7, t6))
				Swap3(t5, t6, t7, t5, t8, t3, t3, t8, t1);
			else
				Swap3(t1, t2, t7, t7, t2, t3, t4, t7, t6);
			return;
		case 6:
			Swap3(t3, t4, t5, t6, t3, t2, t1, t6, t7);
			return;
		case 7:
			Swap3(t6, t5, t8, t2, t1, t4, t8, t5, t4);
			return;
		case 11:
			Swap3(t1, t2, t7, t3, t4, t5, t3, t6, t7);
			return;
		case 12:
			Swap3(t3, t4, t5, t7, t8, t1, t3, t6, t7);
			return;
		case 15:
			Swap3(t3, t4, t5, t3, t6, t7, t8, t3, t2);
			return;
		default:
			printf("Make4OptMove: Internal error");
		}
	}
	void Make5OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
		Node * t5, Node * t6, Node * t7, Node * t8,
		Node * t9, Node * t10, int Case)
	{
		if (Suc(t1) != t2)
			Reversed ^= 1;
		switch (Case) {
		case 1:
			Swap4(t1, t2, t3, t8, t7, t6, t10, t9, t8, t10, t5, t4);
			return;
		case 2:
			if (BETWEEN(t2, t9, t4))
				Swap4(t1, t2, t3, t5, t6, t7, t10, t9, t8, t5, t10, t1);
			else
				Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
			return;
		case 3:
			Swap4(t3, t4, t5, t7, t8, t9, t1, t2, t3, t7, t10, t1);
			return;
		case 4:
			Swap5(t5, t6, t8, t1, t2, t3, t10, t9, t8, t1, t4, t5, t6, t10,
				t1);
			return;
		case 5:
			Swap5(t5, t6, t10, t1, t2, t3, t6, t10, t1, t8, t7, t6, t8, t4,
				t5);
			return;
		case 6:
			Swap4(t1, t2, t3, t9, t10, t1, t7, t8, t9, t6, t5, t4);
			return;
		case 7:
			if (BETWEEN(t3, t9, t7))
				Swap4(t3, t4, t5, t8, t7, t6, t10, t9, t8, t1, t2, t3);
			else if (BETWEEN(t6, t9, t4))
				Swap4(t3, t4, t5, t8, t7, t6, t9, t10, t1, t9, t2, t3);
			else
				Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
			return;
		case 8:
			Swap4(t3, t4, t5, t9, t10, t1, t8, t7, t6, t8, t3, t2);
			return;
		case 9:
			Swap4(t10, t9, t8, t5, t6, t7, t1, t2, t3, t1, t4, t5);
			return;
		case 10:
			if (BETWEEN(t5, t9, t7))
				Swap4(t5, t6, t7, t9, t10, t1, t4, t3, t2, t4, t9, t8);
			else if (BETWEEN(t3, t9, t6))
				Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
			else
				Swap4(t1, t2, t3, t9, t10, t1, t5, t6, t7, t5, t8, t9);
			return;
		case 11:
			if (BETWEEN(t3, t9, t6))
				Swap4(t1, t2, t3, t6, t5, t4, t9, t10, t1, t7, t8, t9);
			else
				Swap4(t5, t6, t7, t10, t9, t8, t2, t1, t10, t4, t3, t2);
			return;
		case 12:
			Swap4(t1, t2, t3, t8, t7, t6, t10, t9, t8, t5, t10, t1);
			return;
		case 13:
			if (BETWEEN(t4, t9, t7))
				Swap5(t7, t8, t10, t5, t6, t7, t1, t2, t3, t5, t9, t1, t9, t1,
				t10);
			else if (BETWEEN(t6, t9, t3))
				Swap5(t10, t9, t1, t7, t8, t9, t3, t4, t5, t3, t6, t7, t3, t1,
				t10);
			else
				Swap5(t10, t9, t1, t4, t3, t2, t5, t6, t7, t5, t8, t10, t9, t1,
				t10);
			return;
		case 14:
			Swap5(t10, t9, t1, t5, t6, t7, t5, t8, t9, t3, t4, t5, t3, t1,
				t10);
			return;
		case 15:
			if (BETWEEN(t6, t9, t3))
				Swap5(t10, t9, t1, t3, t4, t5, t6, t3, t2, t8, t7, t6, t9, t1,
				t10);
			else
				Swap5(t1, t2, t6, t3, t4, t5, t8, t7, t6, t10, t9, t8, t2, t10,
				t1);
			return;
		case 16:
			if (BETWEEN(t4, t9, t7))
				Swap4(t3, t4, t5, t8, t7, t6, t9, t10, t1, t8, t3, t2);
			else if (BETWEEN(t5, t9, t3))
				Swap4(t3, t4, t5, t9, t10, t1, t6, t3, t2, t7, t8, t9);
			else
				Swap4(t3, t4, t5, t1, t2, t3, t7, t8, t9, t7, t10, t1);
			return;
		case 17:
			if (BETWEEN(t7, t9, t3))
				Swap4(t3, t4, t5, t7, t8, t9, t2, t1, t10, t3, t6, t7);
			else
				Swap4(t7, t8, t9, t2, t1, t10, t3, t4, t5, t3, t6, t7);
			return;
		case 18:
			Swap4(t3, t4, t5, t7, t8, t9, t3, t6, t7, t1, t2, t3);
			return;
		case 19:
			Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
			return;
		case 20:
			Swap4(t7, t8, t9, t3, t4, t5, t10, t7, t6, t3, t10, t1);
			return;
		case 21:
			Swap4(t5, t6, t7, t5, t8, t9, t1, t2, t3, t4, t1, t10);
			return;
		case 22:
			Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t1, t9, t10, t1);
			return;
		case 23:
			Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t1, t9, t10, t1);
			return;
		case 24:
			Swap4(t1, t2, t3, t8, t7, t6, t5, t8, t1, t9, t10, t1);
			return;
		case 25:
			Swap4(t1, t2, t3, t8, t7, t6, t5, t8, t1, t9, t10, t1);
			return;
		case 26:
			if (!BETWEEN(t2, t7, t3))
				Swap4(t5, t6, t7, t2, t1, t4, t1, t4, t5, t9, t10, t1);
			else if (BETWEEN(t2, t7, t6))
				Swap4(t5, t6, t7, t5, t8, t3, t3, t8, t1, t9, t10, t1);
			else
				Swap4(t1, t2, t7, t7, t2, t3, t4, t7, t6, t9, t10, t1);
			return;
		case 27:
			Swap4(t3, t4, t5, t6, t3, t2, t1, t6, t7, t9, t10, t1);
			return;
		case 28:
			Swap4(t6, t5, t8, t2, t1, t4, t8, t5, t4, t9, t10, t1);
			return;
		case 29:
			Swap4(t1, t2, t7, t3, t4, t5, t3, t6, t7, t9, t10, t1);
			return;
		case 30:
			if (BETWEEN(t3, t7, t5))
				Swap4(t3, t4, t5, t7, t8, t1, t7, t2, t3, t9, t10, t1);
			else
				Swap4(t3, t4, t5, t3, t6, t7, t1, t2, t3, t9, t10, t1);
			return;
		case 31:
			Swap4(t3, t4, t5, t3, t6, t7, t8, t3, t2, t9, t10, t1);
			return;
		case 32:
			Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
			return;
		case 33:
			if (BETWEEN(t3, t9, t5))
				Swap4(t1, t2, t3, t5, t6, t7, t10, t9, t8, t5, t10, t1);
			else
				Swap4(t1, t2, t3, t7, t8, t9, t7, t10, t1, t5, t6, t7);
			return;
		case 34:
			Swap4(t7, t8, t9, t1, t2, t3, t1, t4, t5, t7, t10, t1);
			return;
		case 35:
			Swap4(t9, t10, t1, t5, t6, t7, t4, t3, t2, t9, t4, t5);
			return;
		case 36:
			Swap4(t9, t10, t1, t7, t8, t9, t3, t4, t5, t6, t3, t2);
			return;
		case 37:
			if (BETWEEN(t6, t9, t4))
				Swap4(t1, t2, t3, t6, t5, t4, t9, t10, t1, t8, t7, t6);
			else
				Swap4(t9, t10, t1, t3, t4, t5, t3, t6, t7, t3, t8, t9);
			return;
		case 38:
			if (BETWEEN(t3, t9, t7))
				Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t6, t1, t10);
			else if (BETWEEN(t6, t9, t4))
				Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
			else
				Swap4(t3, t4, t5, t9, t10, t1, t8, t7, t6, t3, t8, t9);
			return;
		case 39:
			Swap4(t1, t2, t3, t7, t8, t9, t5, t6, t7, t1, t4, t5);
			return;
		case 40:
			Swap4(t9, t10, t1, t4, t3, t2, t5, t6, t7, t5, t8, t9);
			return;
		case 41:
			if (BETWEEN(t5, t9, t7))
				Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
			else if (BETWEEN(t3, t9, t6))
				Swap4(t1, t2, t3, t5, t6, t7, t9, t10, t1, t5, t8, t9);
			else
				Swap4(t5, t6, t7, t9, t10, t1, t2, t9, t8, t3, t4, t5);
			return;
		case 42:
			if (BETWEEN(t3, t9, t6))
				Swap4(t7, t8, t9, t5, t6, t7, t1, t2, t3, t1, t4, t5);
			else
				Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
			return;
		case 43:
			Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
			return;
		case 44:
			if (BETWEEN(t4, t9, t7))
				Swap4(t7, t8, t9, t5, t6, t7, t1, t2, t3, t5, t10, t1);
			else if (BETWEEN(t6, t9, t3))
				Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
			else
				Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
			return;
		case 45:
			Swap4(t9, t10, t1, t3, t4, t5, t7, t8, t9, t3, t6, t7);
			return;
		case 46:
			if (BETWEEN(t6, t9, t3))
				Swap4(t7, t8, t9, t5, t6, t7, t3, t4, t5, t1, t2, t3);
			else
				Swap4(t7, t8, t9, t5, t6, t7, t3, t4, t5, t1, t2, t3);
			return;
		case 47:
			if (BETWEEN(t4, t9, t7))
				Swap4(t5, t6, t7, t1, t2, t3, t9, t10, t1, t5, t8, t9);
			else if (BETWEEN(t5, t9, t3))
				Swap4(t9, t10, t1, t7, t8, t9, t5, t6, t7, t3, t4, t5);
			else
				Swap4(t7, t8, t9, t3, t4, t5, t3, t6, t7, t2, t1, t10);
			return;
		case 48:
			if (BETWEEN(t7, t9, t3))
				Swap4(t3, t4, t5, t8, t7, t6, t2, t1, t10, t8, t3, t2);
			else
				Swap4(t3, t4, t5, t7, t8, t9, t3, t6, t7, t1, t2, t3);
			return;
		case 49:
			Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
			return;
		case 50:
			Swap4(t3, t4, t5, t3, t6, t7, t9, t10, t1, t8, t3, t2);
			return;
		case 51:
			Swap4(t5, t6, t7, t1, t2, t3, t9, t10, t1, t4, t9, t8);
			return;
		case 52:
			Swap4(t5, t6, t7, t3, t4, t5, t9, t10, t1, t3, t8, t9);
			return;
		default:
			printf("Make5OptMove: Internal error");
		}
	}
	void RestoreTour()
	{
		Node *t1, *t2, *t3, *t4;

		/* Loop as long as the stack is not empty */
		while (Swaps > 0) {
			/* Undo topmost 2-opt move */
			Swaps--;
			t1 = SwapStack[Swaps].t1;
			t2 = SwapStack[Swaps].t2;
			t3 = SwapStack[Swaps].t3;
			t4 = SwapStack[Swaps].t4;
			Swap1(t3, t2, t1);
			Swaps--;
			/* Make edges (t1,t2) and (t2,t3) excludable again */
			t1->OldPredExcluded = t1->OldSucExcluded = 0;
			t2->OldPredExcluded = t2->OldSucExcluded = 0;
			t3->OldPredExcluded = t3->OldSucExcluded = 0;
			t4->OldPredExcluded = t4->OldSucExcluded = 0;
		}
	}
	void StoreTour()
	{
		Node *t, *u;
		Candidate *Nt;
		int i;
		while (Swaps > 0) 
		{
			Swaps--;
			for (i = 1; i <= 4; i++) 
			{
				t = i == 1 ? SwapStack[Swaps].t1 :
					i == 2 ? SwapStack[Swaps].t2 :	
					i == 3 ? SwapStack[Swaps].t3 : SwapStack[Swaps].t4;	
			Activate(t);
			t->OldPred = t->Pred;
			t->OldSuc = t->Suc;
			t->OldPredExcluded = t->OldSucExcluded = 0;
			t->Cost = INT_MAX;
			for (Nt = t->CandidateSet; (u = Nt->To); Nt++)
				if (u != t->Pred && u != t->Suc && Nt->Cost < t->Cost)
					t->Cost = Nt->Cost;
			}
		}
	}
	Node* Best2OptMove(Node *t1, Node *t2, long & G0, long & Gain)
	{
		Node *t3;
		Node *t4;
		Node *T3 = 0;
		Node *T4 = 0;
		long G1;
		long G2;
		long BestG2 = LONG_MAX;
		if (t2 != Suc(t1))
			Reversed ^= true;
		Candidate *Nt2;
		for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {

			//i2++;
			//t3 = Nt2.To;
			if (t3 == Pred(t2) || t3 == Suc(t2) || ((G1 = G0 - this->graph->D(t2->Id, t3->Id)) <= 0))
				continue;
			t4 = Pred(t3);
			G2 = G1 + graph->D(t3->Id, t4->Id);
			if (G2 > BestG2 && Swaps < MaxSwaps && Excludable(t3, t4))
			{
				T3 = t3;
				T4 = t4;
				BestG2 = G2;
			}

		}
		Gain = 0;
		if (T4 != 0)
		{
			/* Make the best 2-opt move */
			Swap1(t1, t2, T3);
			Exclude(t1, t2);
			Exclude(T3, T4);
			G0 = BestG2;
		}
		return T4;
	}

	Node* Best5OptMove(Node * t1, Node * t2, long & G0, long & Gain)
	{
		Node *t3, *t4, *t5, *t6 = 0, *t7, *t8 = 0, *t9, *t10 = 0;
		Node *T3 = 0, *T4 = 0, *T5 = 0, *T6 = 0, *T7 = 0, *T8 = 0, *T9 =
			0, *T10 = 0;
		Candidate *Nt2, *Nt4, *Nt6, *Nt8;
		int G1, G2, G3, G4, G5, G6, G7, G8, BestG8 = LONG_MIN;
		int Case6 = 0, Case8 = 0, Case10 = 0, BestCase10 = 0, X4, X6, X8, X10,
			BTW275 = 0, BTW674 = 0,
			BTW571 = 0, BTW376 = 0, BTW574 = 0, BTW671 = 0,
			BTW471 = 0, BTW673 = 0, BTW573 = 0, BTW273 = 0;


		if (t2 != Suc(t1))
			Reversed ^= 1;


		/* 
		* Determine (T3,T4,T5,T6,T7,T8,T9,T10) = (t3,t4,t5,t6,t7,t8,t9,t10)
		* such that
		*
		*     G8 = *G0 - C(t2,T3) + C(T3,T4)
		*              - C(T4,T5) + C(T5,T6)
		*              - C(T6,T7) + C(T7,T8)
		*              - C(T8,T9) + C(T9,T10)
		*
		* is maximum (= BestG8), and (T9,T10) has not previously been included.
		* If during this process a legal move with *Gain > 0 is found, then make
		* the move and exit Best5OptMove immediately. 
		*/

		/* Choose (t2,t3) as a candidate edge emanating from t2 */	
		for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {		
			if (t3 == t2->Pred || t3 == t2->Suc ||
				((G1 = G0 - Nt2->Cost) <= 0))
				continue;
			/* Choose t4 as one of t3's two neighbors on the tour */
			for (X4 = 1; X4 <= 2; X4++) {
				t4 = X4 == 1 ? Pred(t3) : Suc(t3);
				G2 = G1 + graph->D(t3->Id, t4->Id);
				if (X4 == 1 && (Gain = G2 - graph->D(t4->Id, t1->Id)) > 0)
				{
					Make2OptMove(t1, t2, t3, t4);
					return 0;
				}
				/* Choose (t4,t5) as a candidate edge emanating from t4 */
				for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {				
					if (t5 == t4->Pred || t5 == t4->Suc ||
						((G3 = G2 - Nt4->Cost) <= 0))
						continue;

					/* Choose t6 as one of t5's two neighbors on the tour */
					for (X6 = 1; X6 <= 2; X6++) {
						if (X4 == 1) {
							if (X6 == 1) {
								Case6 = 1 + !BETWEEN(t2, t5, t4);
								t6 = Case6 == 1 ? Suc(t5) : Pred(t5);
							} else {
								t6 = t6 == t5->Pred ? t5->Suc : t5->Pred;
								if ((t5 == t1 && t6 == t2) ||
									(t5 == t2 && t6 == t1))
									continue;
								Case6 += 2;
							}
						} else if (BETWEEN(t2, t5, t3)) {
							Case6 = 4 + X6;
							t6 = X6 == 1 ? Suc(t5) : Pred(t5);
							if (t6 == t1)
								continue;
						} else {
							Case6 = 6 + X6;
							t6 = X6 == 1 ? Pred(t5) : Suc(t5);
							if (t6 == t2)
								continue;
						}
						G4 = G3 + graph->D(t5->Id, t6->Id);
						if ((Case6 <= 2 || Case6 == 5 || Case6 == 6) &&
							(Gain = G4 - graph->D(t6->Id, t1->Id)) > 0) {
								Make3OptMove(t1, t2, t3, t4, t5, t6, Case6);
								return 0;
						}
						/* Choose (t6,t7) as a candidate edge emanating from t6 */
						for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
							if (t7 == t6->Pred || t7 == t6->Suc ||
								(t6 == t2 && t7 == t3) ||
								(t6 == t3 && t7 == t2) ||
								((G5 = G4 - Nt6->Cost) <= 0))
								continue;
							/* Choose t8 as one of t7's two neighbors on the tour */
							for (X8 = 1; X8 <= 2; X8++) {
								if (X8 == 1) {
									Case8 = Case6;
									switch (Case6) {
									case 1:
										if ((BTW275 = BETWEEN(t2, t7, t5)))
											t8 = Suc(t7);
										else {
											t8 = Pred(t7);
											BTW674 = BETWEEN(t6, t7, t4);
										}
										break;
									case 2:
										if ((BTW376 = BETWEEN(t3, t7, t6)))
											t8 = Suc(t7);
										else {
											t8 = Pred(t7);
											BTW571 = BETWEEN(t5, t7, t1);
										}
										break;
									case 3:
										t8 = Suc(t7);
										BTW574 = BETWEEN(t5, t7, t4);
										break;
									case 4:
										if ((BTW671 = BETWEEN(t6, t7, t1)))
											t8 = Pred(t7);
										else
											t8 = BETWEEN(t2, t7,
											t4) ? Suc(t7) :
											Pred(t7);
										break;
									case 5:
										t8 = Pred(t7);
										BTW471 = BETWEEN(t4, t7, t1);
										if (!BTW471)
											BTW673 = BETWEEN(t6, t7, t3);
										break;
									case 6:
										if ((BTW471 = BETWEEN(t4, t7, t1)))
											t8 = Pred(t7);
										else {
											t8 = Suc(t7);
											BTW573 = BETWEEN(t5, t7, t3);
										}
										break;
									case 7:
									case 8:
										t8 = Suc(t7);
										BTW273 = BETWEEN(t2, t7, t3);
										break;
									}
								} else {
									t8 = t8 == t7->Pred ? t7->Suc : t7->Pred;
									Case8 += 8;
								}
								if ((t7 == t1 && t8 == t2) ||
									(t7 == t2 && t8 == t1) ||
									(t7 == t3 && t8 == t4) ||
									(t7 == t4 && t8 == t3))
									continue;
								if (Case6 == 3 && !BTW574 &&
									(X8 == 1) == BETWEEN(t3, t7, t1))
									continue;
								if (Case6 == 4 && BTW671 && X8 == 2)
									break;
								if (Case6 == 7 && !BTW273 &&
									(X8 == 1) == BETWEEN(t5, t7, t1))
									continue;
								if (Case6 == 8 && !BTW273
									&& !BETWEEN(t4, t7, t5))
									break;
								G6 = G5 + graph->D(t7->Id, t8->Id);
								if (t8 != t1 &&
									(Case6 == 3 ? BTW574 :
									Case6 == 4 ? !BTW671 :
									Case6 == 7 ? BTW273 :
									Case6 != 8 && X8 == 1) &&
									(Gain = G6 - graph->D(t8->Id, t1->Id)) > 0) {
										Make4OptMove(t1, t2, t3, t4, t5, t6, t7,
											t8, Case8);
										return 0;
								}
								/* Choose (t8,t9) as a candidate edge emanating 
								from t8 */
								for (Nt8 = t8->CandidateSet; (t9 = Nt8->To);
									Nt8++) {
										if (t9 == t8->Pred || t9 == t8->Suc
											|| t9 == t1 || (t8 == t2 && t9 == t3)
											|| (t8 == t3 && t9 == t2) || (t8 == t4
											&& t9 ==
											t5)
											|| (t8 == t5 && t9 == t4)
											|| ((G7 = G6 - Nt8->Cost) <= 0))
											continue;
										/* Choose t10 as one of t9's two neighbors 
										on the tour */
										for (X10 = 1; X10 <= 2; X10++) {
											if (X10 == 1) {
												t10 = 0;
												switch (Case8) {
												case 1:
													t10 = (BTW275 ?
														BETWEEN(t8, t9, t5)
														|| BETWEEN(t3, t9,
														t1) : BTW674
														? BETWEEN(t7, t9,
														t1) :
													BETWEEN(t7, t9,
														t5)) ? Pred(t9)
														: Suc(t9);
													Case10 = 22;
													break;
												case 2:
													t10 = (BTW376 ?
														BETWEEN(t8, t9, t4) :
													BTW571 ?
														BETWEEN(t7, t9, t1)
														|| BETWEEN(t3, t9,
														t6) :
													BETWEEN(t7, t9,
														t1)) ? Pred(t9)
														: Suc(t9);
													Case10 = 23;
													break;
												case 3:
													if (BTW574) {
														t10 = BETWEEN(t5, t9, t1) ?
															Pred(t9) : Suc(t9);
														Case10 = 24;
														break;
													}
													if (!BETWEEN(t5, t9, t4))
														break;
													t10 = Suc(t9);
													Case10 = 1;
													break;
												case 4:
													if (BTW671) {
														if (!BETWEEN(t2, t9, t5))
															break;
														t10 = Suc(t9);
														Case10 = 2;
														break;
													}
													t10 = BETWEEN(t6, t9, t4) ?
														Pred(t9) : Suc(t9);
													Case10 = 25;
													break;
												case 5:
													t10 = (BTW471 ?
														BETWEEN(t7, t9, t1) :
													BTW673 ?
														BETWEEN(t7, t9, t5) :
													BETWEEN(t4, t9, t1)
														|| BETWEEN(t7, t9,
														t5)) ?
														Pred(t9) : Suc(t9);
													Case10 = 26;
													break;
												case 6:
													t10 = (BTW471 ?
														BETWEEN(t7, t9, t3) :
													BTW573 ?
														BETWEEN(t8, t9, t6) :
													BETWEEN(t4, t9, t1)
														|| BETWEEN(t8, t9,
														t6)) ?
														Pred(t9) : Suc(t9);
													Case10 = 27;
													break;
												case 7:
													if (BTW273) {
														t10 = BETWEEN(t5, t9, t3) ?
															Pred(t9) : Suc(t9);
														Case10 = 28;
														break;
													}
													if (!BETWEEN(t2, t9, t3))
														break;
													t10 = Suc(t9);
													Case10 = 3;
													break;
												case 8:
													if (BTW273) {
														if (!BETWEEN(t4, t9, t5))
															break;
														Case10 = 4;
													} else {
														if (!BETWEEN(t2, t9, t3))
															break;
														Case10 = 5;
													}
													t10 = Suc(t9);
													break;
												case 9:
													if (BTW275) {
														if (!BETWEEN(t7, t9, t4))
															break;
														t10 = Suc(t9);
														Case10 = 6;
														break;
													}
													if (!BTW674) {
														if (!BETWEEN(t2, t9, t7))
															break;
														t10 = Suc(t9);
														Case10 = 7;
														break;
													}
													if (!BETWEEN(t6, t9, t7))
														break;
													t10 = Suc(t9);
													Case10 = 8;
													break;
												case 10:
													if (BTW376) {
														if (!BETWEEN(t7, t9, t6))
															break;
														t10 = Suc(t9);
														Case10 = 9;
														break;
													}
													if (BTW571) {
														if (!BETWEEN(t2, t9, t7))
															break;
														t10 = Suc(t9);
														Case10 = 10;
														break;
													}
													if (!BETWEEN(t3, t9, t6) &&
														!BETWEEN(t2, t9, t7))
														break;
													t10 = Suc(t9);
													Case10 = 11;
													break;
												case 11:
													if (BTW574) {
														t10 = BETWEEN(t3, t9, t1) ?
															Pred(t9) : Suc(t9);
														Case10 = 29;
														break;
													}
													if (!BETWEEN(t5, t9, t4))
														break;
													t10 = Suc(t9);
													Case10 = 12;
													break;
												case 12:
													t10 = BETWEEN(t3, t9, t1) ?
														Pred(t9) : Suc(t9);
													Case10 = 30;
													break;
												case 13:
													if (BTW471) {
														if (!BETWEEN(t2, t9, t7))
															break;
														t10 = Suc(t9);
														Case10 = 13;
														break;
													}
													if (BTW673) {
														if (!BETWEEN(t6, t9, t7))
															break;
														t10 = Suc(t9);
														Case10 = 14;
														break;
													}
													if (!BETWEEN(t6, t9, t3) &&
														!BETWEEN(t2, t9, t7))
														break;
													t10 = Suc(t9);
													Case10 = 15;
													break;
												case 14:
													if (BTW471) {
														if (!BETWEEN(t2, t9, t7))
															break;
														t10 = Suc(t9);
														Case10 = 16;
														break;
													}
													if (BTW573) {
														if (!BETWEEN(t7, t9, t3) &&
															!BETWEEN(t2, t9, t6))
															break;
														t10 = Suc(t9);
														Case10 = 17;
														break;
													}
													if (!BETWEEN(t7, t9, t6))
														break;
													t10 = Suc(t9);
													Case10 = 18;
													break;
												case 15:
													if (BTW273) {
														t10 = BETWEEN(t5, t9, t1) ?
															Pred(t9) : Suc(t9);
														Case10 = 31;
														break;
													}
													if (!BETWEEN(t2, t9, t3))
														break;
													t10 = Suc(t9);
													Case10 = 19;
													break;
												case 16:
													if (BTW273) {
														if (!BETWEEN(t4, t9, t5))
															break;
														Case10 = 20;
													} else {
														if (!BETWEEN(t2, t9, t3))
															break;
														Case10 = 21;
													}
													t10 = Suc(t9);
													break;
												}
												if (!t10)
													break;
											} else {
												if (Case10 >= 22)
													continue;
												Case10 += 31;
												t10 =
													t10 ==
													t9->Pred ? t9->Suc : t9->Pred;
											}
											if (t10 == t1 ||
												(t9 == t3 && t10 == t4) ||
												(t9 == t4 && t10 == t3) ||
												(t9 == t5 && t10 == t6) ||
												(t9 == t6 && t10 == t5))
												continue;
											G8 = G7 + graph->D(t9->Id, t10->Id);
											if ((Gain = G8 - graph->D(t10->Id, t1->Id)) > 0) {
												Make5OptMove(t1, t2, t3, t4, t5,
													t6, t7, t8, t9, t10,
													Case10);
												return 0;
											}
											if (G8  < t10->Cost)
												continue;
											if ((G8 > BestG8 ||
												(G8 == BestG8 && Near(T9, T10)
												&& !Near(t9, t10)))
												&& Swaps < MaxSwaps
												&& Excludable(t9, t10)
												&& !InInputTour(t9, t10)) {
													/* Ignore the move if the gain does
													not vary */
													if (G7 == G5 && G5 == G3 && G3 == G1)
														continue;
													T3 = t3;
													T4 = t4;
													T5 = t5;
													T6 = t6;
													T7 = t7;
													T8 = t8;
													T9 = t9;
													T10 = t10;
													BestCase10 = Case10;
													BestG8 = G8;
											}

										}
								}
							}
						}
					}
				}
			}
		}
		Gain = 0;
		if (T10) {
			/* Make the best 5-opt move */
			Make5OptMove(t1, t2, T3, T4, T5, T6, T7, T8, T9, T10, BestCase10);
			Exclude(t1, t2);
			Exclude(T3, T4);
			Exclude(T5, T6);
			Exclude(T7, T8);
			Exclude(T9, T10);
			G0 = BestG8;
		}
		return T10;
	}
	int BridgeGain(Node * s1, Node * s2, Node * s3, Node * s4,
		Node * s5, Node * s6, Node * s7, Node * s8, int Case6, long G)
	{
		Node *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *u2 = 0, *u3 = 0;
		Candidate *Nt2, *Nt4, *Nt6;
		long G0, G1, G2, G3, G4, G5, G6, Gain;
		int X4;

		switch (Case6) {
		case 3:
			if (2 * SegmentSize(s5, s4) <= Dimension) {
				u2 = s5;
				u3 = s4;
			} else {
				u2 = s3;
				u3 = s6;
			}
			break;
		case 4:
			if (2 * SegmentSize(s2, s5) <= Dimension) {
				u2 = s2;
				u3 = s5;
			} else {
				u2 = s6;
				u3 = s1;
			}
			break;
		case 0:
		case 7:
			if (2 * SegmentSize(s2, s3) <= Dimension) {
				u2 = s2;
				u3 = s3;
			} else {
				u2 = s4;
				u3 = s1;
			}
		}

		/* Choose t1 between u2 and u3 */
		for (t1 = u2; t1 != u3; t1 = t2) {
			/* Choose t2 as the successor of t1 */
			t2 = Suc(t1);
			if ((t1 == s1 && t2 == s2) ||
				(t1 == s2 && t2 == s1) ||
				(t1 == s3 && t2 == s4) ||
				(t1 == s4 && t2 == s3) ||
				(t1 == s5 && t2 == s6) ||
				(t1 == s6 && t2 == s5) ||
				(t1 == s7 && t2 == s8) ||
				(t1 == s8 && t2 == s7) )
				continue;
			G0 = G + graph->D(t1->Id, t2->Id);
			/* Choose (t2,t3) as a candidate edge emanating from t2. 
			t3 must not be between u2 and u3 */

			for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
				if (t3 == t2->Pred || t3 == t2->Suc || BETWEEN(u2, t3, u3))
					continue;
				G1 = G0 - Nt2->Cost;
				/* Choose t4 as one of t3's two neighbors on the tour */
				for (X4 = 1; X4 <= 2; X4++) {
					t4 = X4 == 1 ? Suc(t3) : Pred(t3);
					if (t4 == t2 ||
						(t3 == s1 && t4 == s2) ||
						(t3 == s2 && t4 == s1) ||
						(t3 == s3 && t4 == s4) ||
						(t3 == s4 && t4 == s3) ||
						(t3 == s5 && t4 == s6) ||
						(t3 == s6 && t4 == s5) ||
						(t3 == s7 && t4 == s8) ||
						(t3 == s8 && t4 == s7))
						continue;
					G2 = G1 + graph->D(t3->Id, t4->Id);
					/* Test if an improvement can be obtained */
					if ((Gain = G2 - graph->D(t4->Id, t1->Id)) > 0) {
						switch (Case6) {
						case 0:
							if (X4 == 1)
								Swap3(s1, s2, s4, t3, t4, t1, s1, s3, s2);
							else
								Swap2(t1, t2, t3, s1, s2, s3);
							return Gain;
						case 3:
							if ((X4 == 1) ==
								(!BETWEEN(s2, t1, s6) && !BETWEEN(s2, t3, s6)))
								Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
							else
								Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
								t4, t1);
							if (s8)
								Swap1(s7, s8, s1);
							return Gain;
						case 4:
							if ((X4 == 1) ==
								(!BETWEEN(s3, t1, s5) && !BETWEEN(s3, t3, s5)))
								Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
							else
								Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
								t4, t1);
							if (s8)
								Swap1(s7, s8, s1);
							return Gain;
						case 7:
							if ((X4 == 1) ==
								(!BETWEEN(s4, t1, s6) && !BETWEEN(s4, t3, s6)))
								Swap3(s5, s6, s1, t1, t2, t3, s3, s4, s5);
							else
								Swap4(s5, s6, s1, t1, t2, t4, s3, s4, s5, t2,
								t4, t1);
							if (s8)
								Swap1(s7, s8, s1);
							return Gain;
						}
					}
					/* If BridgeGain has been called with a nonfeasible 2-opt move,
					then try to find a 3-opt or 4-opt move which, when composed 
					with the 2-opt move, results in an improvement of the tour */
					if (Case6 != 0)
						continue;
					/* Choose (t4,t5) as a candidate edge emanating from t4 */
					for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
						if (t5 == t4->Pred || t5 == t4->Suc || t5 == t1
							|| t5 == t2)
							continue;
						/* Choose t6 as one of t5's two neighbors on the tour.
						Only one choice! */
						t6 = X4 == 1
							|| BETWEEN(u2, t5, u3) ? Pred(t5) : Suc(t5);
						if ((t5 == s1 && t6 == s2) || (t5 == s2 && t6 == s1)
							|| (t5 == s3 && t6 == s4) || (t5 == s4 && t6 == s3)
							)
							continue;
						G3 = G2 - Nt4->Cost;
						G4 = G3 + graph->D(t5->Id, t6->Id);
						if ((Gain = G4 - graph->D(t6->Id, t1->Id)) > 0) {
							if (X4 == 1)
								Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2, t5,
								t6, t1);
							else
								Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
							return Gain;
						}
						/* Choose (t7,t8) as a candidate edge emanating from t7.
						Only one choice! */
						for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
							if (t7 == t6->Pred || t7 == t6->Suc)
								continue;
							/* Choose t8 as one of t7's two neighbors on the tour.
							Only one choice! */
							if (X4 == 1)
								t8 = (BETWEEN(u2, t5, t1) ? BETWEEN(t5, t7, t1)
								: BETWEEN(t2, t5, u3) ? BETWEEN(u2, t7,
								t1)
								|| BETWEEN(t5, t7, u3) : BETWEEN(Suc(u3),
								t5,
								t3) ?
								BETWEEN(u2, t7, u3)
								|| BETWEEN(t5, t7, t3) : !BETWEEN(t4, t7,
								t6)) ?
								Pred(t7) : Suc(t7);
							else
								t8 = (BETWEEN(u2, t5, t1) ?
								!BETWEEN(u2, t7, t6)
								&& !BETWEEN(t2, t7, u3) : BETWEEN(t2, t5,
								u3) ?
								!BETWEEN(t2, t7, t6) : BETWEEN(Suc(u3),
								t5,
								t4) ?
								!BETWEEN(Suc(u3), t7, t5)
								&& !BETWEEN(t3, t7,
								Pred(u2)) : !BETWEEN(t3, t7,
								t5)) ?
								Pred(t7) : Suc(t7);
							if (t8 == t1
								|| (t7 == t1 && t8 == t2) || (t7 == t3
								&& t8 == t4)
								|| (t7 == t4 && t8 == t3) || (t7 == s1
								&& t8 == s2)
								|| (t7 == s2 && t8 == s1) || (t7 == s3
								&& t8 == s4)
								|| (t7 == s4 && t8 == s3))
								continue;

							G5 = G4 - Nt6->Cost;
							G6 = G5 + graph->D(t7->Id, t8->Id);
							/* Test if an improvement can be achieved */
							if ((Gain = G6 - graph->D(t8->Id, t1->Id)) > 0) {
								if (X4 == 1)
									Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2,
									t5, t6, t1);
								else
									Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
								Swap1(t7, t8, t1);
								return Gain;
							}
						}
					}
				}
			}
		}
		/* No improvement has been found */
		return 0;
	}
	Node *s1;
	int Gain2()
	{    
		Node *s2, *s3, *s4, *s5, *s6 = 0, *s7, *s8 = 0, *s1Stop;
		Candidate *Ns2, *Ns4, *Ns6;
		long G0, G1, G2, G3, G4, G5, G6, Gain, Gain6;
		int X2, X4, X6, X8, Case6 = 0, Case8 = 0;


		if (!s1 || s1->Subproblem != FirstNode->Subproblem)
			s1 = FirstNode;
		s1Stop = s1;
		for (X2 = 1; X2 <= 2; X2++) {
			Reversed = X2 == 1 ? OldReversed : (OldReversed ^= 1);
			do {
				s2 = Suc(s1);

				G0 = graph->D(s1->Id, s2->Id);
				/* Choose (s2,s3) as a candidate edge emanating from s2 */
				for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
					if (s3 == s2->Pred || s3 == s2->Suc)
						continue;
					G1 = G0 - Ns2->Cost;
					for (X4 = 1; X4 <= 2; X4++) {
						s4 = X4 == 1 ? Suc(s3) : Pred(s3);
						G2 = G1 + graph->D(s3->Id, s4->Id);

						/* Try any gainful nonfeasible 2-opt move 
						followed by a 2-, 3- or 4-opt move */
						if (X4 == 1 && s4 != s1 && 2 * SegmentSize(s2, s3) <= Dimension &&
							(G3 = G2 - graph->D(s4->Id, s1->Id)) > 0 &&
							(Gain = BridgeGain(s1, s2, s3, s4, 0, 0, 0, 0, 0,
							G3)) > 0)
							return Gain;
						if (X4 == 2 &&
							(Gain = G2 - graph->D(s4->Id, s1->Id)) > 0) {
								Swap1(s1, s2, s3);
								return Gain;
						}
						if (G2 - s4->Cost <= 0)
							continue;
						/* Try any gainful nonfeasible 3- or 4-opt move 
						folllowed by a 2-opt move */
						/* Choose (s4,s5) as a candidate edge emanating from s4 */
						for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
							if (s5 == s4->Pred || s5 == s4->Suc ||
								(G3 = G2 - Ns4->Cost) <= 0)
								continue;
							/* Choose s6 as one of s5's two neighbors on the tour */
							for (X6 = 1; X6 <= 2; X6++) {
								if (X4 == 2) {
									if (X6 == 1) {
										Case6 = 1 + !BETWEEN(s2, s5, s4);
										s6 = Case6 == 1 ? Suc(s5) : Pred(s5);
									} else {
										s6 = s6 ==
											s5->Pred ? s5->Suc : s5->Pred;
										if (s5 == s1 || s6 == s1)
											continue;
										Case6 += 2;
									}
								} else if (BETWEEN(s2, s5, s3)) {
									Case6 = 4 + X6;
									s6 = X6 == 1 ? Suc(s5) : Pred(s5);
									if (s6 == s1)
										continue;
								} else {
									if (X6 == 2)
										break;
									Case6 = 7;
									s6 = Pred(s5);
								}
								G4 = G3 + graph->D(s5->Id, s6->Id);
								Gain6 = 0;
								if ((Gain6 = G4 - graph->D(s6->Id, s1->Id)) > 0) {
									if (Case6 <= 2 || Case6 == 5 || Case6 == 6) {
										Make3OptMove(s1, s2, s3, s4, s5, s6,
											Case6);
										return Gain6;
									}
									if ((Gain =
										BridgeGain(s1, s2, s3, s4, s5, s6, 0,
										0, Case6, Gain6)) > 0)
										return Gain;
								}

							}
						}
					}
				}
			}
			while ((s1 = s2) != s1Stop);
		}
		return 0;
	}
	int Gain23()
	{    
		Node *s2, *s3, *s4, *s5, *s6 = 0, *s7, *s8 = 0, *s1Stop;
		Candidate *Ns2, *Ns4, *Ns6;
		long G0, G1, G2, G3, G4, G5, G6, Gain, Gain6;
		int X2, X4, X6, X8, Case6 = 0, Case8 = 0;


		if (!s1 || s1->Subproblem != FirstNode->Subproblem)
			s1 = FirstNode;
		s1Stop = s1;
		for (X2 = 1; X2 <= 2; X2++) {
			Reversed = X2 == 1 ? OldReversed : (OldReversed ^= 1);
			do {
				s2 = Suc(s1);

				G0 = graph->D(s1->Id, s2->Id);
				/* Choose (s2,s3) as a candidate edge emanating from s2 */
				for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
					if (s3 == s2->Pred || s3 == s2->Suc)
						continue;
					G1 = G0 - Ns2->Cost;
					for (X4 = 1; X4 <= 2; X4++) {
						s4 = X4 == 1 ? Suc(s3) : Pred(s3);
						G2 = G1 + graph->D(s3->Id, s4->Id);

						/* Try any gainful nonfeasible 2-opt move 
						followed by a 2-, 3- or 4-opt move */
						if (X4 == 1 && s4 != s1 && 2 * SegmentSize(s2, s3) <= Dimension &&
							(G3 = G2 - graph->D(s4->Id, s1->Id)) > 0 &&
							(Gain = BridgeGain(s1, s2, s3, s4, 0, 0, 0, 0, 0,
							G3)) > 0)
							return Gain;
						if (X4 == 2 &&
							(Gain = G2 - graph->D(s4->Id, s1->Id)) > 0) {
								Swap1(s1, s2, s3);
								return Gain;
						}
						if (G2 - s4->Cost <= 0)
							continue;
						/* Try any gainful nonfeasible 3- or 4-opt move 
						folllowed by a 2-opt move */
						/* Choose (s4,s5) as a candidate edge emanating from s4 */
						for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
							if (s5 == s4->Pred || s5 == s4->Suc ||
								(G3 = G2 - Ns4->Cost) <= 0)
								continue;
							/* Choose s6 as one of s5's two neighbors on the tour */
							for (X6 = 1; X6 <= 2; X6++) {
								if (X4 == 2) {
									if (X6 == 1) {
										Case6 = 1 + !BETWEEN(s2, s5, s4);
										s6 = Case6 == 1 ? Suc(s5) : Pred(s5);
									} else {
										s6 = s6 ==
											s5->Pred ? s5->Suc : s5->Pred;
										if (s5 == s1 || s6 == s1)
											continue;
										Case6 += 2;
									}
								} else if (BETWEEN(s2, s5, s3)) {
									Case6 = 4 + X6;
									s6 = X6 == 1 ? Suc(s5) : Pred(s5);
									if (s6 == s1)
										continue;
								} else {
									if (X6 == 2)
										break;
									Case6 = 7;
									s6 = Pred(s5);
								}
								G4 = G3 + graph->D(s5->Id, s6->Id);
								Gain6 = 0;
								if ((Gain6 = G4 - graph->D(s6->Id, s1->Id)) > 0) {
									if (Case6 <= 2 || Case6 == 5 || Case6 == 6) {
										Make3OptMove(s1, s2, s3, s4, s5, s6,
											Case6);
										return Gain6;
									}
									if ((Gain =
										BridgeGain(s1, s2, s3, s4, s5, s6, 0,
										0, Case6, Gain6)) > 0)
										return Gain;
								}
								/* Choose (s6,s7) as a candidate edge
								emanating from s6 */
								for (Ns6 = s6->CandidateSet; (s7 = Ns6->To);
									Ns6++) {
										if (s7 == s6->Pred || s7 == s6->Suc
											|| (s6 == s2 && s7 == s3) || (s6 == s3
											&& s7 ==
											s2)
											|| (G5 = G4 - Ns6->Cost) <= 0)
											continue;
										/* Choose s6 as one of s5's two neighbors 
										on the tour */
										for (X8 = 1; X8 <= 2; X8++) {
											if (X8 == 1) {
												Case8 = Case6;
												switch (Case6) {
												case 1:
													s8 = BETWEEN(s2, s7,
														s5) ? Suc(s7) :
														Pred(s7);
													break;
												case 2:
													s8 = BETWEEN(s3, s7,
														s6) ? Suc(s7) :
														Pred(s7);
													break;
												case 3:
													if (BETWEEN(s5, s7, s4))
														s8 = Suc(s7);
													else {
														s8 = BETWEEN(s3, s7,
															s1) ? Pred(s7)
															: Suc(s7);
														Case8 = 17;
													}
													break;
												case 4:
													if (BETWEEN(s2, s7, s5))
														s8 = BETWEEN(s2, s7,
														s4) ? Suc(s7)
														: Pred(s7);
													else {
														s8 = Pred(s7);
														Case8 = 18;
													}
													break;
												case 5:
													s8 = Pred(s7);
													break;
												case 6:
													s8 = BETWEEN(s2, s7,
														s3) ? Suc(s7) :
														Pred(s7);
													break;
												case 7:
													if (BETWEEN(s2, s7, s3))
														s8 = Suc(s7);
													else {
														s8 = BETWEEN(s5, s7,
															s1) ? Pred(s7)
															: Suc(s7);
														Case8 = 19;
													}
												}
											} else {
												if (Case8 >= 17 ||
													(Case6 != 3 && Case6 != 4
													&& Case6 != 7))
													break;
												s8 = s8 ==
													s7->Pred ? s7->Suc : s7->Pred;
												Case8 += 8;
											}
											if (s8 == s1 ||
												(s7 == s1 && s8 == s2) ||
												(s7 == s3 && s8 == s4) ||
												(s7 == s4 && s8 == s3))
												continue;
											G6 = G5 + graph->D(s7->Id, s8->Id);
											if ((Gain = G6 - graph->D(s8->Id, s1->Id)) > 0) {
												if (Case8 <= 15) {
													Make4OptMove(s1, s2, s3, s4,
														s5, s6, s7, s8,
														Case8);
													return Gain;
												}
												if (Gain > Gain6 &&
													(Gain =
													BridgeGain(s1, s2, s3, s4, s5,
													s6, s7, s8, Case6,
													Gain)) > 0)
													return Gain;
											}
										}
								}
							}
						}
					}
				}
			}
			while ((s1 = s2) != s1Stop);
		}
		return 0;
	}
	int Number_Of_Candidates;	
	int EdgesInFragments;
	long Cost;
	int mark;
	static int compareX(const void *Na, const void *Nb)
	{
		int result;
		double x1 = (*(Node **) Na)->X;
		double y1 = (*(Node **) Na)->Y;
		double x2 = (*(Node **) Nb)->X;
		double y2 = (*(Node **) Nb)->Y;
		result = x1 < x2 ? -1 : x1 > x2 ? 1 : y1 < y2 ? -1 : y1 > y2 ? 1 : 0;
		return result;
	}
	int compareCost(const void *Na, const void *Nb)
	{
		return (*(Node **) Na)->Cost - (*(Node **) Nb)->Cost;
	}
	Node *NearestNeighbor(Node * From)
	{    
		Candidate *NN;
		Node *To, *N, *First = 0, *Last = 0, *Nearest = 0;
		int MaxLevel = Dimension, Min = INT_MAX, d;

 		if (From->Degree == 2)
 			return 0; 
		for (NN = From->CandidateSet; (To = NN->To); NN++) {
			if (FixedOrCommon(From, To) && MayBeAddedToFragments(From, To)) {
				From->Cost = NN->Cost;
				return To;
			}
		}
		From->Level = 0;
		if (++mark == 0)
			mark = 1;
		From->mark = mark;
		/* Insert From into an empty queue */
		First = Last = From;
		From->OldSuc = 0;

		while ((N = First) && N->Level < MaxLevel) {
			/* Remove the first node from the queue */
			if (N == Last)
				First = Last = 0;
			else
				First = N->OldSuc;
			for (NN = N->CandidateSet; (To = NN->To); NN++) {
				if (To->mark != mark) {
					To->mark = mark;
					To->Level = N->Level + 1;
					if (MayBeAddedToFragments(From, To) &&
						(N == From ? (d = NN->Cost) < Min :
						(d = this->graph->D(From->Id, To->Id)) < Min)) {
							Min = From->Cost = d;
							/* Randomization */
							if (!Nearest && rand() % 3 != 0)
								return To;
							Nearest = To;
							MaxLevel = To->Level;
					} else if (To->Level < MaxLevel) {
						/* Insert To as the last element of the queue */
						if (Last)
							Last->OldSuc = To;
						else
							First = To;
						Last = To;
						Last->OldSuc = 0;
					}
				}
			}
		}
		return Nearest;
	}
	Node *NearestInList(Node * From, Node * First)
	{
		Node *To, *Nearest = 0;
		int Min = INT_MAX, d;

		To = First;
		do {
			if (MayBeAddedToFragments(From, To) &&
				(d = graph->D(From->Id, To->Id)) < Min) {
					Min = From->Cost = d;
					Nearest = To;
			}
		}
		while ((To = To->OldSuc) != First);
		return Nearest;
	}
	int MayBeAddedToFragments(Node * From, Node * To)
	{
		return From != To && From->Degree != 2 && To->Degree != 2 &&
			(From->Tail != To || EdgesInFragments == Dimension - 1);
	}
	void AddEdgeToFragments(Node * From, Node * To)
	{
		Node *Temp;

		if (!From->Pred)
			From->Pred = To;
		else
			From->Suc = To;
		if (!To->Pred)
			To->Pred = From;
		else
			To->Suc = From;
		From->Degree++;
		To->Degree++;
		Temp = From->Tail;
		Temp->Tail = To->Tail;
		To->Tail->Tail = Temp;
		EdgesInFragments++;
		Cost += From->Cost;// - From->Pi - To->Pi;
	}
	void RemoveFromList(Node * N, Node ** First)
	{
		if (*First == N)
			*First = N->OldSuc;
		N->OldPred->OldSuc = N->OldSuc;
		N->OldSuc->OldPred = N->OldPred;
	}
	void InitateNodes()
	{
		for (int t = 0; t < Dimension; t++)
		{
			//Activate(NodeSet[t]);
			NodeSet[t]->OldPredExcluded = NodeSet[t]->OldSucExcluded = false;
			NodeSet[t]->OldPred = NodeSet[t]->Pred;
			NodeSet[t]->OldSuc = NodeSet[t]->Suc;
		}
	}
	void InitiateNodesSegments()
	{
		Node* FirstNode = NodeSet[0];
		Node* t1;
		Node* t2;
		int i;
		Segment* S;
		SSegment* SS;
		Reversed = 0;
		S = FirstSegment;
		i = 0;
		do
		{
			S->Size = 0;
			S->Rank = ++i;
			S->Reversed = false;
			S->First = S->Last = 0;
		}
		while ((S = S->Suc) != FirstSegment);
		SS = FirstSSegment;
		i = 0;
		do
		{
			SS->Size = 0;
			SS->Rank = ++i;
			SS->Reversed = false;
			SS->First = SS->Last = 0;
		}
		while ((SS = SS->Suc) != FirstSSegment);
		i = 0;
		t1 = FirstNode;
		do
		{
			t2 = t1->OldSuc = t1->Suc;
			t1->OldPred = t1->Pred;
			t1->Rank = ++i;

			t1->Parent = S;
			S->Size++;
			if (S->Size == 1)
				S->First = t1;
			S->Last = t1;
			if (SS->Size == 0)
				SS->First = S;
			S->Parent = SS;
			SS->Last = S;
			if (S->Size == GroupSize)
			{
				S = S->Suc;
				SS->Size++;
				if (SS->Size == SGroupSize)
					SS = SS->Suc;
			}
			t1->OldPredExcluded = t1->OldSucExcluded = false;

		}
		while ((t1 = t1->Suc) != FirstNode);
		if (S->Size < GroupSize)
			SS->Size++;
	}
	void AllocateSegments()
	{
		Segment *S = 0, *SPrev;
		SSegment *SS = 0, *SSPrev;
		int i;

		GroupSize = (int) sqrt((double) Dimension);
		SegmentSet = new Segment[(Dimension / GroupSize) + 1];

		Groups = 0;int j = 0;
		for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S)
		{
			S = &SegmentSet[j];j++;
			S->Rank = ++Groups;
			if (!SPrev)
				FirstSegment = S;
			else
				SLink(SPrev, S);
		}
		SLink(S, FirstSegment);
		SGroupSize = Dimension;
		SGroups = 0;
		for (i = Groups, SSPrev = 0; i > 0; i -= SGroupSize, SSPrev = SS)
		{
			SS = (SSegment *) malloc(sizeof(SSegment));
			SS->Rank = ++SGroups;
			if (!SSPrev)
				FirstSSegment = SS;
			else
				SLink(SSPrev, SS);
		}
		SLink(SS, FirstSSegment);
		InitateNodes();
		InitiateNodesSegments();
		//free(IsCopied);IsCopied = 0;
	}
	void SetTour(Tour *tour)
	{
		this->s1 = 0;
		int r_node = -1, node = 0, l_node = 0;
		r_node = l_node = -1; 
		for (int i = 0 ; i < Dimension ; i++ )
		{
			NodeSet[i]->Id = i;
			NodeSet[i]->Cost = INT_MAX;
			NodeSet[i]->IsActivated = 0;
			NodeSet[i]->Next = NodeSet[i]->OldPred = NodeSet[i]->OldSuc = NodeSet[i]->Pred = NodeSet[i]->Suc = 0;
			NodeSet[i]->Parent = 0;
		}
		FirstActive = LastActive = 0;
		this->Reversed = 0;	
		for(int i = 0 ; i <= Dimension; i++)
		{
			if (r_node == -1)
			{
				NodeSet[node]->Suc = NodeSet[node];
				NodeSet[node]->Pred = NodeSet[node];
				l_node = r_node = node;
				continue;				
			}
			NodeSet[node]->Pred = NodeSet[r_node];
			NodeSet[node]->Suc = NodeSet[l_node];
			NodeSet[l_node]->Pred = NodeSet[node];
			NodeSet[r_node]->Suc = NodeSet[node];
			r_node = node;
			node = tour->Right(node);
		}
		AllocateSegments();        
	}
	void ReleaseTour()
	{
		delete FirstSSegment;
		delete SegmentSet;
	}
public:
	~Heuristics()
	{
		for(int i = 0; i < this->Dimension; i++) delete NodeSet[i]->CandidateSet;
		for(int i = 0; i < this->Dimension; i++) delete NodeSet[i];
		delete NodeSet;
		delete SwapStack;
	}
	Heuristics(Graph* graph, int NumberOfCandidates)
	{
		this->cand_num = NumberOfCandidates;
		this->Number_Of_Candidates = NumberOfCandidates;
		this->graph = graph;
		MaxSwaps = Dimension = graph->Dimension();
		Swaps = 0;
		this->Reversed = this->OldReversed = 0;
		AllocateStructures();
		CreateNodes();
		CreateNearestNeighborCandidateSet(NumberOfCandidates);
	}	
	Tour* Q_Boruvka()
	{
		Node *Prev = 0, *N = 0;
		for (int i = 0; i < Dimension; i++, Prev = N) 
        {
			N = NodeSet[i];
			if (i == 0)
				FirstNode = N;
			else
				Link(Prev, N);
			N->Id = i; 
		}
		Link(N, FirstNode);

		Node *From, *To, *First, *Last = 0, **Perm;
		int Count, i;
		Heap *heap = new Heap(Dimension); 
		Tour *tour = new Tour(graph);
		mark = 0;
		Cost = 0;
		EdgesInFragments = 0;

		From = FirstNode;
		do {
			From->Degree = 0;
			From->Tail = From;
			From->mark = 0;
			From->Next = From->Suc;
			From->Pred = From->Suc = 0;
		}
		while ((From = From->Next) != FirstNode);
		Count = 0;
		for (;;) {          
			if ((To = NearestNeighbor(From)) && FixedOrCommon(From, To))
				AddEdgeToFragments(From, To);
			else {
				if ((From->Nearest = To)) {                
					Count++;
				}
				if ((From = From->Next) == FirstNode)
					break;
			}
		}
        
		Perm = (Node **) malloc(Count * sizeof(Node *));
		for (From = FirstNode, i = 0; i < Count; From = From->Next)
		{
			if (From->Nearest)
				Perm[i++] = From;
		}

		qsort(Perm, Count, sizeof(Node *), compareX);
		for (i = 0; i < Count; i++) {
			From = Perm[i];
			if ((To = NearestNeighbor(From))) {
				AddEdgeToFragments(From, To);
				i--;
			}
		}

		free(Perm);

		if (EdgesInFragments < Dimension) {
			/* Create a list of the degree-0 and degree-1 nodes */
			First = 0;
			From = FirstNode;
			do {
				if (From->Degree != 2) {
					From->OldSuc = First;
					if (First)
						First->OldPred = From;
					else
						Last = From;
					First = From;
				}
			}
			while ((From = From->Next) != FirstNode);
			First->OldPred = Last;
			Last->OldSuc = First;
			/* Initialize the heap */
			From = First;
			do {
				if ((From->Nearest = NearestInList(From, First))) {
					From->Rank = From->Cost;
					heap->HeapLazyInsert(From);
				}
			}
			while ((From = From->OldSuc) != First);
			heap->Heapify();
			/* Find the remaining fragments */
			while ((From = heap->HeapDeleteMin())) {
				To = From->Nearest;
				if (MayBeAddedToFragments(From, To)) {
					AddEdgeToFragments(From, To);
					if (From->Degree == 2)
						RemoveFromList(From, &First);
					if (To->Degree == 2)
						RemoveFromList(To, &First);
				}
				if (From->Degree != 2
					&& (From->Nearest = NearestInList(From, First))) {
						From->Rank = From->Cost;
						heap->HeapInsert(From);
				}
			}
		}
		/* Orient Pred and Suc so that the list of nodes represents a tour */
		To = FirstNode;
		From = To->Pred;
		do {
			tour->Add(To->Id);
			if (To->Suc == From) {
				To->Suc = To->Pred;
				To->Pred = From;
			}
			From = To;
		}
		while ((To = From->Suc) != FirstNode);
		To->Pred = From;
		delete heap;
		return tour;
	}
	void TwoOpt(Tour* tour, int *ActiveNodes, int NumberOfActiveNodes)
	{

		this->SetTour(tour);
		for (int t = 0; t < NumberOfActiveNodes; t++)
		{
			Activate(NodeSet[ActiveNodes[t]]);
		}

		long Gain = 0;
		long G0;
		int X2;

		Node *t1 = 0;
		Node *SUCt1 = 0;
		Node *t2 = 0;

		Candidate *Nt;

		for (int i = 0; i < Dimension; i++)
		{
			Node *t = NodeSet[i];
			for (Nt = t->CandidateSet; Nt->To; Nt++)
			{
				t2 = Nt->To;
				int cost_candidate = graph->D(t->Id, t2->Id);
				if (t2 != Pred(t) && t2 != Suc(t) && cost_candidate < t->Cost)
					t->Cost = cost_candidate;
			}
		}
		while(1)
		{
			while (t1 = RemoveFirstActive())
			{
				SUCt1 = Suc(t1);
				for (X2 = 1; X2 <= 2; X2++)
				{
					t2 = X2 == 1 ? Pred(t1) : SUCt1;
					if(t2 == t1->BestPred || t2 == t1->BestSuc) continue;
					G0 = graph->D(t1->Id, t2->Id);
					do
					{
						t2 = Best2OptMove(t1, t2, G0, Gain);
					}
					while (t2);
					if (Gain > 0)
					{
						StoreTour();
						Activate(t1);
						break;
					}
					RestoreTour();
				}
			}
			if ((Gain = Gain2()) > 0) {StoreTour();}
			else break;
		}

		tour->Reset();
		int node = 0;
		for(int i = 0 ; i < Dimension; i++)
		{
			tour->Add(node);
			Node *p = Suc(this->NodeSet[node]);
			node = p->Id;
		}
		ReleaseTour();
	}		
	void TwoOpt(Tour* tour)
	{

		this->SetTour(tour);
		for (int t = 0; t < this->Dimension; t++)
		{
			Activate(NodeSet[t]);
		}

		long Gain = 0;
		long G0;
		int X2;

		Node *t1 = 0;
		Node *SUCt1 = 0;
		Node *t2 = 0;

		Candidate *Nt;

		for (int i = 0; i < Dimension; i++)
		{
			Node *t = NodeSet[i];
			for (Nt = t->CandidateSet; Nt->To; Nt++)
			{
				t2 = Nt->To;
				int cost_candidate = graph->D(t->Id, t2->Id);
				if (t2 != Pred(t) && t2 != Suc(t) && cost_candidate < t->Cost)
					t->Cost = cost_candidate;
			}
		}
		while(1)
		{
			while (t1 = RemoveFirstActive())
			{
				SUCt1 = Suc(t1);
				for (X2 = 1; X2 <= 2; X2++)
				{
					t2 = X2 == 1 ? Pred(t1) : SUCt1;
					if(t2 == t1->BestPred || t2 == t1->BestSuc) continue;
					G0 = graph->D(t1->Id, t2->Id);
					do
					{
						t2 = Best2OptMove(t1, t2, G0, Gain);
					}
					while (t2);
					if (Gain > 0)
					{
						StoreTour();
						Activate(t1);
						break;
					}
					RestoreTour();
				}
			}
			if ((Gain = Gain2()) > 0) {StoreTour();}
			else break;
		}

		tour->Reset();
		int node = 0;
		for(int i = 0 ; i < Dimension; i++)
		{
			tour->Add(node);
			Node *p = Suc(this->NodeSet[node]);
			node = p->Id;
		}
		ReleaseTour();
	}		
	void LinKernighan(Tour* tour)
	{

		this->SetTour(tour);
		for (int t = 0; t < this->Dimension; t++)
		{
			Activate(NodeSet[t]);
		}

		long Gain = 0;
		long G0;
		int X2;

		Node *t1 = 0;
		Node *SUCt1 = 0;
		Node *t2 = 0;

		Candidate *Nt;

		for (int i = 0; i < Dimension; i++)
		{
			Node *t = NodeSet[i];
			for (Nt = t->CandidateSet; Nt->To; Nt++)
			{
				t2 = Nt->To;
				int cost_candidate = graph->D(t->Id, t2->Id);
				if (t2 != Pred(t) && t2 != Suc(t) && cost_candidate < t->Cost)
					t->Cost = cost_candidate;
			}
		}
		while(1)
		{
			while (t1 = RemoveFirstActive())
			{
				SUCt1 = Suc(t1);
				for (X2 = 1; X2 <= 2; X2++)
				{
					t2 = X2 == 1 ? Pred(t1) : SUCt1;
					if(t2 == t1->BestPred || t2 == t1->BestSuc) continue;
					G0 = graph->D(t1->Id, t2->Id);
					do
					{
						t2 = Best5OptMove(t1, t2, G0, Gain);
					}
					while (t2);
					if (Gain > 0)
					{
						StoreTour();
						Activate(t1);
						break;
					}
					RestoreTour();
				}
			}		
			if ((Gain = Gain2()) > 0) {StoreTour();}
			else break;
		}

		tour->Reset();
		int node = 0;
		for(int i = 0 ; i < Dimension; i++)
		{
			tour->Add(node);
			Node *p = Suc(this->NodeSet[node]);
			node = p->Id;
		}
		ReleaseTour();
	}		
	void LinKernighan(Tour* tour, int *ActiveNodes, int NumberOfActiveNodes)
	{
		this->SetTour(tour);
		for (int t = 0; t < NumberOfActiveNodes; t++)
		{
			Activate(NodeSet[ActiveNodes[t]]);
		}

		long Gain = 0;
		long G0;
		int X2;

		Node *t1 = 0;
		Node *SUCt1 = 0;
		Node *t2 = 0;

		Candidate *Nt;

		for (int i = 0; i < Dimension; i++)
		{
			Node *t = NodeSet[i];
			for (Nt = t->CandidateSet; Nt->To; Nt++)
			{
				t2 = Nt->To;
				int cost_candidate = graph->D(t->Id, t2->Id);
				if (t2 != Pred(t) && t2 != Suc(t) && cost_candidate < t->Cost)
					t->Cost = cost_candidate;
			}
		}
		while(1)
		{
			while (t1 = RemoveFirstActive())
			{
				SUCt1 = Suc(t1);
				for (X2 = 1; X2 <= 2; X2++)
				{
					t2 = X2 == 1 ? Pred(t1) : SUCt1;
					if(t2 == t1->BestPred || t2 == t1->BestSuc) continue;
					G0 = graph->D(t1->Id, t2->Id);
					do
					{
						t2 = Best5OptMove(t1, t2, G0, Gain);
					}
					while (t2);
					if (Gain > 0)
					{
						StoreTour();
						Activate(t1);
						break;
					}
					RestoreTour();
				}
			}		
			if ((Gain = Gain2()) > 0) {StoreTour();}
			else break;
		}

		tour->Reset();
		int node = 0;
		for(int i = 0 ; i < Dimension; i++)
		{
			tour->Add(node);
			Node *p = Suc(this->NodeSet[node]);
			node = p->Id;
		}
		ReleaseTour();
	}
	void DoubleBridge(Tour* tour, int *WillBeActive)
	{
		//return 8 nodes and these nodes will be activated in LK.
		int Dimension = this->graph->Dimension();
		int* input = new int [Dimension];
		int node = 0;
		for(int i = 0; i < Dimension; i++)
		{
			input[i] = node;
			node = tour->Right(node);
		}
		int A_F = 0;
		int A_E = Random(1, Dimension - 6 - 1);
		int B_F = A_E + 1;
		int B_E = Random(B_F + 1, Dimension - 4 - 1);
		int C_F = B_E + 1;
		int C_E = Random(C_F + 1, Dimension - 2 - 1);
		int D_F = C_E + 1;
		int D_E = Dimension - 1;
		WillBeActive[0] = A_F;
		WillBeActive[1] = A_E;
		WillBeActive[2] = B_F;
		WillBeActive[3] = B_E;
		WillBeActive[4] = C_F;
		WillBeActive[5] = C_E;
		WillBeActive[6] = D_F;
		WillBeActive[7] = D_E;
		int A_LEN = A_E - A_F + 1;
		int B_LEN = B_E - B_F + 1;
		int C_LEN = C_E - C_F + 1;
		int D_LEN = D_E - D_F + 1;
		int* output = (int *) malloc((Dimension) * sizeof(int));

		for(int i = A_F; i <= A_E; i++)
		{
			output[i] = input[i];
		}
		for(int i = 0; i < D_LEN; i++)
		{
			output[A_LEN + i] = input[i + D_F];
		}
		for(int i = 0; i < C_LEN; i++)
		{
			output[i + A_LEN + D_LEN] = input[i + C_F];
		}
		for(int i = 0; i < B_LEN; i++)
		{
			output[i + A_LEN + D_LEN + C_LEN] = input[i + B_F];
		}
		tour->Reset();
		for(int i = 0; i < Dimension; i++)
		{
			tour->Add(output[i]);
		}
		delete output;
		delete input;
	}//DoubleBridge

	void SetCandidates(int node,int*candidates,int NumberOfCandidates)
	{
		if( NumberOfCandidates >= this->cand_num ) 
			NumberOfCandidates = this->cand_num;
		free(NodeSet[node]->CandidateSet);
		NodeSet[node]->CandidateSet = (Candidate *) malloc((NumberOfCandidates + 1) * sizeof(Candidate));
		NodeSet[node]->CandidateSet[NumberOfCandidates].To = 0;
		for(int i = 0; i < NumberOfCandidates; i++)
		{
			NodeSet[node]->CandidateSet[i].To = NodeSet[candidates[i]];
			NodeSet[node]->CandidateSet[i].Cost = (this->graph->D)(NodeSet[node]->Id,NodeSet[candidates[i]]->Id);
		}
	}
	void SetCandidatesCount(int node, int count)
	{
		NodeSet[node]->CandidateSet[count].To = 0;
	}
	int SetCandidates(int node, int candidate, int index)
	{
		Candidate *Ns2;
		int count = 0;
		int location = -1;
		for (Ns2 = NodeSet[node]->CandidateSet; (Ns2->To); Ns2++) 
		{
			if(Ns2->To->Id == candidate)
			{
				location = count;
			}
			count++;
		}
		if(location == -1)
		{
			for(int i = count - 2; i >= index; i--)
			{
				NodeSet[node]->CandidateSet[i + 1].To = NodeSet[node]->CandidateSet[i].To;
				NodeSet[node]->CandidateSet[i + 1].Cost = NodeSet[node]->CandidateSet[i].Cost;				
			}
			NodeSet[node]->CandidateSet[index].To = NodeSet[candidate];
			NodeSet[node]->CandidateSet[index].Cost = (this->graph->D)(node, candidate);

			return 1;
		}
	}
	int SetCandidates_3(int node, int candidate, int index)
	{
		Candidate *Ns2;
		for (Ns2 = NodeSet[node]->CandidateSet; (Ns2->To); Ns2++) 
		{
			if(Ns2->To->Id == candidate) return 0;
		}
		NodeSet[node]->CandidateSet[index].To = NodeSet[candidate];
		NodeSet[node]->CandidateSet[index].Cost = (this->graph->D)(node, candidate);
		return 1;
	}
	int SetCandidates_2(int node, int candidate, int index)
	{
		Candidate *Ns2;
		int count = 0;
		int location = -1;
		for (Ns2 = NodeSet[node]->CandidateSet; (Ns2->To); Ns2++) 
		{
			if(Ns2->To->Id == candidate)
			{
				/**/Node *temp = NodeSet[node]->CandidateSet[index].To;
				int c_temp = NodeSet[node]->CandidateSet[index].Cost;
				NodeSet[node]->CandidateSet[index].To = NodeSet[node]->CandidateSet[count].To;
				NodeSet[node]->CandidateSet[index].Cost = NodeSet[node]->CandidateSet[count].Cost;
				NodeSet[node]->CandidateSet[count].To = temp;
				NodeSet[node]->CandidateSet[count].Cost = c_temp;
				return 0;
			}
			count++;
		}
		for(int i = count - 2; i >= index; i--)
		{
			NodeSet[node]->CandidateSet[i + 1].To = NodeSet[node]->CandidateSet[i].To;
			NodeSet[node]->CandidateSet[i + 1].Cost = NodeSet[node]->CandidateSet[i].Cost;	
		}
		NodeSet[node]->CandidateSet[index].To = NodeSet[candidate];
		NodeSet[node]->CandidateSet[index].Cost = (this->graph->D)(node, candidate);
		return 1;
	}
	int DoesInclude(int node, int other)
	{
		for (Candidate *Ns2 = NodeSet[node]->CandidateSet; (Ns2->To); Ns2++) 
		{
			if(Ns2->To->Id == other)
			{
				return 1;
			}
		}
		return 0;
	}
	void SetBestTour(Tour *best_tour)
	{
		for (int t = 0; t < Dimension; t++)
		{
			NodeSet[t]->BestPred = NodeSet[best_tour->Left(t)];
			NodeSet[t]->BestSuc = NodeSet[best_tour->Right(t)];
		}
	}	
	void ResetBestTour()
	{
		for (int t = 0; t < Dimension; t++)
		{
			NodeSet[t]->BestPred = 0;
			NodeSet[t]->BestSuc = 0;
		}
	}

	int GetCandidate(int node, int index)
	{
		if(NodeSet[node]->CandidateSet[index].To == NULL)
		{
			return -1;
		}
 	    return NodeSet[node]->CandidateSet[index].To->Id;
	}

	int GetNumberOfCandidate()
	{
		return this->cand_num;
	}
};
#endif
