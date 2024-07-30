#ifndef _TOUR_H
#define _TOUR_H

#include "Graph.h"
#include "GainType.h"
#include "Time_Random.h"
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stdio.h"
#include "stdlib.h"
#include <malloc.h>
#include <time.h>
#include <limits.h>
#include <iostream>
using namespace std;

class Tour
{
	struct Node
	{
		short IsExist;
		Node *left;
		Node *right;
		int Id;
	};
	Node *nodes;
	Node *first, *last; // first-> . . . ->last->first
	Graph *graph;
	int Dimension;
	GainType cost;
	int count;
	void CalculateCost()
	{
		int node = 0;
		for(int i = 0; i < Dimension; i++)
		{
			int right = this->Right(node);
			cost = cost + graph->D(node, right);
			node = right;
		}
	}
public:
	Tour(Graph *graph)
	{
		Dimension = graph->Dimension();
		this->graph = graph;
		nodes = new Node[Dimension];
		for(int i = 0; i < Dimension; i++) this->nodes[i].IsExist = 0;
		for(int i = 0; i < Dimension; i++) this->nodes[i].Id = i;
		first = last = 0;
		cost = 0;
		count = 0;
	}
	void Reset()
	{
		for(int i = 0; i < Dimension; i++) this->nodes[i].IsExist = 0;
		first = last = 0;
		cost = 0;
		count = 0;
	}
	~Tour(){ delete nodes;}
	short Add(int node)
	{
		if(node >= Dimension || node < 0 || nodes[node].IsExist) return 0;
		count++;
		nodes[node].IsExist = 1;
		if(last == 0)
		{
			first = last = &nodes[node];
			nodes[node].left = nodes[node].right = &nodes[node];			
			return 1;
		}
		last->right = &nodes[node];
		nodes[node].left = last;
		nodes[node].right = first;
		first->left = &nodes[node];
		last = &nodes[node];
		if(count == Dimension) this->CalculateCost();
		return 1;
	}
	short AddToLeft(int node)
	{
		if(node >= Dimension || node < 0 || nodes[node].IsExist) return 0;
		count++;
		nodes[node].IsExist = 1;
		if(first == 0)
		{
			first = last = &nodes[node];
			nodes[node].left = nodes[node].right = &nodes[node];			
			return 1;
		}
		last->right = &nodes[node];
		nodes[node].left = last;
		nodes[node].right = first;
		first->left = &nodes[node];
		first = &nodes[node];
		if(count == Dimension) this->CalculateCost();
		return 1;
	}
	GainType Cost()	{ return cost;}
	int Right(int node){ return nodes[node].right->Id;}
	int Left(int node){ return nodes[node].left->Id;}
	short IsComplete(){ return (count == Dimension);}
	void InitiateRandomly()
	{
		int *order = RandomArray(Dimension);
		for(int i = 0; i < Dimension; i++)
		{
			Add(order[i]);
		}
		delete order;
	}	
	Tour* Copy()
	{
		Tour *output = new Tour(this->graph);
		int iter = 0;
		for(int i = 0; i < Dimension; i++)
		{
			output->Add(iter);
			iter = Right(iter);
		}
		return output;
	}
	void Write(char* path)
	{
		FILE         *tour_file;
		tour_file = fopen(path, "w");
		int iter = 0;
		for (int i = 0; i < Dimension; i++)
		{
			fprintf(tour_file, "%d\n", iter + 1);
			iter = Right(iter);
		}
		fprintf(tour_file, "\n");
		fclose(tour_file);
	}
	void Print()
	{
		int iter = 0;
		for (int i = 0; i < Dimension; i++)
		{
			cout << iter + 1 << " ";
			iter = Right(iter);
		}
		cout << endl;
	}
	void InitiateOrderly()
	{
		for(int i = 0; i < Dimension; i++)
		{
			Add(i);
		}
	}
	void InitiateNN(int s)
	{
		short *IsAdded = new short[Dimension];
		for(int i = 0; i < Dimension; i++) IsAdded[i] = 0;
		while (!IsComplete())
		{
			IsAdded[s] = 1;
			Add(s);
			int min_node = -1;
			double min_dis = DBL_MAX;
			for(int i = 0; i < Dimension; i++)
			{
				if(!IsAdded[i])
				{
					double dis = graph->D(s, i);
					if(dis < min_dis)
					{
						min_dis = dis;
						min_node = i;
					}
				}
			}
			s = min_node;
		}
		delete IsAdded;
	}
	Graph* getGraph()
	{
		return graph;
	}
};


#endif