#ifndef _IACO_H
#define _IACO_H

#include <conio.h>

#include "Time_Random.h"
#include "Heuristic.h"
#include "Graph.h"
#include "List.h"
#include "Tour.h"
#include "BST.h"
#include "GainType.h"
#include "Result.h"

class IACO
{
public:
	virtual char* Name() = 0;	
	virtual GainType Run(Graph *Func, int number_of_ants, int number_of_iteration, double &TimeB, double &TotT) = 0;
    virtual Result MatlabRun(Graph *graph, int number_of_ants, int number_of_iteration, double alpha, double rho, double beta, int Max_Cl_Count, int NumberOfCandidatesConst, int NumberOfCandidates) = 0;

protected:
	class Ant
	{
		Graph *graph;
		int NumberOfNodes;		

		int LastAdded;
	public:
		short LU_flag;
		short *IsAddedToTour;
		int *PositionOf;
		int *r_nodes;
		int number_of_r_nodes;
		Tour *tour;
		Ant(Graph *graph)
		{
			this->LU_flag = true;
			number_of_r_nodes = NumberOfNodes = graph->Dimension();
			tour = new Tour(graph);
			this->graph = graph;
			IsAddedToTour = new short[NumberOfNodes];
			PositionOf = new int[NumberOfNodes];
			r_nodes = new int[NumberOfNodes];
			for (int i = 0; i < NumberOfNodes; i++) 
			{
				IsAddedToTour[i] = 0;
				PositionOf[i] = i;
				r_nodes[i] = i;
			}
			LastAdded = -1;
		}
		~Ant()
		{
			delete IsAddedToTour;
			delete PositionOf;
			delete r_nodes;
			delete tour;
		}
		void GoTo(int NextCity)
		{
			if (tour->Add(NextCity))
			{
				LastAdded = NextCity;
				IsAddedToTour[NextCity] = true;
				r_nodes[PositionOf[NextCity]] = r_nodes[number_of_r_nodes - 1];
				PositionOf[r_nodes[number_of_r_nodes - 1]] = PositionOf[NextCity];
				PositionOf[NextCity] = -1;
				number_of_r_nodes--;
			}
		}
		int LastAddedCity()
		{
			return LastAdded;
		}
		void Rest()
		{
			this->LU_flag = true;
			number_of_r_nodes = NumberOfNodes;
			this->tour->Reset();
			for (int i = 0; i < NumberOfNodes; i++)
			{
				IsAddedToTour[i] = 0;
				PositionOf[i] = i;
				r_nodes[i] = i;
			}
		}
	};
	enum Bool {False, True};	
};

#endif