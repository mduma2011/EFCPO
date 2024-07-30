#ifndef _ESACO_H
#define _ESACO_H

#include "Result.h"
#include "_IACO.h"
#include "mex.h"

class ESACO : public IACO
{
	enum Bool {False, True};
	class ValuesMatrix
	{
		long double tow0;
		int total_count;
		int Dimension;
		ArrayList<long double>** values;
		ArrayList<int>** cl;
	public:
		ValuesMatrix(long double tow0, int Dimension)
		{
			this->total_count = 0;
			this->Dimension = Dimension;
			this->values = new ArrayList<long double>*[Dimension];
			for (int i = 0; i < Dimension; i++) values[i] = new ArrayList<long double>();
			this->cl = new ArrayList<int>*[Dimension];
			for (int i = 0; i < Dimension; i++) cl[i] = new ArrayList<int>();
			this->tow0 = tow0;
		}
		~ValuesMatrix()
		{
			for(int i = 0; i < Dimension; i++)
			{
				delete values[i];
				delete cl[i];
			}
			delete cl;
			delete values;
		}
		long double Get(int i, int j)
		{
			int index = cl[i]->BinarySearch(j);
			if (index < 0) return tow0;
			return values[i]->Get(index);
		}
		void set(int i, int j, long double value)
		{
			int index = cl[i]->BinarySearch(j);			
			if (index < 0)
			{
				this->total_count++;
				index = ~index;
				cl[i]->Insert(index, j);
				values[i]->Insert(index, value);
				return;
			}
			values[i]->Set(index, value);
		}
		void ReplaceLowest(int i, int j, long double value)
		{
			int index = cl[i]->BinarySearch(j);
			if(index >= 0)
			{
				values[i]->Set(index, value);
				return;
			}

			//find lowest
			long double min = DBL_MAX;
			int ret_index = -1;
			int count = values[i]->Count();
			for (int k = 0; k < count; k++)
			{
				if (min > values[i]->Get(k))
				{
					min = values[i]->Get(k);
					ret_index = k;
				}
			}
			//delete lowest
			values[i]->Delete(ret_index);
			cl[i]->Delete(ret_index);			

			//insert new
			index = cl[i]->BinarySearch(j);
			index = ~index;
			cl[i]->Insert(index, j);
			values[i]->Insert(index, value);
		}
		int GetBestOf(int i)
		{
			long double max = DBL_MIN;
			int ret_index = -1;
			int count = values[i]->Count();
			for (int j = 0; j < count; j++)
			{
				if (max < values[i]->Get(j))
				{
					max = values[i]->Get(j);
					ret_index = j;
				}
			}
			return cl[i]->Get(ret_index);
		}
		int NumberOfCandidatesOf(int i)
		{
			return cl[i]->Count();
		}
		int Count()
		{			
			return total_count;
		}	

		int MaxCC()
		{
			int count = 0;
			for(int i = 0; i < this->Dimension; i++)
			{
				if(cl[i]->Count() > count) count = cl[i]->Count();
			}
			return count;
		}
		int IthCandidateOf(int i, int index)
		{
			if(cl[i] == NULL)
			{
				return -1;
			}
			return cl[i]->Get(index);
		}
		bool IsCandidateOf(int i, int j)
		{
			if (cl[i]->BinarySearch(j) >= 0) return true;
			return false;
		}
	};
	
	
	long double delta(Tour *best_tour, GainType best_cost, int r, int s)
	{
		if (best_tour->Right(r) == s || best_tour->Left(r) == s) return 1.0 / (long double)best_cost;
		return 0.0;
	}
	void global_update(Graph *graph, ValuesMatrix *pheromone, ValuesMatrix *total, int NumberOfNodes, Tour *best_tour, GainType best_cost
		, long double alpha, long double beta)
	{           

		long double d_tau = 1.0 / (long double)best_cost;

		for (int i = 0; i < NumberOfNodes; i++)
		{
			int j = i;
			int h = best_tour->Right(i);

			long double p_value = (1.0 - alpha) * pheromone->Get(j,h) + alpha * d_tau;
			pheromone->set(j, h, p_value);
			long double t_value =pow(pheromone->Get(j, h), alpha) * pow(1.0 / (long double)graph->D(h, j), beta);
			total->set(j, h, t_value);
		}
	}
	void local_update(Graph *graph, ValuesMatrix *pheromone, ValuesMatrix *total, long double alpha, long double beta
		, long double rho, long double tow0, int Max_Cl_Count, int r, int s)
	{			
		long double p_value = (1.0 - rho) * pheromone->Get(r, s) + rho * tow0;
		if(p_value == pheromone->Get(r, s))
			return;

		if(pheromone->NumberOfCandidatesOf(r) < Max_Cl_Count)
			pheromone->set(r, s, p_value);
		else
			pheromone->ReplaceLowest(r, s, p_value);

		/*_t_ change bacic ___*/
		long double t_value = pow(pheromone->Get(r, s), alpha)  * pow(1.0 / (long double)graph->D(r, s), beta);
		/*_t_ change bacic ___*/
		
		if(total->NumberOfCandidatesOf(r) < Max_Cl_Count)
			total->set(r, s, t_value);
		else
			total->ReplaceLowest(r, s, t_value);
	}
	int transition_rule(Graph *graph, ValuesMatrix *pheromone,  ValuesMatrix *total, Heuristics *ConstantHeur
		, long double alpha, long double beta, long double tow0_pow_alpha, int r, Ant *ant)
	{
		int best_city = -1;
		long double max = DBL_MIN;
		ant->LU_flag = true;

		int NumberOfCandidatesOf_r = pheromone->NumberOfCandidatesOf(r);
		for (int i = 0; i < NumberOfCandidatesOf_r; i++)
		{
			if (!ant->IsAddedToTour[pheromone->IthCandidateOf(r, i)])
			{
				for (int index = i; index < NumberOfCandidatesOf_r; index++)
				{
					int u = pheromone->IthCandidateOf(r, index);
					if(u == -1)
					{
						continue;
					}

					if (ant->IsAddedToTour[u]) continue;
					long double value = GetTotal(graph, pheromone, total, alpha, beta, tow0_pow_alpha, r, u);
					if (value > max)
					{
						max = value;
						best_city = u;
					}
				}
				break;
			}
		}

		int NumberOfCandidates = ConstantHeur->GetNumberOfCandidate();
		for (int i = 0; i < NumberOfCandidates; i++)
		{
			if (!ant->IsAddedToTour[ConstantHeur->GetCandidate(r,i)])
			{
				for (int index = i; index < NumberOfCandidates; index++)
				{
					int u = ConstantHeur->GetCandidate(r, index);
					if (ant->IsAddedToTour[u]) continue;
					long double value = GetTotal(graph, pheromone, total, alpha, beta, tow0_pow_alpha, r, u);
					if (value > max)
					{
						max = value;
						best_city = u;
					}
				}
				if (best_city >= 0)
					return best_city;
				break;
			}
		}

		for (int u = 0; u < ant->number_of_r_nodes; u++)
		{
			long double value = GetTotal(graph, pheromone, total, alpha, beta, tow0_pow_alpha, r, ant->r_nodes[u]);
			if (value > max)
			{
				max = value;
				best_city = ant->r_nodes[u];
			}
		}
		return best_city;
	}
	long double GetTotal(Graph *graph, ValuesMatrix *pheromone, ValuesMatrix *total, long double alpha, long double beta
		, long double tow0_pow_alpha, int r, int s)
	{
		long double t_value = total->Get(r, s);
		if (t_value >= 0) return t_value;
		return tow0_pow_alpha * pow(1.0 / (long double)graph->D(r, s), beta);
	}
public:
	char* Name(){return "ESACO";}
	GainType Run(Graph *graph, int number_of_ants, int number_of_iteration, double &TimeB, double &TotT)
	{

        long double alpha = 0.9, rho = 0.6, beta = 2;

		srand(time(0));

		TimeB = TotT = 0;
		double start_time = GetTime();

		int NumberOfNodes = graph->Dimension();	

		GainType best_cost = LLONG_MAX;
		int Max_Cl_Count = 2;  //10
        
        Heuristics *ConstantHeur = new Heuristics(graph, 2);  //2-20

        double t0 = GetTime(), t1;
		Ant **ants = new Ant*[number_of_ants]; for (int i = 0; i < number_of_ants; i++) ants[i] = new Ant(graph);
		int *indexes = new int[NumberOfNodes];for (int i = 0; i < NumberOfNodes; i++) indexes[i] = 0;

		Heuristics *heur = new Heuristics(graph, 2);  //2-10
		Tour *best_tour = heur->Q_Boruvka();
        
        long double tow0 = 1.0 / (long double)(best_tour->Cost() * NumberOfNodes);
		long double tow0_pow_alpha = pow(tow0, alpha);

		ValuesMatrix *pheromone = new ValuesMatrix(tow0, NumberOfNodes);
		ValuesMatrix *total = new ValuesMatrix(-1.1, NumberOfNodes);
        
        for (int step = 1; step <= number_of_iteration; step++)
		{
			/*Initialization phase*/
			for (int k = 0; k < number_of_ants; k++)
			{
				ants[k]->Rest();
				ants[k]->GoTo(Random(0, NumberOfNodes - 1));
			}
            
            /*This is the phase in which ants build their tours.*/
			for (int i = 1; i < NumberOfNodes; i++)
			{	
				for (int k = 0; k < number_of_ants; k++)
				{	
				
					int current = ants[k]->LastAddedCity();
					int next = this->transition_rule(graph, pheromone, total, ConstantHeur, alpha, beta, tow0_pow_alpha, current, ants[k]);	
					ants[k]->GoTo(next);
					this->local_update(graph, pheromone, total, alpha, beta, rho, tow0, Max_Cl_Count, current, next);
					if (i == NumberOfNodes - 1)
					{
						heur->TwoOpt(ants[k]->tour);
					}
				}
			}
            
            /*In this phase best tour is selected*/			
			for (int k = 0; k < number_of_ants; k++)
			{
				GainType cost = ants[k]->tour->Cost();
				if (cost < best_cost)
				{
					delete best_tour;
					best_tour = ants[k]->tour->Copy();
					best_cost = cost;
					TimeB = GetTime() - start_time;
				}
			}
            
            /*In this phase global updating occurs and pheromone is updated.*/
			global_update(graph, pheromone, total, NumberOfNodes, best_tour, best_cost, alpha, beta);

			/*In this phase a one node is injected to candidates of each node according to Pheromones_or_TR*/
			int n_o_c = heur->GetNumberOfCandidate();
			for (int i = 0; i < NumberOfNodes; i++)
			{
				if(pheromone->NumberOfCandidatesOf(i) > 0)
					heur->SetCandidates(i, pheromone->GetBestOf(i), n_o_c - 1);
				heur->SetCandidates(i, best_tour->Right(i), 0);
			}
			heur->SetBestTour(best_tour);
		}
        
       	/*Destroy*/
		for(int i = 0; i < number_of_ants; i++) delete ants[i]; delete ants;
		delete indexes;
		delete pheromone;
		delete total;
		delete heur;
		GainType best_cost_sofar = best_tour->Cost();
        		
		delete best_tour;
		delete ConstantHeur;
		TotT = GetTime() - start_time;
		return best_cost_sofar;        
	}    
    
    Result MatlabRun(Graph *graph, int number_of_ants, int number_of_iteration, double alpha, double rho, double beta, int Max_Cl_Count, int NumberOfCandidatesConst, int NumberOfCandidates)
	{
		srand(time(0));
		int NumberOfNodes = graph->Dimension();
		GainType best_cost = LLONG_MAX;
		        
        Heuristics *ConstantHeur = new Heuristics(graph, NumberOfCandidatesConst);  //2-20

    	Ant **ants = new Ant*[number_of_ants]; for (int i = 0; i < number_of_ants; i++) ants[i] = new Ant(graph);
		int *indexes = new int[NumberOfNodes];for (int i = 0; i < NumberOfNodes; i++) indexes[i] = 0;

		Heuristics *heur = new Heuristics(graph, NumberOfCandidates);  //2-10
		Tour *best_tour = heur->Q_Boruvka();
        
        long double tow0 = 1.0 / (long double)(best_tour->Cost() * NumberOfNodes);
		long double tow0_pow_alpha = pow(tow0, alpha);

		ValuesMatrix *pheromone = new ValuesMatrix(tow0, NumberOfNodes);
		ValuesMatrix *total = new ValuesMatrix(-1.1, NumberOfNodes);
        
        for (int step = 1; step <= number_of_iteration; step++)
		{
			/*Initialization phase*/
			for (int k = 0; k < number_of_ants; k++)
			{
				ants[k]->Rest();
				ants[k]->GoTo(Random(0, NumberOfNodes - 1));
			}
            
            /*This is the phase in which ants build their tours.*/
			for (int i = 1; i < NumberOfNodes; i++)
			{	
				for (int k = 0; k < number_of_ants; k++)
				{	
				
					int current = ants[k]->LastAddedCity();
					int next = this->transition_rule(graph, pheromone, total, ConstantHeur, alpha, beta, tow0_pow_alpha, current, ants[k]);	
					ants[k]->GoTo(next);
					this->local_update(graph, pheromone, total, alpha, beta, rho, tow0, Max_Cl_Count, current, next);
					if (i == NumberOfNodes - 1)
					{
						heur->TwoOpt(ants[k]->tour);
					}
				}
			}
            
            /*In this phase best tour is selected*/			
			for (int k = 0; k < number_of_ants; k++)
			{
				GainType cost = ants[k]->tour->Cost();
				if (cost < best_cost)
				{
					delete best_tour;
					best_tour = ants[k]->tour->Copy();
					best_cost = cost;
				}
			}
            
            /*In this phase global updating occurs and pheromone is updated.*/
			global_update(graph, pheromone, total, NumberOfNodes, best_tour, best_cost, alpha, beta);

			/*In this phase a one node is injected to candidates of each node according to Pheromones_or_TR*/
			int n_o_c = heur->GetNumberOfCandidate();
			for (int i = 0; i < NumberOfNodes; i++)
			{
				if(pheromone->NumberOfCandidatesOf(i) > 0)
					heur->SetCandidates(i, pheromone->GetBestOf(i), n_o_c - 1);
				heur->SetCandidates(i, best_tour->Right(i), 0);
			}
			heur->SetBestTour(best_tour);
		}
        
        GainType best_cost_sofar = best_tour->Cost();
        
        // copy the results
		Result results;		
		int** resants = new int*[number_of_ants];
		for (int k = 0; k < number_of_ants; k++)
		{
			resants[k] = new int[NumberOfNodes];
			for(int j = 0; j < NumberOfNodes; j++)
			{
				resants[k][j] =  ants[k]->tour->Right(j);				
			}			
		}
        results.ants = resants;
        
		double** reseta = new double*[NumberOfNodes];
		for (int r = 0; r < NumberOfNodes; r++)
		{
			reseta[r] = new double[NumberOfNodes];
			for (int c = 0; c < NumberOfNodes; c++)
			{
				reseta[r][c] =  1.0 / (long double)graph->D(r, c);
			}
		}
        results.ETA = reseta;
		results.best_tour = best_tour;
        
       	/*Destroy*/
		for(int i = 0; i < number_of_ants; i++) delete ants[i]; delete ants;
		delete indexes;
		delete pheromone;
		delete total;
		delete heur;		
		delete ConstantHeur;		
        return results;        
	}  
};

#endif

