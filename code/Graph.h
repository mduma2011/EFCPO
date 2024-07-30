#ifndef _GRAPH_H
#define _GRAPH_H

#include <math.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include <float.h>
#include <iostream>
using namespace std;

class Graph
{
	typedef int (Graph::*CostFunction_1)(int a, int b);
	typedef int (Graph::*CostFunction_2)(double x1, double y1, double x2, double y2) ;
	Graph::CostFunction_1 distance_1;
	Graph::CostFunction_2 distance_2;
	struct point {
		double x;
		double y;
	};
#define M_PI 3.14159265358979323846264
#define LINE_BUF_LEN     100
#define M_RRR 6378388.0
	point*nodeptr;
	int n; 
	char *name;
	static double dtrunc (double x)
	{
		int k;
		k = (int) x;
		x = (double) k;
		return x;
	} 
	int round_distance (int i, int j) 
	{
		double xd = nodeptr[i].x - nodeptr[j].x;
		double yd = nodeptr[i].y - nodeptr[j].y;
		double r  = sqrt(xd*xd + yd*yd) + 0.5;
		return (int) r;
	}
	int ceil_distance (int i, int j)
	{
		double xd = nodeptr[i].x - nodeptr[j].x;
		double yd = nodeptr[i].y - nodeptr[j].y;
		double r  = sqrt(xd*xd + yd*yd) + 0.000000001;
		return (int)r;
	}
	int geo_distance (int i, int j) 
	{
		double deg, min;
		double lati, latj, longi, longj;
		double q1, q2, q3;
		int dd;
		double x1 = nodeptr[i].x, x2 = nodeptr[j].x,
			y1 = nodeptr[i].y, y2 = nodeptr[j].y;
		deg = dtrunc (x1);
		min = x1 - deg;
		lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc (x2);
		min = x2 - deg;
		latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		
		deg = dtrunc (y1);
		min = y1 - deg;
		longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc (y2);
		min = y2 - deg;
		longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		
		q1 = cos (longi - longj);
		q2 = cos (lati - latj);
		q3 = cos (lati + latj);
		dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		return dd;
	}
	int att_distance (int i, int j) 
	{
		double xd = nodeptr[i].x - nodeptr[j].x;
		double yd = nodeptr[i].y - nodeptr[j].y;
		double rij = sqrt ((xd * xd + yd * yd) / 10.0);
		double tij = dtrunc (rij);
		int dij;
		if (tij < rij)
			dij = (int) tij + 1;
		else
			dij = (int) tij;
		return dij;
	}
	int geom_distance(int i, int j)
	{
		point Na =  nodeptr[i];
		point Nb =  nodeptr[j];
		double lati = M_PI * (Na.x / 180.0);
		double latj = M_PI * (Nb.x / 180.0);
		double longi = M_PI * (Na.y / 180.0);
		double longj = M_PI * (Nb.y / 180.0);
		double q1 = cos(latj) * sin(longi - longj);
		double q3 = sin((longi - longj) / 2.0);
		double q4 = cos((longi - longj) / 2.0);
		double q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
		double q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
		return (int) (M_RRR * atan2(sqrt(q1 * q1 + q2 * q2), q5) + 1.0);
	}
	int round_distance (double x1, double y1, double x2, double y2) 
	{
		double xd = x1 - x2;
		double yd = y1 - y2;
		double r  = sqrt(xd*xd + yd*yd) + 0.5;
		return (int) r;
	}
	int ceil_distance (double x1, double y1, double x2, double y2) 
	{
		double xd = x1 - x2;
		double yd = y1 - y2;
		double r  = sqrt(xd*xd + yd*yd) + 0.000000001;
		return (int)r;
	}
	int geo_distance (double x1, double y1, double x2, double y2) 
	{
		double deg, min;
		double lati, latj, longi, longj;
		double q1, q2, q3;
		int dd;
		
		deg = dtrunc (x1);
		min = x1 - deg;
		lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc (x2);
		min = x2 - deg;
		latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc (y1);
		min = y1 - deg;
		longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc (y2);
		min = y2 - deg;
		longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		q1 = cos (longi - longj);
		q2 = cos (lati - latj);
		q3 = cos (lati + latj);
		dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		return dd;
	}
	int att_distance (double x1, double y1, double x2, double y2) 
	{
		double xd = x1 - x2;
		double yd = y1 - y2;
		double rij = sqrt ((xd * xd + yd * yd) / 10.0);
		double tij = dtrunc (rij);
		int dij;
		
		if (tij < rij)
			dij = (int) tij + 1;
		else
			dij = (int) tij;
		return dij;
	}
	int geom_distance(double x1, double y1, double x2, double y2) 
	{
		double lati = M_PI * (x1 / 180.0);
		double latj = M_PI * (x2 / 180.0);
		double longi = M_PI * (y1 / 180.0);
		double longj = M_PI * (y2 / 180.0);
		double q1 = cos(latj) * sin(longi - longj);
		double q3 = sin((longi - longj) / 2.0);
		double q4 = cos((longi - longj) / 2.0);
		double q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
		double q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
		return (int) (M_RRR * atan2(sqrt(q1 * q1 + q2 * q2), q5) + 1.0);
	}
	void read_etsp(char *tsp_file_name) 
	{
		FILE         *tsp_file;
		char         buf[LINE_BUF_LEN];
		long int     i, j;
		tsp_file = fopen(tsp_file_name, "r");
		if ( tsp_file == NULL ) {
			fprintf(stderr,"No instance file specified, abort\n");
		}
		assert(tsp_file != NULL);
		fscanf(tsp_file,"%s", buf);
		while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
			if ( strcmp("NAME", buf) == 0 ) 
			{
				fscanf(tsp_file, "%s", buf);
				fscanf(tsp_file, "%s", buf);
				name = new char[100];
				int i = 0;
				for(i = 0; buf[i];i++)
					name[i] = buf[i];
				name[i] = buf[i];
				buf[0]=0;
			}
			else if ( strcmp("NAME:", buf) == 0 ) 
			{
				fscanf(tsp_file, "%s", buf); 
				name = new char[100];
				int i = 0;
				for(i = 0; buf[i];i++)
					name[i] = buf[i];
				name[i] = buf[i];
				buf[0]=0;
			}
			else if ( strcmp("COMMENT", buf) == 0 )
			{
				fgets(buf, LINE_BUF_LEN, tsp_file);
				buf[0]=0;
			}
			else if ( strcmp("COMMENT:", buf) == 0 )
			{
				fgets(buf, LINE_BUF_LEN, tsp_file);
				buf[0]=0;
			}
	else if ( strcmp("TYPE", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    fscanf(tsp_file, "%s", buf);
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
	    }
	    buf[0]=0;
	}
	else if ( strcmp("TYPE:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
	    }
	    buf[0]=0;
	}
	else if( strcmp("DIMENSION", buf) == 0 ){
	    fscanf(tsp_file, "%s", buf);
	    fscanf(tsp_file, "%ld", &n);
	    buf[0]=0;
	}
	else if ( strcmp("DIMENSION:", buf) == 0 ) {
	    fscanf(tsp_file, "%ld", &n);	    
	    buf[0]=0;
	}
	else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    buf[0]=0;
	}
	else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance_1 = &Graph::round_distance;
		distance_2 = &Graph::round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance_1 = &Graph::ceil_distance;
		distance_2 = &Graph::ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance_1 = &Graph::geo_distance;
		distance_2 = &Graph::geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance_1 = &Graph::att_distance;
		distance_2 = &Graph::att_distance;
	    }
		else if ( strcmp("GEOM", buf) == 0 ) {
		distance_1 = &Graph::geom_distance;
		distance_2 = &Graph::geom_distance;
	    }
	    else
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
	       EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance_1 = &Graph::round_distance;
		distance_2 = &Graph::round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance_1 = &Graph::ceil_distance;
		distance_2 = &Graph::ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance_1 = &Graph::geo_distance;
		distance_2 = &Graph::geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance_1 = &Graph::att_distance;
		distance_2 = &Graph::att_distance;
	    }
		else if ( strcmp("GEOM", buf) == 0 ) {
		distance_1 = &Graph::geom_distance;
		distance_2 = &Graph::geom_distance;
	    }
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
	    }
	    buf[0]=0;
	}
	buf[0]=0;
	fscanf(tsp_file,"%s", buf);
    }


    

    if( (nodeptr = (point*)malloc(sizeof(point) * n)) == NULL );
    else {
	for ( i = 0 ; i < n ; i++ ) {
	    fscanf(tsp_file,"%ld %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y );
	}
    }
	fclose(tsp_file);
}
public:
	Graph(char* path)
	{
		read_etsp(path);
	}
	~Graph(){delete nodeptr, name; nodeptr = 0;}
	int D(int node1,int node2)
	{
		return (this->*distance_1)(node1,node2);
	}
	int D(double x1, double y1, double x2, double y2)
	{
		return (this->*distance_2)(x1,y1,x2,y2);
	}
	int compute_tour_length( int *t ) 
	{
		int      i;
		int tour_length = 0;
		for ( i = 0 ; i < n - 1 ; i++ ) 
		{
			tour_length += (this->*distance_1)(t[i] , t[i+1]);
		}
		tour_length += (this->*distance_1)(t[n - 1],t[0]);
		return tour_length;
	}
	int Dimension(){return n;}
	double X(int node){ return nodeptr[node].x; }
	double Y(int node){ return nodeptr[node].y; }
	int CompareTo(double x, double y)
	{
		if(x < y) return -1;
		if(x > y) return 1;
		return 0;
	}
	
	double Max()
	{
		return DBL_MAX;
	}
	
	double Min()
	{
		return 0;
	}
	char* Name()
	{
		return name;
	}
};
#endif
