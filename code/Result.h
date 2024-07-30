/*
 * Result.h
 *
 *  Created on: 22 Oct 2023
 *      Author: MD2020
 */

#ifndef RESULT_H_
#define RESULT_H_

#include "Tour.h"
#include <iostream>
using namespace std;

struct Result
{
   int** ants;
   double** ETA;
   Tour *best_tour;
};

#endif /* RESULT_H_ */
