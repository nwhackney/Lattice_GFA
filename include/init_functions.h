#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "params.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <boost/bind/bind.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//int array_index(int i, int j, int length);
//int mod(int a, int b);
void read_params(xy_params_ *xy_params, mc_params_ *mc_params);
void read_auto_params(xy_params_ *xy_params, auto_params_ *auto_params);
void init_gauge_table(double a_table[][4], double f, int length);
void init_neighbor_table(int neighbor_table[][4], int length);
void init_fixed_BC_neighbor_table(int neighbor_table[][4], int length);
void init_thinning_neighbor_table(int thinning_neighbor_table[][8], int length);

#endif
