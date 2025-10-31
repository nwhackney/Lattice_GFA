#ifndef PA_CLASS_H
#define PA_CLASS_H

#include "init_functions.h"
#include "params.h"
#include "xy_class.h"
#include <iostream>
#include <memory>
#include <list>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include <omp.h>

using namespace std;

struct PA_member {
	XY_conf xy_conf;
	int family;
};

class PA_simulation {
	private:
		std::vector<PA_member> pop_array;
		int nom_pop;
		int max_pop;
		int pop_size;
		int num_threads;
		xy_params_ xy_params;
		gsl_rng *r;
		vector<gsl_rng *> r_arr;
		int unique_families;
		double rho_t;
		void resample(double *delta_beta_F, double *beta, double avg_e, double var_e, gsl_rng *r);
		void energy_calcs(double *avg_e, double *var_e);
		void family_calcs(void);
	public:
		PA_simulation(void);
		PA_simulation(int nom_pop, xy_params_ xy_params, gsl_rng *r);
		void run(void);
		double get_rho_t(void);
};

#endif
