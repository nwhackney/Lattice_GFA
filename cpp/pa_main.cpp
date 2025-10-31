#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "xy_class.h"
#include "pa_class.h"
#include "init_functions.h"
#include "params.h"

using namespace std;

int main(void) {
	gsl_rng *r;
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);

	struct xy_params_ xy_params;
	struct mc_params_ mc_params;
	try{
		read_params(&xy_params, &mc_params);
	}
	catch (const char* e) {
		cerr << e << endl;
		return -1;
	}

	int seed = mc_params.seed;

	gsl_rng_set(r, seed);
	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(r);

	try {
		int pop_size = 200;
		PA_simulation aggregate(pop_size, xy_params, r);
		aggregate.run();

	}
	catch (const char* e) {
		cerr << e << endl;
		return -1;
	}
	return 0;
}
