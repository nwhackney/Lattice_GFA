#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <gsl/gsl_histogram.h>
#include "xy_class.h"
#include "init_functions.h"
#include "params.h"

using namespace std;

void print_histogram(gsl_histogram *h, int bin_flag, FILE *fp) {

	if (bin_flag == 1) {
		for (int i = 0; i < (int) h->n; i++) {
			fprintf(fp, "%g ", h->range[i]);
		}
		fprintf(fp, "\n");
	}
	else {
		double norm = 0.0;
		double dx = h->range[1] - h->range[0];
		for (int i = 0; i < (int) h->n; i++) {
			norm += h->bin[i];
		}
		norm *= dx;
		for (int i = 0; i < (int) h->n; i++) {
			fprintf(fp, "%5.5e ", (h->bin[i]) /norm);
		}
		fprintf(fp, "\n");
	}
	fflush(fp);

}

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
	int num_blocks = mc_params.num_steps;
//	step_size is usually the autocorrelation time
	int step_size = mc_params.step_size;
	
	int stepsize_div_ten = step_size /10;
	if (stepsize_div_ten < 1)
		stepsize_div_ten = 1;

	gsl_rng_set(r, seed);
	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(r);


	double f = xy_params.f;
	int length = xy_params.length;
	int N = xy_params.num_particles;

	int neighbor_table[length *length][4];
	double a_table[length *length][4];
	int thinning_neighbor_table[length *length][8];
	
	init_neighbor_table(neighbor_table, length);
	init_gauge_table(a_table, f, length);
	init_thinning_neighbor_table(thinning_neighbor_table, length);

	XY_conf test_conf;

	FILE *config_fp, *phase_fp, *strain_fp, *phi_n_fp, *w_hist_fp, *energy_fp;
	gsl_histogram *w_hist = gsl_histogram_alloc(150);
	gsl_histogram_set_ranges_uniform(w_hist, 0, 5);

//	If not using an initial file from the autocorrelation simulation, we need to have a burn in step.
	if (mc_params.restart == 0) {
		FILE *fp = fopen("w_bins.dat", "w");
		print_histogram(w_hist, 1, fp);
		fclose(fp);
		energy_fp = fopen("energy.dat", "w");
		config_fp = fopen("configs.dat", "w");
		phase_fp = fopen("phase.dat", "w");
		strain_fp = fopen("strain.dat", "w");
		phi_n_fp = fopen("phi_n_hist.dat", "w");
		w_hist_fp = fopen("w_hist.dat", "w");
		test_conf.initialize(xy_params, r, neighbor_table, a_table);
/*
		for (int t_steps = 0; t_steps < 100; t_steps++) {
			for (int sweep = 0; sweep < step_size; sweep++) {
				for (int i = 0; i < N; i++) {
					int index = gsl_rng_uniform_int(r, N);
					test_conf.local_swap_MC_move(index, r);

					index = gsl_rng_uniform_int(r, N);
					test_conf.global_swap_MC_move(index, r);

					index = gsl_rng_uniform_int(r, N);
					test_conf.rotation_MC_move(index, r);

					index = gsl_rng_uniform_int(r, N);
					test_conf.global_swap_rotation_MC_move(index, r);
				}
			}
		}
*/
	}
	else {
		FILE *fp = fopen("w_bins.dat", "w");
		print_histogram(w_hist, 1, fp);
		fclose(fp);
		energy_fp = fopen("energy.dat", "a");
		config_fp = fopen("configs.dat", "a");
		phase_fp = fopen("phase.dat", "a");
		strain_fp = fopen("strain.dat", "a");
		phi_n_fp = fopen("phi_n_hist.dat", "a");
		w_hist_fp = fopen("w_hist.dat", "a");
		fp = fopen("checkpoint.dat", "r");
		if (fp == NULL)
			throw invalid_argument("No input checkpoint file.");
		test_conf.external_initialize(xy_params, r, neighbor_table, a_table, fp);
		fclose(fp);
	}

	try {
//		Data collection loop. Each block is equal to 
		for (int block = 0; block < num_blocks; block++) {
			test_conf.set_beta(40.0+0.1*float(block));
			gsl_histogram_reset(w_hist);
			double phi_n_hist[xy_params.num_particles];
			for (int i = 0; i < xy_params.num_particles; i++)
				phi_n_hist[i] = 0.0;
			double energy = 0.0;
			int measurement;
//			10 *step_size measurements per block/writeout. Typically this should be equivalent to
//			10 autocorrelation times per block/writeout. Each block can then be treated as independent.
//			(for 10 autocorrelation times per block, there should be 100 measurements -- it is currently set to 1 auto/block)
			for (measurement = 0; measurement < 10; measurement++) {
				for (int sweep = 0; sweep < stepsize_div_ten; sweep++) {
					//int index = gsl_rng_uniform_int(r, N);
					//test_conf.point_reflect_wolff(index, r, f);
					//index = gsl_rng_uniform_int(r, N);
					//test_conf.classic_wolff(index, r, f);
					for (int i = 0; i < N; i++) {
						int index = gsl_rng_uniform_int(r, N);
						test_conf.local_swap_MC_move(index, r);

						index = gsl_rng_uniform_int(r, N);
						test_conf.global_swap_MC_move(index, r);

						index = gsl_rng_uniform_int(r, N);
						test_conf.rotation_MC_move(index, r);

						index = gsl_rng_uniform_int(r, N);
						test_conf.global_swap_rotation_MC_move(index, r);
					}
				}
				test_conf.cluster_calculations_2(phi_n_hist, w_hist);
				energy += test_conf.calc_energy();
			}
			//test_conf.write_config(config_fp, strain_fp, phase_fp);
			print_histogram(w_hist, 0, w_hist_fp);
			for (int i = 0; i < xy_params.num_particles; i++) {
				if (phi_n_hist[i] > 0.0)
					fprintf(phi_n_fp, "%d ", i);
			}
			fprintf(phi_n_fp, "\n");
			for (int i = 0; i < xy_params.num_particles; i++) {
				if (phi_n_hist[i] > 0.0)
					fprintf(phi_n_fp, "%10.10e ", phi_n_hist[i] /(double) measurement);
			}
			fprintf(phi_n_fp, "\n");
			fflush(phi_n_fp);
			fprintf(energy_fp, "%15.15e\n", energy /((double) measurement *N));
			fflush(energy_fp);
			test_conf.write_checkpoint();
		}
		fclose(phi_n_fp);
		fclose(phase_fp);
		fclose(config_fp);
		fclose(strain_fp);
		fclose(energy_fp);
		fclose(w_hist_fp);
		gsl_histogram_free (w_hist);
	}
	catch (const char* e) {
		cerr << e << endl;
		return -1;
	}

	return 0;
}
