#include "xy_class.h"
#include "init_functions.h"
#include "params.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_statistics.h>

using namespace std;

int main(void) {
	gsl_rng *r;
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	r = gsl_rng_alloc(T);

	struct xy_params_ xy_params;
	struct auto_params_ auto_params;
	try{
		read_auto_params(&xy_params, &auto_params);
	}
	catch (const char* e) {
		cerr << e << endl;
		return -1;
	}

	int seed = auto_params.seed;
	int num_steps = auto_params.num_steps;
	int step_size = auto_params.step_size;

	gsl_rng_set(r, seed);
	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(r);

	double f = xy_params.f;
	int length = xy_params.length;
	int N = xy_params.num_particles;

	int neighbor_table[length *length][4];
	double a_table[length *length][4];
	
	init_neighbor_table(neighbor_table, length);
	init_gauge_table(a_table, f, length);

	if (N < auto_params.num_saved_clusters)
		auto_params.num_saved_clusters = N;


	XY_conf test_conf;

	FILE *config_fp, *phase_fp, *strain_fp, *cluster_fp, *cluster_size_fp;
	
	if (auto_params.restart == 0) {
		config_fp = fopen("auto_configs.dat", "w");
		phase_fp = fopen("auto_phase.dat", "w");
		strain_fp = fopen("auto_strain.dat", "w");
		cluster_fp = fopen("clusters.dat", "w");
		cluster_size_fp = fopen("avg_cluster_size.dat", "w");
		test_conf.initialize(xy_params, r, neighbor_table, a_table);
	}
	else {
		config_fp = fopen("auto_configs.dat", "a");
		phase_fp = fopen("auto_phase.dat", "a");
		strain_fp = fopen("auto_strain.dat", "a");
		cluster_fp = fopen("clusters.dat", "a");
		cluster_size_fp = fopen("avg_cluster_size.dat", "a");
		FILE *fp = fopen("checkpoint.dat", "r");
		if (fp == NULL)
			throw invalid_argument("No input checkpoint file.");
		test_conf.external_initialize(xy_params, r, neighbor_table, a_table, fp);
		fclose(fp);
	}

	try {
		for (int step = 1; step < num_steps; step++) {
			if (step % (num_steps /10) == 0) {
				test_conf.write_config(config_fp, strain_fp, phase_fp);
				fflush(config_fp);
			}
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
			
			int particle_cluster_size[N];
			test_conf.cluster_calculations(particle_cluster_size, NULL, NULL);
			for (int i = 0; i < auto_params.num_saved_clusters; i++)
				fprintf(cluster_fp, "%d ", particle_cluster_size[i]);
			fprintf(cluster_fp, "\n");
			fflush(cluster_fp);
			
			fprintf(cluster_size_fp, "%15.15e\n", test_conf.get_avg_cluster_size());
			fflush(cluster_size_fp);
			test_conf.write_checkpoint();
		}
	}
	catch (const char* e) {
		cerr << e << endl;
		return -1;
	}
	fclose(config_fp);
	fclose(phase_fp);
	fclose(strain_fp);
	fclose(cluster_fp);
	return 0;
}
