#include "pa_class.h"
#include "xy_class.h"
#define MAX_BETA 40
#define CULLING_FRAC 0.05
using namespace std;

PA_simulation::PA_simulation(void) {
	throw invalid_argument("Invalid PA simulation initialization parameters.");
}

PA_simulation::PA_simulation(int nom_pop, xy_params_ xy_params, gsl_rng *r) {
	PA_simulation::nom_pop = nom_pop;
	PA_simulation::max_pop = nom_pop + sqrt(nom_pop *10);
	PA_simulation::pop_size = nom_pop;
	PA_simulation::xy_params = xy_params;

	PA_simulation::pop_array.resize(PA_simulation::max_pop);

	PA_simulation::r = r;
	num_threads = omp_get_max_threads();
	num_threads = 2;
	omp_set_num_threads(num_threads);
	r_arr.resize(num_threads);

//	Should probably ensure that no two RNG seeds can be the same (low probability, but still should fix)
	for (int i = 0; i < num_threads; i++) {
		r_arr[i] = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(r_arr[i], gsl_rng_uniform_int(r, 10000000) + 1000000);
		for (int j = 0; j < 10000; j++)
			gsl_rng_uniform(r_arr[i]);
	}
	for (int m = 0; m < nom_pop; m++) {
		int ip[2][4]; double dp[2][4];
		pop_array[m].xy_conf.initialize(xy_params, r, ip, dp);
		pop_array[m].family = m;
	}

//	double dummy_1, dummy_2;
//	energy_calcs(&dummy_1, &dummy_2);
}

/*
void PA_simulation::family_calcs(void) {
	int i;
//	int family_size_dist[nom_pop];
	int family_size[nom_pop];

	for (i = 0; i < nom_pop; i++) {
		family_size[i] = 0;
//		family_size_dist[i] = 0;
	}

//	Count the number of ancestors for each family
	for (i = 0; i < pop_size; i++)
		family_size[pop_array[i].get_family()]++;

//	Count the distribution of family sizes
//	for (i = 0; i < nom_pop; i++)
//		family_size_dist[family_size[i]]++;

	int counter = 0;
	rho_t = 0;

	for (i = 0; i < nom_pop; i++) {
		if (family_size[i] > 0) {
			counter++;
			rho_t += (double) family_size[i] *family_size[i] /pop_size /pop_size;
		}
	}

	unique_families = counter;
	rho_t *= nom_pop;
}
*/

void PA_simulation::resample(double *delta_beta_F, double *beta, double avg_e, double var_e, gsl_rng *r) {
	int new_pop_size = 0;
	double config_weight[pop_size];
	double Q = 0;

	int num_replicas[pop_size];
	for (int i = 0; i < pop_size; i++)
		num_replicas[i] = 0;

	double delta_beta = CULLING_FRAC *sqrt(2 *M_PI /var_e);

	if (*beta + delta_beta < MAX_BETA)
		*beta += delta_beta;
	else
		*beta = MAX_BETA;

//	Calculate configuration weight

	for (int m = 0; m < pop_size; m++) {
		double temp_e = pop_array[m].xy_conf.calc_energy();
		double temp_w = -(delta_beta) *(temp_e);
		config_weight[m] = exp(temp_w);
		Q += config_weight[m];
	}

	Q /= (double) pop_size;
	*delta_beta_F = -log(Q);

//	Calculate the number of replicas (probabilistically)
	for (int m = 0; m < pop_size; m++) {
		double tau = (nom_pop / ((double) pop_size)) *(config_weight[m] /Q);
		int floor = (int) tau;
		int ceiling = floor + 1;

		if (gsl_rng_uniform(r) < ceiling - tau)
			num_replicas[m] = floor;
		else
			num_replicas[m] = ceiling;
		new_pop_size += num_replicas[m];
	}

	if (new_pop_size > max_pop) {
		throw "Maximum population size exceeded.";
	}

//	Fill empty gaps in population array so un-erased members are contiguous
//	This is to avoid deleting vector elements which can be costly.

	int copy_to = 0, copy_from = pop_size - 1;
	while (copy_to < copy_from && copy_to < pop_size && copy_from > 0) {
		while (num_replicas[copy_to] > 0 && copy_to < pop_size)
			copy_to++;
		while (num_replicas[copy_from] <= 0 && copy_from > 0)
			copy_from--;
		if (copy_to < copy_from) {
			pop_array[copy_to] = pop_array[copy_from];
			num_replicas[copy_to] = num_replicas[copy_from];
			num_replicas[copy_from] = 0;
		}
	}

//	Make copies of config's from the beginning to the end of contiguous data block
	copy_to = 0;
	while (num_replicas[copy_to] > 0 && copy_to < pop_size)
		copy_to++;
	int copy_end = copy_to;
	copy_from = 0;

	while (copy_from < copy_end) {
		for (int m = 0; m < num_replicas[copy_from] - 1; m++) {
			pop_array[copy_to] = pop_array[copy_from];
			copy_to++;
		}
		copy_from++;
	}
	PA_simulation::pop_size = new_pop_size;
}

void PA_print(FILE *hist_fp, FILE *w_hist_fp, FILE *config_fp, FILE *energy_fp) {
	
}

void print_histogram(gsl_histogram *h, int bin_flag, FILE *fp) {

	if (bin_flag == 1) {
		for (int i = 0; i < (int) h->n; i++) {
			fprintf(fp, "%g ", h->range[i]);
		}
		fprintf(fp, "\n");
	}
	else {
		for (int i = 0; i < (int) h->n; i++) {
			fprintf(fp, "%g ", h->bin[i]);
		}
		fprintf(fp, "\n");
	}
	fflush(fp);

}

/*
void print_lists(list<double> beta_list, list<double> avg_energy_list, list<double> avg_agg_size_list, list<double> free_energy_list) {
	FILE *fp = fopen("observables.json", "w");

	fprintf(fp, "{\n"}

	fprintf(fp, "\t beta : [");
	for (list<int>::iterator j = beta_list.begin(); j != beta_list.end(); ++j)
		fprintf(

}
*/

void PA_simulation::run(void) {

//	gsl_histogram *w_hist = gsl_histogram_alloc(150);
//	gsl_histogram_set_ranges_uniform(w_hist, 0, 5);

	double avg_e = 0.0, var_e = 0.0;

/*
	list<double> beta_list = {};
	list<double> avg_energy_list = {};
	list<double> avg_agg_size_list = {};
	list<double> free_energy_list = {};
*/

	double beta = 37;
	double beta_F = 0.0;
	int tmp_num_sweeps = 50;
	int num_sweeps = tmp_num_sweeps;

//	FILE *w_hist_fp = fopen("w_hist.dat", "w");
//	FILE *hist_fp = fopen("hist.dat", "w");
	FILE *config_fp = fopen("configs.dat", "w");
	FILE *strain_fp = fopen("strain.dat", "w");
	FILE *phase_fp = fopen("phase.dat", "w");
	FILE *obs_fp = fopen("obs.dat", "w");


	int neighbor_table[num_threads][xy_params.length *xy_params.length][4];
	double a_table[num_threads][xy_params.length *xy_params.length][4];
	for (int i = 0; i < num_threads; i++) {
		init_neighbor_table(neighbor_table[i], xy_params.length);
		init_gauge_table(a_table[i], xy_params.f, xy_params.length);
	}
	
	for (int i = 0; i < pop_size; i++)
		pop_array[i].xy_conf.set_tables(neighbor_table[0], a_table[0]);
		
	#pragma omp parallel for shared(pop_array, r_arr, neighbor_table, a_table)
	for (int m = 0; m < pop_size; m++) {
		int thread = omp_get_thread_num();
		pop_array[m].xy_conf.set_beta(beta);
		pop_array[m].xy_conf.set_tables(neighbor_table[thread], a_table[thread]);
		for (int swp = 0; swp < tmp_num_sweeps *10; swp++) {
			for (int n = 0; n < xy_params.num_particles; n++) {
				int index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
				pop_array[m].xy_conf.local_swap_MC_move(index, r_arr[thread]);

				index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
				pop_array[m].xy_conf.global_swap_MC_move(index, r_arr[thread]);
				
				index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
				pop_array[m].xy_conf.rotation_MC_move(index, r_arr[thread]);

				index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
				pop_array[m].xy_conf.global_swap_rotation_MC_move(index, r_arr[thread]);
			}
		}
	}

	energy_calcs(&avg_e, &var_e);
	while (beta < MAX_BETA) {
		printf("beta = %5.5f, pop_size = %d, avg energy = %5.5e\n", beta, pop_size, avg_e /xy_params.num_particles);
		double delta_beta_F = 0.0;
		resample(&delta_beta_F, &beta, avg_e, var_e, r);
//		beta = beta + 0.01;
//		if (beta > 40)
//			beta = 40;
		beta_F += delta_beta_F;
		xy_params.beta = beta;
//		if (beta > 30)
//			num_sweeps = 5 *tmp_num_sweeps;
		#pragma omp parallel for shared(pop_array, r_arr, neighbor_table, a_table)
		for (int m = 0; m < pop_size; m++) {
			int thread = omp_get_thread_num();
			pop_array[m].xy_conf.set_beta(beta);
			pop_array[m].xy_conf.set_tables(neighbor_table[thread], a_table[thread]);
			for (int swp = 0; swp < num_sweeps; swp++) {
				for (int n = 0; n < xy_params.num_particles; n++) {
					int index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
					pop_array[m].xy_conf.local_swap_MC_move(index, r_arr[thread]);

					index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
					pop_array[m].xy_conf.global_swap_MC_move(index, r_arr[thread]);
					
					index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
					pop_array[m].xy_conf.rotation_MC_move(index, r_arr[thread]);

					index = gsl_rng_uniform_int(r_arr[thread], xy_params.num_particles);
					pop_array[m].xy_conf.global_swap_rotation_MC_move(index, r_arr[thread]);
				}
			}
		}
		energy_calcs(&avg_e, &var_e);
		family_calcs();
//		test_conf.cluster_calculations_4(hist_fp, w_hist, cluster_e_fp);
		fprintf(obs_fp, "%15.15e %15.15e %15.15e %15.15e %15.15e\n", beta, avg_e, var_e, beta_F, rho_t);
		fflush(obs_fp);
//		pop_array[0].xy_conf.write_config(config_fp, strain_fp, phase_fp);
		pop_array[0].xy_conf.write_checkpoint();
	}
	energy_calcs(&avg_e, &var_e);
//	family_calcs();
}

void PA_simulation::energy_calcs(double *avg_e, double *var_e) {
//	int min_index = 0;
//	double temp_min = 0;
	double avg_e2 = 0.0;
	double temp_avg_e = 0.0;
	for (int m = 0; m < pop_size; m++) {
		double temp_e = pop_array[m].xy_conf.calc_energy();
		temp_avg_e += temp_e;
		avg_e2 += temp_e *temp_e;
//		if (temp_e < temp_min) {
//			temp_min = temp_e;
//			min_index = m;
//		}
	}
	avg_e2 /= pop_size;
	temp_avg_e /= pop_size;
	*avg_e = temp_avg_e;
	*var_e = (avg_e2 - temp_avg_e *temp_avg_e);
//	gs = pop_array[min_index];
//	gs_e = gs.get_energy();
}

void PA_simulation::family_calcs() {
	int family_size[nom_pop];

	for (int m = 0; m < nom_pop; m++)
		family_size[m] = 0;

//	Count the number of ancestors for each family
	for (int m = 0; m < pop_size; m++)
		family_size[pop_array[m].family]++;

	int counter = 0;
	rho_t = 0;

	for (int m = 0; m < nom_pop; m++) {
		if (family_size[m] > 0) {
			counter++;
			rho_t += ((double) family_size[m] *family_size[m]) /pop_size /pop_size;
		}
	}

//	unique_families = counter;
	rho_t *= nom_pop;
}

double PA_simulation::get_rho_t(void) {
	return rho_t;
}
