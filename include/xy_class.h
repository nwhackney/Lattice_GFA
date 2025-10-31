#ifndef XY_CLASS_H
#define XY_CLASS_H

#include "params.h"
#include <iostream>
#include <vector>
#include <list>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <math.h>
#include <omp.h>

//using namespace std;
// XY spin struct that contains site occupancy (n) and angle (theta);
struct XY_spin {
	unsigned char n;
	double theta;
};

/*
All unique_ptr variables correspond to arrays of size length**2, to be initialized on creation
num_particles <= length**2
particle_location = vector of coordinates of all particles
lattice = length**2 array of xy spins (note that many will be zero occupancy, n)
*/

class XY_conf {
	private:
		int num_particles;
		int num_clusters;
		double avg_cluster_size;
		int length;
		int (*neighbor_table)[4];
		double (*a_table)[4];
		double beta;
		double K;
		std::vector<int> particle_location;
		std::vector<XY_spin> lattice;
		XY_spin null_spin;

		void swap_particles(int loc_1, int loc_2);
		void move_particle(int particle_index, int vac_loc);
		double calc_local_energy(int loc);
		void calc_bond_energy(int loc, double energy[4]);
//		void wolff_swap(int loc_1, int loc_2, int coord[2], int new_coord[2], int dx, int dy, double f, std::vector<int> particle_location[], int particle_index[]);
		//void find_cluster(int particle_index);
	public:
//		Initialize and instantiate
		XY_conf();
		void ribbon_initialize(xy_params_ xy_params, double ribbon_width, double ribbon_length, gsl_rng *r, int neighbor_table[][4], double a_table[][4]);
		void square_initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4]);
		void initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4]);
		void external_initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4], FILE *init_config_fp);
//		Set observables
		void set_tables(int neighbor_table[][4], double a_table[][4]);
		void set_beta(double beta);
//		Get observables
		double get_avg_cluster_size(void);
//		MC moves
		void classic_wolff(int index, gsl_rng *r, double f);
		void point_reflect_wolff(int index, gsl_rng *r, double f);
		void global_swap_MC_move(int particle_index, gsl_rng *r);
		void local_swap_MC_move(int particle_index, gsl_rng *r);
		void rotation_MC_move(int particle_index, gsl_rng *r);
		void global_swap_rotation_MC_move(int particle_index, gsl_rng *r);
//		void translate_cluster(double f, int dx, int dy);
//		void rotate_cluster(double f, int xc, int yc);
//		void point_reflect_swap(int xc, int yc, int loc, double f, int already_swapped[], int particle_index[]);
//		Cluster calculations
		void cluster_calculations(int particle_cluster_size_ret[], FILE *hist_fp, FILE *cluster_fp);
		void cluster_calculations_2(double phi_n_hist[], gsl_histogram *w_hist);
		void cluster_calculations_4(FILE *hist_fp, gsl_histogram *w_hist, FILE *cluster_e_fp);
//		Energy calculations
		double calc_energy(void);
		void calc_strain_energy(double strain_energy[]);
//		File writing steps
		void write_checkpoint(void);
		void write_config(FILE *config_fp, FILE *strain_fp, FILE *phase_fp);
		void write_strain_energy(FILE *fp, double strain_energy[]);
};

#endif
