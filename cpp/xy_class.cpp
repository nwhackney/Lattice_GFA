#include "xy_class.h"
using namespace std;

XY_conf::XY_conf() {

}

void XY_conf::initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4]) {
	if (xy_params.length <= 0 || xy_params.num_particles <= 0)
		throw invalid_argument("Invalid input parameters.");
	if (xy_params.num_particles > xy_params.length *xy_params.length)
		throw invalid_argument("number of particles too high.");

	XY_conf::num_particles = xy_params.num_particles;
	XY_conf::num_clusters = -1;
	XY_conf::length = xy_params.length;
	XY_conf::neighbor_table = neighbor_table;
	XY_conf::a_table = a_table;
	XY_conf::beta = xy_params.beta;
	K = xy_params.Sigma - 1;
	null_spin.n = 0;
	null_spin.theta = 0;

	lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++) {
		lattice[i] = null_spin;
	}

//	Create dummy index array, shuffle it, then randomly place particles in unique spots with random values of theta.
	int dummy_ind[length *length];
	for (int i = 0; i < length *length; i++) {
		dummy_ind[i] = i;
	}
	
	particle_location.resize(num_particles);
	gsl_ran_shuffle(r, dummy_ind, length *length, sizeof(int));
	for (int s = 0; s < num_particles; s++) {
		int loc = dummy_ind[s];
		particle_location[s] = loc;
		lattice[loc].n = 1;
		lattice[loc].theta = gsl_rng_uniform(r) *2 *M_PI;
	}


}

void XY_conf::square_initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4]) {
	if (xy_params.length <= 0 || xy_params.num_particles <= 0)
		throw invalid_argument("Invalid input parameters.");
	if (xy_params.num_particles > xy_params.length *xy_params.length)
		throw invalid_argument("number of particles too high.");

	XY_conf::num_particles = xy_params.num_particles;
	XY_conf::num_clusters = -1;
	XY_conf::length = xy_params.length;
	XY_conf::neighbor_table = neighbor_table;
	XY_conf::a_table = a_table;
	XY_conf::beta = xy_params.beta;
	K = xy_params.Sigma - 1;
	null_spin.n = 0;
	null_spin.theta = 0;

	lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++) {
		lattice[i] = null_spin;
	}

	particle_location.resize(num_particles);
	int square_length = sqrt(num_particles) + 1;
	int s = 0;
	for (int x = 0; x < square_length; x++) {
		for (int y = 0; y < square_length; y++) {
			int loc = (x + length /2 - square_length /2) + (y + length /2 - square_length /2) *length;
			lattice[loc].n = 1;
			lattice[loc].theta = gsl_rng_uniform(r) *2 *M_PI;
			particle_location[s] = loc;
			s++;
			if (s == num_particles) {
				goto a;
			}
		}
	}
a:;
}

//	For ribbon simulation only. Have to manually set width in main so that it matches the width set here
void XY_conf::ribbon_initialize(xy_params_ xy_params, double ribbon_width, double ribbon_length, gsl_rng *r, int neighbor_table[][4], double a_table[][4]) {
	XY_conf::num_particles = ribbon_length *ribbon_width;
	if (ribbon_length > ribbon_width)
		XY_conf::length = ribbon_length + 100;
	else
		XY_conf::length = ribbon_width + 100;
XY_conf::length = 500;
	XY_conf::neighbor_table = neighbor_table;
	XY_conf::a_table = a_table;
	XY_conf::beta = xy_params.beta;
	K = xy_params.Sigma - 1;
	null_spin.n = 0;
	null_spin.theta = 0;

	lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++) {
		lattice[i] = null_spin;
	}

	particle_location.resize(num_particles);

	int counter = 0;
	for (int j = 0; j < ribbon_width; j++) {
		for (int i = 0; i < ribbon_length; i++) {
			int index = i + j *length;
			particle_location[counter] = index;
			counter++;
			lattice[index].n = 1;
			lattice[index].theta = gsl_rng_uniform(r) *2 *M_PI;
		}
	}
}

void XY_conf::external_initialize(xy_params_ xy_params, gsl_rng *r, int neighbor_table[][4], double a_table[][4], FILE *init_config_fp) {
	if (xy_params.length <= 0 || xy_params.num_particles <= 0)
		throw invalid_argument("Invalid input parameters.");
	if (xy_params.num_particles > xy_params.length *xy_params.length)
		throw invalid_argument("number of particles too high.");
	XY_conf::num_particles = xy_params.num_particles;
	XY_conf::num_clusters = -1;
	XY_conf::length = xy_params.length;
	XY_conf::neighbor_table = neighbor_table;
	XY_conf::a_table = a_table;
	XY_conf::beta = xy_params.beta;
	K = xy_params.Sigma - 1;
	null_spin.n = 0;
	null_spin.theta = 0;

	lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++) {
		lattice[i] = null_spin;
	}

//	Read particle locations and spins from file.
	particle_location.resize(num_particles);
	for (int s = 0; s < num_particles; s++) {
		int loc;
		double phi;
		fscanf(init_config_fp, "%d %lf\n", &loc, &phi);
		particle_location[s] = loc;
		lattice[loc].n = 1;
/*
		while (phi > 2 *M_PI)
			phi -= 2 *M_PI;
		while (phi < 2 *M_PI)
			phi += 2 *M_PI;
*/
		lattice[loc].theta = phi;
	}


}

void XY_conf::cluster_calculations(int particle_cluster_size_in[] = NULL, FILE *hist_fp = NULL, FILE *cluster_fp = NULL) {
//	cluster_size_hist[num_particles], particle_cluster[num_particles], particle_cluster_size[num_particles]
	int particle_cluster[num_particles];
	int cluster_size_hist[num_particles];
	int particle_index[length *length];
	int particle_cluster_size[num_particles];
	for (int l = 0; l < length *length; l++)
		particle_index[l] = -1;
	for (int s = 0; s < num_particles; s++) {
		particle_cluster[s] = -1;
		particle_cluster_size[s] = 0;
		cluster_size_hist[s] = 0;
		particle_index[particle_location[s]] = s;
	}

	int tmp_num_clusters = 0;
	int total_cluster_size = 0;

	for (int s = 0; s < num_particles; s++) {
		if (particle_cluster[s] == -1) {
			list<int> l = { s };
			int cluster_size = 1;
			particle_cluster[s] = tmp_num_clusters;
			for (list<int>::iterator j = l.begin(); j != l.end(); ++j) {
				int loc = particle_location[*j];
				for (int i = 0; i < 4; i++) {
					int n = neighbor_table[loc][i];
					int ss = particle_index[n];
					if (lattice[n].n == 1 && particle_cluster[ss] == -1) {
						l.push_back(ss);
						particle_cluster[ss] = tmp_num_clusters;
						cluster_size++;
					}
				}
			}

			tmp_num_clusters++;
			total_cluster_size += cluster_size;
			cluster_size_hist[cluster_size] += 1;
			for (list<int>::iterator j = l.begin(); j != l.end(); ++j) 
				particle_cluster_size[*j] = cluster_size;
		}
	}

	num_clusters = tmp_num_clusters;
	avg_cluster_size = ((double) total_cluster_size) /(num_clusters);

	if (particle_cluster_size_in != NULL) {
		for (int s = 0; s < num_particles; s++)
			particle_cluster_size_in[s] = particle_cluster_size[s];
	}

	if (hist_fp != NULL) {
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(hist_fp, "%d ", cluster_size_hist[i]);
			}
		}
		fprintf(hist_fp, "\n");
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(hist_fp, "%d ", i);
			}
		}
		fprintf(hist_fp, "\n");
		fflush(hist_fp);
	}

	if (cluster_fp != NULL) {
		for (int s = 0; s < num_particles; s++) {
			fprintf(cluster_fp, "%d ", particle_cluster_size[s]);
		}
		fprintf(cluster_fp, "\n");
	}
}

double euclid_metric(int loc_1, int loc_2, int length) {
	int x_1 = loc_1 % length;
	int x_2 = loc_2 % length;
	int y_1 = loc_1 /length;
	int y_2 = loc_2 /length;

	int delta_x = abs(x_1 - x_2);
	if (abs(delta_x - length) < delta_x)
		delta_x = abs(delta_x - length);
	int delta_y = abs(y_1 - y_2);
	if (abs(delta_y - length) < delta_y)
		delta_y = abs(delta_y - length);

	return sqrt(delta_x *delta_x + delta_y *delta_y);

}

void XY_conf::cluster_calculations_2(double phi_n_hist[] = NULL, gsl_histogram *w_hist = NULL) {
//	cluster_size_hist[num_particles], particle_cluster[num_particles], particle_cluster_size[num_particles]
//	particle_cluster gives the index of the cluster in which each particle resides, -1 means no cluster is assigned (yet).
//	particle_index is lxl array that gives the index of particle at each location. No particle gives an index of -1.
	int particle_cluster[num_particles];
	int cluster_size_hist[num_particles];
	int particle_index[length *length];
	for (int l = 0; l < length *length; l++)
		particle_index[l] = -1;
	for (int s = 0; s < num_particles; s++) {
		particle_cluster[s] = -1;
		cluster_size_hist[s] = 0;
		particle_index[particle_location[s]] = s;
	}

	int tmp_num_clusters = 0;
	for (int s = 0; s < num_particles; s++) {
//		This section builds each cluster
		if (particle_cluster[s] == -1) {
//			If the particle is not already in a cluster (==-1), then build its cluster.
			list<int> cluster = { s };
			int cluster_size = 1;
			particle_cluster[s] = tmp_num_clusters;
			for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
				int loc = particle_location[*j];
				for (int i = 0; i < 4; i++) {
					int n = neighbor_table[loc][i];
					int ss = particle_index[n];
					if (lattice[n].n == 1 && particle_cluster[ss] == -1) {
						cluster.push_back(ss);
						particle_cluster[ss] = tmp_num_clusters;
						cluster_size++;
					}
				}
			}

			tmp_num_clusters++;
			cluster_size_hist[cluster_size] += 1;

//			Measure average distance of each particle in cluster from the edge.
			if (w_hist != NULL) {
				list<int> boundary;
				list<int> bulk;
				double edge_dist = 0.0;
				for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
					int loc = particle_location[*j];
					int lr_sum = 0, ud_sum = 0;
//					i={0,1} = left/right, i = {2,3} = down/up
					for (int i = 0; i < 2; i++) {
						int n = neighbor_table[loc][i];
						lr_sum += lattice[n].n;
					}
					for (int i = 2; i < 4; i++) {
						int n = neighbor_table[loc][i];
						ud_sum += lattice[n].n;
					}
//					If in the bulk, we will calculate the distance to edge in next step
					if (lr_sum + ud_sum == 4) {
						bulk.push_back(*j);
					}
//					If on edge with three neighbors, or an up/down and a left/right neigbor, mean distance to edge = 0.5
					else if (lr_sum && ud_sum) {
						boundary.push_back(*j);
						edge_dist += 0.5;
					}
//					If on edge with up/down or left/right neigbors, mean distance to edge = 0.25 (width = 1)
					else {
						boundary.push_back(*j);
						edge_dist += 0.25;
					}
				}

//				Iterate through bulk spins, find minimum distance to boundary.
				double avg_min_dist = edge_dist;
				for (list<int>::iterator j = bulk.begin(); j != bulk.end(); ++j) {
					double min_dist = num_particles;
					for (list<int>::iterator k = boundary.begin(); k != boundary.end(); ++k) {
						int loc_1 = particle_location[*j];
						int loc_2 = particle_location[*k];
						double tmp = euclid_metric(loc_1, loc_2, length);
						if (tmp < min_dist)
							min_dist = tmp;
					}
					avg_min_dist += min_dist + 0.5;
				}
				avg_min_dist /= cluster.size();
				gsl_histogram_accumulate(w_hist, avg_min_dist, cluster.size());
			}
		}
	}

	num_clusters = tmp_num_clusters;

	if (phi_n_hist != NULL) {
		for (int i = 0; i < num_particles; i++)
			phi_n_hist[i] += (cluster_size_hist[i] *i) /((double) length *length);
	}
}

void XY_conf::cluster_calculations_4(FILE *hist_fp = NULL, gsl_histogram *w_hist = NULL, FILE *cluster_e_fp = NULL) {
//	cluster_size_hist[num_particles], particle_cluster[num_particles], particle_cluster_size[num_particles]
//	particle_cluster gives the index of the cluster in which each particle resides, -1 means no cluster is assigned (yet).
//	particle_index is lxl array that gives the index of particle at each location. No particle gives an index of -1.
	int particle_cluster[num_particles];
	int cluster_size_hist[num_particles];
	int particle_index[length *length];
	for (int l = 0; l < length *length; l++)
		particle_index[l] = -1;
	for (int s = 0; s < num_particles; s++) {
		particle_cluster[s] = -1;
		cluster_size_hist[s] = 0;
		particle_index[particle_location[s]] = s;
	}

	double cluster_energy[num_particles];
	int cluster_count[num_particles];
	for (int i = 0; i < num_particles; i++) {
		cluster_energy[i] = 0.0;
		cluster_count[i] = 0;
	}

	int tmp_num_clusters = 0;

	for (int s = 0; s < num_particles; s++) {
//		This section builds each cluster
		if (particle_cluster[s] == -1) {
//			If the particle is not already in a cluster (==-1), then build its cluster.
			list<int> cluster = { s };
			int cluster_size = 1;
			particle_cluster[s] = tmp_num_clusters;
			for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
				int loc = particle_location[*j];
				for (int i = 0; i < 4; i++) {
					int n = neighbor_table[loc][i];
					int ss = particle_index[n];
					if (lattice[n].n == 1 && particle_cluster[ss] == -1) {
						cluster.push_back(ss);
						particle_cluster[ss] = tmp_num_clusters;
						cluster_size++;
					}
				}
			}

			tmp_num_clusters++;
			cluster_size_hist[cluster_size] += 1;

//			Measure the average energy of each cluster.
			double temp_energy = 0.0;
			for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
				int loc = particle_location[*j];
				temp_energy += calc_local_energy(loc);
			}
			cluster_energy[cluster.size()] += temp_energy /cluster.size();
			cluster_count[cluster.size()] += 1;


//			Measure average distance of each particle in cluster from the edge.
			if (w_hist != NULL) {
				list<int> boundary;
				list<int> bulk;
				double edge_dist = 0.0;
				for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
					int loc = particle_location[*j];
					int lr_sum = 0, ud_sum = 0;
//					i={0,1} = left/right, i = {2,3} = down/up
					for (int i = 0; i < 2; i++) {
						int n = neighbor_table[loc][i];
						lr_sum += lattice[n].n;
					}
					for (int i = 2; i < 4; i++) {
						int n = neighbor_table[loc][i];
						ud_sum += lattice[n].n;
					}
//					If in the bulk, we will calculate the distance to edge in next step
					if (lr_sum + ud_sum == 4) {
						bulk.push_back(*j);
					}
//					If on edge with three neighbors, or an up/down and a left/right neigbor, mean distance to edge = 0.5
					else if (lr_sum && ud_sum) {
						boundary.push_back(*j);
						edge_dist += 0.5;
					}
//					If on edge with up/down or left/right neigbors, mean distance to edge = 0.25 (width = 1)
					else {
						boundary.push_back(*j);
						edge_dist += 0.25;
					}
				}

//				Iterate through bulk spins, find minimum distance to boundary.
				double avg_min_dist = edge_dist;
				for (list<int>::iterator j = bulk.begin(); j != bulk.end(); ++j) {
					double min_dist = num_particles;
					for (list<int>::iterator k = boundary.begin(); k != boundary.end(); ++k) {
						int loc_1 = particle_location[*j];
						int loc_2 = particle_location[*k];
						double tmp = euclid_metric(loc_1, loc_2, length);
						if (tmp < min_dist)
							min_dist = tmp;
					}
					avg_min_dist += min_dist + 0.5;
					double pw=min_dist+0.5; // new
					gsl_histogram_accumulate(w_hist,pw,1.0); // new
				}
//				avg_min_dist *= 4;
				avg_min_dist /= cluster.size();
				//gsl_histogram_accumulate(w_hist, avg_min_dist, cluster.size());
			}
		}
	}

	num_clusters = tmp_num_clusters;

	if (cluster_e_fp != NULL) {
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(cluster_e_fp, "%10.10e ", cluster_energy[i] /cluster_count[i]);
			}
		}
		fprintf(cluster_e_fp, "\n");
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(cluster_e_fp, "%d ", i);
			}
		}
		fprintf(cluster_e_fp, "\n");
		fflush(cluster_e_fp);
	}

	if (hist_fp != NULL) {
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(hist_fp, "%d ", cluster_size_hist[i]);
			}
		}
		fprintf(hist_fp, "\n");
		for (int i = 0; i < num_particles; i++) {
			if (cluster_size_hist[i] > 0) {
				fprintf(hist_fp, "%d ", i);
			}
		}
		fprintf(hist_fp, "\n");
		fflush(hist_fp);
	}
}

void XY_conf::set_tables(int neighbor_table[][4], double a_table[][4]) {
	XY_conf::neighbor_table = neighbor_table;
	XY_conf::a_table = a_table;
}

void XY_conf::set_beta(double beta) {
	XY_conf::beta = beta;
}

void XY_conf::swap_particles(int loc_1, int loc_2) {
	double temp_theta = lattice[loc_1].theta;
	lattice[loc_1].theta = lattice[loc_2].theta;
	lattice[loc_2].theta = temp_theta;
}

void XY_conf::move_particle(int index, int location) {
	int part_loc = particle_location[index];
	lattice[location] = lattice[part_loc];
	particle_location[index] = location;
	lattice[part_loc] =  null_spin;
}

double XY_conf::calc_local_energy(int loc) {
	if (lattice[loc].n == 0)
		return 0.0;
	double energy = 0.0;
	for (int i = 0; i < 4; i++) {
		int n = neighbor_table[loc][i];
		if (lattice[n].n == 1)
			energy += -cos(lattice[loc].theta - lattice[n].theta - a_table[loc][i]) - K;
	}
	return energy;
}

void XY_conf::global_swap_MC_move(int particle_index, gsl_rng *r) {
	int loc = particle_location[particle_index];
	int loc_2 = gsl_rng_uniform_int(r, length *length - 1);
	if (loc_2 >= loc)
		loc_2 += 1;

	while (lattice[loc_2].n == 1) {
		loc_2 = gsl_rng_uniform_int(r, length *length - 1);
			if (loc_2 >= loc)
				loc_2 += 1;
	}

	double energy_i = 0.0, energy_f = 0.0;

	energy_i += calc_local_energy(loc);
	move_particle(particle_index, loc_2);
	energy_f += calc_local_energy(loc_2);
	double weight = exp(-beta *(energy_f - energy_i));
//		If the move is rejected, move particle back
	if (weight < 1 && gsl_rng_uniform(r) > weight)
		move_particle(particle_index, loc);
}


double random_theta(double old_theta, gsl_rng *r, double beta) {
	double width = 0.1*exp(1.5 /beta);
	double new_theta = (gsl_rng_uniform(r) - 0.5) *width + old_theta;
	return new_theta;

//	In principle, we should be taking the modulo as below to avoid a
//	random walk that eventually diverges and causes precision issues.
//	In practice, it appears to be a non-issue since the steps are

/*	
	while (new_theta < 0)
		new_theta += 2 *M_PI;
	while (new_theta > 2 *M_PI)
		new_theta -= 2 *M_PI;
	return new_theta;
*/
}

void XY_conf::rotation_MC_move(int particle_index, gsl_rng *r) {
	int loc = particle_location[particle_index];
	double old_theta = lattice[loc].theta;
	double new_theta = random_theta(lattice[loc].theta, r, beta);

	double energy_i = calc_local_energy(loc);
	lattice[loc].theta = new_theta;
	double energy_f = calc_local_energy(loc);
	double weight = exp(-beta *(energy_f - energy_i));
//	If move is rejected, unrotate the spin
	if (weight < 1 && gsl_rng_uniform(r) > weight)
		lattice[loc].theta = old_theta;
}

void XY_conf::global_swap_rotation_MC_move(int particle_index, gsl_rng *r) {
	int loc = particle_location[particle_index];
	int loc_2 = gsl_rng_uniform_int(r, length *length - 1);
	if (loc_2 >= loc)
		loc_2 += 1;

	while (lattice[loc_2].n == 1) {
		loc_2 = gsl_rng_uniform_int(r, length *length - 1);
			if (loc_2 >= loc)
				loc_2 += 1;
	}

	double energy_i = 0.0, energy_f = 0.0;

	double old_theta = lattice[loc].theta;
	double new_theta = gsl_rng_uniform(r) *2 *M_PI;

	energy_i += calc_local_energy(loc);

	lattice[loc].theta = new_theta;
	move_particle(particle_index, loc_2);

	energy_f += calc_local_energy(loc_2);
	double weight = exp(-beta *(energy_f - energy_i));
//		If the move is rejected, move particle back and unrotate
	if (weight < 1 && gsl_rng_uniform(r) > weight) {
		move_particle(particle_index, loc);
		lattice[loc].theta = old_theta;
	}
}

void XY_conf::write_config(FILE *config_fp = NULL, FILE *strain_fp = NULL, FILE *phase_fp = NULL) {
	if (config_fp == NULL) {
		for (int i = 0; i < num_particles; i++) {
			int loc = particle_location[i];
			printf("%d ", loc);
			printf("%5.5e ", lattice[loc].theta);
			printf("%5.5e ", calc_local_energy(loc) /4);
		}
		printf("\n");
		fflush(stdout);
	}
	else {
		for (int i = 0; i < num_particles; i++) {
			int loc = particle_location[i];
			fprintf(config_fp, "%d ", particle_location[i]);
			fprintf(strain_fp, "%5.5e ", calc_local_energy(loc) /4);
			fprintf(phase_fp, "%5.5e ", lattice[loc].theta);
		}
		fprintf(config_fp, "\n");
		fprintf(strain_fp, "\n");
		fprintf(phase_fp, "\n");
		fflush(config_fp);
		fflush(strain_fp);
		fflush(phase_fp);
	}
}

void XY_conf::write_strain_energy(FILE *fp, double strain_energy[]) {
	if (fp == NULL) {
		for (int i = 0; i < num_particles; i++)
			printf("%15.15e ", strain_energy[i]);
		printf("\n");
		fflush(stdout);
	}
	else {
		for (int i = 0; i < num_particles; i++)
			fprintf(fp, "%15.15e ", strain_energy[i]);
		fprintf(fp, "\n");
		fflush(fp);
	}	
}


double XY_conf::calc_energy(void) {

	double energy = 0.0;
	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		energy += calc_local_energy(loc);
	}
	return energy /2.0;
}

void XY_conf::local_swap_MC_move(int particle_index, gsl_rng *r) {
	int loc = particle_location[particle_index];
	int direction = gsl_rng_uniform_int(r, 4);
	int loc_2 = neighbor_table[loc][direction];

	double energy_i = 0.0, energy_f = 0.0;

	if (lattice[loc_2].n == 0 && loc_2 != length *length) {
		energy_i += calc_local_energy(loc);
		move_particle(particle_index, loc_2);
		energy_f += calc_local_energy(loc_2);
		double weight = exp(-beta *(energy_f - energy_i));
//		If the move is rejected, move particle back
		if (weight < 1 && gsl_rng_uniform(r) > weight)
			move_particle(particle_index, loc);
	}
}

void XY_conf::calc_strain_energy(double strain_energy[]) {
	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		double tmp_E = 0.0;
		int nn = 0;
		for (int i = 0; i < 4; i++) {
			int n = neighbor_table[loc][i];
			if (lattice[n].n == 1) {
				tmp_E += -cos(lattice[loc].theta - lattice[n].theta - a_table[loc][i]) - K;
				nn++;
			}
		}
		strain_energy[s] = tmp_E /nn;
	}	
}

void XY_conf::write_checkpoint(void) {
	rename("checkpoint.dat", "checkpoint.bak");
	FILE *fp = fopen("checkpoint.dat", "w");
	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		fprintf(fp, "%d %20.20e\n", loc, lattice[loc].theta);
	}
	fclose(fp);
}

double XY_conf::get_avg_cluster_size(void) {
	return avg_cluster_size;
}

/*
void XY_conf::find_cluster(int particle_index) {
	int cluster_lattice[length *length];
	for (int i = 0; i < length *length; i++)
		cluster_lattice[i] = 0;
	
	int s = particle_index;
	int loc = particle_location[s];
	cluster_lattice[loc] = 1;
	list<int> cluster = { loc };

	for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) {
		int loc = *j;
		for (int i = 0; i < 4; i++) {
			int n = neighbor_table[loc][i];
			if (lattice[n].n == 1 && cluster_lattice[n] == 0) {
				cluster.push_back(ss);
				particle_cluster[ss] = 0;
			}
		}
	}
	
	int min_x = length, max_x, min_y, max_y;
	for (list<int>::iterator j = cluster.begin(); j != cluster.end(); ++j) { 
		
	}

}
*/

void get_coords(int coord[2], int ind, int length) {
	coord[0] = ind % length;
	coord[1] = (ind - coord[0]) /length;
}

double phase_shift(int dx, int dy, int i, int j, double f, double theta) {
	double tmp = theta + M_PI *f *(-dy *i + dx *j);
	while (tmp >= 2 *M_PI)
		tmp -= 2 *M_PI;
	while (tmp < 0)
		tmp += 2 *M_PI;

	return tmp;
}

int loc_mod(int a, int b) {
	int ret;

	if (b < 0)
		return loc_mod(-a, -b);
	ret = a % b;
	if (ret < 0)
		ret += b;
	return ret;
}

typedef struct {
	int loc;
	double delta_bond_E;
	int dir;
} Swap;

void XY_conf::calc_bond_energy(int loc, double energy[4]) {
	if (lattice[loc].n == 0)
		return;
	for (int i = 0; i < 4; i++) {
		int n = neighbor_table[loc][i];
		if (lattice[n].n == 1)
			energy[i] += -cos(lattice[loc].theta - lattice[n].theta - a_table[loc][i]) - K;
	}
	return;
}

void XY_conf::classic_wolff(int index, gsl_rng *r, double f) {

	int already_rotated[length *length];
	for (int i = 0; i < length *length; i++) {
		already_rotated[i] = 0;
	}
	
	int loc = particle_location[index];
	list<Swap> check_for_rotation = { };
	
	double dphi = (gsl_rng_uniform(r) - 0.5) *2 *M_PI;
	
	double E_bond_i[4], E_bond_f[4];
	for (int i = 0; i < 4; i++) {
		E_bond_i[i] = 0.0;
		E_bond_f[i] = 0.0;
	}
	calc_bond_energy(loc, E_bond_i);
	lattice[loc].theta += dphi;
	calc_bond_energy(loc, E_bond_f);
	
	for (int i = 0; i < 4; i++) {
		int n = neighbor_table[loc][i];
		if (lattice[n].n == 1) {
			check_for_rotation.push_back({n, E_bond_f[i] - E_bond_i[i], 0});
		}
	}
	already_rotated[loc] = 1;
	
	for (list<Swap>::iterator j = check_for_rotation.begin(); j != check_for_rotation.end(); ++j) {
		loc = j->loc;
		double weight = 1 - exp(-beta *(j->delta_bond_E));
		if (already_rotated[loc] == 0 && (weight >= 1 || gsl_rng_uniform(r) < weight)) {
			for (int i = 0; i < 4; i++) {
				E_bond_i[i] = 0.0;
				E_bond_f[i] = 0.0;
			}
			calc_bond_energy(loc, E_bond_i);
			lattice[loc].theta += dphi;
			calc_bond_energy(loc, E_bond_f);
			
			for (int i = 0; i < 4; i++) {
				int n = neighbor_table[loc][i];
				if (lattice[n].n == 1 && already_rotated[n] == 0) {
					check_for_rotation.push_back({n, E_bond_f[i] - E_bond_i[i], 0});
				}
			}
			already_rotated[loc] = 1;
		}
	}


}

/*
void XY_conf::translate_cluster(double f, int dx, int dy) {
	std::vector<XY_spin> tmp_lattice;
	tmp_lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++)
		tmp_lattice[i] = null_spin;
	
	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		XY_spin tmp_spin = lattice[loc];
		int coord[2];
		get_coords(coord, loc, length);
		tmp_spin.theta = phase_shift(dx, dy, coord[0], coord[1], f, tmp_spin.theta);
		coord[0] = loc_mod(coord[0] + dx, length);
		coord[1] = loc_mod(coord[1] + dy, length);
		int new_loc = coord[0] + coord[1]*length;
		tmp_lattice[new_loc] = tmp_spin;
		particle_location[s] = new_loc;
	}
	for (int i = 0; i < length *length + 1; i++)
		lattice[i] = tmp_lattice[i];
}
*/

/*
void XY_conf::rotate_cluster(double f, int xc, int yc) {
	std::vector<XY_spin> tmp_lattice;
	tmp_lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++)
		tmp_lattice[i] = null_spin;

	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		XY_spin tmp_spin = lattice[loc];
		int coord[2];
		get_coords(coord, loc, length);
		int new_coord[2];
		new_coord[0] = loc_mod(2 *xc - coord[0], length);
		new_coord[1] = loc_mod(2 *yc - coord[1], length);
		int dx = new_coord[0] - coord[0];
		int dy = new_coord[1] - coord[1];
		
		tmp_spin.theta = phase_shift(dx, dy, coord[0], coord[1], f, tmp_spin.theta);
		int new_loc = new_coord[0] + new_coord[1]*length;
		tmp_lattice[new_loc] = tmp_spin;
		particle_location[s] = new_loc;
	}
	for (int i = 0; i < length *length + 1; i++)
		lattice[i] = tmp_lattice[i];
}
*/

/*
void XY_conf::rotate_cluster(double f, int xc, int yc) {
	std::vector<XY_spin> tmp_lattice;
	tmp_lattice.resize(length *length + 1);
	for (int i = 0; i < length *length + 1; i++)
		tmp_lattice[i] = null_spin;

	for (int s = 0; s < num_particles; s++) {
		int loc = particle_location[s];
		XY_spin tmp_spin = lattice[loc];
		int coord[2];
		get_coords(coord, loc, length);
		int new_coord[2];
		new_coord[0] = loc_mod(2 *xc - coord[0], length);
		new_coord[1] = loc_mod(2 *yc - coord[1], length);
		int dx = 2 *(coord[0] - xc);
		int dy = 2 *(coord[1] - yc);
		
		tmp_spin.theta = phase_shift(dx, dy, coord[0], coord[1], f, tmp_spin.theta);
		int new_loc = new_coord[0] + new_coord[1]*length;
		tmp_lattice[new_loc] = tmp_spin;
		particle_location[s] = new_loc;
	}
	for (int i = 0; i < length *length + 1; i++)
		lattice[i] = tmp_lattice[i];
}

void XY_conf::point_reflect_swap(int xc, int yc, int loc, double f, int already_swapped[], int particle_index[]) {
	int coord[2], new_coord[2];
	get_coords(coord, loc, length);
	new_coord[0] = loc_mod(2 *xc - coord[0], length);
	new_coord[1] = loc_mod(2 *yc - coord[1], length);
	int new_loc = new_coord[0] + new_coord[1] *length;
	int dx = 2 *(coord[0] - xc);
	int dy = 2 *(coord[1] - yc);
	lattice[loc].theta = phase_shift(dx, dy, coord[0], coord[1], f, lattice[loc].theta);
	already_swapped[loc] = 1;
	int index = particle_index[loc];
	if (lattice[new_loc].n == 1) {
		int new_index = particle_index[new_loc];
		lattice[new_loc].theta =  phase_shift(-dx, -dy, new_coord[0], new_coord[1], f, lattice[new_loc].theta);
		XY_spin tmp_spin = lattice[new_loc];
		lattice[new_loc] = lattice[loc];
		lattice[loc] = tmp_spin;
		already_swapped[new_loc] = 1;
		particle_location[index] = new_loc;
		particle_location[new_index] = loc;
		particle_index[loc] = new_index;
		particle_index[new_loc] = index;
	}
	else {
		lattice[new_loc] = lattice[loc];
		lattice[loc] = null_spin;
		particle_location[index] = new_loc;
		particle_index[loc] = -1;
		particle_index[new_loc] = index;
	}
}
*/

void XY_conf::point_reflect_wolff(int index, gsl_rng *r, double f) {
	int already_swapped[length *length], particle_index[length *length];
	for (int i = 0; i < length *length; i++) {
		already_swapped[i] = 0;
		particle_index[i] = -1;
	}
	
	for (int s = 0; s < num_particles; s++) {
		particle_index[particle_location[s]] = s;
	}
	
	int xc = gsl_rng_uniform_int(r, length);
	int yc = gsl_rng_uniform_int(r, length);

	double dphi = (gsl_rng_uniform(r) - 0.5) *2 *M_PI;
	dphi = 0.0;

	int loc = particle_location[index];
	int coord[2];
	get_coords(coord, loc, length);
	int new_coord[2];
	new_coord[0] = loc_mod(2 *xc - coord[0], length);
	new_coord[1] = loc_mod(2 *yc - coord[1], length);
	int new_loc = new_coord[0] + new_coord[1]*length;
	int dx = 2 *(coord[0] - xc);
	int dy = 2 *(coord[1] - yc);
	
	double E_bond_loc_i[4], E_bond_loc_f[4];
	double E_bond_new_loc_i[4], E_bond_new_loc_f[4];
	
	list<Swap> check_for_swap = { };
	for (int i = 0; i < 4; i++) {
		E_bond_loc_i[i] = 0.0;
		E_bond_loc_f[i] = 0.0;
		E_bond_new_loc_i[i] = 0.0;
		E_bond_new_loc_f[i] = 0.0;
	}
	calc_bond_energy(new_loc, E_bond_new_loc_i);
	calc_bond_energy(loc, E_bond_loc_i);

	lattice[loc].theta = phase_shift(dx, dy, coord[0], coord[1], f, lattice[loc].theta) + dphi;
	if (lattice[new_loc].n == 1) {
		lattice[new_loc].theta = phase_shift(-dx, -dy, new_coord[0], new_coord[1], f, lattice[new_loc].theta) + dphi;
		XY_spin tmp_spin = lattice[new_loc];
		lattice[new_loc] = lattice[loc];
		lattice[loc] = tmp_spin;
		already_swapped[loc] = 1;
		already_swapped[new_loc] = 1;
	}
	else {
		lattice[new_loc] = lattice[loc];
		lattice[loc] = null_spin;
		already_swapped[new_loc] = 1;
	}
	calc_bond_energy(new_loc, E_bond_new_loc_f);
	calc_bond_energy(loc, E_bond_loc_f);

	for (int i = 0; i < 4; i++) {
		int n = neighbor_table[loc][i];
		if (lattice[n].n == 1) {
			check_for_swap.push_back({n, E_bond_loc_f[i] - E_bond_loc_i[i], 0});
		}
		n = neighbor_table[new_loc][i];
		if (lattice[n].n == 1) {
			check_for_swap.push_back({n, E_bond_new_loc_f[i] - E_bond_new_loc_i[i], 0});
		}
	}
	
	for (list<Swap>::iterator j = check_for_swap.begin(); j != check_for_swap.end(); ++j) {
		loc = j->loc;
		double weight = 1 - exp(-beta *(j->delta_bond_E));
		if (weight < 0)
			weight = 0;
		if (lattice[loc].n == 1 && already_swapped[loc] == 0 && (weight >= 1 || gsl_rng_uniform(r) < weight)) {
			get_coords(coord, loc, length);
			new_coord[0] = loc_mod(2 *xc - coord[0], length);
			new_coord[1] = loc_mod(2 *yc - coord[1], length);
			new_loc = new_coord[0] + new_coord[1]*length;
			dx = 2 *(coord[0] - xc);
			dy = 2 *(coord[1] - yc);
			
			for (int i = 0; i < 4; i++) {
				E_bond_loc_i[i] = 0.0;
				E_bond_loc_f[i] = 0.0;
				E_bond_new_loc_i[i] = 0.0;
				E_bond_new_loc_f[i] = 0.0;
			}
			calc_bond_energy(new_loc, E_bond_new_loc_i);
			calc_bond_energy(loc, E_bond_loc_i);

			lattice[loc].theta = phase_shift(dx, dy, coord[0], coord[1], f, lattice[loc].theta) + dphi;
			if (lattice[new_loc].n == 1) {
				lattice[new_loc].theta = phase_shift(-dx, -dy, new_coord[0], new_coord[1], f, lattice[new_loc].theta) + dphi;
				XY_spin tmp_spin = lattice[new_loc];
				lattice[new_loc] = lattice[loc];
				lattice[loc] = tmp_spin;
				already_swapped[loc] = 1;
				already_swapped[new_loc] = 1;
			}
			else {
				lattice[new_loc] = lattice[loc];
				lattice[loc] = null_spin;
				already_swapped[new_loc] = 1;
			}
			calc_bond_energy(new_loc, E_bond_new_loc_f);
			calc_bond_energy(loc, E_bond_loc_f);

			for (int i = 0; i < 4; i++) {
				int n = neighbor_table[loc][i];
				if (lattice[n].n == 1) {
					check_for_swap.push_back({n, E_bond_loc_f[i] - E_bond_loc_i[i], 0});
				}
				n = neighbor_table[new_loc][i];
				if (lattice[n].n == 1) {
					check_for_swap.push_back({n, E_bond_new_loc_f[i] - E_bond_new_loc_i[i], 0});
				}
			}
		}
	}

	int c = 0;
	for (int i = 0; i < length *length; i++) {
		if (lattice[i].n == 1) {
			particle_location[c++] = i;
		}
	}
	if (c != num_particles) {
		printf("%d\n", c);
		printf("error!\n");
		exit(1);
	}

}
















