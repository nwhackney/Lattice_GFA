#include "init_functions.h"

/* local functions */
int array_index(int j, int i, int length) {
	return (i + j *length);
}

int mod(int a, int b) {
	int ret;

	if (b < 0)
		return mod(-a, -b);
	ret = a % b;
	if (ret < 0)
		ret += b;
	return ret;
}

/* public functions */
void read_params(xy_params_ *xy_params, mc_params_ *mc_params) {
	boost::property_tree::ptree root;
	boost::property_tree::read_json("params.json", root);

	mc_params->step_size = root.get<int>("mc_params.step_size");
	mc_params->num_steps = root.get<int>("mc_params.num_steps");
	mc_params->seed = root.get<long>("mc_params.seed");
	mc_params->restart = root.get<int>("mc_params.restart", 0);
	
	xy_params->length = root.get<int>("xy_params.length");
	xy_params->f = root.get<double>("xy_params.f");
	xy_params->Phi = root.get<double>("xy_params.Phi");
	xy_params->Sigma = root.get<double>("xy_params.Sigma");
	xy_params->beta = root.get<double>("xy_params.beta");

	xy_params->num_particles = (int) (xy_params->length *xy_params->length *xy_params->Phi);
}

void read_auto_params(xy_params_ *xy_params, auto_params_ *auto_params) {
	boost::property_tree::ptree root;
	boost::property_tree::read_json("params.json", root);

	auto_params->step_size = root.get<int>("auto_params.step_size");
	auto_params->num_steps = root.get<int>("auto_params.num_steps");
	auto_params->seed = root.get<long>("auto_params.seed");
	auto_params->restart = root.get<int>("auto_params.restart");
	auto_params->num_saved_clusters = root.get<int>("auto_params.num_saved_clusters");

	xy_params->length = root.get<int>("xy_params.length");
	xy_params->f = root.get<double>("xy_params.f");
	xy_params->Phi = root.get<double>("xy_params.Phi");
	xy_params->Sigma = root.get<double>("xy_params.Sigma");
	xy_params->beta = root.get<double>("xy_params.beta");

	xy_params->num_particles = (int) (xy_params->length *xy_params->length *xy_params->Phi);
}

void init_gauge_table(double a_table[][4], double f, int length) {
	for (int j = 0; j < length; j++) {
		for (int i = 0; i < length; i++) {
			a_table[array_index(j, i, length)][0] = -j *M_PI *f;
			a_table[array_index(j, i, length)][1] = j *M_PI *f;
			a_table[array_index(j, i, length)][2] = i *M_PI *f;
			a_table[array_index(j, i, length)][3] = -i *M_PI *f;
		}
	}
}

void init_neighbor_table(int neighbor_table[][4], int length) {
	for (int j = 0; j < length; j++) {
		for (int i = 0; i < length; i++) {
			neighbor_table[array_index(j, i, length)][0] = array_index(j, mod(i - 1, length), length);
			neighbor_table[array_index(j, i, length)][1] = array_index(j, mod(i + 1, length), length);
			neighbor_table[array_index(j, i, length)][2] = array_index(mod(j - 1, length), i, length);
			neighbor_table[array_index(j, i, length)][3] = array_index(mod(j + 1, length), i, length);
		}
	}
}

void init_fixed_BC_neighbor_table(int neighbor_table[][4], int length) {
	for (int j = 0; j < length; j++) {
		for (int i = 0; i < length; i++) {
			if (i - 1 >= 0)
				neighbor_table[array_index(j, i, length)][0] = array_index(j, i - 1, length);
			else
				neighbor_table[array_index(j, i, length)][0] = length *length;
			if (i + 1 < length)
				neighbor_table[array_index(j, i, length)][1] = array_index(j, i + 1, length);
			else
				neighbor_table[array_index(j, i, length)][1] = length *length;
			if (j - 1 >= 0)
				neighbor_table[array_index(j, i, length)][2] = array_index(j - 1, i, length);
			else
				neighbor_table[array_index(j, i, length)][2] = length *length;
			if (j + 1 < length)
				neighbor_table[array_index(j, i, length)][3] = array_index(j + 1, i, length);
			else
				neighbor_table[array_index(j, i, length)][3] = length *length;
		}
	}
}

void init_thinning_neighbor_table(int thinning_neighbor_table[][8], int length) {
	for (int j = 0; j < length; j++) {
		for (int i = 0; i < length; i++) {
			thinning_neighbor_table[array_index(j, i, length)][0] = array_index(mod(j + 1, length), i, length);
			thinning_neighbor_table[array_index(j, i, length)][1] = array_index(mod(j + 1, length), mod(i + 1, length), length);
			thinning_neighbor_table[array_index(j, i, length)][2] = array_index(j, mod(i + 1, length), length);
			thinning_neighbor_table[array_index(j, i, length)][3] = array_index(mod(j - 1, length), mod(i + 1, length), length);
			thinning_neighbor_table[array_index(j, i, length)][4] = array_index(mod(j - 1, length), i, length);
			thinning_neighbor_table[array_index(j, i, length)][5] = array_index(mod(j - 1, length), mod(i - 1, length), length);
			thinning_neighbor_table[array_index(j, i, length)][6] = array_index(j, mod(i - 1, length), length);
			thinning_neighbor_table[array_index(j, i, length)][7] = array_index(mod(j + 1, length), mod(i - 1, length), length);
		}
	}
}
