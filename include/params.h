#ifndef PARAMS_H
#define PARAMS_H

struct xy_params_ {
	int length;
	int num_particles;
	double f;
	double Phi;
	double Sigma;
	double beta;
};

struct mc_params_ {
	int num_steps;
	int step_size;
	long seed;
	int restart;
};

struct auto_params_ {
	int num_steps;
	int step_size;
	long seed;
	int num_saved_clusters;
	int restart;
};
#endif
