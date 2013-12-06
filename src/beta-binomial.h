#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdbool.h>

/*******************************************************************************
 * The Run Pamareters
 *
 * See init_run_params
 *
*******************************************************************************/

typedef struct Run_Params_s {
  unsigned int iterations;
  double equilibrium_ratio;
  double step_size;
  double max_alpha;
  gsl_rng* rng;
} Run_Params;


/*******************************************************************************
 * Sample passed function on model parameters given the number of
 * wins, games and team number.
 * 
 * Required parameters: fxn is a four argument function taking alpha,
 * P, the number of teams, and the run_parameters. It will be called
 * on each realization of a model parameter set. Run_Params rp should
 * contain information on how many iterations, convergance criteria,
 * etc.
 *
 * Optional parameters: alpha and P. Pass null to initialize to 1 and 0.5, respectively
 *
 * Returns: acceptance ratio of alpha moves, useful for tuning step size
 *
*******************************************************************************/
double sample_model(const unsigned int *wins, const unsigned int *games, 
		  double *alpha, double *P,
		    void (*fxn)(const double*, 
				const double*, 
				unsigned int, 
				unsigned int, 
				const Run_Params*,
				void*),
		    void *fxn_args,
		  const unsigned int team_number, const Run_Params *rp);


/*******************************************************************************
 * Generate a win matrix given the game number matrix, the alpha
 * model parameters. 
 *
 * Required parameters: alpha, and games. Wins is where the result
 * of the calculation will be stored. May be null.
 *
 * Optional parameters: fxn contains an optional function which will
 * be executed on each realiziation of the win matrix. The elements
 * passed to fxn in order are: wins, games, team_number, run_parameters
 *
 * Returns: void
 *******************************************************************************/
void generate_wins(unsigned int *wins, const unsigned int *games,
		   const double *alpha,
		   void (*fxn)(const unsigned int*, 
			       const unsigned int*, 
			       unsigned int, 
			       unsigned int,
			       const Run_Params*,
			       void*),
                   		    void *fxn_args,
		   const unsigned int team_number, const Run_Params *rp);


/*******************************************************************************
 * Create and initialize a Run_Param structure.
 *
 * Required parameters: Number of iterations to run and the ratio of
 * those iterations that are ignored due to equilibration of the model. Step_size
 * is the size of steps to take for the alpha monte carlo moves.
 *
 * Returns: Run_Params 
 *
 *******************************************************************************/
Run_Params* init_run_params(unsigned int iterations, double equilibrium_ratio, double step_size, double max_alpha);


/*******************************************************************************
 * A function which prints the wins and satisfies the requirements for
 * being passed to generate_wins.
 *
 *******************************************************************************/
void print_wins_fxn(const unsigned int *wins, const unsigned int *games,
		    unsigned int team_number, 
		    unsigned int iterations, 
		    const Run_Params *rp, 
		    void *arg);

/*******************************************************************************
 * A function which histograms and prints alpha and satisfies the
 * requirements for being passed to sample_model.
 *
 *******************************************************************************/
void accum_sample_fxn(const double* alpha, const double* P,
		      unsigned int team_number, 
		      unsigned int iterations, const Run_Params *rp,
		      void* arg);
/*******************************************************************************
 * This function can be passed to sample_model and will generate a
 * matrix of wins then histogram the result.
 *
 * *******************************************************************************/
void sample_wins_wrapper(const double* alpha, const double* P,
			 unsigned int team_number, 
			 unsigned int iterations, const Run_Params *rp,
			 void* win_params);

typedef struct Matrix_Histogram_s {
  unsigned int dimension;
  unsigned int nbins;
  double width;
  unsigned long int *hist; 
} Matrix_Histogram;

Matrix_Histogram* build_matrix_histogram(unsigned int dimension, 
					 unsigned int nbins);
					 

typedef struct Wins_Parameters_s {
  Matrix_Histogram* mhist;
  unsigned int* win_sample;
  unsigned int* games;
} Win_Parameters;

unsigned int* load_uint_matrix(char* filename, unsigned int nrow, unsigned int ncol);
