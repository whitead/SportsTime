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
Run_Params* init_run_params(unsigned int iterations, double equilibrium_ratio, double step_size);


/*******************************************************************************
 * A function which prints the wins and satisfies the requirements for
 * being passed to generate_wins.
 *
 *******************************************************************************/
void print_wins_fxn(const unsigned int *wins, const unsigned int *games,
		    unsigned int team_number, unsigned int iterations, const Run_Params *rp, 
		    void *arg);

/*******************************************************************************
 * A function which prints the wins and satisfies the requirements for
 * being passed to sample_model.
 *
 *******************************************************************************/
void print_sample_fxn(const double* alpha, const double* P,
		      unsigned int team_number, unsigned int iterations, const Run_Params *rp,
		      void* arg);

