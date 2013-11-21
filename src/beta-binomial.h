/*******************************************************************************
 * Sample passed function on model parameters given the number of
 * wins, games and team number.
 * 
 * Required parameters: fxn is a three argument function taking alpha,
 * P, and the number of teams. It will be called on each realization
 * of a model parameter set. Run_Params rp should contain information on 
 * how many iterations, convergance criteria, etc.
 *
 * Optional parameters: alpha and P. Pass null to initialize to 1 and 0.5, respectively
 *
 * Returns: acceptance ratio of alpha moves, useful for tuning step size
 *
*******************************************************************************/
double sample_model(const unsigned int *wins, const unsigned int *games, 
		  double *alpha, double *P,
		  void *fxn,
		  const unsigned int team_number, const Run_Params *rp);


/*******************************************************************************
 * Generate a win matrix given the game number matrix, the alpha
 * model parameters. 
 *
 * Required parameters: alpha, and games. Wins is where the result
 * of the calculation will be stored. May be null.
 *
 * Optional parameters: fxn contains an optional function which will
 * be executed on each realiziation of the win matrix.
 *
 * Returns: void
 *******************************************************************************/
void generate_wins(unsigned int *wins, unsigned int *games,
		   const double *alpha,
		   void *fxn,
		   const unsigned int team_number, const Run_Params *rc);


/*******************************************************************************
 * Create and initialize a Run_Param structure.
 *
 * Required parameters: 
 *
 * Returns: Run_Params 
 *
 *******************************************************************************/
Run_Params* init_run_params(unsigned int iterations, double discard_ratio);


