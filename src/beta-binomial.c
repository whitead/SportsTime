#include "beta-binomial.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>


inline
double alpha_lprob(double alpha_trial, double alpha, double alpha_other, double p,
		   unsigned int wins, unsigned int games) {
  return (gsl_sf_lngamma(alpha_other + alpha_trial) - gsl_sf_lngamma(alpha_other + alpha)) + 
    (gsl_sf_lngamma(alpha) - gsl_sf_lngamma(alpha_trial)) + 
    (alpha_trial - alpha) * log(p);
}

inline
double alpha_lprob_stirling(double alpha_trial, double alpha, double alpha_other, double p,
		   unsigned int wins, unsigned int games) {
  return (alpha_trial - alpha) * (log(alpha_other) + log(p)) + alpha * (log(alpha) - 1) - alpha_trial * (log(alpha_trial) - 1);
}

inline
unsigned int symm_index(unsigned int i, unsigned int j, unsigned int n) {
  //       size of pairs        * number of pairs     +
  //((n - 1) + (n - 1 - i + 1)) * (i) / 2 +
  //   depth into last row - 1 because sum formula above starts at 1
  // j - i - 1
  //simplify that polynomial
  return (2 * n * i - i *i - i) / 2 + j - i - 1;
}

inline
bool lMH(double lp, gsl_rng* rng){
  if(lp >= 0)
    return true;

  return gsl_rng_uniform(rng) < exp(lp);
}

int test() {
  
  unsigned int team_number = 5;
  unsigned int i;
  double *alpha = (double*) malloc(sizeof(double) * team_number);
  for(i = 0; i < team_number; i++)
    alpha[i] = 10 - i * 2;
  Run_Params* rp = init_run_params(1, 0, 0.2);
  unsigned int *games = (unsigned int*) malloc(sizeof(unsigned int) * team_number * (team_number - 1) / 2);
  for(i = 0; i < team_number * (team_number - 1) / 2; i++)
    games[i] = i; 

  unsigned int *wins = (unsigned int*) calloc(team_number * (team_number - 1) / 2, sizeof(unsigned int));
  double *alpha_fit = NULL;
  double *P = NULL;
  
  generate_wins(wins, games, alpha, print_wins_fxn, NULL, team_number, rp);

  free(rp);
  rp = init_run_params(100000, 0.5, 20);
  for(i = 0; i < team_number; i++)
    alpha[i] = 0;

  printf("Acceptance: %g\n", sample_model(wins, games, alpha_fit, P, print_sample_fxn, alpha, team_number, rp));
  
  return 0;
}

int main(int argc, char* argv[]) {


  unsigned int team_number;
  sscanf(argv[1], "%ud", &team_number);
  unsigned int* square_wins = load_uint_matrix(argv[2], team_number, team_number);

  unsigned int i,j;
  unsigned int *games = (unsigned int*) malloc(sizeof(unsigned int) * team_number * (team_number - 1) / 2);
  unsigned int *wins = (unsigned int*) calloc(team_number * (team_number - 1) / 2, sizeof(unsigned int));

  for(i = 0; i < team_number; i++) {
    for(j = i + 1; j < team_number; j++) {
      games[symm_index(i,j,team_number)] = square_wins[i * team_number + j];
      games[symm_index(i,j,team_number)] += square_wins[j * team_number + i];
      wins[symm_index(i,j,team_number)] = square_wins[i * team_number + j];
    }      
  }
    

  Run_Params* rp = init_run_params(200000, 0.5, 0.3);
  //  Run_Params* rp = init_run_params(2, 0, 0.3);
  double *alpha = (double*) calloc(team_number, sizeof(double));

  //  print_wins_fxn(wins, games, team_number, 0, NULL, NULL);

  printf("Acceptance: %g\n", sample_model(wins, games, NULL, NULL, print_sample_fxn, alpha, team_number, rp));

  
  return 0;
}

double sample_model(const unsigned int *wins, const unsigned int *games, 
		  double *alpha, double *P,
		    void (*fxn)(const double*, 
				const double*, 
				unsigned int, 
				unsigned int, 
				const Run_Params*,
				void*),
		    void *fxn_args,
		    const unsigned int team_number, const Run_Params *rp) {
  
  unsigned int i,j,k, iters, acceptance, index;
  double alpha_trial, rel_lprob, ratio;
  const unsigned int max_alpha = 5;


  //initalize parameters as needed
  if(!alpha) {
    alpha = (double*) malloc(sizeof(double) * team_number);
    for(i = 0; i < team_number; i++)
      alpha[i] = 1;
  }
  
  if(!P)
    P = (double*) malloc(sizeof(double) * team_number * (team_number - 1) / 2);

  for(iters = 0 ;iters < rp->iterations; iters++) {

    //sample new pij
    for(i = index = 0; i < team_number; i++) {
      for(j = i + 1; j < team_number; j++) {	
	P[index] = gsl_ran_beta(rp->rng, alpha[i] + wins[index], alpha[j] + games[index] - wins[index]);
	index++;
      }
    }
    
    //sample new alpha_ij with replacement
    acceptance = 0;
    for(k = 0; k < team_number; k++) {
      //choose index
      i = gsl_rng_uniform_int(rp->rng, team_number);
      //perturb alpha
      alpha_trial = alpha[i] + gsl_rng_uniform(rp->rng) * 2 * rp->step_size - rp->step_size;
      //project
      if(alpha_trial < 1)
	alpha_trial = 1.;
      if(alpha_trial > max_alpha)
	alpha_trial = max_alpha;
      //accumate terms
      rel_lprob = 0;
      for(j = i + 1; j < team_number; j++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[i], 
				 alpha[j], 
				 P[symm_index(i,j,team_number)],
				 wins[symm_index(i,j,team_number)], 
				 games[symm_index(i,j,team_number)]);

      for(j = 0; j < i; j++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[i], 
				 alpha[j], 
				 1 - P[symm_index(i,j,team_number)],
				 games[symm_index(i,j,team_number)] - wins[symm_index(i,j,team_number)], 
				 games[symm_index(i,j,team_number)]);
      //Metrpolis-Hastings 
      if(lMH(rel_lprob, rp->rng)) {
	alpha[i] = alpha_trial;	
	acceptance += 1;
      }      
    }

    //calculate acceptance ratio
    ratio = ratio * (iters) / (iters + 1) + ((double) acceptance / team_number) / (iters + 1);
    
    //run the function
    if(fxn && iters > rp->iterations * rp->equilibrium_ratio)
      fxn(alpha, P, team_number, iters - rp->iterations * rp->equilibrium_ratio, rp, fxn_args);

  }
  return ratio;
}

void generate_wins(unsigned int *wins, const unsigned int *games,
		   const double *alpha,
		   void (*fxn)(const unsigned int*, 
			       const unsigned int*, 
			       unsigned int, 
			       unsigned int,
			       const Run_Params*,
			       void*),
		    void *fxn_args,
		   const unsigned int team_number, const Run_Params *rp){

  unsigned int iters, i, j, index;
  double pij;
  
  //initialize parameters as needed
  if(!wins)
    wins = (unsigned int*) calloc(team_number * (team_number - 1) / 2, sizeof(unsigned int));
  
  //sample realization
  for(iters = 0 ;iters < rp->iterations; iters++) {
    //sample
    for(i = index = 0; i < team_number; i++) {
      for(j = i + 1; j < team_number; j++) {
	pij = gsl_ran_beta(rp->rng, alpha[i],alpha[j]);
	wins[index] = gsl_ran_binomial(rp->rng, pij, games[index]);
	index++;
      }
    }
    fxn(wins, games, team_number, iters, rp, fxn_args);
  }
}

Run_Params* init_run_params(unsigned int iterations, double equilibrium_ratio, double step_size) {
  Run_Params* rp = (Run_Params*) malloc(sizeof(Run_Params));
  rp->iterations = iterations;
  rp->equilibrium_ratio = equilibrium_ratio;
  rp->step_size = step_size;

  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rp->rng = gsl_rng_alloc(T);

  return rp;
}


void log_array(FILE* file, double* array, unsigned int n_cols, unsigned int n_rows, bool do_sum) {
  
  if(file == NULL)
    return;

  unsigned int i,j;
  double sum = 0;

  for(i = 0; i < n_rows; i++) {
    if(do_sum)
      sum = 0;
    for(j = 0; j < n_cols; j++) {
      if(do_sum)
	sum += array[i * n_cols + j];
      fprintf(file, "%12g ", array[i * n_cols + j]); 
    }
    if(do_sum)
      fprintf(file, "%12g\n", sum);
    else
      fprintf(file, "\n");
  }

  fflush(file);
  
}


void print_wins_fxn(const unsigned int *wins, const unsigned int *games,
		    unsigned int team_number, unsigned int iterations, const Run_Params *rp, void* arg) {

  unsigned int i,j;
  unsigned int sum = 0;  

  for(i = 0; i < team_number; i++) {
    sum = 0;
    for(j = 0; j < i; j++) {
      sum += (games[symm_index(j,i,team_number)] - wins[symm_index(j,i,team_number)]);
      fprintf(stdout, "%12u ", games[symm_index(j,i,team_number)] - wins[symm_index(j,i,team_number)]); 
    }
    fprintf(stdout, "%12u ", 0);
    for(j = i + 1; j < team_number; j++) {
      sum += wins[symm_index(i,j,team_number)];
      fprintf(stdout, "%12u ",wins[symm_index(i,j,team_number)]);
    }
    fprintf(stdout, "%12u\n", sum);
  }
  fprintf(stdout, "\n");
}

void print_sample_fxn(const double* alpha, const double* P,
		      unsigned int team_number, unsigned int iteration, const Run_Params *rp, void* arg) {

  double* running_alpha = (double*) arg;
  unsigned int i,j;

  for(i = 0; i < team_number; i++)
    running_alpha[i] = running_alpha[i] * (iteration - 1) / iteration +
      alpha[i] / iteration;


  if(iteration % 1000 != 0)
    return;

  printf("P:\n");
  for(i = 0; i < 80; i++)
    printf("-");
  printf("\n");
  for(i = 0; i < team_number; i++) {
    for(j = 0; j < i; j++) {
      fprintf(stdout, "%12g ", 1 - P[symm_index(j,i,team_number)]); 
    }
    fprintf(stdout, "%12g ", 0.);
    for(j = i + 1; j < team_number; j++) {
      fprintf(stdout, "%12g ", P[symm_index(i,j,team_number)]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
  printf("Alpha:\n");
  for(i = 0; i < 80; i++)
    printf("-");
  printf("\n");
  for(i = 0; i < team_number; i++)
    printf("%12g ", running_alpha[i]);
  printf("\n");

}

unsigned int* load_uint_matrix(char* filename, unsigned int nrow, unsigned int ncol) {

  FILE *mfile;
  mfile = fopen(filename, "r");
  if(mfile != NULL) {

    unsigned int i, j;

    unsigned int *matrix =  (unsigned int*) malloc(sizeof(unsigned int) * nrow * ncol);
    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%ud", &(matrix[i*ncol + j])) == 0) {
	  fprintf(stderr, "Incorrect number of rows or columns"
		  "at i = %d, and j=%d, nrow=%d, ncol=%d\n", i, j, nrow, ncol);
	  exit(1);
	}
      }  
    }

    fclose(mfile);
    return matrix;
    
  }else {
    perror("Could not open file\n");
  }

  return NULL;

}
