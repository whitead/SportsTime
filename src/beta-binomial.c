#include "beta-binomial.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#define FILENAME_BUFFER_LENGTH 64
#define NBINS 1024
#define tsize(x) (x * (x - 1) / 2)
#define MAX_WINS 25

inline
double alpha_lprob(double alpha_trial, double alpha, double alpha_other, double p) {
  return (gsl_sf_lngamma(alpha_other + alpha_trial) - gsl_sf_lngamma(alpha_other + alpha)) + (gsl_sf_lngamma(alpha) - gsl_sf_lngamma(alpha_trial)) +  (alpha_trial - alpha) * log(p);
}

inline
double lprob(double *alpha, double *P, unsigned int team_number) {
  double lprob = 0;
  unsigned int i, j, index;
  for(i = index = 0; i < team_number; i++) {
    for(j = i + 1; j < team_number; j++) {
      lprob += gsl_sf_lngamma(alpha[i] + alpha[j]) + (alpha[i] - 1) * log(P[index]) + (alpha[j] - 1) * log(1 - P[index]) - gsl_sf_lngamma(alpha[i]) - gsl_sf_lngamma(alpha[j]);
      index++;
    }
  }
  return lprob;
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

/* 
 * Calculate the log-likelihood of a single new alpha 
 */
inline
double single_alpha_move_lprob(unsigned int trial_index, 
		  double alpha_trial, 
		  double *alpha, double *P, 
		  unsigned int team_number) {
      //accumate log probability terms
      double rel_lprob = 0;
      unsigned int j;
      for(j = trial_index + 1; j < team_number; j++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[trial_index], 
				 alpha[j], 
				 P[symm_index(trial_index,j,team_number)]);

      for(j = 0; j < trial_index; j++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[trial_index], 
				 alpha[j], 
				 1 - P[symm_index(trial_index,j,team_number)]);
      return rel_lprob;
}
/* 
 * Calculate the log-likelihood for an array of new alphas
 */
inline
double multi_alpha_move_lprob(double *alpha_trial_array, 
		  double *alpha, double *P, 
		  unsigned int team_number) {

  return lprob(alpha_trial_array, P, team_number) - lprob(alpha, P, team_number);
}

/* L_2 measure with a skipped value
 *
 */
inline
double l2(double* array, unsigned int length, unsigned int skip_index) {
  unsigned int i;
  double sum = 0;

  for(i = 0; i < skip_index; i++)
    sum += array[i] * array[i];
  for(i = skip_index + 1; i < length; i++)
    sum += array[i] * array[i];

  return sum;
}

int test() {


  //test my equations
  //alpha_lprob 
  //should not improve
  double test = alpha_lprob(10, 1, 1, 0.5);
  printf("test 1: %g < 0 ?\n", test);
  //should improve
  test = alpha_lprob(2, 1, 1, 0.75);
  printf("test 1: %g > 0 ?\n", test);

  test = alpha_lprob(1, 1, 1, 0.5);
  printf("test 1: %g == 0 ?\n", test);


  //set up a test system
  //5 teams, alpha is some simple linear pattern  
  unsigned int team_number = 5;
  unsigned int i,j;
  double *alpha = (double*) malloc(sizeof(double) * team_number);
  for(i = 0; i < team_number; i++) {
    alpha[i] = 21 - i * 4;
  }
  double  max_alpha = 15 * sqrt(team_number);

  //run parameters. Just one win matrix.
  Run_Params* rp = init_run_params(1, 0, 0.2, 10);

  //make a simple schedule
  unsigned int *games = (unsigned int*) malloc(sizeof(unsigned int) * tsize(team_number));

  for(i = 0; i < tsize(team_number); i++)
    games[i] = 3; 

  unsigned int *wins = (unsigned int*) calloc(tsize(team_number), sizeof(unsigned int));
  double *alpha_fit = NULL;
  double *P = NULL;

  //generate a set of wins given those alphas.   
  generate_wins(wins, games, alpha, print_wins_fxn, NULL, team_number, rp);

  free(rp);  

  //now we do a large run to try to recover our alphas from the win
  //matrix
  rp = init_run_params(1000000, 0.2, 2, max_alpha);
  //max alpha should be 
  Array_Histogram* alpha_hist = build_array_histogram(team_number, NBINS, max_alpha / NBINS);    

  printf("Acceptance: %g\n", sample_model(wins, games, alpha_fit, P, accum_sample_fxn, alpha_hist, team_number, rp));

  //now we check to see how close the mode of alpha is to our original parameters

  FILE *mfile = fopen("alpha.txt", "w");
  log_histogram(mfile, alpha_hist);
  fclose(mfile);

  unsigned long int sample_sum;
  double alpha_median;
  printf("TRUE | ESTIMATOR\n");
  for(i = 0; i < team_number; i++) {
    //for each alpha
    sample_sum = 0;    
    for(j = 0; j < alpha_hist->nbins; j++) {
      sample_sum += alpha_hist->hist[i * alpha_hist->nbins + j];
      if(sample_sum > rp->sample_number / 2) {
	alpha_median = j * alpha_hist->width;
	break;
      }
    }
    //compare with true value
    printf("%g | %g\n", alpha[i], alpha_median);
  }
  
  free(alpha_hist);
  free(alpha);

  return 0;
}

void print_help() {

  printf("Usage: beta-binomial -w [filename] -p [filename] -m [team number] -n [iterations]\n");
  printf("w (wins filename) - file containing matrix of wins between M teams\n");
  printf("p (predict filename) - optional file which contains a matrix of games whose outcome should be predicted\n");
  printf("m (integer) - number of teams\n");
  printf("n (integer) - number of model samples to use\n");

}

int main(int argc, char* argv[]) {
  
  if(argc < 6) {
    print_help();
    return 0;
  }

  unsigned int team_number, iterations;
  team_number = iterations = 0;
  char win_filename[FILENAME_BUFFER_LENGTH];
  char* predict_filename = NULL;
  extern char *optarg;
  char c;
  while((c = getopt(argc, argv, "w:p:m:n:")) != -1) {
    
    switch(c) {      
    case 'w':
      strcpy(win_filename, optarg);
      break;
    case 'p':
      predict_filename = (char*) malloc(sizeof(char) * FILENAME_BUFFER_LENGTH);
      strcpy(predict_filename, optarg);
      break;
    case 'm':
      team_number = atoi(optarg);
      break;
    case 'n':
      iterations = atoi(optarg);
      break;
    default:
      print_help();
      return 0;
    }
  }
 
#ifdef DEBUG
  printf("Parsed arguments: \n");
  printf("win_file: %s\n", win_filename);
  printf("predict_file: %s\n", predict_filename ? predict_filename : "N/A");
  printf("team_number: %u\n", team_number);
  printf("iterations: %u\n", iterations);
#endif //DEBUG
     
  unsigned int* square_wins = load_uint_matrix(win_filename, team_number, team_number);  

  unsigned int i,j;
  unsigned int *games = (unsigned int*) malloc(sizeof(unsigned int) * tsize(team_number));
  unsigned int *wins = (unsigned int*) calloc(tsize(team_number), sizeof(unsigned int));
  unsigned int max_games = 0;
  double max_alpha;

  for(i = 0; i < team_number; i++) {
    for(j = i + 1; j < team_number; j++) {
      games[symm_index(i,j,team_number)] = square_wins[i * team_number + j];
      games[symm_index(i,j,team_number)] += square_wins[j * team_number + i];
      if(games[symm_index(i,j,team_number)] > max_games)
	max_games = games[symm_index(i,j,team_number)];
      wins[symm_index(i,j,team_number)] = square_wins[i * team_number + j];
    }      
  }

#ifdef DEBUG

  unsigned int w, g, k;
  //print record of all the teams
  for(i = 0; i < team_number; i++) {
    w = g = 0;
    for(j = 0; j < i; j++) {
      k = symm_index(j,i,team_number);
      g += games[k];
      w += games[k] - wins[k];
    }
    for(j = i + 1; j < team_number; j++) {
      k = symm_index(i,j,team_number);
      w += wins[k];
      g += games[k];
    }
    printf("Record %d %d - %d\n", i, w, g - w);
  }

#endif //DEBUG
  

  max_alpha = fmin(max_games, tsize(team_number)) * sqrt(team_number) / 2;
  Run_Params* rp = init_run_params(iterations, 0.5, max_alpha / 25, max_alpha);

  if(predict_filename == NULL) {
    Array_Histogram* alpha_hist = build_array_histogram(team_number, NBINS, max_alpha / NBINS);    
    double acceptance = sample_model(wins, games, NULL, NULL, accum_sample_fxn, alpha_hist, team_number, rp);
    printf("Acceptance: %g\n", acceptance);
    //write histogram
    FILE *mfile = fopen("alpha.txt", "w");
    log_histogram(mfile, alpha_hist);
    fclose(mfile);
  } else {
    //sample wins while sampleing the model
    //setup win parameters
    Wins_Parameters* wp = (Wins_Parameters*) malloc(sizeof(Wins_Parameters));
    wp->mhist = build_array_histogram(tsize(team_number), MAX_WINS, 1.0); 
    wp->win_sample = (unsigned int*) malloc(sizeof(unsigned int) * tsize(team_number));
    //load the games to predict
    wp->games = diagonalize(load_uint_matrix(predict_filename, team_number, team_number), team_number);

    //sample from loaded games
    sample_model(wins, games, 
		 NULL, NULL, 
		 sample_wins_wrapper, wp, 
		 team_number, rp);
    //write histogram
    FILE *mfile = fopen("predicted_wins.txt", "w");
    log_histogram(mfile, wp->mhist);
    fclose(mfile);
    
  }    
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
  double alpha_trial, rel_lprob, radius, ratio = 0;
  double max_alpha = rp->max_alpha;
  double *alpha_trial_array = (double*) malloc(sizeof(double)  * team_number);


  //initalize parameters as needed
  if(!alpha) {
    alpha = (double*) malloc(sizeof(double) * team_number);
    for(i = 0; i < team_number; i++)
      alpha[i] = 1;
  }
  
  if(!P)
    P = (double*) malloc(sizeof(double) * tsize(team_number));

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
      //perturb alpha. Reject samples less than 1.
      alpha_trial = alpha[i] + gsl_rng_uniform(rp->rng) * 2 * rp->step_size - rp->step_size;
      while(alpha_trial < 1)
	alpha_trial = alpha[i] + gsl_rng_uniform(rp->rng) * 2 * rp->step_size - rp->step_size;
      //check if it's in feasible set      
      radius =  l2(alpha, team_number, i) + alpha_trial * alpha_trial;
      //if it is, we only have to evaluate the trial
      if(radius <= max_alpha*max_alpha) {
	rel_lprob = single_alpha_move_lprob(i, alpha_trial, alpha, P, team_number);
	//Metrpolis-Hastings 
	if(lMH(rel_lprob, rp->rng)) {
	  alpha[i] = alpha_trial;	
	  acceptance += 1;
	}      
      }
      else {
	//alpha_trial is now a matrix
	radius = sqrt(radius);
	for(j = 0; j < team_number; j++)
	  alpha_trial_array[j] = max_alpha * alpha[j] / radius;
	alpha_trial_array[i] = max_alpha * alpha_trial / radius;
	rel_lprob = multi_alpha_move_lprob(alpha_trial_array, alpha, P, team_number);
	//Metrpolis-Hastings 
	if(lMH(rel_lprob, rp->rng)) {
	  for(j = 0; j < team_number; j++)
	    alpha[j] = alpha_trial_array[j];
	  acceptance += 1;
	} 
      }
    }

    //calculate acceptance ratio
    ratio = ratio * (iters) / (iters + 1) + ((double) acceptance / team_number) / (iters + 1);
    
    //run the function
    if(fxn && iters > rp->iterations * rp->equilibrium_ratio)
      fxn(alpha, P, team_number, iters, rp, fxn_args);

  }
  
  free(alpha_trial_array);

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
    wins = (unsigned int*) calloc(tsize(team_number), sizeof(unsigned int));
  
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
    if(fxn)
      fxn(wins, games, team_number, iters, rp, fxn_args);
  }
}

Run_Params* init_run_params(unsigned int iterations, double equilibrium_ratio, double step_size, double max_alpha) {
  Run_Params* rp = (Run_Params*) malloc(sizeof(Run_Params));
  rp->iterations = iterations;
  rp->equilibrium_ratio = equilibrium_ratio;
  rp->sample_number = (int) floor(iterations * (1 - equilibrium_ratio));
  rp->step_size = step_size;
  rp->max_alpha = max_alpha;

  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rp->rng = gsl_rng_alloc(T);

  return rp;
}


void log_array(FILE* file, double* array, unsigned int n_cols, unsigned int n_rows) {
  
  if(file == NULL)
    return;

  unsigned int i,j;
  for(i = 0; i < n_rows; i++) {
    for(j = 0; j < n_cols; j++) {
      fprintf(file, "%5g ", array[i * n_cols + j]); 
    }
    fprintf(file, "\n");
  }
  fflush(file);  
}

void log_histogram(FILE* file, Array_Histogram* hist) {
  
  if(file == NULL)
    return;

  fprintf(file, "# dimension %d\n", hist->dimension);
  fprintf(file, "# nbins %d\n", hist->nbins);
  fprintf(file, "# width %g\n", hist->width);
  
  unsigned int i,j;
  fprintf(file, "Value ");
  for(i = 0; i < hist->dimension; i++)
    fprintf(file, "%5u ", i);
  fprintf(file, "\n");

  for(i = 0; i < hist->nbins; i++) {
    fprintf(file, "%5g ", hist->width * i + hist->width / 2.);
    for(j = 0; j < hist->dimension; j++) {
      fprintf(file, "%5lu ", hist->hist[j * hist->nbins + i]); 
    }
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

void accum_sample_fxn(const double* alpha, const double* P,		      
		      unsigned int team_number, unsigned int iteration, 
		      const Run_Params *rp, void* arg) {

  Array_Histogram *alpha_hist = (Array_Histogram*) arg;
  unsigned int i,j;
  
  //histogram current alpha values
  for(i = 0; i < team_number; i++) {
    j = (unsigned int) fmin(alpha[i] / alpha_hist->width, 
			    alpha_hist->nbins - 1);
    alpha_hist->hist[i * alpha_hist->nbins + j]++;
  }
    
}

unsigned int* load_uint_matrix(char* filename, unsigned int nrow, unsigned int ncol) {

  FILE *mfile;
  mfile = fopen(filename, "r");
  if(mfile != NULL) {

    unsigned int i, j;

    unsigned int *matrix =  (unsigned int*) malloc(sizeof(unsigned int) * nrow * ncol);
    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%u", &(matrix[i*ncol + j])) == 0) {
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

void sample_wins_wrapper(const double* alpha, const double* P,
			 unsigned int team_number, 
			 unsigned int iterations, const Run_Params *rp,
			 void* win_params) {


  Wins_Parameters* wp = (Wins_Parameters*) win_params;
  Array_Histogram *mhist = wp->mhist;
  
  unsigned int i,j,k,l;
  
  //sample wins
  generate_wins(wp->win_sample, wp->games, alpha,
		NULL, NULL, team_number, rp);
  
  //histogram value
  for(i = 0; i < team_number; i++){ 
    for(j = i + 1; j < team_number; j++) {
      k = symm_index(i,j,team_number);      
      l = (unsigned int) fmin(wp->win_sample[k] / mhist->width, 
			      mhist->nbins - 1);
      mhist->hist[k * mhist->nbins + l]++;
    }
  }
}

Array_Histogram* build_array_histogram(unsigned int dimension, 
					 unsigned int nbins,
					 double width) {

  Array_Histogram* mhist = (Array_Histogram*) malloc(sizeof(Array_Histogram));
  mhist->dimension = dimension;
  mhist->nbins = nbins;
  mhist->width = width;
  mhist->hist = (unsigned long int*) malloc(sizeof(unsigned long int) * 
					    dimension * nbins);
  return mhist;

}

unsigned int* diagonalize(unsigned int* matrix, unsigned int dim) {
  unsigned int* new_matrix = (unsigned int*) malloc(sizeof(unsigned int) * tsize(dim));
  unsigned int i, j, index;
  for(i = index = 0; i < dim; i++)
    for(j = i+1; j < dim; j++)
      new_matrix[index++] = matrix[i * dim + j];
  free(matrix);
  return new_matrix;
}

