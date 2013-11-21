#include "beta-binomial.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <bool.h>

inline
double alpha_lprob(const double alpha_trial, const double alpha, const double alpha_other, const double p,
		   const unsigned int wins, const unsigned int games) {
  return (gsl_sf_lngamma(games + alpha_other + alpha_trial) - gsl_sf_lngamma(games + alpha_other + alpha)) + 
    (gsl_sf_lngamma(games + alpha) - gsl_sf_lngamma(games + alpha_trial)) + 
    (alpha_trial - alpha) * log(p);
}

inline
double alpha_lprob_stirling(const double alpha_trial, const double alpha, const double alpha_other, const double p,
		   const unsigned int wins, const unsigned int games) {
  return (alpha_trial - alpha) * (log(alpha_other + games) + log(games) + log(p));
}

inline
bool lMH(double lp, gsl_rng* rng){
  if(lp >= 0)
    return true;

  return gsl_rng_uniform(rng) < exp(lp);
}


double sample_model(const unsigned int *wins, const unsigned int *games, 
		  double *alpha, double *P,
		  void *fxn,
		  const unsigned int team_number, const Run_Params *rp) {
  
  unsigned int i,j,k, iters, acceptance;
  double alpha_trial, rel_lprob, ratio;

  //initalize parameters as needed
  if(!alpha) {
    alpha = (double*) malloc(sizeof(double) * team_number);
    for(i = 0; i < team_number; i++)
      alpha[i] = 1;
  }
  
  if(!P) {
    P = (double*) malloc(sizeof(double) * team_number * (team_number - 1) / 2);
    for(i = 0; i < team_number * (team_number - 1) / 2; i++)
      P[i] = 0.5;
  }


  for(iters = 0 ;iters < rp->iter_max; iters++) {

    //sample new pij
    for(i = 0; i < team_number; i++) {
      for(j = i + 1; j < team_number; j++) {
	P[i * team_number + j] = gsl_ran_beta(rp->rng, alpha[i] + wins[i * team_number + j], alpha[j] + games[i * team_number + j] - wins[i * team_number + j]);
      }
    }
    
    //sample new alpha_ij with replacement
    acceptance = 0;
    for(k = 0; k < team_number; k++) {
      //choose index
      i = gsl_rng_uniform_int(rp->rng, team_number);
      //perturb alpha
      alpha_trial = alpha[i] + gsl_rng_unfiform(rp->rng, 2 * rp->step_size) - rp->step_size;
      //accumate terms
      rel_lprob = 0;
      for(j = k + 1; j < team_number; k++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[k], alpha[j], 
				 wins[k * team_number + j], 
				 games[k * team_number + j],
				 P[k * team_number + j]);
      for(j = 0; j < k; j++)
	rel_lprob += alpha_lprob(alpha_trial, 
				 alpha[k], 
				 alpha[j], 
				 wins[j * team_number + k], 
				 games[j * team_number + k],
				 1 - P[j * team_number + k]);
      //Metrpolis-Hastings 
      if(lMH(rel_lprob)) {
	alpha[k] = alpha_trial;	
	acceptance += 1;
      }      
    }

    //calculate acceptance ratio
    ratio = ratio * (iters) / (iters + 1) + ((double) acceptance / team_number) / (iters + 1);
    
    //run the function
    if(fxn)
      fxn(alpha, P, team_number);

  }
  return ratio
}

void generate_wins(unsigned int *wins, unsigned int *games,
		   const double *alpha,
		   void *fxn,
		   const unsigned int team_number, const Run_Params *rc) {

  unsigned int iters;
  doouble pij;
  
  //initialize parameters as needed
  if(!wins)
    wins = (unsigned int*) calloc(sizeof(unsigned int) * team_number * (team_nubmer - 1) / 2);
  
  //sample realization
  for(iters = 0 ;iters < rp->iter_max; iters++) {
    //sample
    for(i = 0; i < team_number; i++) {
      for(j = i + 1; j < team_number; j++) {
	pij = gsl_ran_beta(rp->rng, alpha[i],alpha[j]);
	wins[i * team_number + j] = gsl_ran_binomial(rp->rng, pij, games[i * team_number + j]);
      }
    }
    fxn(wins);
}
