#include "mendel.h"
#ifdef WINDOWS
#define MAX_DEL_MUTN 20000
#endif

/* This routine eliminates the least fit individuals in a new        *
 * generation to reduce the population size to a level not to exceed *
 * pop_size.  If the population is recovering from a bottlenecking   *
 * event, let half the excess reproduction be used to increase       *
 * population size and the other half be used for selection.         */

void selection(int *** dmutn, int *** fmutn, int **** lb_mutn_count, 
               double *** linkage_block_fitness, double * fitness, 
               double *pheno_fitness, 
               double *work_fitness, double *sorted_score,
	       float *initial_allele_effects, int max_size, 
	       int total_offspring, int gen)
{
  int remaining;
  int i, j, k, l, lb, mutn, m, *zygous; //, zygous[num_linkage_subunits];
  double homozygous_fitness_loss=0, noise;
  double homozygous_fitness_gain, fitness_loss;
  double covariance;
  double max_work_fitness, score_cutoff;
  double geno_fitness_variance, pheno_fitness_variance;
  double mean_pheno_fitness;
  float w, x, count, effect, factor;
  
  int my_max_size;
  my_max_size = max(max_del_mutn_per_indiv,max_size);

 /* items is a buffer used for quick_sort.  For some reason,
  * there is a problem if its memory is dynamically allocated,
  * And Windows Visual C++ compiler gives an error if you try
  * to allocated the array using the variable max_del_mutn_per_indiv. */
#ifdef WINDOWS
  double items[MAX_DEL_MUTN];
#else
  double items[my_max_size];
  //double items[max_del_mutn_per_indiv];
#endif
  
  zygous = malloc(num_linkage_subunits*sizeof(int));

  sprintf(version_selection,
          "$Id: selection.c,v 1.37 2009/07/23 18:56:29 wes Exp $");

  w = multiplicative_weighting;
  
  /* If the population is recovering from a bottlenecking event, *
   * compute the new population size that accounts for selection *
   * as well as population growth.  As a place holder here, let  *
   * half the excess reproduction be used to increase population * 
   * size and the other half be used for selection.              */
  
  if(gen > bottleneck_generation + num_bottleneck_generations 
     && current_pop_size < pop_size) 
    current_pop_size = min(pop_size, (int)(1. + current_pop_size*(1. + 0.25*
		(offspring_per_female*(1. - fraction_random_death) - 2.0))));
  
  /* Compute the fitness of each member of the new generation. */
  
  for (i = 0; i < total_offspring; i++)
    fitness[i] = 1.0;
  
  if(fitness_distrib_type == 0) {
    homozygous_fitness_loss = uniform_fitness_effect;
    homozygous_fitness_gain = uniform_fitness_effect*max_fav_fitness_gain;
  }
  
  for (i=0; i < total_offspring; i++) {
    for (lb=0; lb < num_linkage_subunits; lb++) {
      
      fitness[i] = (fitness[i] - (1. - w)*((double)2.0
					   - linkage_block_fitness[i][0][lb]
					   - linkage_block_fitness[i][1][lb]))
	*((double)1. - ((double)1. - linkage_block_fitness[i][0][lb])*w)
	*((double)1. - ((double)1. - linkage_block_fitness[i][1][lb])*w);
    }

    /* Apply the appropriate fitness degradation adjustment for  *
     * homozygous deleterious mutations.  Skip this step for the *
     * cases of clonal reproduction and co-dominance.            */
    
    if(!clonal_reproduction && dominant_hetero_expression != 0.5) {
      
      j = 1;
      
      for (k=1; k < dmutn[i][0][0]+1; k++) {
	while(abs(dmutn[i][0][k]) >
	      abs(dmutn[i][1][j]) && j <= dmutn[i][1][0]) j++;  
	
	if(dmutn[i][0][k] == dmutn[i][1][j]) {
	  
	  if(dmutn[i][0][k] == num_linkage_subunits*lb_modulo + 1)
	    printf("ERROR: dmutn[%d] range invalid\n",i);
	  
	  mutn = abs(dmutn[i][0][k])%lb_modulo;
	  
	  if(fitness_distrib_type > 0)  
	    homozygous_fitness_loss = 
	      exp(-alpha_del*pow((float)mutn*del_scale,gamma_del));
	  
	  /* Apply the proper fitness decrease associated with a *
	   * homozygous mutation, giving it 100% of the nominal  * 
	   * mutation effect.                                    */
	  
	  fitness[i] = (fitness[i] 
			- (1. - w)*homozygous_fitness_loss)
	    *((double)1. - w*homozygous_fitness_loss);
	  if(dmutn[i][0][k] < 0)  homozygous_fitness_loss =
				    recessive_hetero_expression
				    *homozygous_fitness_loss;
	  if(dmutn[i][0][k] > 0)  homozygous_fitness_loss =
				    dominant_hetero_expression
				    *homozygous_fitness_loss;
	  
	  /* Remove the fitness decreases that were applied elsewhere 
	   * when it was assumed the mutation was heterozygous.  
	   * Remove the heterozygous effect by adding it back twice.
	   * since it was carried out on both haplotypes.           */
	  fitness[i] = fitness[i] / 
	    pow((double)1. - w*homozygous_fitness_loss,2)
	    + (1. - w) *homozygous_fitness_loss*2.;
	}
      }
      
      /* Apply the appropriate fitness enhancement adjustment for
       * homozygous favorable mutations.  */
      
      j = 1;
      
      for (k=1; k < fmutn[i][0][0]+1; k++) {
	
	while(abs(fmutn[i][0][k]) > abs(fmutn[i][1][j]) &&
	      j <= fmutn[i][1][0]) j++;
	if(fmutn[i][0][k] == fmutn[i][1][j]) {
	  
	  if(fmutn[i][0][k] == num_linkage_subunits*lb_modulo + 1)
	    printf("ERROR: fmutn[%d] range invalid\n",i);
	  
	  mutn = abs(fmutn[i][0][k])%lb_modulo;
	  if(fitness_distrib_type > 0) { 
	    homozygous_fitness_gain = max_fav_fitness_gain
	      *exp(-alpha_fav*pow((float)mutn*fav_scale,gamma_fav));
	    fitness[i] = (fitness[i] 
			  + (1. - w)*homozygous_fitness_gain)
	      *((double)1. + w*homozygous_fitness_gain);
	    if(fmutn[i][0][k] < 0) homozygous_fitness_gain =
				     recessive_hetero_expression
				     *homozygous_fitness_gain;
	    if(fmutn[i][0][k] > 0) homozygous_fitness_gain =
				     dominant_hetero_expression
				     *homozygous_fitness_gain;
	    fitness[i] = (fitness[i] 
			  - (1. - w) *homozygous_fitness_gain*2.)
	      / pow(((double)1. + w*homozygous_fitness_gain),2);
	  }
	}
      }
    }
  }
  
  if(synergistic_epistasis && !clonal_reproduction) {
    
    /* In our synergistic epistasis (SE) treatment, we break its 
     * effect into two parts, one involving interactions between 
     * mutations occurring on the same linkage block (linked 
     * interactions) and the other part involving interactions of 
     * mutations on different linkage blocks (nonlinked interactions). 
     * SE effects from linked interactions are inherited, while those
     * from nonlinked interactions are treated as temporary and act
     * in effect as a type of noise as far as the selection process 
     * is concerned.  The main input parameter governing the SE 
     * treatment is the one that specifies the fraction of total SE 
     * contributions from linked interactions relative to those from 
     * nonlinked interactions, namely, linked_mutn_se_fraction.  From
     * this parameter we can readily compute the ratio of the
     * nonlinked SE contribution to the linked contribution as the
     * ratio (1. - linked_mutn_se_fraction)/linked_mutn_se_fraction. 
     
     * Here we treat the nonlinked SE interactions.  The treatment
     * involves a mean SE effect that applies to all the possible
     * nonlinked interactions.  This mean effect for the nonlinked
     * interactions is proportional to (1) the nonlinked to linked
     * SE contribution ratio, (2) the mean SE effect used for the
     * linked interactions, and (3) the ratio of the number of linked
     * to nonlinked interactions.  If M is the total number of 
     * mutations in the genome and n is the number of linkage blocks,
     * then the total number of pairwise interactions between 
     * mutations is M factorial = M(M-1)/2, the mean number of 
     * mutations per linkage block is M/n, and the approximate number
     * of linked interactions is n(M/n)[(M/n)-1]/2.  Since SE 
     * contributions become significant only when M becomes 
     * moderately large, we approximate M-1 by M and (M/n)-1 by M/n.
     * With these approximations, the number of linked interactions 
     * becomes M**2/(2n) and the number of nonlinked interactions
     * becomes (1 - 1/n)*M**2/2.  We used 0.5*n/haploid_genome_size 	
     * for the factor we applied for the mean SE contribution of the
     * linked interactions.  When we combine these component parts
     * together, we get a cancellation of a factor n and obtain a 
     * multiplying factor for the mean nonlinked SE effect given by 
     * 
     *  0.5*(1 - linked_mutn_se_fraction)/(linked_mutn_se_fraction
     *     *haploid_genome_size*(1. - 1./num_linkage_subunits))
     *
     * Note that when linked_mutn_se_fraction = 0.5, that is, when
     * half the total SE effect is from linked mutations and the
     * other half is from nonlinked mutations, then the mean effect
     * for the nonlinked mutations is about 1/n of that for the
     * linked mutations [valid When the factor (1 - 1/n) ~ 1].

     * The approach we apply for incorporating SE contributions for
     * nonlinked interactions, given that these interactions are 
     * predominantly transient and almost entirely disrupted by the
     * genomic shuffling of the deck that occurs in gamete formation,
     * is to generate 100 temporary mutation values.  The SE factor
     * is multiplied by the number of total interactions and then
     * divided by 100 to perfectly compensate for the smaller number
     * in a such a way that the mean SE value is achieved.     */
    
    linked_mutn_se_fraction = max(0.1, linked_mutn_se_fraction);
    
    factor = 0.5*(1. - linked_mutn_se_fraction)
           /(linked_mutn_se_fraction*haploid_genome_size
           *(1. - 1./num_linkage_subunits)*100.);

    for (i=0; i < total_offspring; i++) {
      
      /* Compute the number of deleterious mutations carried by *
       * this offspring.                                        */
      
      count = 0;
      
      for(lb=0; lb < num_linkage_subunits; lb++)
	count += lb_mutn_count[i][0][0][lb] 
	  + lb_mutn_count[i][1][0][lb];
      
      /* In routine offspring we accounted for linked SE inter- *
       * actions.  Here we are dealing with the nonlinked ones. */
      
      for (k=0; k <100; k++) {
	
	if(fitness_distrib_type > 0) {
          x = randomnum();
	  mutn = min(lb_modulo-2, ((int) lb_modulo*x));
	  fitness_loss = exp(-alpha_del
			     *pow((float)mutn*del_scale,gamma_del));
	} else {
	  fitness_loss = uniform_fitness_effect;
	}
	
	fitness_loss = fitness_loss*factor*pow(count,2)/2.;
	
	fitness[i] = (fitness[i] - (1. - w) *fitness_loss)
	  *((double)1. - w*fitness_loss);
	
      }
      
    }
    
  }
  
  /* Account for possible homozygosity in initial contrasting alleles. */
  
  if(num_contrasting_alleles > 0) {
    
    for(i=0; i < current_pop_size; i++) {
      
      for (j = 0; j < num_linkage_subunits; j++)
	zygous[j] = 0;
      
      for(m=1; m < dmutn[i][0][0]+1;m++) {	
	if(dmutn[m][0][i]%lb_modulo == lb_modulo-1) {
	  lb = dmutn[i][0][m]/lb_modulo;
	  zygous[lb]++;
	}
      }
      
      for(m=1; m < dmutn[i][1][0] + 1; m++) {
	if(dmutn[i][1][m]%lb_modulo == lb_modulo-1) {
	  lb = dmutn[i][1][m]/lb_modulo;
	  zygous[lb]++;
	}
      }
      
      for(lb=0; lb < num_linkage_subunits; lb++) {
	if(zygous[lb] == 2) {
	  effect = initial_allele_effects[lb];
	  fitness[i] = (fitness[i] - (1. - w)*effect)
	    *((double)1. - w*effect);
	  effect = recessive_hetero_expression*effect;
	  fitness[i] = fitness[i] / pow((double)1. - w*effect,2)
	    + (1. - w) *effect*2.;
	}
      }
      
      for (j = 0; j < num_linkage_subunits; j++)
	zygous[j] = 0;
      
      for(m=1; m < fmutn[i][0][0]+1; m++) {
	if(fmutn[i][0][m]%lb_modulo == lb_modulo-1) {
	  lb = fmutn[i][0][m]/lb_modulo; 
	  zygous[lb]++;
	}
      }
      
      for(m=1; m < fmutn[i][1][0]+1; m++) {
	if(fmutn[i][1][m]%lb_modulo == lb_modulo-1) {
	  lb = fmutn[i][1][m]/lb_modulo; 
	  zygous[lb]++;
	}
      }
      
      for(lb=0; lb < num_linkage_subunits; lb++) {
	if(zygous[lb] == 2) {
	  effect = initial_allele_effects[lb];
	  fitness[i] = (fitness[i] + (1. - w)*effect)
	    *((double)1. + w*effect);
	  effect =  dominant_hetero_expression*effect;
	  fitness[i] = (fitness[i] - (1. - w) *effect*2.)
	    / pow(((double)1. + w*effect),2);
	}
      }
      
    }
    
  }

  /* Compute the mean genotypic fitness of the new generation. */
  
  pre_sel_fitness = 0.;
  
  for(i=0; i < total_offspring; i++) 
    pre_sel_fitness += fitness[i];
  
  pre_sel_fitness /= total_offspring;
  
  /* Compute the genotypic fitness variance of the new generation. */
  
  geno_fitness_variance = 0.0;
  
  for(i=0; i < total_offspring; i++) 
    geno_fitness_variance += pow(fitness[i] - pre_sel_fitness,2);
  
  geno_fitness_variance /= total_offspring;
  
  pre_sel_geno_sd = sqrt(geno_fitness_variance);
  
  /* If population has collapsed to a single individual, skip the *
   * selection process and return.                                */
  
  if(total_offspring == 1) {
    current_pop_size = 1;
    return;
  }
  
  /* Compute the noise variance required to yield the specified  *
   * heritability.  Add to this fitness-dependent noise a noise  *
   * component that is fitness independent. Take the square root *
   * to obtain the standard deviation.                           */
  
  noise = sqrt(geno_fitness_variance*(1. - heritability)
	       /heritability + pow(non_scaling_noise,2));
  
  /* Add noise to the fitness to create a phenotypic fitness score. *
   * Add a tiny variable positive increment to eliminate identical  * 
   * fitness values when the noise is zero.                         */
  
  
  for(i=0; i < total_offspring; i++) {
    x = random_normal(0,1);
    pheno_fitness[i] = fitness[i] + x*noise + 1.e-15*i;
  }
  
  /* Compute the mean phenotypic fitness of offspring. */
  
  mean_pheno_fitness = 0.;
  
  for(i=0; i < total_offspring; i++)
    mean_pheno_fitness += pheno_fitness[i];
  
  mean_pheno_fitness /= total_offspring;
  
  /* Compute the phenotypic fitness variance, the covariance of   *
   * genotypic and phenotypic fitness, and the genotype-phenotype *
   * correlation.                                                 */
  
  pheno_fitness_variance = 0.;
  covariance = 0.;
  
  for (i=0; i < total_offspring; i++) {
    pheno_fitness_variance += 
      pow(pheno_fitness[i] - mean_pheno_fitness,2);
    covariance += fitness[i]*pheno_fitness[i];
  }
  
  pheno_fitness_variance /= total_offspring;
  
  pre_sel_pheno_sd = sqrt(pheno_fitness_variance);
  
  covariance = covariance/total_offspring 
    - pre_sel_fitness*mean_pheno_fitness;
  
  pre_sel_corr = 0.;
  effect = sqrt(geno_fitness_variance*pheno_fitness_variance);
  if(effect > 0.) pre_sel_corr = covariance/effect;
  
  /* Move, in effect, those offspring whose phenotypic fitness is      *
   * negative to the end of the list of offspring, and then, in effect,*
   * truncate the list so that these individuals cannot reproduce and  *
   * do not even participate in the subsequent selection process.      */
  
  remaining = total_offspring - 1;
  //remaining = total_offspring;
  for(i=0; i < total_offspring; i++) {
    if(pheno_fitness[i] < (double) 0.) {
      while(pheno_fitness[remaining] < (double)0.0) remaining--;
      if(remaining < 0) {
	fprintf(stderr,"ERROR: remaining less than zero, line number:476\n");
	exit(EXIT_FAILURE);
      }
      k = dmutn[remaining][0][0] + 1;
      for (j = 0; j < k; j++)
	dmutn[i][0][j] = dmutn[remaining][0][j];
      k = dmutn[remaining][1][0] + 1;
      for (j = 0; j < k; j++)
	dmutn[i][1][j] = dmutn[remaining][1][j];
      k = fmutn[remaining][0][0] + 1;
      for (j = 0; j < k; j++)
	fmutn[i][0][j] = fmutn[remaining][0][j];
      k = fmutn[remaining][1][0] + 1;
      for (j = 0; j < k; j++)
	fmutn[i][1][j] = fmutn[remaining][1][j];
      
      for (l = 0; l < 2; l++)
	for (j = 0; j < 2; j++)
	  for (m = 0; m < num_linkage_subunits; m++)
	    lb_mutn_count[i][l][j][m] = 
	      lb_mutn_count[remaining][l][j][m];
      
      for (j = 0; j < 2; j++)
	for (m = 0; m < num_linkage_subunits; m++)
	  linkage_block_fitness[i][j][m] = 
	    linkage_block_fitness[remaining][j][m]; 
      
      fitness[i] = fitness[remaining];
      pheno_fitness[i] = pheno_fitness[remaining];
      remaining--;
    }
  }
 
  remaining++; 
  total_offspring = remaining;
  
  /* Adjust the population size for the next generation such that it *
   * does not exceed the number of offspring after removal of those  *
   * with negative phenotypic fitness.                               */
  
  current_pop_size = min(current_pop_size, remaining);
  
  /* Allow the population size for the next generation potentially *
   * to rebound from an earlier reduction in previous generations  * 
   * because of individuals with negative phenotypic fitness.      */
  
  if(!bottleneck_yes && remaining > current_pop_size)
    current_pop_size = min(remaining, pop_size);
  
  /* Copy the phenotypic fitnesses into array work_fitness. */
  
  for (i=0; i < total_offspring; i++) 
    work_fitness[i] = pheno_fitness[i];
  
  if (selection_scheme == 2) {
    
    /* For unrestricted probability selection, divide the phenotypic   *
     * fitness by a uniformly distributed random number prior to       *
     * ranking and truncation.  This procedure allows the probability  *
     * of surviving and reproducing in the next generation to be       *
     * directly related to phenotypic fitness and also for the correct *
     * number of individuals to be eliminated to maintain a constant   *
     * population size.                                                */
    
    for (i=0; i < total_offspring; i++) { 
      work_fitness[i] /= randomnum() + (double)1.e-15;
    }
    
  }
  
  if (selection_scheme == 3) {
    
    /* For strict proportionality probability selection, rescale the
     * phenotypic fitness values such that the maximum value is one.
     * Then divide the scaled phenotypic fitness by a uniformly
     * distributed random number prior to ranking and truncation.  
     * Allow only those individuals to reproduce whose resulting 
     * ratio of scaled phenotypic fitness to the random number value
     * exceeds one.  This approach ensures that no individual 
     * automatically survives to reproduce regardless of the value
     * of the random number.  But it restricts the fraction of the 
     * offspring that can survive.  Therefore, when the reproduction
     * rate is low, the number of surviving offspring may not be
     * large enough to sustain a constant population size.  */
    
    max_work_fitness = 0.;
    for(i=0; i < total_offspring; i++) 
      max_work_fitness = max(max_work_fitness, work_fitness[i]);
    
    for(i=0; i < total_offspring; i++) {
      work_fitness[i] /= max_work_fitness + (double)1.e-15;
      work_fitness[i] /= randomnum() + (double)1.e-15;
    }
    
  }
  
  if (selection_scheme == 4) {
    
    /* For partial truncation selection, divide the phenotypic  
     * fitness by the sum of theta and (1. - theta) times a random 
     * number distributed uniformly between 0.0 and 1.0 prior to 
     * ranking and truncation, where theta is the parameter 
     * partial_truncation_value.  This selection scheme is 
     * intermediate between truncation selection and unrestricted 
     * probability selection.  The procedure allows for the correct 
     * number of individuals to be eliminated to maintain a constant 
     * population size. */
    
    for (i=0; i < total_offspring; i++) 
      work_fitness[i] /= partial_truncation_value 
	+ (1. - partial_truncation_value)
	*randomnum();
  }
  
  /* Sort the resulting work fitnesses in ascending order. */
  
  for (i=0; i < total_offspring; i++)
    sorted_score[i] = work_fitness[i];
  
  for(i=0; i< max_size; i++){
    items[i] = (double) sorted_score[i];
  }
  
  //quick_sort_engine(items, 0, sizeof(items) / sizeof(items[0])-1);
  //myheapsort(sorted_score, total_offspring);
  myheapsort(items, total_offspring);
  
  for(i=0; i< max_size; i++) {
    sorted_score[i] = items[i];
  }
  
  if (selection_scheme <= 4) {
    
    /* Apply truncation selection to reduce the population size to
       current_pop_size.  */
    
    /* Compute the score cutoff value. */
    
    if(total_offspring > current_pop_size) 
      score_cutoff = sorted_score[total_offspring - current_pop_size - 1];
    else 
      score_cutoff = -1000.0;
    
    if(selection_scheme == 3) 
      score_cutoff = max((double) 1.0, score_cutoff);
    
    /* Copy pheno_fitness into array sorted_score for diagnostics *
     * purposes.                                                  */
    
    for(i=0; i< total_offspring; i++) 
      sorted_score[i] = pheno_fitness[i];
    
    /* Remove those individuals whose score lies below the cutoff
       value to reduce the population size to its appropriate value.  */
    
    //printf("SELECTION: current_pop_size = %d total_offspring=%d\n",
    //        current_pop_size, total_offspring);
    current_pop_size = min(current_pop_size, total_offspring);
    remaining = total_offspring;
    remaining--;
 
    for (i=0; i < current_pop_size; i++) {
      
      /* If the work fitness if individual i is below the cutoff *
       * value, find another individual in the pool of excess    *
       * offspring whose work fitness is equal to or above the   *
       * cutoff value and replace the first individual with the  *
       * second in the list of reproducing individuals for that  *
       * generation.                                             */
      
      if(work_fitness[i] < score_cutoff && i < remaining) {
	while(work_fitness[remaining] < score_cutoff) remaining--;
	if(remaining < 0) {
	  fprintf(stderr,"ERROR: remaining less than zero, line number:668 \n");
	  exit(EXIT_FAILURE);
	}
	if(i < remaining) {
	  k = dmutn[remaining][0][0] + 1;
	  for (j = 0; j < k; j++)
	    dmutn[i][0][j] = dmutn[remaining][0][j];
	  k = dmutn[remaining][1][0] + 1;
	  for (j = 0; j < k; j++)
	    dmutn[i][1][j] = dmutn[remaining][1][j];
	  k = fmutn[remaining][0][0] + 1;
	  for (j = 0; j < k; j++)
	    fmutn[i][0][j] = fmutn[remaining][0][j];
	  k = fmutn[remaining][1][0] + 1;
	  for (j = 0; j < k; j++)
	    fmutn[i][1][j] = fmutn[remaining][1][j];
	  
	  for (l = 0; l < 2; l++)
	    for (j = 0; j < 2; j++)
	      for (m = 0; m < num_linkage_subunits; m++)
		lb_mutn_count[i][l][j][m] = 
		  lb_mutn_count[remaining][l][j][m];
	  
	  for (j = 0; j < 2; j++)
	    for (m = 0; m < num_linkage_subunits; m++)
	      linkage_block_fitness[i][j][m] = 
		linkage_block_fitness[remaining][j][m];
	  fitness[i] = fitness[remaining];
	  pheno_fitness[i] = pheno_fitness[remaining];
	  remaining--;
	}
      }
    }
    
    current_pop_size = min(current_pop_size, remaining);
    
  } else {
    
    fprintf(stderr,"ERROR: invalid selection scheme %d\n", selection_scheme);
    fprintf(fp9,"ERROR: invalid selection scheme %d\n", selection_scheme);
    exit(EXIT_FAILURE); 
    
  }
  
  /* Compute the mean genotypic and phenotypic fitnesses of the new *
   * generation after selection.                                    */
  
  post_sel_fitness   = 0.0;
  mean_pheno_fitness = 0.0;
  
  for (i=0; i < current_pop_size; i++) {
    post_sel_fitness   += fitness[i];
    mean_pheno_fitness += pheno_fitness[i];
  }
  
  post_sel_fitness   /= current_pop_size;
  mean_pheno_fitness /= current_pop_size;
  
  /* Compute the genotypic and phenotypic fitness variances, the *
   * covariance of genotypic and phenotypic fitness, and the     *
   * genotype-phenotype correlation of the new generation.       */
  
  geno_fitness_variance  = 0.0;
  pheno_fitness_variance = 0.0;
  covariance = 0.0;
  
  for (i=0; i < current_pop_size; i++) {
    geno_fitness_variance  += pow(fitness[i] - post_sel_fitness,2);
    pheno_fitness_variance += pow(pheno_fitness[i] - mean_pheno_fitness,2);
    covariance += fitness[i]*pheno_fitness[i];
  }
  
  geno_fitness_variance  /= current_pop_size;
  pheno_fitness_variance /= current_pop_size;
  
  post_sel_geno_sd  = sqrt(geno_fitness_variance);
  post_sel_pheno_sd = sqrt(pheno_fitness_variance);
  
  covariance = covariance/current_pop_size 
    - post_sel_fitness*mean_pheno_fitness;
  
  post_sel_corr = 0.;
  effect = sqrt(geno_fitness_variance*pheno_fitness_variance);
  if(effect > 0.) post_sel_corr = covariance/effect;
  
  post_sel_fitness = max(0., post_sel_fitness);
  
  for (i = current_pop_size; i < max_size; i++) 
    fitness[i] = 0.0;
  
}

