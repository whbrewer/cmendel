#include "mendel.h"

/* This routine creates an offspring from a pair of individuals,      *
 * indexed 'dad' and 'mom', in the current population.  This          *
 * offspring inherits one set of mutations from each linkage block    *
 * pair from each parent.  This mutation information is loaded into   *
 * array 'offsprng'.  Also, a set of new deleterious mutations,       *
 * equal in number to new_mutn_per_offspring*(1. - frac_fav_mutn) and *
 * chosen randomly, are generated for this offspring.  The fitness    *
 * associated with each of these new mutations is applied to update   *
 * the fitness of the linkage block in which it occurs.               */

void offspring(int*** dmutn_offsprng, int*** fmutn_offsprng, 
               int**** offsprng_lb_mutn_count, double*** offsprng_lb_fitness, 
               int*** dmutn, int*** fmutn, int**** lb_mutn_count, 
               double*** linkage_block_fitness, int dad, int mom, 
               int child)
{
  double fitness_effect, se_effect;
  float w, x;
  int chr_length, ch, ls0, ls1, ls2, iseg, iseg_max;
  int md1, md2, mf1, mf2, mdd_off, mfd_off, mdm_off, mfm_off;
  int lb, hap_id=0, new_mutn, mut, mutn, mutn_indx, num_mutn, i, j, k;
  bool fav;
  
  sprintf(version_offspring,
          "$Id: offspring.c,v 1.40 2008/11/25 05:21:59 wes Exp $");
  
  w = multiplicative_weighting;
  
  for (i = 0; i < 2; i++) {
    dmutn_offsprng[child][i][0] = 0;
    fmutn_offsprng[child][i][0] = 0;
  }
  
  if(!clonal_reproduction) {
    
    if(dynamic_linkage) iseg_max = 3;
    else iseg_max = 1;
    
    chr_length = num_linkage_subunits/haploid_chromosome_number;
    //printf("chr_length is %d %d %d\n",chr_length, num_linkage_subunits, haploid_chromosome_number);
     
    md1     = 1;
    md2     = 1;
    mf1     = 1;
    mf2     = 1;
    mdd_off = 1;
    mfd_off = 1;
    
    for(ch=1; ch <= haploid_chromosome_number; ch++) {
      ls0 = (ch - 1)*chr_length + 1;
      x = randomnum();
      ls1 = min(chr_length-1,abs((int)(chr_length*x))) + ls0;
      x = randomnum();
      ls2 = min(chr_length-1,abs((int)(chr_length*x))) + ls0;
      //printf("ls1 = %d ls2 = %d\n",ls1,ls2);
      
      if(ls1 > ls2) {
	lb  = ls1;
	ls1 = ls2;
	ls2 = lb;
      }
      
      if(!dynamic_linkage) ls1 = num_linkage_subunits;
      
      for(iseg=1; iseg <= iseg_max; iseg++) {
	
	if(iseg == 2) {
	  ls0 = ls1 + 1;
	  ls1 = ls2;
	} else if(iseg == 3) {
	  ls0 = ls2 + 1;
	  ls1 = ch*chr_length;
	}
	
	if(dynamic_linkage){
           x = randomnum(); 
           hap_id = min(1,(int)(2*x));
	}
        
	for(lb=ls0 - 1; lb < ls1; lb++) {
	  
	  /* Copy the mutation list from the randomly selected haplotype *
	   * from the father to form gamete one for the offspring.  Also *
	   * copy the corresponding linkage block mutation count and     *
	   * fitness.                                                    */
	  
	  if(!dynamic_linkage){
             x = randomnum(); 
             hap_id = min(1,(int)(2*x));
          }
	  
	  if(tracking_threshold != 1.) {
	    
	    while(abs(dmutn[dad][0][md1]) < (lb+1)*lb_modulo &&
		  md1 <= dmutn[dad][0][0]) {
	      if(hap_id == 0) {
		dmutn_offsprng[child][0][mdd_off] = 
		  dmutn[dad][0][md1];
		mdd_off++;
	      }
	      md1++;
	    }
	    
	    while(abs(dmutn[dad][1][md2]) < (lb+1)*lb_modulo &&
		  md2 <= dmutn[dad][1][0]) {
	      if(hap_id == 1) {
		dmutn_offsprng[child][0][mdd_off] = 
		  dmutn[dad][1][md2];
		mdd_off++;
	      }
	      md2++;
	    }
	    
	    while(abs(fmutn[dad][0][mf1]) < (lb+1)*lb_modulo &&
		  mf1 <= fmutn[dad][0][0]) {
	      if(hap_id == 0) {
		fmutn_offsprng[child][0][mfd_off] =
		  fmutn[dad][0][mf1];
		mfd_off++;
	      }
	      mf1++;
	    }
	    
	    while(abs(fmutn[dad][1][mf2]) < (lb+1)*lb_modulo &&
		  mf2 <= fmutn[dad][1][0]) {
	      if(hap_id == 1) {
		fmutn_offsprng[child][0][mfd_off] =
		  fmutn[dad][1][mf2];
		mfd_off++;
	      }
	      mf2++;
	    }
	    
	  }
	  
	  offsprng_lb_mutn_count[child][0][0][lb] = 
	    lb_mutn_count[dad][0][hap_id][lb];
	  offsprng_lb_mutn_count[child][1][0][lb] = 
	    lb_mutn_count[dad][1][hap_id][lb];
	  offsprng_lb_fitness[child][0][lb] = 
	    linkage_block_fitness[dad][hap_id][lb];
	  
	}
      }
    }

    md1     = 1;
    md2     = 1;
    mf1     = 1;
    mf2     = 1;
    mdm_off = 1;
    mfm_off = 1;
    
    for(ch=1; ch <= haploid_chromosome_number; ch++) {
      
      ls0 = (ch - 1)*chr_length + 1;
      x = randomnum();
      ls1 = min(chr_length-1,abs((int)(chr_length*x))) + ls0;
      x = randomnum();
      ls2 = min(chr_length-1,abs((int)(chr_length*x))) + ls0;
      
      if(ls1 > ls2) {
	lb  = ls1;
	ls1 = ls2;
	ls2 = lb;
      }
      
      if(!dynamic_linkage) ls1 = num_linkage_subunits;
      
      for(iseg=1; iseg <= iseg_max; iseg++) {
	
	if(iseg == 2) {
	  ls0 = ls1 + 1;
	  ls1 = ls2;
	} else if(iseg == 3) {
	  ls0 = ls2 + 1;
	  ls1 = ch*chr_length;
	}
	
	if(dynamic_linkage){
          x = randomnum();
          hap_id = min(1,(int)(2.*x));
        }
		
	for (lb=ls0 - 1; lb < ls1; lb++) {
	  
	  /* Copy the mutation list from the randomly selected haplotype *
	   * from the mother to form gamete two for the offspring.  Also *
	   * copy the corresponding linkage block mutation count and     *
	   * fitness.                                                    */
	  
	  if(!dynamic_linkage){
            x = randomnum();
            hap_id = min(1,(int)(2.*x));
          }
	  
	  if(tracking_threshold != 1.) {
	    
	    while(abs(dmutn[mom][0][md1]) < (lb+1)*lb_modulo &&
		  md1 <= dmutn[mom][0][0]) {
	      if(hap_id == 0) {
		dmutn_offsprng[child][1][mdm_off] = 
		  dmutn[mom][0][md1];
		mdm_off++;
	      }
	      md1++;
	    }
	    
	    while(abs(dmutn[mom][1][md2]) < (lb+1)*lb_modulo &&
		  md2 <= dmutn[mom][1][0]) {
	      if(hap_id == 1) {
		dmutn_offsprng[child][1][mdm_off] = 
		  dmutn[mom][1][md2]; 
		mdm_off++;
	      }
	      md2++;
	    }
	    
	    while(abs(fmutn[mom][0][mf1]) < (lb+1)*lb_modulo &&
		  mf1 <= fmutn[mom][0][0]) {
	      if(hap_id == 0) {
		fmutn_offsprng[child][1][mfm_off] = 
		  fmutn[mom][0][mf1];
		mfm_off++;
	      }
	      mf1++;
	    }
	    
	    while(abs(fmutn[mom][1][mf2]) < (lb+1)*lb_modulo &&
		  mf2 <= fmutn[mom][1][0]) {
	      if(hap_id == 1) { 
		fmutn_offsprng[child][1][mfm_off] = 
		  fmutn[mom][1][mf2];
		mfm_off++;
	      }
	      mf2++;
	    }
	    
	  }
	  
	  offsprng_lb_mutn_count[child][0][1][lb] = 
	    lb_mutn_count[mom][0][hap_id][lb];
	  offsprng_lb_mutn_count[child][1][1][lb] = 
	    lb_mutn_count[mom][1][hap_id][lb];
	  offsprng_lb_fitness[child][1][lb] = 
	    linkage_block_fitness[mom][hap_id][lb];
	  
	}
	
      }
      
    }
    
    /* Load the mutation count into the first location in the first *
     * index of each of the mutation arrays.                        */
    
    dmutn_offsprng[child][0][0] = mdd_off - 1;
    dmutn_offsprng[child][1][0] = mdm_off - 1;
    fmutn_offsprng[child][0][0] = mfd_off - 1;
    fmutn_offsprng[child][1][0] = mfm_off - 1;
    
  } else { /* Treat case of clonal reproduction. */
    
    mdd_off = dmutn[dad][0][0];
    mdm_off = dmutn[dad][1][0];
    mfd_off = fmutn[dad][0][0];
    mfm_off = fmutn[dad][1][0];
    
    for (i = 0; i <= mdd_off; i++) 
      dmutn_offsprng[child][0][i] = dmutn[dad][0][i];
    for (i = 0; i <= mdm_off; i++) 
      dmutn_offsprng[child][1][i] = dmutn[dad][1][i];
    for (i = 0; i <= mfd_off; i++) 
      fmutn_offsprng[child][0][i] = fmutn[dad][0][i];
    for (i = 0; i <= mfm_off; i++) 
      fmutn_offsprng[child][1][i] = fmutn[dad][1][i];
    
    for (i = 0; i < 2; i++)
      for (j = 0; j < 2; j++)
	for (k = 0; k < num_linkage_subunits; k++)
	  offsprng_lb_mutn_count[child][i][j][k] = 
	    lb_mutn_count[dad][i][j][k];
    
    for (i = 0; i < 2; i++)
      for (j = 0; j < num_linkage_subunits; j++)
	offsprng_lb_fitness[child][i][j] = 
	  linkage_block_fitness[dad][i][j]; 
    
  } /* End clonal_reproduction */
  
    /* Generate new_mutn new random mutations in randomly chosen     *
     * linkage block locations, where new_mutn is a random deviate   *
     * drawn from a Poisson distribution with a mean value given by  *
     * the parameter poisson_mean, which under most circumstances is *
     * very close to the value of new_mutn_per_offspring.            */
  
  //new_mutn = random_poisson(poisson_mean);
  new_mutn = poisson(poisson_mean);
  // printf("NEW_MUTN is %d\n",new_mutn);
  
  new_mutn_count += new_mutn;
  
  for(mut=1; mut <= new_mutn; mut++) {
    
    /* Select a random linkage block for the new mutation. */
    x = randomnum();
    lb = min(num_linkage_subunits,
	     1 + abs((int)(num_linkage_subunits*x)));
    
    x = randomnum();
    hap_id = min(1,(int)(2*x));
    
    /* Determine whether new mutation is deleterious or favorable. */
    
    if(randomnum() < frac_fav_mutn) fav = true;
    else fav = false;
 
    /* Compute the mutation index that is used to specify its fitness. */
    
    x = randomnum();
    
    if(tracking_threshold != 1.0) {
      if(fav) mutn = min(lb_modulo-2, (int)(x/fav_scale));
      else mutn = min(lb_modulo-2, (int)(x/del_scale)); 
    } else {
      mutn = 1;
    }
    
    mutn_indx = (lb - 1)*lb_modulo + mutn;
    
    /* Identify the appropriate fraction of new mutations as   *
     * recessive.  To distinguish recessive mutations from the *
     * dominant ones, label the recessives by assigning them a *
     * negative mutation index.                                */
    
    if(fraction_recessive > 0.) 
      if(randomnum() < fraction_recessive) mutn_indx = -mutn_indx;
    
    /*
     * Compute the fitness effect associated with this new mutation.
     *
     * When parameter fitness_distrib_type is 1, the fitness
     * effect e is obtained from the mutation index mutn using a
     * distribution function of the form
     *
     * e = exp(-alpha_del*x**gamma_del) ,
     *
     * where alpha_del is log(haploid_genome_size) and x is a
     * random number uniformly distributed between zero and one.
     *
     * When parameter fitness_distrib_type is 0, the fitness
     * effect is constant and given by the expression
     *
     * e = uniform_fitness_effect
     *
     * For favorable mutations, this fitness effect is further
     * multiplied by the input paramenter max_fav_fitness_gain,
     * which specifies the upper limit on the favorable fitness
     * effect.                                                     */
    
    if(fitness_distrib_type > 0) {
       if(fav) fitness_effect = max_fav_fitness_gain*
                                exp(-alpha_fav*pow(x,gamma_fav));
       else fitness_effect = exp(-alpha_del*pow(x,gamma_del));
    } else {
       if(fav) fitness_effect = max_fav_fitness_gain*
                                uniform_fitness_effect;
       else    fitness_effect = uniform_fitness_effect;
    }
    
    /* Track this mutation if its fitness effect exceeds the value * 
     * of tracking_threshold.                                      */

    if(fitness_effect > tracking_threshold) {
      
      if(fav) {
	
	/* Test to see if the storage limit of array offspring.fmutn *
	 * has been exceeded.  (Note that we are using the first     *
	 * slot to hold the actual mutation count.)                  */
	
	num_mutn = fmutn_offsprng[child][hap_id][0] + 1;
	
	if(num_mutn+1 >= max_fav_mutn_per_indiv/2) {
	  fprintf(stderr,"ERROR: Favorable mutation count exceeds limit\n");
	  fprintf(fp9,"ERROR: Favorable mutation count exceeds limit\n");
	  exit(EXIT_FAILURE); 
	}
	
	fmutn_offsprng[child][hap_id][0] = num_mutn;
	
	/* Insert new mutation such that mutations are maintained *
	 * in ascending order of their absolute value.            */
	
	i = num_mutn - 1;
	
	while(abs(fmutn_offsprng[child][hap_id][i]) > abs(mutn_indx) 
	      && i > 0) {
	  fmutn_offsprng[child][hap_id][i+1] = 
	    fmutn_offsprng[child][hap_id][i];
	  i--;
	}
	
	fmutn_offsprng[child][hap_id][i+1] = mutn_indx;

      } else {
	
	/* Test to see if the storage limit of array offspring.dmutn *
	 * has been exceeded.  (Note that we are using the first     *
	 * slot to hold the actual mutation count.)                  */
	
	num_mutn = dmutn_offsprng[child][hap_id][0] + 1;
	
	if(num_mutn+1 >= max_del_mutn_per_indiv/2) {
	  fprintf(stderr,"ERROR: Deleterious mutation count exceeds limit\n");
	  fprintf(stderr,"ERROR: num_mutn = %d,max_del_mutn_perindiv/2 = %d\n", 
		 num_mutn,max_del_mutn_per_indiv/2);
	  fprintf(fp9,"ERROR: Deleterious mutation count exceeds limit\n");
	  exit(EXIT_FAILURE);
	}
	
	dmutn_offsprng[child][hap_id][0] = num_mutn;
	
	/* Insert new mutation such that mutations are maintained *
	 * in ascending order of their absolute value.            */
	
	i = num_mutn - 1;
	
	while(abs(dmutn_offsprng[child][hap_id][i]) > abs(mutn_indx) 
	      && i > 0) {
	  dmutn_offsprng[child][hap_id][i+1] = 
	    dmutn_offsprng[child][hap_id][i];
	  i--;
	}
	
	dmutn_offsprng[child][hap_id][i+1] = mutn_indx;
      }
    }
    
    
    /* If synergistic epistasis (SE) treatment is to be included, 
     * we must account for the effect of the interactions of this 
     * new mutation with the other mutations already present in the
     * linkage block.  We apply the following considerations to 
     * compute the mean SE effect resulting from these interactions. 
     * First, we require amplitude of the mean effect to be
     * _proportional_ to the mutation's effect, apart from SE, 
     * already computed above.  This means that if the new 
     * mutation's effect on the non-mutant genome is small, then the
     * SE contributions from interactions with other mutations is
     * assumed likewise to be small, and vice versa.  Next we require
     * that if _every_ site on the linkage block is mutant, then the
     * resulting total additional SE effect be equal to the new
     * mutation's effect computed above on the total fitness apart
     * from SE interactions.  In addition, we assume co-dominance, 
     * that is, 50% expression of the mutations base value.  These
     * considerations then imply that the mean SE effect the new
     * mutation adds through its interaction with each of the other
     * mutations already present in the linkage block is given by
     *
     * 0.5*fitness_effect/(haploid_genome_size/num_linkage_subunits)
     *
     * Note that we here account only for SE effects from linked
     * mutations, that is, those which occur on the same linkage 
     * block.  We treat the SE for the nonlinked mutation inter- 
     * actions elsewhere.  Also note that the SE effects from linked 
     * mutations are inherited, that is, passed from parent to
     * offspring.  Thus we add the SE effects to the linkage block
     * fitness.  Also note that we apply this SE treatment only to 
     * deleterious mutations. */
    
    if(synergistic_epistasis && !fav)
      se_effect = (0.5*fitness_effect*num_linkage_subunits
		   /haploid_genome_size)*
	offsprng_lb_mutn_count[child][0][hap_id][lb-1];
    else 
      se_effect = 0.;
    
    /* Recessive mutations (identified as such with a negative
     * mutation index) here incur only recessive_hetero_expression
     * times of their fitness effect, while dominant mutations incur
     * only dominant_hetero_expression times their fitness effect,
     * relative to the case of heterozygous expression.  The full 
     * fitness effect is realized only when a mutation occurs in  
     * both instances of its linkage block, that is, is homozygous.  */
    
    if(mutn_indx < 0) {
      fitness_effect *= recessive_hetero_expression;
    } else {
      fitness_effect *= dominant_hetero_expression;
    }
    
    fitness_effect += se_effect;
    
    /* Update linkage subunit fitness. */
    
    if(fav) {
      
      offsprng_lb_fitness[child][hap_id][lb-1] = 
	(offsprng_lb_fitness[child][hap_id][lb-1] + 
	 (1. - w)*fitness_effect)*((double)1.0 + w *fitness_effect);
      
      /* Increment the mutation count for the linkage subunit *
       * in which the mutation occurs.                        */
      
      offsprng_lb_mutn_count[child][1][hap_id][lb-1]++;
      
    } else {
      
      offsprng_lb_fitness[child][hap_id][lb-1] =
	(offsprng_lb_fitness[child][hap_id][lb-1] - 
	 (1. - w)*fitness_effect)*((double)1.0 - w *fitness_effect);
      
      /* Increment the mutation count for the linkage subunit *
       * in which the mutation occurs.                        */
      
      offsprng_lb_mutn_count[child][0][hap_id][lb - 1]++;
      
    }
    
  }
}

