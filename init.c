#include "mendel.h"
void output_fitness_effect_distribution(FILE*,float,float,float,float);

/* This routine generates a random mutation in a randomly chosen  *
 * linkage block with a randomly chosen haploid identity in a     *
 * randomly chosen individual in the population and modifies the  *
 * linkage block fitness to reflect the resulting fitness change. */
void favorable_mutn(int *** fmutn , int **** lb_mutn_count, 
                    double *** linkage_block_fitness)
{
  double fitness_gain;
  float w, x;
  int i, lb, hap_id, mutn, mutn_indx, num_mutn, j;
  
  /* Specify the new random favorable mutation. */
  
  /* Generate the index of the random individual. */
  x = randomnum();
  i  = min(current_pop_size-1,  abs((int)(current_pop_size*x)));
  
  /* Generate the linkage block index. */
  
  x = randomnum();
  lb = min(num_linkage_subunits-1,  abs((int)(num_linkage_subunits*x)));
  
  /* Generate the haploid identity. */
  
  x = randomnum();
  //hap_id = min(1, abs((int)(2.*x)));
  hap_id = 1 - min(1, abs((int)(2.*x)));
  
  /* Generate a random index mutn to specify the fitness effect
     associated with the mutation. */
  
  x = randomnum();
  
  if(tracking_threshold != 1.0) {
    mutn = min(lb_modulo-2, (int)(x/fav_scale));
  } else {
    mutn = 1;
  }
  
  /* Add an offset to assign it to the appropriate linkage block. */
  
  //mutn_indx = mutn + (lb - 1)*lb_modulo;
  mutn_indx = mutn + lb*lb_modulo;
  
  /* Specify whether the mutation is dominant or recessive.
     (Recessives have a negative mutation index.)  */
  
  if(fraction_recessive > 0.) {
    if(randomnum() < fraction_recessive) mutn_indx = -mutn_indx;
  }
  
  /* Increment the favorable mutation count for the appropriate
     individual, linkage block, and haploid index. */
  
  lb_mutn_count[i][1][hap_id][lb]++;
  
  /* Compute the fitness factor associated with this new mutation.
   * Incorporate this fitness contribution into the fitness of the
   * the appropriate linkage block. 
   *
   * When parameter fitness_distrib_type is 1, the fitness
   * factor f is obtained from the mutation index mutn using a
   * distribution function of the form
   * 
   * f = (1. + max_fav_fitness_gain*exp(-alpha_fav*pow(x,gamma_fav))
   * 
   * where max_fav_fitness_gain is an input parameter, alpha_fav is 
   * is log(haploid_genome_size*max_fav_fitness_gain) and x is a 
   * random number uniformly distributed between zero and one.
   * 
   * When parameter fitness_distrib_type is 0, the fitness
   * factor is constant and given by the expression    
   * 
   * f = 1. + max_fav_fitness_gain*uniform_fitness_effect             */
  
  if(fitness_distrib_type > 0) {
    //fitness_gain = max_fav_fitness_gain*exp(pow(-alpha_fav*x,(int)gamma_fav));
    fitness_gain = max_fav_fitness_gain*exp(-alpha_fav*pow(x,(int)gamma_fav));
  } else {
    fitness_gain = max_fav_fitness_gain*uniform_fitness_effect;
  }
  
  /* Track this mutation if its fitness gain exceeds the value of
   * tracking_threshold. */
  
  if(fitness_gain > tracking_threshold) {
    
    /* Test to see if the storage limit of array fmutn has been
     * exceeded.  (Note that we are using the first slot to hold the
     * actual mutation count.) */
    
    num_mutn = fmutn[i][hap_id][0] + 1;
    printf("num_mutn is %d\n",num_mutn);
    
    if(num_mutn + 1 >= max_fav_mutn_per_indiv/2) {
      fprintf(stderr,"Favorable mutations exceed the storage limit\n");
      fprintf(fp9,"Favorable mutations exceed the storage limit\n");
      exit(EXIT_FAILURE);
    }
    
    fmutn[i][hap_id][0] = num_mutn;
    
    /* Insert new mutation such that mutations are maintained
     * in ascending order of their absolute value.  */
    
    j = num_mutn;
    
    while(abs(fmutn[i][hap_id][j]) > abs(mutn_indx) && j > 1) {
      fmutn[i][hap_id][j+1] = fmutn[i][hap_id][j];
      j--;
    }
    
    fmutn[i][hap_id][j+1] = mutn_indx;
    
  }
  
  /* Recessive mutations (identified as such with a negative
   * mutation index) here incur only recessive_hetero_expression
   * times of their fitness gain, while dominant mutations incur
   * only dominant_hetero_expression times their fitness gain.
   * The full fitness gain is realized only when a mutation occurs 
   * in both instances of its linkage block, that is, is homozygous.  */
  
  if(mutn_indx < 0) {
    fitness_gain = recessive_hetero_expression*fitness_gain;
  } else {
    fitness_gain =  dominant_hetero_expression*fitness_gain;
  }
  
  w = multiplicative_weighting;
  
  linkage_block_fitness[i][hap_id][lb] =
    (linkage_block_fitness[i][hap_id][lb] + (1. - w)*fitness_gain)
    * ((double)1 + w *fitness_gain);
}


/* This routine generates a small number (no larger than the number
 * of linkage subunits) of paired alleles, with a random fitness
 * effect on one haplotype set and an effect with the same magnitude
 * but the opposite sign on the other haplotype set.  Variation of 
 * of fitness effect is according to a uniform random distribution
 * with a user-specified mean value.  */
void gen_initial_contrasting_alleles(int *** dmutn, int *** fmutn,
				     double *** linkage_block_fitness, 
				     float * initial_allele_effects, 
				     int max_size)
{
  
  float w, x, effect, expressed;
  int h1_id, h2_id, lb, mutn, mutn_indx;
  int i, m, n, nskp;
  
  w = multiplicative_weighting;
  
  if(num_contrasting_alleles > 0) {
    num_contrasting_alleles = min(num_linkage_subunits,
                                  num_contrasting_alleles);
    nskp = num_linkage_subunits/num_contrasting_alleles;
  } else return; 
  
  for(n=0; n < num_linkage_subunits ; n++) 
    initial_allele_effects[n] = 0.;
  
  for(n=0; n < num_contrasting_alleles; n++) {
    
    lb = 1 + n*nskp;
    
    /* Use the same mutation effect index for all of these paired 
     * alleles. This index, lb_modulo-1, is reserved exclusively for
     * these alleles.  When treated as an ordinary mutation, the 
     * fitness effect it would imply is the smallest effect possible.
     * The fitness effects associated with these alleles, however, 
     * are handled via the linkage_block_fitness array.  We tag these
     * alleles with a mutation index to be able to track them over
     * successive generations and obtain statistics on them at the 
     * end of a run.  */
    
    mutn = lb_modulo - 1;
    
    /* Add an offset to assign it to the appropriate linkage block. */
    
    mutn_indx = mutn + (lb - 1)*lb_modulo;

    /* Generate random haplotype identities. */
    
    x = randomnum();
    h1_id = min(1,(int)(2*x));
    h2_id = 1 - h1_id;
    
    m = dmutn[0][h1_id][0] + 1;
    for (i = 0; i < max_size; i++) {
      dmutn[i][h1_id][m+1] = mutn_indx;
      dmutn[i][h1_id][0] = m;
    }

    m = fmutn[0][h2_id][0] + 1;
    for (i = 0; i < max_size; i++) {
      fmutn[i][h2_id][m+1] = mutn_indx;
      fmutn[i][h2_id][0] = m; 
      //printf("fmutn[%d][%d][%d] = %d\n",i,h2_id,m+1,mutn_indx);
    }
    
    /* Generate the uniformly distributed random fitness effect
     * associated with the allele pair.  */
    
    effect = 2.*initial_alleles_mean_effect*randomnum();
    if(num_contrasting_alleles < 11) 
      effect = initial_alleles_mean_effect;
    
    /* Store the value of the fitness effect for each allele pair in
     * array initial_allele_effects. */
    
    initial_allele_effects[lb] = effect;
    
    /* We assume deleterious alleles behave in a recessive manner
     * and when heterozygous have an effect given by the allele
     * fitness effect multiplied by recessive_hetero_expression.
     * Similarly, we assume favorable alleles behave in a dominant
     * manner and when heterozygous have an effect given by the 
     * allele fitness effect times dominant_hetero_expression.  The
     * full allele fitness effect is realized only when the same
     * version of the allele occurs on both instances of its linkage 
     * block, that is, is homozygous.
     * Apply the appropriate fitness effects to the appropriate
     * linkage blocks.  */
    
    expressed = recessive_hetero_expression*effect;
    for (i = 0; i < pop_size; i++) {
      linkage_block_fitness[i][h1_id][lb] = 
	((double)1. - (1. - w)*expressed)
	* ((double)1. - w *expressed);
    }
    
    expressed =  dominant_hetero_expression*effect;
    for (i = 0; i < pop_size; i++) {
      linkage_block_fitness[i][h2_id][lb] = 
	((double)1. + (1. - w)*expressed)
	* ((double)1. + w *expressed);
    }
    
  }
  
}

void initialize(char myid_str[])
{
  FILE *myfp;
  int i, k;
  float d1, d2, sum, del_mean, fav_mean, alpha, mygamma;
  char zero_str[3];
  char filename[80];
  
  time_t curtime;
  struct tm *loctime;
  static long iseed1, iseed2;
  
  /* Get the current time. */
  curtime = time (NULL);
  /* Convert it to local time representation. */
  loctime = localtime (&curtime);

  if(dynamic_linkage && num_linkage_subunits < 50) {
     fprintf(stderr,"ERROR: dynamic_linkage requires more than 50 linkage\n");
     fprintf(stderr,"       subunits.  Either turn dynamic_linkage off, \n");
     fprintf(stderr,"       Or increase the number of linkage subunits.\n");
     exit(EXIT_FAILURE);
  }
  
  if(is_parallel) {
    sprintf(myid_str,"%.3d",myid+1);
    
    if (!homogenous_tribes) {
      fp5 = fopen("mendel.in", "r");
      if (fp5 == NULL) {
	fprintf(stderr,"Cannot read mendel.in\n");
	exit(EXIT_FAILURE);
      }
      read_parameters(fp5);
      fclose(fp5);
    }
    
    if (myid==0) printf("subpopulation size is %d\n", pop_size);
    
    if(num_indiv_exchanged > pop_size) {
      fprintf(fp6,"ERROR: num_indiv_exchanged >= tribal pop_size\n");
      fprintf(fp6,"ERROR: decrease num_indiv_exchanged.\n");
      mpi_myabort();
    } 
  }else {
    sprintf(myid_str,"%.3d",0);
    myid = 0;
  }
  
  /* Output version information.  RCS will automatically update *
   * the following $Id string on check-in                      */
  
  sprintf(version_init,"$Id: init.c,v 1.60 2009/05/28 21:23:05 wes Exp $");

  sprintf(filename,"%s.%s.hap",case_id,myid_str);
  fp4 = fopen(filename, "w");
  
  sprintf(filename,"%s.%s.hst",case_id,myid_str);
  if (restart_case) {
    fp7 = fopen(filename, "a");
  } else {
    fp7 = fopen(filename, "w");
    fprintf(fp7,"# generation      fitness       fitness_sd");
    fprintf(fp7," num_dmutns     num_fmutns     pop_size\n");
  }
  
  sprintf(filename,"%s.%s.dst",case_id,myid_str);
  fp8 = fopen(filename, "w");
  
  sprintf(filename,"%s.%s.out",case_id,myid_str);
  fp9 = fopen(filename, "w");
  
  /* This file is monitored during each main generation         *
   * loop.  If the value is 0, the program keeps running.       *
   * If the value gets changed to 1 (by user or user interface) *
   * then the program will write a restart dump and exit.       *
   * This is for nicely shutting down the program.              */

  sprintf(filename,"%s.%s.st8",case_id,myid_str);
  myfp = fopen(filename, "w");
  if (myfp == NULL) { 
    fprintf(stderr,"ERROR: cannot open st8 file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(myfp,"%d\n",0); 
  fclose(myfp);
  
  sprintf(filename,"%s.%s.plm",case_id,myid_str);
  fp11 = fopen(filename, "w");

  // Maintain a separate file for windows which just has a snapshot
  // of the most recent polymorphism information.
  sprintf(filename,"%s.%s.pls",case_id,myid_str);
  fp13 = fopen(filename, "w");

  sprintf(filename,"%s.%s.pmd",case_id,myid_str);
  fp19 = fopen(filename, "w"); 

#ifdef PPM
  sprintf(filename,"%s.%s.ppm",case_id,myid_str);
  fp15 = fopen(filename, "w");
  
  // Write header for PPM file
  fprintf(fp15,"P3\n");
  fprintf(fp15,"%d %d\n", pop_size, num_generations); 
  fprintf(fp15,"%d\n",255);
#endif
  
  sprintf(filename,"%s.%s.fit",case_id,myid_str);
  fp16 = fopen(filename, "w");
  
  sprintf(filename,"%s.%s.tim",case_id,myid_str);
  fp22 = fopen(filename, "w");
  
  fprintf(fp22,"# generation  time_per_gen(s) "); 
  fprintf(fp22,"time_offspring(s)  time_selection(s)\n");
  
  sprintf(filename,"%s.%s.sel",case_id,myid_str);
  fp24 = fopen(filename, "w");

  sprintf(filename,"%s.%s.thr",case_id,myid_str);
  fp25 = fopen(filename, "w");
  fprintf(fp25,"# generation\tselection threshold\n");
  fprintf(fp25,"#  generation\tselection thresholds\n");
  fprintf(fp25,"#\tdel_dom_thres  del_rec_thres  fav_dom_thres  fav_rec_thres\n");
  
  sprintf(filename,"%s.%s.acc",case_id,myid_str);
  fp26 = fopen(filename, "w");

  /* If parallel, write additional average files with name-tag 000. */
  
  if (is_parallel) {
    
    sprintf(zero_str,"%.3d",0);
    
    sprintf(filename,"%s.%s.hap",case_id,zero_str);
    fp14 = fopen(filename, "w");
    
    sprintf(filename,"%s.%s.hst",case_id,zero_str);
    fp17 = fopen(filename, "w");
    
    sprintf(filename,"%s.%s.dst",case_id,zero_str);
    fp18 = fopen(filename, "w");
    
    sprintf(filename,"%s.%s.plm",case_id,zero_str);
    fp21 = fopen(filename, "w");
    
    sprintf(filename,"%s.%s.tim",case_id,zero_str);
    fp23 = fopen(filename, "w");
    fprintf(fp23,"# generation  time_per_gen(s)  ");
    fprintf(fp23,"time_offspring(s)  time_selection(s)\n");

    sprintf(filename,"%s.%s.sel",case_id,zero_str);
    fp34 = fopen(filename, "w");

    sprintf(filename,"%s.%s.thr",case_id,zero_str);
    fp35 = fopen(filename, "w");

    fprintf(fp35,"#  generation                   selection thresholds\n");
    fprintf(fp35,"#             del_dom_thres  del_rec_thres  ");
    fprintf(fp35,              "fav_dom_thres  fav_rec_thres\n");

  }
  
  if(myid==0) printf("Case started: %s\n",asctime(loctime));
  
  fprintf(fp9,"Case started: %s\n",asctime(loctime));
  
  /* Check to ensure parameter fraction_random_death does not 
   * reduce actual fertility below two per female.  */
  
  if(offspring_per_female*(1. - fraction_random_death) < 2.) {
    fprintf(stderr,"ERROR: Input value of fraction_random_death \n");
    fprintf(stderr,"implies a fertility of less than 2 per female.\n");
    fprintf(fp9,"ERROR: Input value of fraction_random_death \n");
    fprintf(fp9,"implies a fertility of less than 2 per female.\n");
    exit(EXIT_FAILURE); 
  }
  
  /* Limit the minimum value of heritability to be 10**-20. */
  
  heritability = max(1.e-20, heritability);
  
  if(!dynamic_linkage) haploid_chromosome_number = 1;
  
  /* Initialize some quantities from the input values. */
     
  /* For clonal reproduction, set the number of chromosomes to one
   * and the number of linkage blocks to one.  */
  
  if(clonal_reproduction) {
    haploid_chromosome_number  = 1;
    num_linkage_subunits       = 1;
    dominant_hetero_expression = 1.;
    fraction_recessive         = 0.;
  }
  
  /* For the case of dynamic linkage, ensure that the number of linkage 
     subunits is an integer times the haploid chromosome number.  */
  
  if(dynamic_linkage) num_linkage_subunits = 
		     (num_linkage_subunits/haploid_chromosome_number)*
		      haploid_chromosome_number;
  
  /* Echo the parameters for this case to the output file. */
  
  if(is_parallel) {
     if(myid==0) write_parameters(stdout);
  } else {
     write_parameters(stdout);
  }
  write_parameters(fp9);
  
  /* When tracking_threshold is input as zero, this means that all      *
   * mutations are to be tracked.  In this case, set tracking_threshold *
   * to be equal to the minimum fitness value to prevent various        *
   * numerical overflow problems that would arise otherwise.            */
  
  //if (tracking_threshold == 0.) tracking_threshold = 1.e-8;
  ////tracking_threshold = 1./haploid_genome_size;
  tracking_threshold = max(1./haploid_genome_size,tracking_threshold); 
  
  lb_modulo  = (pow(2,30)-2)/num_linkage_subunits;
  
  alpha_del  = log(haploid_genome_size);
  if(max_fav_fitness_gain > 0.) {
     alpha_fav  = log(haploid_genome_size*max_fav_fitness_gain);
  } else {
     alpha_fav = alpha_del;
  }

  gamma_del = log(-log((double)high_impact_mutn_threshold)/alpha_del) 
              /log((double)high_impact_mutn_fraction);

  gamma_fav = log(-log((double)high_impact_mutn_threshold)/alpha_fav) 
    /log((double)high_impact_mutn_fraction);
  
  if(tracking_threshold != 1. && tracking_threshold > 0.) {
    del_scale = exp(log(-log(tracking_threshold)/alpha_del)
                /gamma_del)/(lb_modulo-2);
    if(max_fav_fitness_gain > 0.)
       fav_scale = exp(log(-log(tracking_threshold
                   /max_fav_fitness_gain)/alpha_fav)
                   /gamma_fav)/(lb_modulo-2);
    else fav_scale = 0.;

  } else if(tracking_threshold == 1.0) {
      del_scale = 0.;
      fav_scale = 0.;
  } else {
      del_scale = 1./lb_modulo;
      fav_scale = 1./lb_modulo;
  }

  /* Compute mean absolute fitness effect for deleterious mutations.*/
  sum = 0.;
  d2  = 1.;
  
  for (i=1; i <= 1000000; i++) {
    d1 = d2;
    d2 = exp(-alpha_del*pow(0.000001*i,gamma_del));
    sum += d1 + d2;
  }

  del_mean = 0.0000005*sum;

  // Compute mean absolute fitness effect for favorable mutations.
 
  sum = 0.;
  d2  = 1.;
 
  for(i=1; i<=1000000; i++) {
    d1 = d2;
    d2 = exp(-alpha_fav*pow(0.000001*i,gamma_fav));
    sum += d1 + d2;
  }
  
  fav_mean = 0.0000005*sum*max_fav_fitness_gain;
 
  if(max_fav_fitness_gain > 0.) {
     alpha = alpha_fav;
     mygamma = gamma_fav;
  } else {
     alpha = 0.;
     mygamma = 0.;
  }
  
  if(myid==0) output_fitness_effect_distribution(stdout,del_mean,fav_mean,
                                                 alpha,mygamma);
  output_fitness_effect_distribution(fp9,del_mean,fav_mean,alpha,mygamma);
  
  if (myid == 0) {
     printf(" Tracking threshold = %g\n", tracking_threshold);
     printf(" Fraction deleterious mutations tracked = %g\n", 
			del_scale*(lb_modulo-2));
     printf(" Fraction favorable   mutations tracked = %g\n",
                        fav_scale*(lb_modulo-2));
  }
     
  fprintf(fp9," Tracking threshold = %g\n", tracking_threshold);
  fprintf(fp9," Fraction deleterious mutations tracked = %g\n",  
	        del_scale*(lb_modulo-2));
  fprintf(fp9," Fraction favorable   mutations tracked = %g\n",  
	        fav_scale*(lb_modulo-2));
  
  /* Impose a reasonable limit of the number of tracked mutations 
   * based on the number of new mutations per offspring, the number
   * of generations, and the fraction of deleterious mutations tracked.
   * If the run is a restart run, double the number again.  Limit the 
   * number by the input value for max_tracted_mutn_per_indiv.  */
  
  k = 1.8*new_mutn_per_offspring*num_generations*del_scale*lb_modulo;
  if (restart_case) k = 2*k;

  // Temporarily remove this constraint imposed by k to allow larger
  // numbers of favorable mutations.
  //max_tracked_mutn_per_indiv = min(k ,max_tracked_mutn_per_indiv );
  
  max_del_mutn_per_indiv = max(100, (int)(max_tracked_mutn_per_indiv
					  *(1. - frac_fav_mutn)));
  max_fav_mutn_per_indiv = max(100, (int)(max_tracked_mutn_per_indiv
					  *frac_fav_mutn));
  
  max_del_mutn_per_indiv = max_del_mutn_per_indiv 
    + 2*num_contrasting_alleles;
  max_fav_mutn_per_indiv = max_fav_mutn_per_indiv 
    + 2*num_contrasting_alleles;

  /* Prevent array overflow for cases with large numbers of initial
   * favorable mutations.  */
  
  max_fav_mutn_per_indiv = max(max_fav_mutn_per_indiv,
			       30*num_initial_fav_mutn/pop_size);
  
  /* When the parameter tracking_threshold is set to one, it signals
   * that no mutations are to be tracked.  In this case, the size of
   * the dmutn and fmutn arrays can be reduced to a minimum.  */
  
  if(tracking_threshold == 1.0) {
    max_del_mutn_per_indiv = 4;
    max_fav_mutn_per_indiv = 4;
    fraction_recessive     = 0.;
    dominant_hetero_expression  = 0.5;
    recessive_hetero_expression = 0.5;
  }
  
  if(myid == 0) {
     printf(" Maximum  deleterious mutations tracked = %d\n", 
                        max_del_mutn_per_indiv);
     printf(" Maximum  beneficial mutations tracked = %d\n", 
                        max_fav_mutn_per_indiv);
  }
  fprintf(fp9," Maximum  deleterious mutations tracked = %d\n",
               max_del_mutn_per_indiv);
  fprintf(fp9," Maximum  beneficial mutations tracked = %d\n",
               max_fav_mutn_per_indiv);
  
  /* Initialize random number generator. */
  iseed1 = abs(random_number_seed)+myid;
  iseed2 = iseed1;
  setall(iseed1,iseed2);

  //printf("RANDOM NUMBER SEED IS: %ld %d\n",iseed1,myid);
  
  //c-built-in random number generator
  srand(iseed1);
  
  poisson_mean = new_mutn_per_offspring;
  
}

/* Print some properties of the fitness effect distribution. */
void output_fitness_effect_distribution(FILE *fp, float del_mean, 
                                        float fav_mean, float alpha, 
                                        float mygamma) 
{
  
  fprintf(fp,"  Properties of the Weibull fitness effect distribution function:\n\n");

  fprintf(fp,"              e(x) = exp(-alpha*x**gamma), 0 < x < 1\n\n");
  
  fprintf(fp,"  genome_size       = %8.2e\n",haploid_genome_size);
  fprintf(fp,"  e_high_impact     = %8.2e",high_impact_mutn_threshold);
  fprintf(fp," defining value of *high impact* mutation\n");
  fprintf(fp,"  frac_high_impact  = %8.2e",high_impact_mutn_fraction);
  fprintf(fp," fraction *high impact* mutations of total\n\n");
  
  fprintf(fp,"  mutation:             deleterious         favorable\n");
  fprintf(fp,"  alpha             = %12.5f\t%12.5f log(genome_size)\n",alpha_del,alpha);
  fprintf(fp,"\t\t\t\t\t\t\t for deleterious case\n");
  fprintf(fp,"  gamma             = %12.6f\t%12.6f\n\n",gamma_del,mygamma);
  
  fprintf(fp,"  mean   effect     = %12.2e\t%12.2e\n", -del_mean, fav_mean);
  fprintf(fp,"  median effect     = %12.2e\t%12.2e (x = 0.5)\n",
	  -expf(-alpha_del*pow(0.5,gamma_del)), 
	   expf(-alpha_fav*pow(0.5,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  10th   percentile = %12.2e\t%12.2e (x = 0.9)\n",
	  -expf(-alpha_del*pow(0.9,gamma_del)), 
	   expf(-alpha_fav*pow(0.9,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  20th   percentile = %12.2e\t%12.2e (x = 0.8)\n",
	  -expf(-alpha_del*pow(0.8,gamma_del)), 
	   expf(-alpha_fav*pow(0.8,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  30th   percentile = %12.2e\t%12.2e (x = 0.7)\n",
	  -expf(-alpha_del*pow(0.7,gamma_del)), 
	   expf(-alpha_fav*pow(0.7,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  40th   percentile = %12.2e\t%12.2e (x = 0.6)\n",
	  -expf(-alpha_del*pow(0.6,gamma_del)), 
	   expf(-alpha_fav*pow(0.6,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  50th   percentile = %12.2e\t%12.2e (x = 0.5)\n",
	  -expf(-alpha_del*pow(0.5,gamma_del)), 
	   expf(-alpha_fav*pow(0.5,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  60th   percentile = %12.2e\t%12.2e (x = 0.4)\n",
	  -expf(-alpha_del*pow(0.4,gamma_del)), 
	   expf(-alpha_fav*pow(0.4,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  70th   percentile = %12.2e\t%12.2e (x = 0.3)\n",
	  -expf(-alpha_del*pow(0.3,gamma_del)), 
	   expf(-alpha_fav*pow(0.3,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  80th   percentile = %12.2e\t%12.2e (x = 0.2)\n",
	  -expf(-alpha_del*pow(0.2,gamma_del)), 
	   expf(-alpha_fav*pow(0.2,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  90th   percentile = %12.2e\t%12.2e (x = 0.1)\n",
	  -expf(-alpha_del*pow(0.1,gamma_del)), 
	   expf(-alpha_fav*pow(0.1,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  99th   percentile = %12.2e\t%12.2e (x = 0.01)\n",
	  -expf(-alpha_del*pow(0.01,gamma_del)), 
	   expf(-alpha_fav*pow(0.01,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  99.9   percentile = %12.2e\t%12.2e (x = 0.001)\n",
	  -expf(-alpha_del*pow(0.001,gamma_del)), 
	   expf(-alpha_fav*pow(0.001,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  99.99  percentile = %12.2e\t%12.2e (x = 0.0001)\n",
	  -expf(-alpha_del*pow(0.0001,gamma_del)), 
	   expf(-alpha_fav*pow(0.0001,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  99.999 percentile = %12.2e\t%12.2e (x = 0.00001)\n\n",
	  -expf(-alpha_del*pow(0.00001,gamma_del)), 
	   expf(-alpha_fav*pow(0.00001,gamma_fav))*max_fav_fitness_gain); 
  fprintf(fp,"  Notes:\n");
  fprintf(fp,"  (1) The e(x) values above are for a homozygous pair of mutations.\n");
  fprintf(fp,"  (2) For favorables, e(x) also includes the factor ");
  fprintf(fp,"max_fav_fitness_gain.\n\n");
}

