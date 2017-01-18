#include "mendel.h"

int main(int argc, char* argv[]) {
  int i, j, k, l,l1,l2, n, gen, dad, mom, child, slot;
  int max_size, this_size, gen_0, shutdown_gen=0, run_status=0;
  int total_offspring, actual_offspring;
  int offspring_count, empty, replace, parent;
  int tracked_fav_mutn=0, ica_count[2], cumulative_offspring=0;
  float num_offspring, d, x;
  float fitness_adjusted_offspring;
#ifdef DYNAMIC_POP_SIZE
  int max_pop_size;
  double accum[50]={0}; 
  float carrying_capacity;
#endif
  TIME_T start_gen, stop_gen;
  TIME_T start, stop, start_mendel, stop_mendel;
  double time_offspring, time_selection, time_per_gen;
  const float extinction_fitness_threshold = 0.1;
  bool print_flag;
  char myid_str[4] = "000\0";
  char filename[40];
  int *** dmutn = NULL;
  int *** fmutn = NULL;
  int *** dmutn_offsprng = NULL;
  int *** fmutn_offsprng = NULL;
  int **** lb_mutn_count = NULL;
  int **** offsprng_lb_mutn_count = NULL;
  double *** linkage_block_fitness = NULL;
  double *** offsprng_lb_fitness = NULL;
  float  * initial_allele_effects = NULL;
  double * pheno_fitness = NULL;
  double * fitness = NULL;
  double * work_fitness = NULL; 
  double * sorted_score = NULL;
#ifdef PPM
  int red, green, blue;
#endif
#ifdef MPICH
  double par_time_per_gen;
  double par_time_offspring, par_time_selection;
#endif
#ifdef WRITE_SAMPLE
  float ***defect;
  float ***improve;
  float *effect;
#endif
  FILE *fp;
  
  TIME(&start_mendel);
  
  /* Read input parameters from input file mendel.in */
  
  fp = fopen("mendel.in", "r");
  if (fp == NULL) {
    printf("Cannot read mendel.in\n");
    return 0;
  }
  
  read_parameters(fp);
  fclose(fp);
  
  /* Perform certain initializations including opening files */
  
  if(is_parallel) {
    if(num_indiv_exchanged > pop_size) {
      fprintf(stderr,"ERROR: num_indiv_exchanged > pop_size.\n");
      exit(EXIT_FAILURE);
    }
    myid = mpi_myinit(argc,argv);
  }
  
  TIME(&start);
  initialize(myid_str);
  TIME(&stop);
  sec[3] += DIFFTIME(stop,start);
 
  sprintf(version_mendel,
          "$Id: mendel.c,v 1.66 2009/06/04 22:32:11 wes Exp $");
  
  if(pop_growth_model > 0) { 
     /* Pass in max_pop_size through num_generations input parameter */
     max_pop_size = num_generations;
     max_size = (int) (0.55*offspring_per_female*max_pop_size
   		       *(1. - fraction_random_death));
  } else {
     max_size = (int) (0.55*offspring_per_female*pop_size 
		    *(1. - fraction_random_death));
  }
  //printf("max_size is %d\n",max_size);
  
  n = offspring_per_female*(1. - fraction_random_death) + 0.999;
  
  /* Allocate memory for large arrays */
  dmutn = (int ***)malloc_int3d(max_size,2,max_del_mutn_per_indiv/2);
  fmutn = (int ***)malloc_int3d(max_size,2,max_fav_mutn_per_indiv/2);
  dmutn_offsprng = (int ***)malloc_int3d(n,2,max_del_mutn_per_indiv/2);
  fmutn_offsprng = (int ***)malloc_int3d(n,2,max_fav_mutn_per_indiv/2);  
  lb_mutn_count = (int ****)malloc_int4d(max_size,2,2,num_linkage_subunits);
  offsprng_lb_mutn_count = (int ****)malloc_int4d(n,2,2,num_linkage_subunits);
  linkage_block_fitness = (double ***)malloc_double3d(max_size, 2, 
						      num_linkage_subunits);
  offsprng_lb_fitness = (double ***)malloc_double3d(n, 2, 
						    num_linkage_subunits);
  initial_allele_effects = (float *)malloc_float1d(num_linkage_subunits);
  pheno_fitness = (double *)malloc_double1d(max_size);
  fitness = (double *)malloc_double1d(max_size);
  work_fitness = (double *)malloc_double1d(max_size);
  
  sorted_score = (double *)malloc_double1d(max_size);
  available = (bool *)malloc_int1d(pop_size);
  replaced_by_offspring = (bool *)malloc_int1d(pop_size);
  
  /* If this is a restart case, read the restart dump file and       *
   * set the current dump number to the restart dump number.         *
   * Otherwise, set it to zero.  The variable gen_0 is the initial   *
   * generation number, retrieved from the restart dump in a restart *
   * case and zero otherwise.                                        */
  
  if(restart_case) {
    read_restart_dump(dmutn,fmutn,lb_mutn_count,
		      linkage_block_fitness, initial_allele_effects,
		      &gen_0,max_size,myid_str);
    dump_number = restart_dump_number;
  } else {
    gen_0 = 0;
    dump_number = 0;
  }
  
  
  
  /* If the bottleneck flag, bottleneck_yes, is false, set the value
   * of bottleneck_generation beyond the generation range for this run. */
  
  if(!bottleneck_yes) bottleneck_generation = 1 + gen_0 + num_generations;
  
  /* Initialize the population size to be equal to the parameter
   * pop_size unless the parameter bottleneck_generation has the
   * value zero.  In the latter case, initialize the population size 
   * to bottleneck_pop_size. */
  
  if(bottleneck_generation > 0) current_pop_size = pop_size;
  else current_pop_size = bottleneck_pop_size;
  
  /* If not a restart case, initialize entire population to have no *
   * initial mutations.                                             *
   *                                                                *
   * Initialize the linkage block fitness such that all individuals *
   * in the population have identical haplotypes.  If initial       *
   * contrasting alleles are to be included, generate them here.    */
  
  if(!restart_case) {
    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++) 
	for (k = 0; k < max_del_mutn_per_indiv/2; k++)
	  dmutn[i][j][k] = num_linkage_subunits*lb_modulo+1;

    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++)
	for (k = 0; k < max_fav_mutn_per_indiv/2; k++)
	  fmutn[i][j][k] = num_linkage_subunits*lb_modulo+1;

    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++)
	dmutn[i][j][0] = 0;
    
    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++)
	fmutn[i][j][0] = 0;
    
    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++)
	for (k = 0; k < 2; k++)
	  for (l = 0; l < num_linkage_subunits; l++)
	    lb_mutn_count[i][j][k][l] = 0;
    
    for (i = 0; i < max_size; i++)
      for (j = 0; j < 2; j++)
	for (k = 0; k < num_linkage_subunits; k++)
	  linkage_block_fitness[i][j][k] = (double)1; 
    
    if(num_contrasting_alleles > 0)
      gen_initial_contrasting_alleles(dmutn,fmutn,linkage_block_fitness, 
     				      initial_allele_effects, max_size);
  }  

  /* Generate num_initial_fav_mutn random initial favorable mutations */
  
  for (k = 1; k <= num_initial_fav_mutn; k++)
    favorable_mutn(fmutn,lb_mutn_count,linkage_block_fitness);
  
  post_sel_fitness = (double) 1;
  ica_count[0] = 0;
  ica_count[1] = 0;
  
 /*******************************************************
  * Step population through num_generations generations *
  *******************************************************/ 
  TIME(&stop);
  sec[0] += DIFFTIME(stop,start);
  
  for(gen = gen_0+1; gen <= gen_0+num_generations; gen++) {
    
    msg_num = 1;
    TIME(&start_gen);
    
    /* If the generation number lies within the bottleneck interval, *
     * set the current population size to bottleneck_pop_size.       */
    
    if (gen >= bottleneck_generation && gen < 
	(bottleneck_generation + num_bottleneck_generations)) 
      current_pop_size = bottleneck_pop_size;
    
#ifdef MPICH
    /* Move individuals between tribes/processors */
    
    //printf("myid is: %d gen is: %d\n",myid,gen);
    if (is_parallel && mod(gen,migration_generations)==0) {
      TIME(&start); 
      mpi_migration(dmutn,fmutn,linkage_block_fitness, 
		    lb_mutn_count,gen,myid,msg_num,ierr,available);
      TIME(&stop);
      sec[2] += DIFFTIME(stop,start);
    }
    //printf("done migrating\n");
#endif
    
    for (i = 0; i < pop_size; i++) { 
      available[i] = true;
      replaced_by_offspring[i] = false;
    }
    
    total_offspring = current_pop_size;
    offspring_count = 0;
    num_offspring   = offspring_per_female *(1. - fraction_random_death);
    
    if((gen - gen_0)%10 == 1 || gen - gen_0 <= 10) {
      cumulative_offspring = 0;
      new_mutn_count       = 0;
    }
    
    TIME(&start);
    for (i = 0; i < current_pop_size/2; i++) {

      /* Randomly mate one half of the population with members *
       * from the other half.                                  */
      
      x = randomnum();
      dad = min(current_pop_size-1, 
		abs((int)(current_pop_size*x)));
      
      while(!available[dad]){
        if (dad%(current_pop_size - 1) != 0 || dad == 0){
	  dad = dad%(current_pop_size - 1) + 1;
	} else {
	  dad = 0;
	}
      } 
      available[dad] = false;
      
      x = randomnum();
      mom = min(current_pop_size-1, 
		abs((int)(current_pop_size*x)));
      
      while(!available[mom]){
	if (mom%(current_pop_size -1) != 0 || mom == 0){ 
	  mom = mom%(current_pop_size - 1) + 1;
	} else {
	  mom = 0;
	}
      }
      available[mom] = false;
      
      /* Generate an appropriate number offspring from *
       * the two parents                               */
      
      x = randomnum();
      if(fitness_dependent_fertility) {
	fitness_adjusted_offspring = 
	  num_offspring*sqrt(min((double)1.0, post_sel_fitness));
	actual_offspring = (int)fitness_adjusted_offspring;
	if(fitness_adjusted_offspring - actual_offspring > x) 
	  actual_offspring++;
      } else {
	actual_offspring = (int)num_offspring;
	if((num_offspring - (int)num_offspring) > x)
	  actual_offspring++;
      }
      
      if(i == 1) actual_offspring = max(1, actual_offspring);
      
      actual_offspring = min((int)num_offspring + 1, actual_offspring);
      
      /* If the parameter fraction_self_fertilization is non-zero *
       * (as can be the case for many types of plants), implement *
       * self-fertilization for the appropriate portion of the    *
       * population by calling routine offspring using 'dad' for  *
       * both parents for half the offspring and 'mom' for both   *
       * parents for the other half.                              */
      
      for(child=0; child < actual_offspring; child++) {
        x = randomnum();
	if(fraction_self_fertilization < x && !clonal_reproduction) {
	  
	  /* This call generates an offspring from a sexual *
	   * union of the two individuals dad and mom.      */
	  
	  offspring(dmutn_offsprng, fmutn_offsprng, 
		    offsprng_lb_mutn_count, 
		    offsprng_lb_fitness, dmutn, fmutn, lb_mutn_count, 
		    linkage_block_fitness, dad, mom, child);
	  
	} else {
	  if(x < 0.5) parent = dad;
	  else parent = mom;
	  
	  /* This call generates an offspring exclusively    *
	   * from the genetic makeup of individual parent.   *
	   * If the parameter clonal_reproduction is .true., *
	   * the offspring is a clone of parent except for   *
	   * possible new mutations.  If not, the genotype   *
	   * of the offspring is the product of gametic      * 
	   * shuffling of the parent's chromosomes via       *
	   * self-fertilization.                             */
	  offspring(dmutn_offsprng, fmutn_offsprng, 
		    offsprng_lb_mutn_count, 
		    offsprng_lb_fitness, dmutn, fmutn, lb_mutn_count, 
		    linkage_block_fitness, parent, parent, child);
	}
	
      }

      offspring_count += actual_offspring;
      
      /* Copy mutation list arrays for each of the first two     *
       * offspring into locations of the two parents.  Update    *
       * the linkage block mutation count and the linkage block  *
       * fitness for these two offspring.                        */

      if(actual_offspring >= 1) {
	k = dmutn_offsprng[0][0][0] + 1;
	for(j=0; j < k; j++)
	  dmutn[dad][0][j] = dmutn_offsprng[0][0][j];
	k = dmutn_offsprng[0][1][0] + 1;
	for(j=0; j < k; j++)
	  dmutn[dad][1][j] = dmutn_offsprng[0][1][j];
	k = fmutn_offsprng[0][0][0] + 1;
	for(j=0; j < k; j++)
	  fmutn[dad][0][j] = fmutn_offsprng[0][0][j];
	k = fmutn_offsprng[0][1][0] + 1;
	for(j=0; j < k; j++)
	  fmutn[dad][1][j] = fmutn_offsprng[0][1][j];
	
	for (j = 0; j < 2; j++)
	  for (k = 0; k < 2; k++)
	    for (l = 0; l < num_linkage_subunits; l++)
	      lb_mutn_count[dad][j][k][l] = 
		offsprng_lb_mutn_count[0][j][k][l];
	
	for (j = 0; j < 2; j++)
	  for (k = 0; k < num_linkage_subunits; k++)
	    linkage_block_fitness[dad][j][k] = 
	      offsprng_lb_fitness[0][j][k];
	
	replaced_by_offspring[dad] = true;
      }
      
      if(actual_offspring >= 2) {
	k = dmutn_offsprng[1][0][0] + 1;
	for (j = 0; j < k; j++)
	  dmutn[mom][0][j] = dmutn_offsprng[1][0][j];
	k = dmutn_offsprng[1][1][0] + 1;
	for (j = 0; j < k; j++)
	  dmutn[mom][1][j] = dmutn_offsprng[1][1][j];
	k = fmutn_offsprng[1][0][0] + 1;
	for (j = 0; j < k; j++)
	  fmutn[mom][0][j] = fmutn_offsprng[1][0][j];
	k = fmutn_offsprng[1][1][0] + 1;
	for (j = 0; j < k; j++)
	  fmutn[mom][1][j] = fmutn_offsprng[1][1][j];
	for (j = 0; j < 2; j++)
	  for (k = 0; k < 2; k++)
	    for (l = 0; l < num_linkage_subunits; l++)
	      lb_mutn_count[mom][j][k][l] =  
		offsprng_lb_mutn_count[1][j][k][l];
	for (j = 0; j < 2; j++)
	  for (k = 0; k < num_linkage_subunits; k++)
	    linkage_block_fitness[mom][j][k] = 
	      offsprng_lb_fitness[1][j][k];
      }		
      
      /* Copy the mutation list for any other offspring into arrays   *
       * dmutn and fmutn with an index greater than current_pop_size. *
       * Update array linkage_block_fitness appropriately.            */
      
      for (slot=2; slot < actual_offspring; slot++) {
	total_offspring++;
	j = total_offspring;
	k = dmutn_offsprng[slot][0][0] + 1;
	for (l = 0; l < k; l++)
	  dmutn[j][0][l] = dmutn_offsprng[slot][0][l];
	k = dmutn_offsprng[slot][1][0] + 1;
	for (l = 0; l < k; l++)
	  dmutn[j][1][l] = dmutn_offsprng[slot][1][l];
	k = fmutn_offsprng[slot][0][0] + 1;
	for (l = 0; l < k; l++)
	  fmutn[j][0][l] = fmutn_offsprng[slot][0][l];
	k = fmutn_offsprng[slot][1][0] + 1;
	for (l = 0; l < k; l++)
	  fmutn[j][1][l] = fmutn_offsprng[slot][1][l];
	for (l = 0; l < 2; l++)
	  for (l1 = 0; l1 < 2; l1++)
	    for (l2 = 0; l2 < num_linkage_subunits; l2++)
	      lb_mutn_count[j][l][l1][l2] = 
		offsprng_lb_mutn_count[slot][l][l1][l2];
	for (l = 0; l < 2; l++)
	  for (l1 = 0; l1 < num_linkage_subunits; l1++)
	    linkage_block_fitness[j][l][l1] =  
	      offsprng_lb_fitness[slot][l][l1];
      }

    }

    /* For slots not overwritten by new offspring, move data from   *
     * offspring with higher index numbers to populate these slots. */
    
    if(fitness_dependent_fertility) {
      
      replace = total_offspring-1;
      empty   = 0;
      
      while(empty < offspring_count && replace >= offspring_count) {
        while(replaced_by_offspring[empty] && empty < current_pop_size - 1)
	  empty++;
	while(!replaced_by_offspring[replace] && replace > 0) 
	  replace--;
	
	if(empty < offspring_count && replace >= offspring_count) {
	  k = dmutn_offsprng[replace][0][0] + 1;
	  for (i = 0; i < k; i++)
	    dmutn[empty][0][i] = dmutn_offsprng[replace][0][i];
	  k = dmutn_offsprng[replace][1][0] + 1;
	  for (i = 0; i < k; i++)
	    dmutn[empty][1][i] = dmutn_offsprng[replace][1][i];
	  k = fmutn_offsprng[replace][0][0] + 1;
	  for (i = 0; i < k; i++)
	    fmutn[empty][0][i] = fmutn_offsprng[replace][0][i];
	  k = fmutn_offsprng[replace][1][0] + 1;
	  for (i = 0; i < k; i++)
	    fmutn[empty][1][i] = fmutn_offsprng[replace][1][i];
	  for (i = 0; i < 2; i++)
	    for (j = 0; j < 2; j++)
	      for (l1 = 0; l1 < num_linkage_subunits; l1++)
		lb_mutn_count[empty][i][j][l1] = 
		  lb_mutn_count[replace][i][j][l1];
	  for (i = 0; i < 2; i++)
	    for (j = 0; j < num_linkage_subunits; j++)
	      linkage_block_fitness[empty][i][j] = 
		linkage_block_fitness[replace][i][j];
	  replaced_by_offspring[empty] = true;
	}
      }	
    }

    /* Because the Poisson random number generator does not yield
     * the specified mean number of new mutations to sufficient
     * accuracy, to improve accuracy make an adjustment to the value
     * fed to the generator. */
    
    cumulative_offspring += offspring_count;
    
    if(((gen - gen_0)%10 == 0 || gen - gen_0 < 10) &&
       new_mutn_per_offspring >= 1.) { 
      d = (float)new_mutn_count/(float)cumulative_offspring;
      poisson_mean += 0.3*(new_mutn_per_offspring - d);
    }
    
    TIME(&stop);
    time_offspring = DIFFTIME(stop,start);
    sec[4] += time_offspring;
    
    current_pop_size = min(current_pop_size, offspring_count);
    
    /* Impose selection based on fitness to reduce the population *
     * size to a value not to exceed the parameter pop_size.      */
    
    TIME(&start);
    selection(dmutn, fmutn, lb_mutn_count, 
	      linkage_block_fitness, fitness, pheno_fitness, 
	      work_fitness, sorted_score, initial_allele_effects,
	      max_size, total_offspring, gen);
    TIME(&stop);
    sec[5] += DIFFTIME(stop,start);
    
    /* Write diagnostic information to output files. */
    
    current_pop_size = max(1, current_pop_size);
    
    if(current_pop_size < 0) 
      printf("ERROR: current_pop_size less than zero.\n");
    
    TIME(&start);
    
    if(mod(gen, 10) == 0 
       || (!bottleneck_yes && current_pop_size <= pop_size/20) 
       || gen < 4) print_flag = true;
    else           print_flag = false;
    
    if(num_contrasting_alleles > 0) 
      diagnostics_contrasting_alleles(dmutn, fmutn,
    				      offsprng_lb_mutn_count, work_fitness, 
    				      initial_allele_effects, ica_count, 
                                      max_size, false);

    diagnostics_history_plot1(dmutn, fmutn, lb_mutn_count, 
			      tracked_fav_mutn, ica_count, gen, print_flag);
    
    if(gen <= 3 || gen%20 == 0) {
      
      /* Output fitness of each individual in population */
      
      rewind(fp16);
      fprintf(fp16,"# individual  fitness\n");
      for(i=0; i < current_pop_size; i++) 
	fprintf(fp16,"%d %f\n",i, fitness[i]);
      fflush(fp16);
    }
    
    if(fitness_distrib_type > 0 && gen%20 == 0) {
      if(tracking_threshold != 1.0) 
	diagnostics_mutn_bins_plot2(dmutn, fmutn, accum, gen);
      
      diagnostics_mutn_bins_plot4(dmutn, fmutn, 
				  linkage_block_fitness, lb_mutn_count, gen);
      
      diagnostics_selection(sorted_score,pheno_fitness,
			    total_offspring,gen);
    }
    
    TIME(&stop);
    time_selection = DIFFTIME(stop, start);
    sec[6] += time_selection;

    TIME(&start);
    if((gen == 200 || gen%1000 == 0) &&
       gen != num_generations && tracking_threshold != 1.0)
       diagnostics_polymorphisms_plot7(dmutn, fmutn, max_size, gen);
    TIME(&stop);
    sec[7] += DIFFTIME(stop,start);
    
#ifdef WRITE_SAMPLE
    if(gen%100 == 0) {
      
      /* Output user-friendly details of the mutations carried by
       * a few representative individuals in the population. */
      
      /* Reuse the following arrays as buffer space, since *
       * they are not being used right now                 */
      defect   = (float ***)dmutn_offsprng;
      improve  = (float ***)fmutn_offsprng;
      effect   = (float *)   pheno_fitness;

      write_sample(dmutn, fmutn, lb_mutn_count,
		   linkage_block_fitness, fitness,
		   defect, improve, effect, gen);
    }	
#endif
    
    /* If the population size or the mean fitness has collapsed, *
     * print message and shutdown all processors.                */
    
    if(post_sel_fitness < extinction_fitness_threshold) {
      diagnostics_history_plot1(dmutn, fmutn, lb_mutn_count, 
				tracked_fav_mutn, ica_count, gen, true);
      
      fprintf(stdout, "Shutdown due to extinction\n");
      fprintf(fp9, "Shutdown due to extinction\n");
      if(is_parallel) run_status = -1;
      goto shutdown;	
    }
    
    if((is_parallel && current_pop_size < 0.1 * pop_size) 
       || current_pop_size <= 1) {
      fprintf(stdout, "Shutdown due to mutational meltdown\n");
      fprintf(fp9, "Shutdown due to mutational meltdown\n");
      if(is_parallel) run_status = -1;
      goto shutdown;
    }
    
    TIME(&stop_gen);
    time_per_gen = DIFFTIME(stop_gen,start_gen);
    
    fprintf(fp22,"%12d %17.7lf %19.7f %19.7f\n",
	    gen, time_per_gen, time_offspring, time_selection);
    fflush(fp22);
    
#ifdef MPICH
    if(is_parallel) {
      mpi_davg(&time_per_gen,&par_time_per_gen,1);
      mpi_davg(&time_offspring,&par_time_offspring, 3);
      mpi_davg(&time_selection,&par_time_selection,1);
      if (myid==0) {
	fprintf(fp23,"%12d %17.lf %19.7lf %19.7lf\n",gen, par_time_per_gen, 
		par_time_offspring, par_time_selection);
        fflush(fp23);
	//if (myid==0 && (mod(gen,10)==0 || gen<4)) 
	//printf("iteration time:%6d  milliseconds\n",(int)(1000*time_per_gen));
      }
    }
#endif
    
    /* Monitor state file for shutdown flag */
    if (run_status >= 0) {
      sprintf(filename,"%s.%s.st8",case_id,myid_str);
      fp = fopen(filename, "r");
      fscanf(fp,"%d",&run_status);
      fclose(fp);
      
      /* Premature shutdown */
      shutdown_gen = gen;
      if (run_status == 1) {
	printf("STATE: WRITING RESTART FILE & EXITING RUN\n");
	write_dump = true;
	restart_dump_number = 8; //shutdown_gen;
	dump_number = restart_dump_number;
	goto shutdown;
      }
    }
    
#ifdef PPM
    /* Write PPM file.  This is used to create PPM file which 
       can be converted into GIF file which shows an image of the
       fitness of the entire simulation representing pixels as
       individuals, with white representing fitness=1, and
       black representing fitness=0.  Must use ppm2gif to convert
       the PPM file to a viewable GIF file.                       */
    
    for(i=0; i<pop_size; i++) {
      if (fitness(i) > 1) {
	red = 255;
	green = 0;
	blue = 0;
      } else {
	red = (int)fitness[i]*255;
	if (red < 0) red = 0;
	green = red;
	blue = red;
      }
      fprintf(fp15,"%4d %4d %4d") red,green,blue;
    }
    fprintf(fp15,"\n");
#endif
    
#ifdef DYNAMIC_POP_SIZE
    if(pop_growth_model > 0) {
      /* For dynamic population sizes compute new pop_size */
      if(pop_growth_model==1) {
         pop_size = ceil(pop_growth_rate*pop_size);
      } else if (pop_growth_model==2) {
         /* pass carrying_capacity in through num_generations */
         carrying_capacity = num_generations;
         pop_size = ceil(pop_size*(1. + pop_growth_rate*
                     (1. - pop_size/carrying_capacity)));
      } else {
         fprintf(stderr,"ERROR: pop_growth_model = %d not supported\n\n",
            pop_growth_model);
         exit(EXIT_FAILURE); 
      }

      this_size = (int) (0.55*offspring_per_female*pop_size 
		    *(1. - fraction_random_death));
      if(this_size > max_size) {
         printf("\n\n\tOUT OF MEMORY! SHUTTING DOWN! POP SIZE = %d\n",pop_size);
         goto shutdown; 
      }
    }
#endif

    fflush(stdout);
    fflush(fp9);
    
  } /* End generation loop */
  
  
 shutdown: /* Shutdown procedures */
  
  /* Perform diagnositics on initial contrasting alleles and
     polymorphisms and create file for plotting. */
  
  TIME(&start);
  
  if(num_contrasting_alleles > 0) 
    diagnostics_contrasting_alleles(dmutn, fmutn,
  				    offsprng_lb_mutn_count, work_fitness, 
  				    initial_allele_effects, ica_count, 
                                    max_size, true);
  
  if(tracking_threshold != 1.0)
    diagnostics_polymorphisms_plot7(dmutn, fmutn, max_size, gen-1);
  
  TIME(&stop);
  sec[7] += DIFFTIME(stop, start);
  
  /* Write an output dump that contains the current set of parameter
     values, the stored mutation array mutn, and the linkage block
     fitness array linkage_block_fitness. */
  
  dump_number++;
  
  if(write_dump) write_output_dump(dmutn,fmutn,lb_mutn_count,
				   linkage_block_fitness,initial_allele_effects,
				   shutdown_gen,myid_str);
  
#ifdef MPICH
  if (is_parallel){
    if (run_status == -1) mpi_myabort();
    mpi_myfinalize(ierr);
  }
#endif
  
  /* Report timing statistics */
  TIME(&stop_mendel);
  sec[1] = DIFFTIME(stop_mendel, start_mendel);
  profile(stdout);
  profile(fp9);

  /* Close files */
  
  fclose(fp4);
  fclose(fp7);
  fclose(fp8);
  fclose(fp9);
  fclose(fp11);
  fclose(fp13);
#ifdef PPM
  fclose(fp15);
#endif
  fclose(fp16);
  fclose(fp19);
  fclose(fp24);
  fclose(fp25);
  fclose(fp26);
  
  if (is_parallel) {
    fclose(fp14);
    fclose(fp17);
    fclose(fp18);
    fclose(fp21);
    fclose(fp23);
    fclose(fp34);
    fclose(fp35);
  }
  
  /* Free memory */
  free(dmutn);
  free(fmutn);
  free(dmutn_offsprng);
  free(fmutn_offsprng);
  free(lb_mutn_count);
  free(offsprng_lb_mutn_count);
  free(linkage_block_fitness);
  free(offsprng_lb_fitness);
  //for some reason following statement causes core dump
  //free(initial_allele_effects);
  free(pheno_fitness);
  free(fitness);
  free(work_fitness);
  free(sorted_score);
  free(available);
  free(replaced_by_offspring);
  dmutn = NULL;
  fmutn = NULL;
  dmutn_offsprng = NULL;
  fmutn_offsprng = NULL;
  lb_mutn_count = NULL;
  offsprng_lb_mutn_count = NULL;
  linkage_block_fitness =  NULL;
  offsprng_lb_fitness = NULL;
  initial_allele_effects = NULL;
  pheno_fitness = NULL;
  fitness =  NULL;
  work_fitness = NULL;
  sorted_score = NULL;
  available = NULL;
  replaced_by_offspring = NULL;

  /* Print version information for various source files */
  sprintf(filename,"%s.ver",case_id);
  fp3 = fopen(filename, "a");
  fprintf(fp3,"%30s\n",version_mendel);
  fprintf(fp3,"%30s\n",version_offspring);
  fprintf(fp3,"%s\n",version_selection);
  fprintf(fp3,"%s\n",version_diagnostics);
  fprintf(fp3,"%s\n",version_fileio);
  fprintf(fp3,"%s\n",version_init);
  fprintf(fp3,"%s\n",version_qsort);
  fprintf(fp3,"%s\n",version_ranlib);
  fprintf(fp3,"%s\n",version_mem);
#ifdef MPICH
  fprintf(fp3,"%s\n",version_mpi);
#endif
  fclose(fp3);
  
  exit(EXIT_SUCCESS);
}
