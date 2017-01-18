#include "mendel.h"

/* This routine reads a dump file containing the mutation arrays *
 * dmutn and fmutn and the linkage block mutation count array    *
 * lb_mutn_count and the linkage block fitness array             *
 * linkage_block_fitness for purposes of restart.  The argument  *
 * generation_number is the generation number of the dump file   *
 * being read.                                                   */
void read_restart_dump(int *** dmutn, int *** fmutn, int **** lb_mutn_count, 
                       double *** linkage_block_fitness, 
                       float * initial_allele_effects, 
                       int *generation_number, int max_size, char myid_str[]) 
{
  FILE *fp;
  int i, j, k, l, n1, nm;
  float f1;
  bool l1;
  TIME_T start, stop;
  char char12[12], char20[20], path[80];
  char filename[40];
  
  TIME(&start);
  
  sprintf(filename,"%s.%s.dmp.%d",case_id,myid_str,restart_dump_number);

  fp = fopen(filename, "r");
  if (fp == NULL) {
    printf("Cannot read %s\n",filename);
    return;
  }
  
  fscanf(fp,"%12d %*s", &n1); //pop_size
  fscanf(fp,"%12d %*s", &n1); //num_generations
  fscanf(fp,"%12d %*s", &n1); //fitness_distrib_type
  fscanf(fp,"%12d %*s", &n1); //selection_scheme
  fscanf(fp,"%12d %*s", &n1); //haploid_chromosome_number
  fscanf(fp,"%12d %*s", &n1); //num_linkage_subunits
#ifdef DYNAMIC_POP_SIZE
  fscanf(fp,"%f   %*s", &f1); //pop_growth_rate
  fscanf(fp,"%12d %*s", &n1); //pop_growth_model
#endif
  fscanf(fp,"%e %*s",   &f1); //haploid_genome_size

  //parameters offspring_per_female to partial_truncation_value
  for(i=0; i<=16; i++)  
    fscanf(fp,"%f %*s",   &f1);

  fscanf(fp,"%d %*s",     &num_contrasting_alleles);
  
  fscanf(fp,"%f %*s",     &f1); //initial_alleles_mean_effect
  fscanf(fp,"%f %*s",     &f1); //linked_mutn_se_fraction
  fscanf(fp,"%f %*s",     &f1); //se_scaling_factor
  fscanf(fp,"%1d %*s",    &l1); //synergistic_epistasis
  fscanf(fp,"%1d %*s",    &l1); //clonal_reporduction
  fscanf(fp,"%1d %*s",    &l1); //clonal_haploid
  fscanf(fp,"%1d %*s",    &l1); //dynamic_linkage
  fscanf(fp,"%1d %*s",    &l1); //fitness_dependent_fertility
  fscanf(fp,"%1d %*s",    &l1); //is_parallel
  fscanf(fp,"%1d %*s",    &l1); //bottleneck_yes
  fscanf(fp,"%12d %*s",   &n1); //bottleneck_generation
  fscanf(fp,"%12d %*s",   &n1); //bottleneck_pop_size
  fscanf(fp,"%12d %*s",   &n1); //num_bottleneck_generations
  fscanf(fp,"%12d %*s",   &n1); //num_initial_fav_mutn
  fscanf(fp,"%12d %*s",   &n1); //num_indiv_exchanged
  fscanf(fp,"%12d %*s",   &n1); //migration_generations
  fscanf(fp,"%12d %*s",   &n1); //migration_model
  fscanf(fp,"%1d %*s",    &l1); //homogenous_tribes
  fscanf(fp,"%12d %*s",   &n1); //max_tracked_mutn_per_indiv
  fscanf(fp,"%12d %*s",   &n1); //random_number_seed
  fscanf(fp,"%1d %*s",    &l1); //write_dump
  fscanf(fp,"%1d %*s",    &l1); //restart_case
  fscanf(fp,"%12d %*s",   &n1); //restart_dump_number
  fscanf(fp,"%12s %*s",   char12); //case_id
  fscanf(fp,"%20s %80s",  char20, path );

  fscanf(fp,"%12d %*s", generation_number);

  /* Initialize arrays with max values */
  for (i = 0; i < max_size; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < max_del_mutn_per_indiv/2; k++)
        dmutn[i][j][k] = num_linkage_subunits*lb_modulo + 1;

  for (i = 0; i < max_size; i++)
    for (j = 0; j < 2; j++)
      for (k = 0; k < max_fav_mutn_per_indiv/2; k++)
        fmutn[i][j][k] = num_linkage_subunits*lb_modulo + 1;
  
  /* Read in restart data */
  for(i = 0; i < pop_size; i++) {

    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        for (l = 0; l < num_linkage_subunits; l++ )
           fscanf(fp,"%6d",&lb_mutn_count[i][j][k][l]);

    for (j = 0; j < 2; j++)
       for (l = 0; l < num_linkage_subunits; l++ ) 
          fscanf(fp, "%lf",&linkage_block_fitness[i][j][l]);
    
    fscanf(fp,"%12d",&dmutn[i][0][0]);
    nm = min(max_del_mutn_per_indiv/2, dmutn[i][0][0]+1);

    for (j = 1; j <= dmutn[i][0][0]; j++) 
       fscanf(fp,"%12d",&dmutn[i][0][j]);

    fscanf(fp,"%12d", &dmutn[i][1][0]);
    nm = min(max_del_mutn_per_indiv/2, dmutn[i][1][0]+1);

    for (j = 1; j <= dmutn[i][1][0]; j++) 
       fscanf(fp,"%12d",&dmutn[i][1][j]);

    fscanf(fp," %12d",&fmutn[i][0][0]);
    nm = min(max_fav_mutn_per_indiv/2,fmutn[i][0][0]+1);

    for (j = 1; j <= fmutn[i][0][0]; j++) 
       fscanf(fp,"%12d",&fmutn[i][0][j]);

    fscanf(fp," %12d",&fmutn[i][1][0]);
    nm = min(max_fav_mutn_per_indiv/2, fmutn[i][1][0]+1);

    if (fmutn[i][0][0] > 0) 
    for (j = 1; j <= fmutn[i][1][0]; j++) 
       fscanf(fp,"%12d",&fmutn[i][1][j]);
  }
  
  
  if(num_contrasting_alleles > 0) {
    for (l = 0; l < num_linkage_subunits; l++ )
       fscanf(fp,"%f",&initial_allele_effects[l]);

    if(!is_parallel)  {
      printf("Restart run will use the previous value for");
      printf("parameter num_contrasting_alleles = %10d", 
	     num_contrasting_alleles);
    }

    fprintf(fp9,"Restart run will use the previous value for");
    fprintf(fp9,"parameter num_contrasting_alleles = %10d",
	    num_contrasting_alleles);
  }

  fclose(fp);
  
  TIME(&stop);
  sec[8] += DIFFTIME(stop,start);
}

/* This routine printfs an output file containing the current    *
 * parameter values, the stored mutation arrays dmutn and fmutn, *
 * the linkage block mutation count array lb_mutn_count and      *
 * the linkage block fitness array linkage_block_fitness for     *
 * the current generation specified by generation_number.        */
void write_output_dump(int *** dmutn, int *** fmutn, 
                       int **** lb_mutn_count, 
                       double *** linkage_block_fitness,
		       float *initial_allele_effects,
		       int generation_number,char myid_str[])
{
  int i, j, k, l;
  char filename[40];
  TIME_T start, stop;
  FILE *fp;
  
  TIME(&start);
  
  sprintf(filename,"%s.%s.dmp.%d",case_id,myid_str,restart_dump_number);
  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Cannot read %s\n",filename);
    return;
  }
  
  write_parameters(fp);
  
  fprintf(fp,"%12d\tgeneration_number\n",generation_number);
  
  for (i = 0; i < pop_size; i++){
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        for (l = 0; l < num_linkage_subunits; l++) 
           fprintf(fp,"%6d",lb_mutn_count[i][j][k][l]);
    fprintf(fp,"\n");

    for (j = 0; j < 2; j++)
      for (l = 0; l < num_linkage_subunits; l++) 
           fprintf(fp,"%12.8f",linkage_block_fitness[i][j][l]);

    fprintf(fp,"\n%12d\n",dmutn[i][0][0]);
    
    for (j = 1; j <= dmutn[i][0][0]; j++) 
       fprintf(fp,"%12d",dmutn[i][0][j]);
    
    fprintf(fp,"\n%12d\n",dmutn[i][1][0]);
    
    for (j = 1; j <= dmutn[i][1][0]; j++) 
       fprintf(fp,"%12d",dmutn[i][1][j]);

    fprintf(fp,"\n%12d\n",fmutn[i][0][0]);

    for (j = 1; j <= fmutn[i][0][0]; j++) 
       fprintf(fp,"%12d",fmutn[i][0][j]);

    fprintf(fp,"\n%12d\n",fmutn[i][1][0]);

    for (j = 1; j <= fmutn[i][1][0]; j++) 
       fprintf(fp,"%12d",fmutn[i][1][j]);
    
    fprintf(fp,"\n");
  }
  
  if(num_contrasting_alleles > 0) 
    for(i=0; i<num_linkage_subunits; i++) 
      fprintf(fp,"%12.9f", initial_allele_effects[i]);
    fprintf(fp,"\n");
  
  fclose(fp);
  
  TIME(&stop);
  sec[9] += DIFFTIME(stop,start);
}


#ifdef WRITE_SAMPLE
/* This routine prints an output file containing details concerning *
 * the mutations carried by five members of the total population in *
 * a format that is (hopefully) readily understandable.             */
void write_sample(int ***dmutn, int ***fmutn, int ****lb_mutn_count, 
                  double ***linkage_block_fitness, double *fitness, 
                  float ***defect, float ***improve, float *effect, 
                  int generation_number)
{
  int i, j, lb, m, mutn;
  TIME_T start, stop;
  float d;
  int del_mutn[2][num_linkage_subunits];
  int fav_mutn[2][num_linkage_subunits];
  FILE *myfp;
  char filename[20];

  TIME(&start);
  
  sprintf(filename,"%s.%.3d.sam",case_id,myid);
  myfp = fopen(filename, "w");

  write_parameters(myfp);
  
  fprintf(myfp,"                       generation number = %6d\n",
	  generation_number);
  
  for (i=current_pop_size/10; i < current_pop_size/3; i+=current_pop_size/5) { 
    
    fprintf(myfp,"\nindividual number = %6d  fitness = %86.f\n", i, fitness[i]);
    
    for (lb=0; lb<num_linkage_subunits; lb++) {
       del_mutn[0][lb] = 0;
       del_mutn[1][lb] = 0;
       fav_mutn[0][lb] = 0;
       fav_mutn[1][lb] = 0;
    }
    
    for (m=1; m <= dmutn[i][0][0]; m++) {
      lb = abs(dmutn[i][0][m])/lb_modulo;
      del_mutn[0][lb]++;
      mutn = abs(dmutn[i][0][m])%lb_modulo;
      d    = exp(-alpha_del*pow((float)mutn*del_scale,gamma_del));
      if(dmutn[i][0][m] < 0) d = -d;
      defect[del_mutn[0][lb]][0][lb] = d;
    }
    
    for (m=1; m <= dmutn[i][1][0]; m++) {
      lb = abs(dmutn[i][1][m])/lb_modulo;
      del_mutn[1][lb]++;
      mutn = abs(dmutn[i][1][m])%lb_modulo;
      d    = exp(-alpha_del*pow((float)mutn*del_scale,gamma_del));
      if(dmutn[i][1][m] < 0) d = -d;
      defect[del_mutn[1][lb]][1][lb] = d;
    }
    
    for (m=1; m <= fmutn[i][0][0]; m++) {
      lb = abs(fmutn[i][0][m])/lb_modulo;
      fav_mutn[0][lb]++;
      mutn = abs(fmutn[i][0][m])%lb_modulo;
      d    = exp(-alpha_fav*pow((float)mutn*fav_scale,gamma_fav)
		 *max_fav_fitness_gain);
      if(fmutn[i][0][m] < 0) d = -d;
      improve[fav_mutn[0][lb]][0][lb] = d;
    }
    
    for (m=1; m <= fmutn[i][1][0]; m++) {
      lb = abs(fmutn[i][1][m])/lb_modulo;
      fav_mutn[1][lb]++;
      mutn = abs(fmutn[i][1][m])%lb_modulo;
      d    = exp(-alpha_fav*pow((float)mutn*fav_scale,gamma_fav)
		 *max_fav_fitness_gain);
      if(fmutn[i][1][m] < 0) d = -d;
      improve[fav_mutn[1][lb]][1][lb] = d;
    }
    
    for(lb=0; lb < num_linkage_subunits; lb++) {
      
      fprintf(myfp,"\n                           lb number = %6d\n",lb);
      
      fprintf(myfp,"Haplotype 1: total deleterious mutn count = %6d  composite fitness = %8.6f\n",
	      lb_mutn_count[i][0][0][lb], linkage_block_fitness[i][0][lb]);
      
      if(tracking_threshold != 1.0) {
	
	j = 0;
	
	for(m=0; m<del_mutn[0][lb]; m++) {
	  if(defect[m][0][lb] < 0) {
	    effect[j] = -defect[m][0][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness degradations of tracked deleterious recessive mutations: ");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<del_mutn[0][lb]; m++) {
	  if(defect[m][0][lb] > 0) {
	    effect[j] = defect[m][0][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness degradations of tracked deleterious dominant mutations: ");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<fav_mutn[0][lb]; m++) {
	  if(improve[m][0][lb] < 0) {
	    effect[j] = -improve[m][0][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness improvements of tracked favorable recessive mutations:\n");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<fav_mutn[0][lb]; m++) {
	  if(improve[m][0][lb] > 0) {
	    effect[j] = improve[m][0][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness improvements of tracked favorable dominant mutations:\n");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
      } 
      
      fprintf(myfp,"Haplotype 2: total deleterious mutn count = %6d  composite fitness = %8.6f\n", 
	      lb_mutn_count[i][0][1][lb], linkage_block_fitness[i][1][lb]);
      
      if(tracking_threshold != 1.0) {
	
	j = 0;
	for(m=0; m<del_mutn[1][lb]; m++) {
	  if(defect[m][1][lb] < 0) {
	    effect[j] = -defect[m][1][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness degradations of tracked deleterious recessive mutations: ");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<del_mutn[1][lb]; m++) {
	  if(defect[m][1][lb] > 0) {
	    effect[j] = defect[m][1][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness degradations of tracked deleterious dominant mutations: ");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<fav_mutn[1][lb]; m++) {
	  if(improve[m][1][lb] < 0) {
	    effect[j] = -improve[m][1][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness improvements of tracked favorable recessive mutations:");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
	j = 0;
	for(m=0; m<fav_mutn[1][lb]; m++) {
	  if(improve[m][1][lb] > 0) {
	    effect[j] = improve[m][1][lb];
	    j++;
	  }
	}
	if(j > 0) {
	  fprintf(myfp,"\nFitness improvements of tracked favorable dominant mutations:");
	  for(m=0; m<j; m++) fprintf(myfp,"%9.6f ",effect[m]);
	  fprintf(myfp,"\n");
	} 
	
      } 
      
    }
    
  }
  
  fclose(myfp);

  TIME(&stop);
  sec[10] += DIFFTIME(stop,start);

}
#endif


void read_parameters(FILE *fp) {
  TIME_T start, stop;

  sprintf(version_fileio,
          "$Id: fileio.c,v 1.38 2008/12/01 05:51:42 wes Exp $");
  
  TIME(&start);
  
  fscanf(fp,"%12d %*s", &pop_size);
  fscanf(fp,"%12d %*s", &num_generations);
  fscanf(fp,"%12d %*s", &fitness_distrib_type);
  fscanf(fp,"%12d %*s", &selection_scheme);
  fscanf(fp,"%12d %*s", &haploid_chromosome_number);
  fscanf(fp,"%12d %*s", &num_linkage_subunits);
#ifdef DYNAMIC_POP_SIZE
    fscanf(fp,"%12f %*s", &pop_growth_rate);
    fscanf(fp,"%12d %*s", &pop_growth_model);
#endif
  fscanf(fp,"%12e %*s", &haploid_genome_size);
  fscanf(fp,"%12f %*s", &offspring_per_female);
  fscanf(fp,"%12f %*s", &fraction_random_death);
  fscanf(fp,"%12f %*s", &fraction_self_fertilization);
  fscanf(fp,"%12f %*s", &new_mutn_per_offspring);
  fscanf(fp,"%12f %*s", &high_impact_mutn_fraction); 
  fscanf(fp,"%12f %*s", &high_impact_mutn_threshold);
  fscanf(fp,"%12f %*s", &uniform_fitness_effect);
  fscanf(fp,"%12f %*s", &multiplicative_weighting);
  fscanf(fp,"%12f %*s", &tracking_threshold);
  fscanf(fp,"%12f %*s", &fraction_recessive);
  fscanf(fp,"%12f %*s", &recessive_hetero_expression);
  fscanf(fp,"%12f %*s", &dominant_hetero_expression);
  fscanf(fp,"%12f %*s", &frac_fav_mutn);
  fscanf(fp,"%12f %*s", &max_fav_fitness_gain);
  fscanf(fp,"%12f %*s", &heritability);
  fscanf(fp,"%12f %*s", &non_scaling_noise);
  fscanf(fp,"%12f %*s", &partial_truncation_value);
  fscanf(fp,"%12d %*s", &num_contrasting_alleles);
  fscanf(fp,"%12f %*s", &initial_alleles_mean_effect);
  fscanf(fp,"%12f %*s", &linked_mutn_se_fraction);
  fscanf(fp,"%12f %*s", &se_scaling_factor);
  fscanf(fp,"%12d %*s", &synergistic_epistasis);
  fscanf(fp,"%12d %*s", &clonal_reproduction);
  fscanf(fp,"%12d %*s", &clonal_haploid);
  fscanf(fp,"%12d %*s", &dynamic_linkage);
  fscanf(fp,"%12d %*s", &fitness_dependent_fertility);
  fscanf(fp,"%12d %*s", &is_parallel);
  fscanf(fp,"%12d %*s", &bottleneck_yes);
  fscanf(fp,"%12d %*s", &bottleneck_generation);
  fscanf(fp,"%12d %*s", &bottleneck_pop_size);
  fscanf(fp,"%12d %*s", &num_bottleneck_generations);
  fscanf(fp,"%12d %*s", &num_initial_fav_mutn);
  fscanf(fp,"%12d %*s", &num_indiv_exchanged);
  fscanf(fp,"%12d %*s", &migration_generations);
  fscanf(fp,"%12d %*s", &migration_model);
  fscanf(fp,"%12d %*s", &homogenous_tribes);
  fscanf(fp,"%12d %*s", &max_tracked_mutn_per_indiv);
  fscanf(fp,"%12d %*s", &random_number_seed);
  fscanf(fp,"%12d %*s", &write_dump);
  fscanf(fp,"%12d %*s", &restart_case);
  fscanf(fp,"%12d %*s", &restart_dump_number);
  fscanf(fp,"%s %*s", case_id);
  fscanf(fp,"%s %*s", data_file_path);
  
  TIME(&stop);
  sec[11] += DIFFTIME(stop,start);
}

/* This routine fprintfs the current parameter values to unit fp. */
void write_parameters(FILE *fp) {
  TIME_T start, stop;
  
  TIME(&start);
  
  fprintf(fp,"%12d\tpop_size\n", pop_size);
  fprintf(fp,"%12d\tnum_generations\n", num_generations);
  fprintf(fp,"%12d\tfitness_distrib_type\n", fitness_distrib_type); 
  fprintf(fp,"%12d\tselection_scheme\n", selection_scheme);
  fprintf(fp,"%12d\thaploid_chromosome_number\n", haploid_chromosome_number);
  fprintf(fp,"%12d\tnum_linkage_subunits\n", num_linkage_subunits); 
#ifdef DYNAMIC_POP_SIZE
  fprintf(fp,"%12.7f\tpop_growth_rate\n", pop_growth_rate);
  fprintf(fp,"%12d\tpop_growth_model\n",  pop_growth_model);
#endif
  fprintf(fp,"%12.3e\thaploid_genome_size\n", haploid_genome_size);
  fprintf(fp,"%12.7f\toffspring_per_female\n", offspring_per_female); 
  fprintf(fp,"%12.7f\tfraction_random_death\n", fraction_random_death);
  fprintf(fp,"%12.7f\tfraction_self_fertilization\n", fraction_self_fertilization); 
  fprintf(fp,"%12.7f\tnew_mutn_per_offspring\n", new_mutn_per_offspring); 
  fprintf(fp,"%12.7f\thigh_impact_mutn_fraction\n", high_impact_mutn_fraction);
  fprintf(fp,"%12.7f\thigh_impact_mutn_threshold\n", high_impact_mutn_threshold);
  fprintf(fp,"%12.7f\tuniform_fitness_effect\n", uniform_fitness_effect);
  fprintf(fp,"%12.7f\tmultiplicative_weighting\n", multiplicative_weighting);
  fprintf(fp,"%12.3e\ttracking_threshold\n", tracking_threshold); 
  fprintf(fp,"%12.7f\tfraction_recessive\n", fraction_recessive); 
  fprintf(fp,"%12.7f\trecessive_hetero_expression\n", recessive_hetero_expression);
  fprintf(fp,"%12.7f\tdominant_hetero_expression\n", dominant_hetero_expression);
  fprintf(fp,"%12.7f\tfrac_fav_mutn\n", frac_fav_mutn);
  fprintf(fp,"%12.7f\tmax_fav_fitness_gain\n", max_fav_fitness_gain);
  fprintf(fp,"%12.7f\theritability\n", heritability); 
  fprintf(fp,"%12.7f\tnon_scaling_noise\n", non_scaling_noise);
  fprintf(fp,"%12.7f\tpartial_truncation_value\n", partial_truncation_value);
  fprintf(fp,"%12d\tnum_contrasting_alleles\n", num_contrasting_alleles);
  fprintf(fp,"%12.7f\tinitial_alleles_mean_effect\n", initial_alleles_mean_effect);
  fprintf(fp,"%12.7f\tlinked_mutn_se_fraction\n", linked_mutn_se_fraction);
  fprintf(fp,"%12.7f\tse_scaling_factor\n", se_scaling_factor);
  fprintf(fp,"%12d\tsynergistic_epistasis\n", synergistic_epistasis);
  fprintf(fp,"%12d\tclonal_reproduction\n", clonal_reproduction);
  fprintf(fp,"%12d\tclonal_haploid\n", clonal_haploid); 
  fprintf(fp,"%12d\tdynamic_linkage\n", dynamic_linkage); 
  fprintf(fp,"%12d\tfitness_dependent_fertility\n", fitness_dependent_fertility);
  fprintf(fp,"%12d\tis_parallel\n", is_parallel); 
  fprintf(fp,"%12d\tbottleneck_yes\n", bottleneck_yes); 
  fprintf(fp,"%12d\tbottleneck_generation\n", bottleneck_generation); 
  fprintf(fp,"%12d\tbottleneck_pop_size\n", bottleneck_pop_size); 
  fprintf(fp,"%12d\tnum_bottleneck_generations\n", num_bottleneck_generations);
  fprintf(fp,"%12d\tnum_initial_fav_mutn\n", num_initial_fav_mutn);
  fprintf(fp,"%12d\tnum_indiv_exchanged\n", num_indiv_exchanged); 
  fprintf(fp,"%12d\tmigration_generations\n", migration_generations);
  fprintf(fp,"%12d\tmigration_model\n", migration_model);
  fprintf(fp,"%12d\thomogenous_tribes\n", homogenous_tribes);
  fprintf(fp,"%12d\tmax_tracked_mutn_per_indiv\n", max_tracked_mutn_per_indiv); 
  fprintf(fp,"%12d\trandom_number_seed\n", random_number_seed); 
  fprintf(fp,"%12d\twrite_dump\n", write_dump);
  fprintf(fp,"%12d\trestart_case\n", restart_case);
  fprintf(fp,"%12d\trestart_dump_number\n", restart_dump_number); 
  fprintf(fp,"      %.6s\tcase_id\n", case_id);
  fprintf(fp,"data_file_path: %s\n\n", data_file_path);
  
  TIME(&stop);
  sec[12] += DIFFTIME(stop,start);
}


void write_status(FILE *unit, int gen, int current_pop_size,
		  float frac_recessive, double total_del_mutn, 
		  double tracked_del_mutn, double total_fav_mutn, 
		  double pre_sel_fitness, double pre_sel_geno_sd, 
		  double pre_sel_pheno_sd, double pre_sel_corr, 
		  double post_sel_fitness, double post_sel_geno_sd, 
		  double post_sel_pheno_sd, double post_sel_corr)
{
  fprintf(unit,"\ngeneration     = %6d  population size = %6d", 
	  gen, current_pop_size);
  fprintf(unit,"  frac recessive =%7.4f\n", frac_recessive); 
  fprintf(unit,"before sel: geno fitness = %9.5f",pre_sel_fitness);
  fprintf(unit,"  geno s.d. =%8.5f  pheno s.d. =%8.5f\n",
	  pre_sel_geno_sd, pre_sel_pheno_sd);
  fprintf(unit,"after  sel:                 %8.5f             %8.5f              %8.5f\n",
	  post_sel_fitness, post_sel_geno_sd,post_sel_pheno_sd);
  fprintf(unit,"before sel geno-pheno corr =%7.4f", pre_sel_corr);
  fprintf(unit,"     after sel geno-pheno corr =%7.4f\n",post_sel_corr);
  fprintf(unit," del mutn/indiv =%6d", (int)(total_del_mutn/current_pop_size));
  fprintf(unit,"  tracked del/ind = %6d  fav mutn/indiv =%7.4f\n",
	  (int)(tracked_del_mutn/current_pop_size),
	  total_fav_mutn/current_pop_size);
  fflush(unit);
}
