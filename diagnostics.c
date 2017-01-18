#include "mendel.h"
#define NBINS 100

void diagnostics_history_plot1(int *** dmutn, int *** fmutn, 
                               int ****lb_mutn_count, int tracked_fav_mutn,
                               int ica_count[], int gen, bool print_flag)
{
  int i, j, lb, num_recessive;
  double total_del_mutn, total_fav_mutn, tracked_del_mutn;
  double frac_recessive;
  double par_total_del_mutn, par_total_fav_mutn;
  double par_tracked_del_mutn, frac_accum, total_mutn, st;
  double par_pre_sel_fitness, par_post_sel_fitness,
    par_pre_sel_geno_sd, par_pre_sel_pheno_sd, 
    par_pre_sel_corr, par_post_sel_geno_sd, 
    par_post_sel_pheno_sd, par_post_sel_corr;
  
  total_del_mutn   = 0;
  total_fav_mutn   = 0;
  tracked_del_mutn = 0;
  tracked_fav_mutn = 0;
  
  for(i=0; i < current_pop_size; i++) {
    for(lb=0; lb < num_linkage_subunits; lb++) {
      total_del_mutn += lb_mutn_count[i][0][0][lb] + lb_mutn_count[i][0][1][lb];
      total_fav_mutn += lb_mutn_count[i][1][0][lb] + lb_mutn_count[i][1][1][lb];
    }
    tracked_del_mutn += dmutn[i][0][0] + dmutn[i][1][0];
    tracked_fav_mutn += fmutn[i][0][0] + fmutn[i][1][0];
  }
  
  tracked_del_mutn -= ica_count[0];
  tracked_fav_mutn -= ica_count[1];
  
  /* Compute averages across processors. */
  
#ifdef MPICH
  if (is_parallel) {
    mpi_davg(&pre_sel_fitness,&par_pre_sel_fitness,1);
    mpi_davg(&post_sel_fitness,&par_post_sel_fitness,1);
    mpi_mybcast(par_post_sel_fitness,1);
    mpi_davg(&pre_sel_geno_sd,&par_pre_sel_geno_sd,1);
    mpi_davg(&pre_sel_pheno_sd,&par_pre_sel_pheno_sd,1);
    mpi_davg(&pre_sel_corr,&par_pre_sel_corr,1);
    mpi_davg(&post_sel_geno_sd,&par_post_sel_geno_sd,1);
    mpi_davg(&post_sel_pheno_sd,&par_post_sel_pheno_sd,1);
    mpi_davg(&post_sel_corr,&par_post_sel_corr,1);
    mpi_davg(&total_del_mutn,&par_total_del_mutn,1);
    mpi_davg(&tracked_del_mutn,&par_tracked_del_mutn,1);
    mpi_davg(&total_fav_mutn,&par_total_fav_mutn,1);
  }
#endif
  
  /*  Output to a file for plotting the generation number, the mean *
   *  fitness, the average number of mutations per individual, and  *
   *  the number of fixed favorable mutations.                      */
  
  fprintf(fp7,"%12d ",gen);
  fprintf(fp7,"%12.4e ",post_sel_fitness);
  fprintf(fp7,"%12.4e ",post_sel_geno_sd);
  fprintf(fp7,"%14.4e ",total_del_mutn/current_pop_size);
  fprintf(fp7,"%14.4e ",total_fav_mutn/current_pop_size);
  fprintf(fp7,"%12d\n",current_pop_size);
  fflush(fp7);
  
  if(is_parallel && myid==0) {
    fprintf(fp17,"%d %lf %16.4e %16.4e %14.4e %12d\n", 
	    gen, par_post_sel_fitness,
	    par_post_sel_geno_sd, 
	    par_total_del_mutn/current_pop_size,
	    par_total_fav_mutn/current_pop_size,
            current_pop_size);
    fflush(fp17);
  }
  
  if (!print_flag) return;
  
  num_recessive  = 0;
  frac_recessive = 0.;
  
  for(i=0; i < current_pop_size; i++) {
    for(j=1; j <= dmutn[i][0][0]; j++)
      if(dmutn[i][0][j] < 0) num_recessive++;
    for(j=1; j <= dmutn[i][1][0]; j++)
      if(dmutn[i][1][j] < 0) num_recessive++;
  }
  
  if(tracked_del_mutn > 0) 
    frac_recessive = (float)num_recessive/tracked_del_mutn;
  
  write_status(fp9, gen, current_pop_size, 
	       frac_recessive, total_del_mutn, tracked_del_mutn, 
	       total_fav_mutn, pre_sel_fitness, pre_sel_geno_sd, 
	       pre_sel_pheno_sd, pre_sel_corr, post_sel_fitness,
	       post_sel_geno_sd, post_sel_pheno_sd, post_sel_corr);
  
  if(is_parallel) {
    
    if (myid == 0) { 
      write_status(stdout, gen, current_pop_size*num_procs, 
		   frac_recessive, par_total_del_mutn*num_procs, 
		   par_tracked_del_mutn*num_procs,
		   par_total_fav_mutn*num_procs, 
		   par_pre_sel_fitness, par_pre_sel_geno_sd, 
		   par_pre_sel_pheno_sd, par_pre_sel_corr,
		   par_post_sel_fitness, par_post_sel_geno_sd, 
		   par_post_sel_pheno_sd, par_post_sel_corr);
    }
    
  } else {
    
    write_status(stdout, gen, current_pop_size, 
		 frac_recessive, total_del_mutn, tracked_del_mutn, 
		 total_fav_mutn, pre_sel_fitness, pre_sel_geno_sd, 
		 pre_sel_pheno_sd, pre_sel_corr, post_sel_fitness,
		 post_sel_geno_sd, post_sel_pheno_sd, post_sel_corr);
    
  }
  
  if(tracking_threshold == 1.0 && gen >= 200 && gen%20==0) {
    total_mutn = gen*current_pop_size*new_mutn_per_offspring;
    if(frac_fav_mutn != 1.) {
      frac_accum = total_del_mutn/(total_mutn*(1.-frac_fav_mutn));
      st = 1.3*exp(-alpha_del*pow(1. - frac_accum,gamma_del));
      fprintf(fp9, "deleterious selection threshold  =%10.3e\n",st);
      fprintf(fp9, "deleterious fraction accumulated =%7.4f\n",frac_accum);
      if(!is_parallel) {
	printf("deleterious selection threshold  =%10.3e\n",st);
	printf("deleterious fraction accumulated =%7.4f\n",frac_accum);
      }
      fprintf(fp25, "%10d%15.3e NaN NaN\n", gen, st);
      if(is_parallel && myid==0) 
	fprintf(fp35, "%10d%15.3e NaN NaN\n", gen, st);
    }
  }
  
  if(num_contrasting_alleles > 0) {
    fprintf(fp9,"\n                  Statistics for initial contrasting alleles\n");
    fprintf(fp9,"mean fitness contrib = %7.4f  ", ica_mean_effect);
    fprintf(fp9,"fav mean freq = %7.4f  ", fav_mean_freq);
    fprintf(fp9,"fixed = %4d  ", fav_fixed);
    fprintf(fp9,"lost = %4d\n", fav_lost);
  }
  
  if(num_contrasting_alleles > 0 && !is_parallel) {
    printf("\n                  Statistics for initial contrasting alleles\n");
    printf("mean fitness contrib = %7.4f  ", ica_mean_effect);
    printf("fav mean freq = %7.4f  ", fav_mean_freq);
    printf("fixed = %4d  ", fav_fixed);
    printf("lost = %4d\n", fav_lost);
  }
  
}

void diagnostics_mutn_bins_plot2(int ***dmutn, int ***fmutn, double *accum, 
                                 int gen) {
  
  int i, j, k, k0, accum_gen;
  int oneortwo;
  FILE *fid;
  double fitness_bins[2][100] = {{0}};
  double par_fitness_bins[2][100] = {{0}};
  double bin_fitness_boxwidth[101], bin_fitness_midpoint[101], work[2][100];
  double refr_bins[100], bin_fitness[101], del_bin_width, b0, b1;
  double d, x0, x1, y0, y1, s, mutn_sum, fav_bin_width;
#ifdef MPICH
  double par_del_dom_thres, par_fav_dom_thres;
#endif
  double av1, av2, fm1, fm2, sum, del, del_no_sel, ratio;
  double del_dom_thres, del_rec_thres, fav_dom_thres;
  
  // Compute the total number of mutations and the bin widths.
  
  mutn_sum  = current_pop_size*gen*new_mutn_per_offspring;
  del_bin_width = -log(tracking_threshold)/50;
  if(max_fav_fitness_gain > 0.) 
    fav_bin_width = -log(tracking_threshold/max_fav_fitness_gain)/50.;
  else
    fav_bin_width = del_bin_width;
  
  x0 = 0.;
  y0 = 0.;
  
  for(k=0; k<50; k++) {
    x1 = pow(del_bin_width*(k+1)/alpha_del,1./gamma_del);
    refr_bins[k] = (1. - frac_fav_mutn)*mutn_sum*(x1 - x0);
    y1 = pow(fav_bin_width*(k+1)/alpha_fav,1./gamma_fav);
    refr_bins[50+k] = frac_fav_mutn*mutn_sum*(y1 - y0);
    x0 = x1;
    y0 = y1;
  }
  
  //  Compute statistics on favorable and recessive mutations.
  
  for(i=0; i<current_pop_size; i++) {
    for(j=1; j <= dmutn[i][0][0]; j++) {
      d = alpha_del*pow(mod(abs(dmutn[i][0][j]),lb_modulo)*del_scale,gamma_del);
      k = (int)(d/del_bin_width);
      //printf("d: %lf,   del_bin_width: %lf,   k: %d\n",d,del_bin_width,k);
      if(dmutn[i][0][j] < 0) {
      	if(k < 50) fitness_bins[0][k]++;
      } else {
	      if(k < 50) fitness_bins[1][k]++;
      }
    }
    
    for(j=1; j <= dmutn[i][1][0]; j++) {
      d = alpha_del*pow(mod(abs(dmutn[i][1][j]),lb_modulo)*del_scale,gamma_del);
      k = (int)(d/del_bin_width);
      //printf("d is %lf k is %d del_bin_width is %lf\n",d,k,del_bin_width);
      if(dmutn[i][0][j] < 0) {
	      if(k < 50) fitness_bins[0][k]++;
      } else {
	      if(k < 50) fitness_bins[1][k]++;
      }
    }
    
    //printf("fmutn[%d][0][0] is: %d\n",i,fmutn[i][0][0]);
    for(j=1; j <= fmutn[i][0][0]; j++) {
      d = alpha_fav*pow(mod(abs(fmutn[i][0][j]),lb_modulo)*fav_scale,gamma_fav);
      k = 50 + (int)(d/fav_bin_width);
      //printf("k1, d is: %d %lf\n",k,d);
      if(fmutn[i][0][j] < 0) {
	      if(k < 100) fitness_bins[0][k]++;
      } else {
	      if(k < 100) fitness_bins[1][k]++;
      }
    }
    
    //printf("fmutn[%d][1][0] is: %d\n",i,fmutn[i][1][0]);
    for(j=1; j <= fmutn[i][1][0]; j++) {
      d = alpha_fav*pow(mod(abs(fmutn[i][1][j]),lb_modulo)*fav_scale,gamma_fav);
      k = 50 + (int)(d/fav_bin_width);
      //printf("k2, d is: %d %lf\n",k,d);
      if(fmutn[i][0][j] < 0) {
	      if(k < 100) fitness_bins[0][k]++;
      } else {
	      if(k < 100) fitness_bins[1][k]++;
      }
    }
  }
  
  //  Compute fitness values for bin boundaries and bin centers.
  
  for(k=0; k < 51; k++) {
    bin_fitness[k] = exp(-del_bin_width*k);
    if (k > 0) {
      bin_fitness_boxwidth[k-1] = fabs(bin_fitness[k] - bin_fitness[k-1]);
      bin_fitness_midpoint[k-1] = (bin_fitness[k] + bin_fitness[k-1])/2.;
    }
  }

  for(k=50; k < 101; k++) {
    bin_fitness[k] = max_fav_fitness_gain
      *exp(-fav_bin_width*(k - 50));
    if (k > 50) {
      bin_fitness_boxwidth[k-1] = fabs(bin_fitness[k] - bin_fitness[k-1]);
      bin_fitness_midpoint[k-1] = (bin_fitness[k] + bin_fitness[k-1])/2.;
    }
  }
  

   // Compute the frequency of output of the mutation accumulation
   // statistics.

  accum_gen = 1000000;
    if(gen%100 == 0) {
     if(gen <= 500) {
       accum_gen = 100;
     } else if(gen <=  1000) {
       accum_gen = 500;
     } else if(gen <=  5000) {
       accum_gen = 1000;
     } else if(gen <= 10000) {
       accum_gen = 5000;
     } else {
       accum_gen = 10000;
     }
     if(gen == 100) for(k=0; k<50; k++) accum[k] = 0.;
   }
   
   /* Output the number of accumulated deleterious dominant mutations *
    * in each of 50 bins during the previous accum_gen generations at *
    * appropriate intervals, along with other relevant information.   */
   
   if(gen%accum_gen == 0) {
     
     fprintf(fp26,"#\n#           Generation    Accumulation Interval\n");
     fprintf(fp26,"%19d %19d\n#\n", gen, accum_gen);
     fprintf(fp26,"#           Accumulation Over Previous Interval\n");
     fprintf(fp26,"#\n# bin  fitness effect    actual      expected       ");
     
     if(gen == 100) fprintf(fp26,"ratio  expected fraction\n");
     else           fprintf(fp26,"ratio  total accum\n");
     
     for(k=0; k<50; k++) {
       del = fitness_bins[1][k] - accum[k];
       del_no_sel = refr_bins[k]*(1. - fraction_recessive)            
	 *accum_gen/gen;
       ratio = del/del_no_sel;
       if(gen == 100) {
	 x0    = refr_bins[k]*(1. - fraction_recessive)/mutn_sum;
	 fprintf(fp26,"%4d %14.3e %12d %14.2f %14.7f %14.3e\n",
		 k, bin_fitness_midpoint[k], (int)del, del_no_sel, ratio, x0);
       } else 
	 fprintf(fp26,"%4d %14.3e %12d %14.2f %14.7f %14d\n",
		 k, bin_fitness_midpoint[k], (int)del, del_no_sel, ratio, 
		 (int)fitness_bins[1][k]);
     }
     
     for(k=0; k<50; k++) accum[k] = fitness_bins[1][k];

     fflush(fp26);
     
   }
   
   // Normalize the binned mutations by the reciprocal of the expected
   // number of mutations per bin in the absence of selection.
   
   for(k=0; k<100; k++) {
     
     if(refr_bins[k] > 0. && fraction_recessive > 0.) 
       fitness_bins[0][k] /= fraction_recessive*refr_bins[k];
     else
       fitness_bins[0][k] = 0.;
     
     if(refr_bins[k] > 0. && fraction_recessive < 1.) 
       fitness_bins[1][k] /= (1. - fraction_recessive)*refr_bins[k];
     else
       fitness_bins[1][k] = 0.;
     
   }

   // Perform an iteration of smoothing on the fitness_bin values
   // using a three-point average.  Iterate three times.
  
   for(i=0; i<3; i++) {
     fm1 = fitness_bins[0][0];
     fm2 = fitness_bins[1][0];
     for(k=1; k<49; k++) {
       av1 = fitness_bins[0][k] + 0.5*(fm1 + fitness_bins[0][k+1]);
       fm1 = fitness_bins[0][k];
       work[0][k] = 0.5*av1;
       av2 = fitness_bins[1][k] + 0.5*(fm2 + fitness_bins[1][k+1]);
       fm2 = fitness_bins[1][k];
       work[1][k] = 0.5*av2;
     }
     fitness_bins[0][49] = 0.5*(fitness_bins[0][48] + fitness_bins[0][49]);
     fitness_bins[1][49] = 0.5*(fitness_bins[1][48] + fitness_bins[1][49]);
     for(k=1; k<49; k++) {
       fitness_bins[0][k] = work[0][k];
       fitness_bins[1][k] = work[1][k];
     }
   }
   
   // For favorable distribution, limit maximum to a value of 100.
   // To increase the smoothness, iterate the smoothing two times.
   
   for(k=50; k<100; k++) {
     fitness_bins[0][k] = min(100., fitness_bins[0][k]);
     fitness_bins[1][k] = min(100., fitness_bins[1][k]);
   }
   
   for(i=0; i<2; i++) {
     fm1 = fitness_bins[0][50];
     fm2 = fitness_bins[1][50];
     for(k=51; k<99; k++) {
       av1 = fitness_bins[0][k] + 0.5*(fm1 + fitness_bins[0][k+1]);
       fm1 = fitness_bins[0][k];
       work[0][k] = 0.5*av1;
       av2 = fitness_bins[1][k] + 0.5*(fm2 + fitness_bins[1][k+1]);
       fm2 = fitness_bins[1][k];
       work[1][k] = 0.5*av2;
       fitness_bins[0][k] = work[0][k];
       fitness_bins[1][k] = work[1][k];
     }
   }

   //  Write the fitness bin information for deleterious mutations
   //  to standard output.

   if (!is_parallel) {
      printf("\n              Fraction of mutations retained versus ");
      printf("fitness effect");
      printf("\n effect:");
      for(k=2; k<50; k+=5) printf("%7.4f",bin_fitness_midpoint[k]);
      printf("\n recess:");
      for(k=2; k<50; k+=5) printf("%7.4f",fitness_bins[0][k]);
      printf("\n domint:");
      for(k=2; k<50; k+=5) printf("%7.4f",fitness_bins[1][k]);
   }

   fprintf(fp9,"              Fraction of mutations retained versus fitness effect\n");
   fprintf(fp9,"\n effect:");
   for(k=2; k<50; k+=5) fprintf(fp9,"%7.4f",bin_fitness_midpoint[k]);
   fprintf(fp9,"\n recess:");
   for(k=2; k<50; k+=5) fprintf(fp9,"%7.4f",fitness_bins[0][k]);
   fprintf(fp9,"\n domint:");
   for(k=2; k<50; k+=5) fprintf(fp9,"%7.4f",fitness_bins[1][k]);

   rewind(fp8);

   fprintf(fp8,"# generation = %d\n",gen);
   fprintf(fp8,"# deleterious mutations\n");
   fprintf(fp8,"# bin_fitness   recessive  dominant   box_width\n");

   for(k=0; k<50; k++) fprintf(fp8, "%13.5e %11.5f %11.5f %13.5e\n",
         bin_fitness_midpoint[k], fitness_bins[0][k], fitness_bins[1][k],
         bin_fitness_boxwidth[k]);

      fprintf(fp8,"# favorable mutations\n");
      fprintf(fp8,"# bin_fitness   recessive  dominant   box_width\n");
      for(k=50; k<100; k++) fprintf(fp8, "%13.5e %11.5f %11.5f %13.5e %13.5e\n",
           bin_fitness_midpoint[k], fitness_bins[0][k], fitness_bins[1][k],
           bin_fitness_boxwidth[k], refr_bins[k]);
      fflush(fp8);

     if(is_parallel) {
        mpi_davg((double *)fitness_bins,(double *)par_fitness_bins,200);
        // For some reason in this specific case, the above call only
        // does a summation
        for(k=0; k<100; k++) {
           par_fitness_bins[0][k] /= num_procs;
           par_fitness_bins[1][k] /= num_procs;
        }
        rewind(fp18);
        if(myid==0) {
           fprintf(fp18,"# generation = %d\n",gen);
           fprintf(fp18,"# deleterious mutations\n");
           fprintf(fp18,"# bin_fitness   recessive  dominant   box_width\n");
           for(k=0; k<50; k++) fprintf(fp18, "%13.5e %11.5f %11.5f %13.5f\n",
              bin_fitness[k], par_fitness_bins[0][k], par_fitness_bins[1][k],
              bin_fitness_boxwidth[k]);

           fprintf(fp18,"# favorable mutations\n");
           fprintf(fp18,"# bin_fitness   recessive  dominant   box_width\n");
           for(k=50; k<100; k++) fprintf(fp18, "%13.5e %11.5f %11.5f %13.5f\n", 
              bin_fitness[k], par_fitness_bins[0][k], par_fitness_bins[1][k],
              bin_fitness_boxwidth[k]);
           fflush(fp18);
        }
     }
     
     // Compute the current values of the selection thresholds for both
     // dominant and recessive deleterious mutations.
     
     //  Compute estimate for the current deleterious dominant threshold.
     
     k   = 0;
     k0  = 0;
     sum = (1. - fraction_recessive)*refr_bins[49];
     del_dom_thres = 0.;
     
     while(k < 50 && sum > 5000) { 
       if(fitness_bins[1][k] > 0.25 && k0 == 0) k0 = k;
       if(fitness_bins[1][k] > 0.75 && k > k0 + 1) {
	 x0 = 0.;
	 y0 = 0.;
	 for(i=k0-1; i<k-1; i++) {
	   x0 = x0 + i - 0.5 + 1;
	   y0 = y0 + fitness_bins[1][i];
	 }
	 x0 = x0/(k - k0);
	 y0 = y0/(k - k0);
	 s = 0.;
	 d = 0.;
	 for(i=k0-1; i<k-1; i++) {
	   s += (i - 0.5 - x0 + 1)*(fitness_bins[1][i] - y0);
	   d += pow(i - 0.5 - x0 + 1,2);
	 }
	 b1 = s/d;
	 b0 =   y0 - b1*x0;
	 x1 = (0.5 - b0)/b1;
	 del_dom_thres = exp(-x1*del_bin_width);
	 k = 50;
       }
       k++;
     }
     
     // Compute estimate for the current deleterious recessive threshold.
     
     k   = 0;
     k0  = 0;
     sum = fraction_recessive*refr_bins[49];
     del_rec_thres = 0.;
     
     while(k < 50 && sum > 5000) {
       if(fitness_bins[0][k] > 0.25 && k0 == 0) k0 = k;
       if(fitness_bins[0][k] > 0.75 && k > k0+1) {
	 x0 = 0.;
	 y0 = 0.;
	 for(i=k0; i<k; i++) {
           x0 += i - 0.5;
           y0 += fitness_bins[0][i-1];
	 }
	 x0 /= k - k0;
	 y0 /= k - k0;
	 s = 0.;
	 d = 0.;
	 for(i=k0; i<k; i++) {
           s += (i - 0.5 - x0)*(fitness_bins[0][i] - y0);
           d += pow(i - 0.5 - x0,2);
	 }
	 b1 = s/d;
	 b0 = y0 - b1*x0;
	 x1 = (0.5 - b0)/b1;
	 del_rec_thres = exp(-x1*del_bin_width);
	 k = 50;
       }
       k++;
     }
     
     // Compute estimate for the current favorable dominant threshold.
     
     // First find the bin with the maximum ratio of actual to expected
     // mutations if there were no selection, but restricted to ratios
     // less than 3.
     
     y0 = fitness_bins[1][95];
     k0 = 4;
     k  = 5;
     while(k < 50) {
       if(fitness_bins[1][100-k] > y0) {
	 y0 = fitness_bins[1][100-k];
	 k0 = k;
       }
       if(fitness_bins[1][100-k] >= 3.) k = 50;
       k++;
     }
     
     // Now find the first bin with k < k0 that is bracketed by ratios
     // below and above 2.0.
     
     j = k0 - 2;
     for(k=k0-1; k < k0-5; k--) 
       if(fitness_bins[1][99-k] > 2.0 && fitness_bins[1][100-k] <= 2.0) 
	 j = k;
     
     sum = (1. - fraction_recessive)*refr_bins[99];
     fav_dom_thres = 0.;
     
     if(sum > 2000 && y0 > 2.5) {
       
       // Use simple linear interpolation to find the fitness effect
       // value corresponding to the ratio of 2.0.
       
       s  = (fitness_bins[1][99-j] - fitness_bins[1][100-j])/
	 (bin_fitness_midpoint[99-j] - bin_fitness_midpoint[100-j]);
       y0 =  fitness_bins[1][100-j] - s*bin_fitness_midpoint[101-j];
       fav_dom_thres = max(0., (2.0 - y0)/s);
       
     }
     
     if(!is_parallel) {
       if(del_dom_thres > 0.) {
	 x0 = pow(-log(del_dom_thres)/alpha_del,1/gamma_del);
	 printf("\ndeleterious selection threshold   =%10.3e\n", del_dom_thres);
	 printf("deleterious fraction unselectable =%6.3f\n", 1. - x0);
       }
       if(fav_dom_thres > 0.) {
	 x0 = pow(-log(fav_dom_thres/max_fav_fitness_gain)/alpha_fav,
		  1/gamma_fav);
	 printf("  favorable selection threshold   =%10.3e\n", fav_dom_thres);
	 printf("  favorable fraction unselectable =%6.3f\n", 1. - x0);
       }
     } else {
#ifdef MPICH
       mpi_davg(&del_dom_thres,&par_del_dom_thres,1);
       if(myid == 0. && par_del_dom_thres > 0.) {
	 x0 = pow(-log(par_del_dom_thres)/alpha_del,1/gamma_del);
	 printf("\ndeleterious selection threshold   =%10.3e\n",
		par_del_dom_thres);
	 printf("deleterious fraction unselectable =%6.3f\n", 1. - x0);
       }
       mpi_davg(&fav_dom_thres,&par_fav_dom_thres,1);
       if(myid == 0 && par_fav_dom_thres > 0.) {
	 x0 = pow(-log(par_fav_dom_thres/max_fav_fitness_gain)
		  /alpha_fav,1/gamma_fav);
	 printf("  favorable selection threshold   =%10.3e\n",
		par_fav_dom_thres);
	 printf("  favorable fraction unselectable =%6.3f\n", 1. - x0);
       }
#endif
     }
     
     if(del_dom_thres > 0.) {
       x0 = pow(-log(del_dom_thres)/alpha_del,1/gamma_del);
       fprintf(fp9, "deleterious selection threshold   =%10.3e\n", 
	       del_dom_thres);
       fprintf(fp9, "deleterious fraction unselectable =%6.3f\n", 
	       1. - x0);
     }
     if(fav_dom_thres > 0.) {
       x0 = pow(-log(fav_dom_thres/max_fav_fitness_gain)/alpha_fav,
		1/gamma_fav);
       fprintf(fp9,"  favorable selection threshold   =%10.3e\n", fav_dom_thres);
       fprintf(fp9,"  favorable fraction unselectable =%6.3f\n", 1. - x0);
     }
     
     if(is_parallel && myid==0) oneortwo=2;
     else oneortwo=1;
     
     fid=fp25;
     
     for(i=1; i<=2; i++) {
       if(oneortwo==2) {
	 fid=fp35;
#ifdef MPICH
	 del_dom_thres = par_del_dom_thres;
	 fav_dom_thres = par_fav_dom_thres;
#endif
       }
       
       // Do not write out zero threshold values, instead use NaN's
       if (del_dom_thres == 0 && del_rec_thres == 0 && fav_dom_thres == 0) 
	 fprintf(fid, "%10d   NaN   NaN   NaN\n", gen);
       else if (del_dom_thres == 0 && del_rec_thres == 0) 
	 fprintf(fid, "%10d   NaN   NaN %15.3e\n", gen, fav_dom_thres);
       else if (del_dom_thres == 0 && fav_dom_thres == 0) 
	 fprintf(fid, "%10d   NaN  %15.3e  NaN\n", gen, del_rec_thres);
       else if (del_rec_thres == 0 && fav_dom_thres == 0) 
	 fprintf(fid, "%10d  %15.3e  NaN   NaN\n", gen, del_dom_thres);
       else if(del_dom_thres == 0) 
	 fprintf(fid, "%10d NaN %15.3e %15.3e\n", 
		 gen, del_rec_thres, fav_dom_thres);
       else if(del_rec_thres == 0) 
	 fprintf(fid, "%10d %15.3e NaN %15.3e\n",
		 gen, del_dom_thres, fav_dom_thres);
       else if (fav_dom_thres == 0) 
	 fprintf(fid, "%10d %15.3e %15.3e NaN\n", 
		 gen, del_dom_thres, del_rec_thres);
       else
	 fprintf(fid, "%10d %15.3e %15.3e %15.3e\n",
		 gen, del_dom_thres, del_rec_thres, fav_dom_thres);
     }
     
     fflush(fid);
     
}

void diagnostics_mutn_bins_plot4(int *** dmutn, int *** fmutn, double *** linkage_block_fitness, int **** lb_mutn_count, int gen) {
  
  int i, j, k, lb;
  double bins_mutns[2][200] = {{0}}, expn_bins[200] = {0};
  double haplotype_bins[200] = {0}, haplotype_bin_width, avg_lb_effect;
  double lb_fitness_frac_positive, num_del_lb, num_fav_lb;
  double y0;
  double d, x0, x1, z0, z1, del_bin_width, fav_bin_width;
#ifdef MPICH
  double par_haplotype_bins[200] = {0}, par_avg_lb_effect;
  double par_lb_fitness_frac_pos;
#endif
  
  //  Generate the theoretical distribution curves for plotting.
  
  haplotype_bin_width = 1.e-04;
  del_bin_width = haplotype_bin_width;
  fav_bin_width = 0.01*min(0.01, max_fav_fitness_gain);
  
  x0 = 1.;
  for(k=99; k >= 0; k--) {
    x1 = pow(-log(del_bin_width*(100-k))/alpha_del,1./gamma_del);
    expn_bins[k] = (1. - frac_fav_mutn)*(x0 - x1);
    x0 = x1;
  }
  
  z0 = 1.;

  for(k=99; k >= 0; k--) {
    x0 = fav_bin_width*(100-k)/max_fav_fitness_gain;
    x0 = min(1., x0);
    z1 = pow(-log(x0)/alpha_fav,1./gamma_fav);
    expn_bins[199-k] = frac_fav_mutn*(z0 - z1);
    if(x0 > fav_scale*lb_modulo) expn_bins[199-k] = 0.;
    z0 = z1;
  }
  
  if(expn_bins[99] > 0.) 
     for(i=0;   i<100; i++) expn_bins[i] /= expn_bins[99];
  if(expn_bins[100] > 0.) 
     for(i=100; i<200; i++) expn_bins[i] /= expn_bins[100];
  
  //  Compute statistics on favorable and recessive mutations.
  
  for(i=0; i<200; i++) haplotype_bins[i] = 0;
  
  for(i=0; i < current_pop_size; i++) {
    
    for(lb=0; lb < num_linkage_subunits; lb++) {
      
      y0 = (linkage_block_fitness[i][0][lb] - (double)1)
	/haplotype_bin_width;
      
      if(y0 < 1.e-20) {
	k = 99 + (int)y0;
	if(k >= 0)    haplotype_bins[k]++;
      } else {
	k = 100 + (int)y0;
	if(k < 200) haplotype_bins[k]++;
      }
      
      y0 = (linkage_block_fitness[i][1][lb] - (double)1)
	/haplotype_bin_width;
      
      if(y0 < 1.e-20) {
	k = 99 + (int)y0;
	if(k >= 0)    haplotype_bins[k]++;
      } else {
	k = 100 + (int)y0;
	if(k < 200) haplotype_bins[k]++;
      }
      
    }
    
    //printf("\nDEBUG1: dmutn %d %d\n",i,dmutn[i][1][1]);
    //for(j=2; j <= dmutn[i][1][1]+1; j++) printf("%d ",dmutn[i][1][j]);
    //printf("\nDEBUG2: dmutn %d %d\n",i,dmutn[i][2][1]);
    //for(j=2; j <= dmutn[i][2][1]+1; j++) printf("%d ",dmutn[i][2][j]);
    
    for(j=1; j <= dmutn[i][0][0]; j++) {
      d = exp(-alpha_del*pow(mod(abs(dmutn[i][0][j]),lb_modulo)
			      *del_scale,gamma_del));
      k = 99 - (int)(d/del_bin_width);
      if(dmutn[i][0][j] < 0) {
	if(k >= 0 && k < 100) bins_mutns[0][k]++;
      } else {
	if(k >= 0 && k < 100) bins_mutns[1][k]++;
      }
    }
    
    for(j=1; j <= dmutn[i][1][0]; j++) {
      d = exp(-alpha_del*pow(mod(abs(dmutn[i][1][j]),lb_modulo)
			      *del_scale,gamma_del));
      k = 99 - (int)(d/del_bin_width);
      if(dmutn[i][1][j] < 0) {
	if(k >= 0 && k < 100) bins_mutns[0][k]++;
      } else {
	if(k >= 0 && k < 100) bins_mutns[1][k]++;
      }
    }
    
    for(j=1; j <= fmutn[i][0][0]; j++) {
      d = exp(-alpha_fav*pow(mod(abs(fmutn[i][0][j]),lb_modulo)
			      *fav_scale,gamma_fav));
      k = 100 + (int)(d*max_fav_fitness_gain/fav_bin_width);
      
      if(fmutn[i][0][j] < 0) {
	if(k >= 100 && k < 200) bins_mutns[0][k]++;
      } else {
	if(k >= 100 && k < 200) bins_mutns[1][k]++;
      }
    }
    
    for(j=1; j <= fmutn[i][1][0]; j++) {
      d = exp(-alpha_fav*pow(mod(abs(fmutn[i][1][j]),lb_modulo)
			      *fav_scale,gamma_fav));
      k = 100 + (int)(d*max_fav_fitness_gain/fav_bin_width);
      
      if(fmutn[i][1][j] < 0) {
	if(k >= 100 && k < 200) bins_mutns[0][k]++;
      } else {
	if(k >= 100 && k < 200) bins_mutns[1][k]++;
      }
    }
    
  }
  
  x0 = 1.e-10;
  for(k=0; k < 200; k++) x0 = max(x0, haplotype_bins[k]);
  for(k=0; k < 200; k++) haplotype_bins[k] /= x0;
  
  if(tracking_threshold != 1.0) {
    for (k=0; k<100; k++) {
       bins_mutns[0][k] *= expn_bins[98]/(bins_mutns[0][98] + 1.e-10);
       bins_mutns[1][k] *= expn_bins[98]/(bins_mutns[1][98] + 1.e-10);
    }
    for (k=100; k<200; k++) {
       bins_mutns[0][k] = expn_bins[101]/(bins_mutns[0][101] + 1.e-10);
       bins_mutns[1][k] = expn_bins[101]/(bins_mutns[1][101] + 1.e-10);
    }
  }
  
  rewind(fp4);
  fprintf(fp4,"# generation = %8d\n", gen);
  fprintf(fp4,"# effect-bin  theory(red) lb-fitns(gr) dominants  recessives\n");
  
  for(k=0; k < 200; k++) 
    fprintf(fp4, "%12.4e %12.4e %12.4e %12.4e %12.4e\n",
	    (k - 99.5)*haplotype_bin_width, expn_bins[k], 
	    haplotype_bins[k], 
	    bins_mutns[1][k], bins_mutns[0][k]);
  
  avg_lb_effect = (post_sel_fitness - (double)1)/
    (2*num_linkage_subunits);
  
  num_del_lb = 0;
  num_fav_lb = 0;
  
  for(k=0; k < 100; k++) 
    num_del_lb += haplotype_bins[k];
  
  for(k=100; k < 200; k++) 
    num_fav_lb += haplotype_bins[k];
  
  lb_fitness_frac_positive = num_fav_lb/(num_del_lb + num_fav_lb);
  
  fprintf(fp4,"# favorable x-axis scaling = %12.4e\n",
          fav_bin_width/del_bin_width);
  fprintf(fp4,"# avg_linkage_block_effect = %12.4e\n", avg_lb_effect);
  fprintf(fp4,"# lb_fitness_percent_positive = %12.4f\n",
	  lb_fitness_frac_positive*100.);
  fflush(fp4);
  
#ifdef MPICH
  if (is_parallel) {
    mpi_davg(&avg_lb_effect,&par_avg_lb_effect,1);
    mpi_davg(&lb_fitness_frac_positive,&par_lb_fitness_frac_pos,1);
    mpi_davg(haplotype_bins,par_haplotype_bins,200);
    
    //     Note: currently only averaging haplotype_bins--
    //     all the other values are currently coming from processor 0
    if(myid==0) {
      rewind(fp14);
      fprintf(fp14,"# generation = %8d\n", gen);
      fprintf(fp14,"# effect-bin  theory(red) lb-fitns(g)  dominants  recessives\n");

      for(k=0; k < 200; k++) 
	fprintf(fp14,"%12.4e %12.4e %12.4e %12.4e %12.4e\n",
		(k - 99.5)*haplotype_bin_width,
		expn_bins[k], par_haplotype_bins[k],
		bins_mutns[1][k], bins_mutns[0][k]);
      fprintf(fp14,"# avg_linkage_block_effect = %12.4e\n", avg_lb_effect);
      fprintf(fp14,"# lb_fitness_percent_positive = %12.4f\n",
	      par_lb_fitness_frac_pos*100.);
      fflush(fp14);
    }
  }
#endif
}

/* This routine analyzes the distribution of the initial contrasting *
 *  alleles and their effects on overall fitness.  When bool         *
 *  variable list is TRUE, routine outputs a list of allele          *
 *  frequencies.                                                     */
void diagnostics_contrasting_alleles(int *** dmutn, int *** fmutn, 
            int **** count, double *cum_effect, float *initial_allele_effects, 
            int ica_count[], int max_size, bool list) 
{
  float w, effect, freq;
  int i, j, lb, m, indx;
  int *zygous;

  zygous = malloc(num_linkage_subunits*sizeof(int));

  w = multiplicative_weighting;
  
  for (i = 0; i < 2; i++)
    for (lb = 0; lb < num_linkage_subunits; lb++) 
      count[0][0][i][lb] = 0;

  for(i=0; i < pop_size; i++) cum_effect[i] = 1.0;
  
  for(i=0; i < current_pop_size; i++) {
    
    for(j=0; j < num_linkage_subunits; j++) zygous[j] = 0;
    
    for(m=1; m <= dmutn[i][0][0]; m++) {
      if(mod(dmutn[i][0][m], lb_modulo) == lb_modulo - 1) { 
	lb = dmutn[i][0][m]/lb_modulo;
	zygous[lb]++;
      }
    }
    
    for(m=1; m <= dmutn[i][1][0]; m++) {
       if(mod(dmutn[i][1][m],lb_modulo) == lb_modulo - 1) { 
 	  lb = dmutn[i][1][m]/lb_modulo;
	  zygous[lb]++;
       }
    }
    
    for(lb=0; lb < num_linkage_subunits; lb++) {
      if(zygous[lb] == 1) {
	effect = initial_allele_effects[lb]
	  *recessive_hetero_expression;
	count[0][0][0][lb]++;
	cum_effect[i] = (cum_effect[i] - 
			 ((double)1.0 - w)*effect)
	  * ((double)1.0 - w *effect);
      } else if(zygous[lb] == 2) {
	effect = initial_allele_effects[lb];
	count[0][0][0][lb] += 2; 
	cum_effect[i] = (cum_effect[i] - 
			 ((double)1.0 - w)*effect)
	  * ((double)1.0 - w *effect);
      }
    }

    
    for(j=0; j < num_linkage_subunits; j++) zygous[j] = 0;
    
    for(m=1; m <= fmutn[i][0][0]; m++) {
      if(mod(fmutn[i][0][m], lb_modulo) == lb_modulo - 1) { 
	lb = fmutn[i][0][m]/lb_modulo;
	zygous[lb]++;
      }
    }
    
    for(m=1; m <= fmutn[i][1][0]; m++) {
      if(mod(fmutn[i][1][m], lb_modulo) == lb_modulo - 1) { 
	lb = fmutn[i][1][m]/lb_modulo;
	zygous[lb]++;
      }
    }
    
    //printf("\n\nindividual %d\n",i);
    for(lb=0; lb < num_linkage_subunits; lb++) {
      if(zygous[lb] == 1) {
	effect = initial_allele_effects[lb]
	  *dominant_hetero_expression;
	count[0][0][1][lb]++;
	cum_effect[i] = (cum_effect[i] + ((double)1.0 - w)*effect)
	  * ((double)1.0 + w *effect);
      } else if(zygous[lb] == 2) {
	effect = initial_allele_effects[lb];
	count[0][0][1][lb] += 2;
	cum_effect[i] = (cum_effect[i] + ((double)1.0 - w)*effect)
	  * ((double)1.0 + w *effect);
      }
      //if(count[0][0][1][lb]>0) 
      //   printf("lb. count is %d.\t%d\n",lb,count[0][0][1][lb]);
    }
    
 }

  ica_count[0] = 0;
  ica_count[1] = 0;
  
  for(lb=0; lb < num_linkage_subunits; lb++) {
    ica_count[0] += count[0][0][0][lb];
    ica_count[1] += count[0][0][1][lb];
  }
  
  if(!list) {
    
    // Compute average effect of initial contrasting alleles.
    
    ica_mean_effect = 0.;
    
    for(i=0; i < current_pop_size; i++) {
      ica_mean_effect += cum_effect[i] - (double)1.;
    }
    
    ica_mean_effect /= current_pop_size;
    
    //     Compute mean frequency of positive initial contrasting alleles
    //     and the number of positive initial contrasting alleles fixed 
    //     and lost.
    
    fav_lost  = 0;
    fav_fixed = 0;
    fav_mean_freq = 0;
    
    for(lb=0; lb < num_linkage_subunits; lb++) {
      fav_mean_freq += count[0][0][1][lb];
      if(fabs(initial_allele_effects[lb]) > 0.) {
	if(count[0][0][1][lb] == 2*current_pop_size) fav_fixed++;
	if(count[0][0][1][lb] == 0) fav_lost++;
      }
    }
    
    //printf("*** fav_mean_freq is %f\n",fav_mean_freq);
    fav_mean_freq /= 2*current_pop_size*num_contrasting_alleles;
    
  } else {
    
    printf("   List of initial contrasting allele freqencies\n");
    printf("   and fitness effect values at end of run:\n");
    printf("      allele   linkage  favorable  homozygous\n");
    printf("               subunit  frequency    effect\n");

    if(is_parallel && myid==0) {
       fprintf(fp14,"   List of initial contrasting allele freqencies\n");
       fprintf(fp14,"   and fitness effect values at end of run:\n");
       fprintf(fp14,"      allele   linkage  favorable  homozygous\n");
       fprintf(fp14,"               subunit  frequency    effect\n");
    }

    indx = 0;
    
    for(lb=0; lb < num_linkage_subunits; lb++) {
      if(fabs(initial_allele_effects[lb]) > 0.) {
	indx++;
	freq = 0.5*(float)count[0][0][1][lb]/(float)current_pop_size;
	
	if(!is_parallel) 
	  printf("%10d %10d %12.4f %11.4f", 
		 indx, lb, freq, fabsf(initial_allele_effects[lb]));
	fprintf(fp9,"%10d %10d %12.4f %11.4f", 
		indx, lb, freq, fabsf(initial_allele_effects[lb]));
	
	if(!is_parallel) printf("\n");
	printf("\n");
	
      }
    }
  }
}


void diagnostics_polymorphisms_plot7(int *** dmutn, int *** fmutn,
                                     int max_size, int gen)
{
  int *mutn_count, *mutn_list;
  int dcount, fcount;
  int i, j, k, m, lb, it, it0, ie, list_count;
  int mymax;
  int num_falleles[3], num_dalleles[3], dwarn, fwarn;
  double dpbin[NBINS] = {0}, dpbin_count[NBINS], dpbin_max, pbin_width;
  double fpbin[NBINS] = {0}, fpbin_count[NBINS], fpbin_max;
  double dsum, fsum, fe_bin_width, febin[100][10] = {{0}};
  float bin_fitness[11], bin_center[11];
  double **mfirst;
#ifdef MPICH
  int par_dcount, par_fcount;
  int par_num_falleles[3], par_num_dalleles[3];
  double par_dpbin[NBINS], par_fpbin[NBINS];
  double par_dpbin_count[NBINS], par_fpbin_count[NBINS];
#endif
  int  mutn_limit;
  bool new_mutn;

  mymax = max(25000,max_del_mutn_per_indiv);
  mutn_count = malloc(mymax*sizeof(int));
  mutn_list = malloc(mymax*sizeof(int));
  mfirst = (double **)malloc_double2d(mymax, 2); 

  pbin_width = 2.*current_pop_size/NBINS;
  fe_bin_width = -log(tracking_threshold)/10.;
  
  if((int)pbin_width == 0) {
    printf("Polymorphism analysis skipped because population");
    printf(" size is too small.\n");
    return;
  }
  
  // Compute statistics on deleterious polymorphisms.  To keep the
  // required computer time reasonable, limit the maximum number of
  // polymorphisms to be considered to 25000.  This may result in
  // the analysis being performed only over a portion of the total
  // total population.

  dsum   = 0.;
  for(i=0; i<NBINS; i++) dpbin[i] = 0.;

  for(i=0; i<mymax; i++)
    for(j=0; j<2; j++)  
      mfirst[i][j] = 2;
  dwarn  = 0;

  // Loop over recessives first and dominants second.
  
  it0 = 1;

  if(fraction_recessive == 0.) it0 = 2;

  //printf("num_linkage_subunits == %d",num_linkage_subunits);

  for(it=it0; it<=2; it++) {

    for(lb=1; lb <= num_linkage_subunits; lb++) {

      if(it == 1) 
        mutn_limit = -lb_modulo*(num_linkage_subunits - lb);
      else
        mutn_limit =  lb_modulo*lb;

// mutn_list is the list of mutation indices being analyzed
// for their frequency.
// mutn_count is the number of occurrences of each mutation
// in mutn_list.
// list_count is the total number of deleterious alleles
// being analyzed.

      for(i=0; i<mymax; i++) mutn_list[i] = 0;
      for(i=0; i<mymax; i++) mutn_count[i] = 0;
      list_count = 0;

      for(i=0; i<current_pop_size; i++) {
        for(j=0; j<2; j++) {
          m = mfirst[i][j];
          while(dmutn[i][j][m-1] < mutn_limit &&
            m <= dmutn[i][j][0]+1 && list_count < 25000) {
            if(abs(dmutn[i][j][m-1])%lb_modulo != lb_modulo-1) {
              new_mutn = true;
              for(k=0; k<list_count; k++) {
                if(dmutn[i][j][m-1] == mutn_list[k]) {
                  mutn_count[k]++;
                  new_mutn = false;
                }
              }
              if(new_mutn) {
                list_count = min(25000, list_count+1);
                mutn_list[list_count-1] = dmutn[i][j][m-1];
                mutn_count[list_count-1] = 1;
              }
            }
            m++;
          }
          mfirst[i][j] = m;
        }
      }

// Load the allele hit counts into the proper bins and count
// the total number of hits accumulated.

      for(k=0; k<list_count; k++) {
        j = min(99, (int)(mutn_count[k]/pbin_width));
        dpbin[j]++;
        dsum += mutn_count[k];
        ie = (int)(pow(alpha_del*(abs(mutn_list[k])%lb_modulo)*del_scale,
                       gamma_del)/fe_bin_width);
        febin[j][ie]++;
      }
      if(list_count == 25000) dwarn = 1;
    }
  }

// Compute statistics on favorable polymorphisms.  To keep the
// required computer time reasonable, limit the maximum number of
// polymorphisms to be considered to 25000.  This may result in
// the analysis being performed only over a portion of the total
// total population.

  fsum   = 0.;
  for(i=0; i<NBINS; i++) fpbin[i] = 0.;
  for(i=0; i<mymax; i++)
    for(j=0; j<2; j++) 
      mfirst[i][j] = 2;
  fwarn  = 0;

  if(frac_fav_mutn > 0.) {

    // Loop over recessives first and dominants second.

    it0 = 1;
    if(fraction_recessive == 0.) it0 = 2;

    for(it=it0; it<2; it++) {

      for(lb=1; lb <= num_linkage_subunits; lb++) {

        if(it == 1) 
          mutn_limit = -lb_modulo*(num_linkage_subunits - lb);
        else
          mutn_limit =  lb_modulo*lb;

        for(i=0; i<mymax; i++) mutn_list[i]  = 0;
        for(i=0; i<mymax; i++) mutn_count[i]  = 0;
        list_count = 0;

        for(i=0; i < current_pop_size; i++) {
          for(j=0; j<2; j++) {
            m = mfirst[i][j];
            while(fmutn[i][j][m-1] < mutn_limit &&
              m <=  fmutn[i][j][0]+1 && list_count < 25000) {
              if(abs(fmutn[i][j][m-1])%lb_modulo != lb_modulo-1) {
                new_mutn = true;
                for(k=0; k<list_count; k++) {
                  if(fmutn[i][j][m-1] == mutn_list[k]) {
                    mutn_count[k]++;
                    new_mutn = false;
                  }
                }
                if(new_mutn) {
                  list_count = min(25000, list_count + 1);
                  mutn_list [list_count-1] = fmutn[i][j][m-1];
                  mutn_count[list_count-1] = 1;
                }
              }
              m++;
            }
            mfirst[i][j] = m;
          }
        }

        for(k=0; k<list_count; k++) {
          j = min(99, (int)(mutn_count[k]/pbin_width));
          fpbin[j]++; 
          fsum += mutn_count[k];
        }

        if(list_count == 25000) fwarn = 1;

      }

    }

  }

  dpbin_max = 1;
  fpbin_max = 1;

  num_dalleles[1] = 0;
  num_falleles[1] = 0;

  for(j=0; j<NBINS; j++) {
    dpbin_max = max(dpbin_max, dpbin[j]);
    fpbin_max = max(fpbin_max, fpbin[j]);
    if(j >= 1 && j < NBINS-1) {
      num_dalleles[1] += dpbin[j];
      num_falleles[1] += fpbin[j];
    }
  }

  num_dalleles[0] = dpbin[0];
  num_falleles[0] = fpbin[0];
  num_dalleles[2] = dpbin[NBINS-1];
  num_falleles[2] = fpbin[NBINS-1];

  dcount = dpbin[0] + num_dalleles[1] + dpbin[99];
  fcount = fpbin[0] + num_falleles[1] + fpbin[99];
  
  for(i=0; i<NBINS; i++) {
    dpbin_count[i] = dpbin[i];
    dpbin[i] /= dpbin_max;
    
    fpbin_count[i] = fpbin[i];
    fpbin[i] /= fpbin_max;
  }

  // Compute values for fitness effect bin centers.
  for(k=1; k<=11; k++) {
    bin_fitness[k] = exp(-fe_bin_width*(k - 1));
    if (k > 0) bin_center[k-1] = 0.5*(bin_fitness[k] + bin_fitness[k-1]);
  }
  
  fprintf(fp11,"# generation = %8d\n", gen);
  fprintf(fp11,"# frequency del_normalized fav_normalized del_count fav_count\n");
  for(k=0; k<NBINS; k++) 
    fprintf(fp11,"%11d %15.11f %15.11f %11.0f %11.0f\n", 
	    k+1, dpbin[k], fpbin[k], dpbin_count[k], fpbin_count[k]);

  fprintf(fp19,"# generation = %d\n",gen);
  fprintf(fp19,"#         Table of polymorphism frequency vs. fitness effect category\n");
  fprintf(fp19,"#\n#freq                 fitness effect category center value\n");
  for(k=1; k<5; k++)  fprintf(fp19,"%7.0e",bin_center[k]);
  for(k=5; k<11; k++) fprintf(fp19,"%8.1e",bin_center[k]);
  fprintf(fp19,"\n");
  for(k=0; k<100; k++) { 
    fprintf(fp19,"%3d",k);
    for(j=0; j<10; j++) fprintf(fp19,"%7d",(int)febin[k][j]);
    fprintf(fp19,"\n");
  }
  fprintf(fp19,"\n");
  for(j=1; j<5;  j++) fprintf(fp19,"%7.0e",bin_center[j]);
  for(j=5; j<11; j++) fprintf(fp19,"%8.1e",bin_center[j]);
  fprintf(fp19,"\n");
  fflush(fp19);
  
  // keep a second "snapshot" file of latest polymorphism data
  rewind(fp13);
  fprintf(fp13,"# generation = %d\n",gen);
  fprintf(fp13,"# frequency del_normalized fav_normalized   del_count fav_count\n");
  fprintf(fp13,"#           Allele summary statistics (tracked mutations only):\n");
  for(k=0; k<100; k++)
    fprintf(fp13,"%11d %15.11f %15.11f %11.0f %11.0f\n", 
                  k, dpbin[k], fpbin[k], dpbin_count[k], fpbin_count[k]);

  fflush(fp13);

  fprintf(fp11,"#           Allele summary statistics (tracked mutations only):\n");
  fprintf(fp11,"#   (Statistics are based on %8.1f tracked deleterious mutations.)\n",dsum);
  fprintf(fp11,"#         (Results are scaled to account for all linkage blocks.)\n");
  fprintf(fp11,"#    Very rare     Polymorphic        Fixed          Total\n");
  fprintf(fp11,"#      (0-1%%)        (1-99%%)         (100%%)\n");
  fprintf(fp11,"#%12d%12d%12d%12d\tdeleterious\n", num_dalleles[0],
                num_dalleles[1], num_dalleles[2], dcount);
  fprintf(fp11,"#%12d%12d%12d%12d\tfavorable\n", num_falleles[0],
                num_falleles[1], num_falleles[2], fcount);

  fflush(fp11);

  if(dwarn == 1) fprintf(fp11,"# Warning: Number of deleterious polymorhisms exceeded the linkage block limit of 25000\n");
  if(fwarn == 1) fprintf(fp11,"# Warning: Number of   favorable polymorhisms exceeded the linkage block limit of 25000");
  
  if(myid == 0) {
    printf("\n\n           Allele summary statistics (tracked mutations only):\n");
    printf("   (Statistics are based on %8.1f tracked deleterious mutations.)\n",dsum);
    printf("         (Results are scaled to account for all linkage blocks.)\n");
    printf("    Very rare     Polymorphic        Fixed          Total\n");
    printf("      (0-1%%)        (1-99%%)         (100%%)\n");
    printf("%12d%12d%12d%12d\tdeleterious\n", num_dalleles[0],
                num_dalleles[1], num_dalleles[2], dcount);
    printf("%12d%12d%12d%12d\tfavorable\n", num_falleles[0],
                num_falleles[1], num_falleles[2], fcount);
  }
  
  fprintf(fp9,"\n\n           Allele summary statistics (tracked mutations only):\n");
  fprintf(fp9,"   (Statistics are based on %8.1f tracked deleterious mutations.)\n",dsum);
  fprintf(fp9,"         (Results are scaled to account for all linkage blocks.)\n");
  fprintf(fp9,"    Very rare     Polymorphic        Fixed          Total\n");
  fprintf(fp9,"      (0-1%%)        (1-99%%)         (100%%)\n");
  fprintf(fp9,"%12d%12d%12d%12d\tdeleterious\n", num_dalleles[0],
                num_dalleles[1], num_dalleles[2], dcount);
  fprintf(fp9,"%12d%12d%12d%12d\tfavorable\n", num_falleles[0],
                num_falleles[1], num_falleles[2], fcount);
  fflush(fp9);
  
#ifdef MPICH
  if (is_parallel) {
    mpi_dsum(dpbin_count,par_dpbin_count,NBINS);
    mpi_dsum(fpbin_count,par_fpbin_count,NBINS);

    // Reanalyze the parallel data to get global statistics

    if (myid==0) {

       for(i=0; i<NBINS; i++) {
          dpbin[i] = 0;
          fpbin[i] = 0;
       }

       for(j=0; j<NBINS; j++) {
          // note: the following statement assumes all tribes
          //       have equal population sizes.  Otherwise, 
          //       would need to make an MPI_REDUCE call to 
          //       sum all the populations from the various tribes
          k = (int) j/num_procs;
          dpbin[k] += par_dpbin_count[j];
          fpbin[k] += par_fpbin_count[j];
          printf("%d %d %f %f\n",j,k,dpbin[k],par_dpbin[j]);
       }

       dpbin_max = 1;
       fpbin_max = 1;
      
       num_dalleles[1] = 0;
       num_falleles[1] = 0;
      
       for(j=0; j<NBINS; j++) {
          dpbin_max = max(dpbin_max, dpbin[j]);
          fpbin_max = max(fpbin_max, fpbin[j]);
          if(j >= 1 && j < NBINS-1) {
             num_dalleles[1] += dpbin[j];
             num_falleles[1] += fpbin[j];
          }
       }
       
      num_dalleles[0] = dpbin[0];
      num_falleles[0] = fpbin[0];
      num_dalleles[2] = dpbin[NBINS-1];
      num_falleles[2] = fpbin[NBINS-1];
      
      dcount = dpbin[0] + num_dalleles[1] + dpbin[99];
      fcount = fpbin[0] + num_falleles[1] + fpbin[99];
      
      for(i=0; i<NBINS; i++) {
	dpbin_count[i] = dpbin[i];
	dpbin[i] /= dpbin_max;
	
	fpbin_count[i] = fpbin[i];
	fpbin[i] /= fpbin_max;
      }

      rewind(fp21);
      fprintf(fp21,"# Note: data averaged over %4d tribes.\n", num_procs);
      fprintf(fp21,"# generation = %8d\n", gen);
      fprintf(fp21,"# frequency del_normalized fav_normalized  del_count fav_count\n");
      for(k=0; k<NBINS; k++) 
	fprintf(fp21,"%11d %15.11f %15.11f %11.2e %11.2e\n", 
		k, dpbin[k], fpbin[k], dpbin_count[k], fpbin_count[k]);
      fprintf(fp21,"# Allele summary statistics:\n");
      fprintf(fp21,"#    Very rare      Polymorphic        Fixed           Total\n");
      fprintf(fp21,"#     (0-1%%)         (1-99%%)           (100%%)\n");
      fprintf(fp21,"# %10d %10d %10d %10d deleterious\n",
	      num_dalleles[0], num_dalleles[1], num_dalleles[2], dcount);
      fprintf(fp21,"# %10d %10d %10d %10d favorable\n",
	      num_falleles[0], num_falleles[1], num_falleles[2], fcount);
      fflush(fp21);
    }
  }
#endif
}

void diagnostics_selection(double *fitness_pre_sel, double *fitness_post_sel, 
                           int total_offspring, int gen) 
{
  int i, j, jmax, max_bin;
  int sel_bins[2][150] = {{0}};
  int oneortwo;
#ifdef MPICH
  int par_sel_bins[2][150] = {{0}};
#endif
  float select_ratio, srm, srp;
  //char filename[80];
  //FILE *myfp;
  FILE *fid;
  double par_pre_sel_fitness, par_post_sel_fitness,
    par_pre_sel_geno_sd, par_pre_sel_pheno_sd,
    par_pre_sel_corr, par_post_sel_geno_sd,
    par_post_sel_pheno_sd, par_post_sel_corr;

  for(i=0; i<total_offspring; i++) {
    j = min(150, (int)(100.*fitness_pre_sel[i]));
    sel_bins[0][j]++;
  }
  
  for(i=0; i<current_pop_size; i++) {
    j = min(150, (int)(100.*fitness_post_sel[i]));
    sel_bins[1][j]++;
  }
  
  jmax    = 1;
  max_bin = 0;
  
  for(j=0; j<150; j++) {
    if(sel_bins[1][j] > max_bin) {
      max_bin = sel_bins[1][j];
      jmax = j;
    }
  }
  
#ifdef MPICH
  // Compute averages across processors.
  if(is_parallel) 
     mpi_isum((int*)sel_bins,(int*)par_sel_bins,300);

  if(is_parallel) {
     //sprintf(filename,"%s.%.3d.sel",case_id,myid+1);
     mpi_davg(&pre_sel_fitness,&par_pre_sel_fitness,1);
     mpi_davg(&post_sel_fitness,&par_post_sel_fitness,1);
     mpi_mybcast(par_post_sel_fitness,1);
     mpi_davg(&pre_sel_geno_sd, &par_pre_sel_geno_sd,1);
     mpi_davg(&pre_sel_pheno_sd, &par_pre_sel_pheno_sd,1);
     mpi_davg(&pre_sel_corr, &par_pre_sel_corr,1);
     mpi_davg(&post_sel_geno_sd, &par_post_sel_geno_sd,1);
     mpi_davg(&post_sel_pheno_sd, &par_post_sel_pheno_sd,1);
     mpi_davg(&post_sel_corr, &par_post_sel_corr,1);
  } else {
     //sprintf(filename,"%s.%.3d.sel",case_id,0);
  }
  //myfp = fopen(filename, "w");

  // If parallel is turned on, write two files
  // caseid.001.sel which contains results for one tribe
  // and caseid.000.sel which contains data averaged from
  // all tribes.
  if (is_parallel && myid==0) oneortwo = 2;
  else oneortwo = 1;
#endif

  fid = fp24;

  for(i = 1; i<=2; i++) {
#ifdef MPICH
    if(oneortwo==2) {
      fid = fp34;
      for(j = 0; j<150; j++) {
	sel_bins[0][j] = par_sel_bins[0][j];
	sel_bins[1][j] = par_sel_bins[1][j];
      }
      pre_sel_fitness  = par_pre_sel_fitness;
      post_sel_fitness = par_post_sel_fitness;
      pre_sel_geno_sd  = par_pre_sel_geno_sd;
      post_sel_geno_sd = par_post_sel_geno_sd;
      pre_sel_pheno_sd = par_pre_sel_pheno_sd;
    }
#endif
    
    rewind(fid);  
    
    fprintf(fid,"# generation = %8d\n",gen);
    fprintf(fid,"# heritability            = %9.5f\n",heritability);
    fprintf(fid,"# non scaling noise       = %9.5f\n",non_scaling_noise);
    fprintf(fid,"# pre  selection fitness  = %9.5f\n",pre_sel_fitness);
    fprintf(fid,"# post selection fitness  = %9.5f\n",post_sel_fitness);
    fprintf(fid,"# pre  selection geno  sd = %9.5f\n",pre_sel_geno_sd);
    fprintf(fid,"# post selection geno  sd = %9.5f\n",post_sel_geno_sd);
    fprintf(fid,"# pre  selection pheno sd = %9.5f\n",pre_sel_pheno_sd);
    fprintf(fid,"# post selection pheno sd = %9.5f\n",post_sel_pheno_sd);
    fprintf(fid,"# Effect of Selection on Phenotypic Fitness Distribution\n");
    fprintf(fid,"# Note: Individuals with phenotypic fitness > 1.5 are all put into last bin\n");
    fprintf(fid,"# max bin fitness     before selection     after selection     ratio\n");
    
    for(j=0; j<149; j++) {
      srm = sel_bins[1][j]/((float)sel_bins[0][j] + 0.000001);
      if(sel_bins[0][j] == 0) srm = -2.;
      
      srp = sel_bins[1][j+1]/((float)sel_bins[0][j+1] + 0.000001);
      if(sel_bins[0][j+1] == 0) srp = -2.;
      
      select_ratio = 0.5*(srm + srp);
      if(sel_bins[1][j] + sel_bins[1][j+1] <= max(4, max_bin/50)) 
	select_ratio = -1.;
      
      if(select_ratio > 0.) {
	fprintf(fid,"%10.3f %20d %20d %18.3f\n",
		j*0.01, sel_bins[0][j], sel_bins[1][j], select_ratio);
      } else {
	fprintf(fid,"%10.3f %20d %20d\t\t\?\n",
		j*0.01, sel_bins[0][j], sel_bins[1][j]);
      }
    }
    fflush(fid);
  }
}

/* Routine to dump timing information at end of run */

void profile(FILE *fp)
{ 
  sprintf(version_diagnostics,
          "$Id: diagnostics.c,v 1.73 2009/11/06 19:20:58 wes Exp $");

  if (myid == 0) {
    fprintf(fp,"\n          CPU SECONDS USED PER PROCESSOR:\n\n" );
    
    fprintf(fp,"          TOTAL     %10.3lf\n\n",sec[1]);
    
    sec[2] /= num_generations*migration_generations;
    sec[3] /= num_generations; //*fav_mutn+num_initial_fav_mutn
    sec[4] /= num_generations; //*actual_offspring
    sec[5] /= num_generations;
    sec[6] /= num_generations; // every 20 generations
    sec[7] /= num_generations; // every 100 generations

    fprintf(fp,"PRIMARY SUBROUTINES (MILLISECONDS/GEN/PROC):\n\n");
    
    fprintf(fp,"          MIGRATION        %lf\n",1000*sec[2]);
    fprintf(fp,"          FAVORABLE_MUTN   %lf\n",1000*sec[3]);
    fprintf(fp,"          OFFSPRING        %lf\n",1000*sec[4]);
    fprintf(fp,"          SELECTION        %lf\n\n",1000*sec[5]);
    
    fprintf(fp,"          AUXILIARY ROUTINES (MILLISEC/GEN/PROC):\n\n");
    
    fprintf(fp,"          MUTN STATISTICS  	%lf\n",1000*sec[6]);
    fprintf(fp,"          ALLELE STATISTICS	%lf\n",1000*sec[7]);
    fprintf(fp,"          READ_RESTART_DUMP	%lf\n",1000*sec[8]);
    fprintf(fp,"          WRITE_OUTPUT_DUMP	%lf\n",1000*sec[9]);
    fprintf(fp,"          WRITE_SAMPLE     	%lf\n",1000*sec[10]);
    fprintf(fp,"          READ_PARAMETERS  	%lf\n",1000*sec[11]);
    fprintf(fp,"          WRITE_PARAMETERS 	%lf\n",1000*sec[12]);
  }
}
