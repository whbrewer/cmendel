/* Set following options before compiling */
#include "config.h"
 
//#define MPICH
//#define PTHREADS

#define UNIX
//#define WINDOWS

/* Uncomment PPM line to output large file that 
 * provides graphical depiction of simulation representing
 * each individual as pixel with color white representing fitness=1 
 * and color black representing fitness=0 and various shades
 * of gray in between */
//#define PPM

//#define DEBUG

/*  Uncomment WRITE_SAMPLE line to print an output file 
 *  containing details concerning the mutations carried by 
 *  five members of the total population in a format that 
 *  is (hopefully) readily understandable.             */
//#define WRITE_SAMPLE

#define DYNAMIC_POP_SIZE

/* Should not need to change anything after this line */

typedef int bool;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ranlib.h"
#include <sys/time.h>

/* Macro functions */
#define random_poisson ignpoi
#define randomnum() (float)(rand()/(float)RAND_MAX)
//#define randomnum() ranf()
#define random_normal gennor

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

//#ifdef UNIX
#ifdef HAVE_SYS_TIME_H
   #define TIME(a) gettimeofday(a,NULL)
   #define DIFFTIME(a,b) difftime(a.tv_sec,b.tv_sec) + difftime(a.tv_usec,b.tv_usec)/1e6
   #define TIME_T struct timeval 
#else
   #define TIME(a) time(a)
   #define DIFFTIME(a,b) difftime(a,b)
   #define TIME_T time_t
#endif

#define true  1
#define false 0

#ifdef WINDOWS
   #define fabsf(x)   ((float)fabs((double)(x)))
   #define logf(x)    ((float)log((double)(x)))
   #define expf(x)    ((float)exp((double)(x)))
#endif

#define mod(a, b) a % b

#if 0
typedef struct {
   int num_dmutn; 
   int num_fmutn;
   int *dmutn;
   int *fmutn;
   double *linkage_block_fitness;
   int **lb_mutn_count;   
} diploid;

typedef struct {
   diploid hap[2];
} individual;
//individual *id;
#endif

/* Function prototype declarations */

void read_parameters(FILE *);//ok
void write_parameters(FILE *);//ok
void favorable_mutn(int ***, int ****, double ***);//ok
void initialize(char *);

void gen_initial_contrasting_alleles(int ***, int ***, double ***, 
                                     float *, int );

void mpi_migration(int ***dmutn, int ***fmutn, double ***linkage_block_fitness,
                   int ****lb_mutn_count, int gen, int ownid, int msg_num, 
                   int ierr, bool *available);

void selection(int ***dmutn, int ***fmutn, int ****lb_mutn_count,
               double ***linkage_block_fitness, double *fitness,
               double *pheno_fitness,
               double *work_fitness, double *sorted_score,
               float *initial_allele_effects, int max_size,
               int total_offspring, int gen);

void offspring(int*** dmutn_offsprng, int*** fmutn_offsprng,
               int**** offsprng_lb_mutn_count, double*** offsprng_lb_fitness,
               int*** dmutn, int*** fmutn, int**** lb_mutn_count,
               double*** linkage_block_fitness, int dad, int mom,
               int child);

void profile(FILE *);

void mpi_ravg(float *,float *,int);
void mpi_myfinalize();
void mpi_myabort(void);
void mpi_davg(double*, double*, int);
void mpi_mybcast(double, int);
void mpi_isum(int *,int *,int);
int  mpi_myinit(int, char **);

int**     malloc_int2d(int m,int n);
int***    malloc_int3d(int m,int n,int o);
int****   malloc_int4d(int m,int n,int o,int p);
double**  malloc_double2d(int m,int n);
double*** malloc_double3d(int m,int n,int o);
float*    malloc_float1d(int m);
double*   malloc_double1d(int m);
int*      malloc_int1d(int m);

void quick_sort(double []);
void quick_sort_engine(double [], int, int);
void myheapsort(double a[], int n);

int poisson(float);

void read_restart_dump(int ***dmutn, int ***fmutn, int ****lb_mutn_count, 
                       double ***linkage_block_fitness, 
                       float *initial_allele_effects, int *generation_number, 
                       int max_size, char myid_str[]);

void diagnostics_contrasting_alleles(int ***dmutn, int ***fmutn, 
                                     int ****count, double *cum_effect, 
                                     float *initial_allele_effects, 
                                     int *ica_count, int max_size, bool list);

void diagnostics_history_plot1(int ***dmutn, int ***fmutn, 
                               int ****lb_mutn_count, int tracked_fav_mutn, 
                               int ica_count[], int gen, bool print_flag);

void diagnostics_mutn_bins_plot2(int ***dmutn, int ***fmutn, 
                                 double *accum, int gen);

void diagnostics_mutn_bins_plot4(int ***dmutn, int ***fmutn, 
                                 double ***linkage_block_fitness, 
                                 int ****lb_mutn_count, int gen); 

void diagnostics_polymorphisms_plot7(int *** dmutn, int *** fmutn,
                                     int max_size, int gen);

void diagnostics_selection(double *fitness_pre_sel, double *fitness_post_sel,
                           int total_offspring, int gen);

void write_status(FILE *unit, int gen, int current_pop_size,
                  float frac_recessive, double total_del_mutn,
                  double tracked_del_mutn, double total_fav_mutn,
                  double pre_sel_fitness, double pre_sel_geno_sd,
                  double pre_sel_pheno_sd, double pre_sel_corr,
                  double post_sel_fitness, double post_sel_geno_sd,
                  double post_sel_pheno_sd, double post_sel_corr);

void write_output_dump(int ***dmutn, int ***fmutn, int ****lb_mutn_count, 
                       double ***linkage_block_fitness, 
                       float *initial_allele_effects, 
                       int generation_number,char myid_str[]);

#ifdef WRITE_SAMPLE
void write_sample(int ***dmutn, int ***fmutn, int ****lb_mutn_count,
                  double ***linkage_block_fitness, double *fitness,
                  float ***defect, float ***improve, float *effect,
                  int generation_number);
#endif

//for debug
void showInt_3dArray(int d1,int d2,int d3,int *** arr);

//for debug
void showInt_4dArray(int d1,int d2,int d3,int d4,int **** arr);

FILE *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8,
     *fp9, *fp11, *fp13, *fp14, *fp16, *fp17, *fp18, 
     *fp19, *fp21, *fp22, *fp23, *fp24, *fp25, *fp26, 
     *fp34, *fp35;

int fh24, fh26;

#ifdef WRITE_STATUS
FILE *fp12;
#endif

#ifdef PPM
FILE *fp15;
#endif

/* Global variable declarations */

char   version_mendel[80], version_offspring[80], version_selection[80], 
       version_mpi[80], version_diagnostics[80], version_fileio[80],
       version_init[80], version_qsort[80], version_ranlib[80], 
       version_mem[80];

int    pop_size, num_generations, num_linkage_subunits, 
       bottleneck_generation, bottleneck_pop_size,
       num_bottleneck_generations, fitness_distrib_type,
       max_tracked_mutn_per_indiv, new_mutn_count, 
       max_del_mutn_per_indiv, max_fav_mutn_per_indiv,
       num_initial_fav_mutn, num_indiv_exchanged,
       random_number_seed, restart_dump_number, 
       dump_number, lb_modulo, current_pop_size,
       mutn_per_indiv, haploid_chromosome_number,
       selection_scheme, migration_generations,
       migration_model, num_contrasting_alleles,
       myid, num_procs, ierr, msg_num,
       par_mutn_per_indiv, pop_growth_model, 
       fav_fixed, fav_lost;

float  gamma_del, gamma_fav; 

float  offspring_per_female, new_mutn_per_offspring, 
       haploid_genome_size, high_impact_mutn_fraction,
       high_impact_mutn_threshold, fraction_recessive,
       dominant_hetero_expression, max_fav_fitness_gain,
       recessive_hetero_expression, frac_fav_mutn, 
       heritability, uniform_fitness_effect, 
       multiplicative_weighting, tracking_threshold,
       initial_alleles_mean_effect, non_scaling_noise,
       del_scale, fav_scale, 
       fraction_random_death, fraction_self_fertilization,
       linked_mutn_se_fraction, poisson_mean,
       pop_growth_rate, ica_mean_effect, fav_mean_freq,
       partial_truncation_value, se_scaling_factor;

double sec[14];

double alpha_del, alpha_fav, fitness_cutoff, synergistic_factor,
       pre_sel_fitness, post_sel_fitness,  
       pre_sel_geno_sd, pre_sel_pheno_sd, pre_sel_corr,
       post_sel_geno_sd, post_sel_pheno_sd, post_sel_corr;

bool   fitness_dependent_fertility, dynamic_linkage,
       synergistic_epistasis, is_parallel, bottleneck_yes, 
       restart_case, write_dump, homogenous_tribes,
       clonal_reproduction, clonal_haploid;

char case_id[6+1], data_file_path[80+1];

int *available, *replaced_by_offspring;

