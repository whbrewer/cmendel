#include <stdlib.h>
#include <stdio.h>

/* experimental function to remove dependency on type */
void malloc3d(int size,int m, int n, int o)
{
  void ***p = NULL;
  int i,j;

  p = (void ***)malloc(size*m);
  if (p != NULL){
    for (i = 0;i < m;i++){
      p[i] = (void **)malloc(size*n);
      if (p[i] != NULL){ 
        for (j = 0;j < n;j++) p[i][j] = (void *)malloc(size*o);
      }else{
        printf("ERROR: not enough memory\n");
        exit(EXIT_FAILURE);
      }
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
//  return p;
}

/* Allocate memory for three dimensional int array and 
 * return the pointer which point to the array */
int*** malloc_int3d(int m,int n,int o)
{
  int ***p = NULL;
  int i,j,k;
  extern char version_mem[80];

  sprintf(version_mem,
          "$Id: mem.c,v 1.10 2008/07/11 07:38:20 wes Exp $");
  
  p = (int ***)malloc(sizeof(int **)*m);
  if (p != NULL){
    for (i = 0;i < m;i++){
      p[i] = (int **)malloc(sizeof(int *)*n);
      if (p[i] != NULL){ 
        for (j = 0;j < n;j++){
          p[i][j] = (int *)malloc(sizeof(int)*o);
          if (p[i][j] != NULL){
            for (k = 0;k < o;k++) 
              p[i][j][k] = 0;
          }else{
            printf("ERROR: not enough memory\n");
            exit(EXIT_FAILURE);
          }
        }
      }else{
        printf("ERROR: not enough memory\n");
        exit(EXIT_FAILURE);
      }
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
  return p;
}

/* Allocate memory for four dimensional int array and 
 * return the pointer which point to the array */
int**** malloc_int4d(int m,int n,int o,int p)
{
  int ****q = NULL;
  int i,j,k,l;
  
  q = (int ****)malloc(sizeof(int ***)*m);
  if (q != NULL){
    for (i = 0;i < m;i++){
      q[i] = (int ***)malloc(sizeof(int **)*n);
      if (q[i] != NULL){ 
        for (j = 0;j < n;j++){
          q[i][j] = (int **)malloc(sizeof(int *)*o);
          if (q[i][j] != NULL){
            for (k = 0;k < o;k++){
              q[i][j][k] = (int *)malloc(sizeof(int)*p);
              if (q[i][j][k] != NULL){
                for (l = 0;l < p;l++) 
                  q[i][j][k][l] = 0;
              }else{
                printf("ERROR: not enough memory\n");
                exit(EXIT_FAILURE);
              }
            }
          }else{
            printf("ERROR: not enough memory\n");
            exit(EXIT_FAILURE);
          }
        }
      }else{
        printf("ERROR: not enough memory\n");
        exit(EXIT_FAILURE);
      }
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  return q;
}

/* Allocate memory for two dimensional double array and *
 * return a pointer which point to the array            */
double** malloc_double2d(int m,int n)
{
  double **a = NULL;
  int i,j;
  a = (double **)malloc(sizeof(double *)*m);
  if (a != NULL){
    for (i = 0; i < m;i++){
      a[i] = (double *)malloc(sizeof(double)*n);
      if (a[i]!=NULL){
        for (j = 0;j < n;j++)
          a[i][j] = 0; 
      }else{
        printf("ERROR: not enough memory\n");
        exit(EXIT_FAILURE);
      }
    } 
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  } 
  
  return a;
}

/* Allocate memory for three dimensional double array and * 
 * return a pointer which point to the array              */
double*** malloc_double3d(int m,int n,int o)
{
  double ***p = NULL;
  int i,j,k;
  
  p = (double ***)malloc(sizeof(double **)*m);
  if (p != NULL){
    for (i = 0;i < m;i++){
      p[i] = (double **)malloc(sizeof(double *)*n);
      if (p[i] != NULL){ 
        for (j = 0;j < n;j++){
          p[i][j] = (double *)malloc(sizeof(double)*o);
          if (p[i][j] != NULL){
            for (k = 0;k < o;k++)
              p[i][j][k] = 0; 
          }else{
            printf("ERROR: not enough memory\n");
            exit(EXIT_FAILURE);
          }
        }
      }else{
        printf("ERROR: not enough memory\n");
        exit(EXIT_FAILURE);
      }
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
  return p;
}

/* Allocate memory for one dimensional float array and *
 * return a pointer which point to the array           */
float * malloc_float1d(int m)
{
  
  float * p = NULL;
  int i;
  
  p = (float *)malloc(sizeof(float)*m);
  if (p != NULL){
    for (i = 0; i < m; i++){
      p[i] = 0;
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
  return p;
} 

/* Allocate memory for one dimensional double array and * 
 * return a pointer which point to the array            */
double * malloc_double1d(int m)
{
  double * p = NULL;
  int i;
  
  p = (double *)malloc(sizeof(double)*m);
  if (p != NULL){
    for (i = 0; i < m; i++){
      p[i] = 0;
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
  return p;
}

/* Allocate memory for one dimensional int array and * 
 * return a pointer which point to the array         */
int * malloc_int1d(int m)
{
  int * p = NULL;
  int i;
  
  p = (int *)malloc(sizeof(int)*m);
  if (p != NULL){
    for (i = 0; i < m; i++){
      p[i] = 0;
    }
  }else{
    printf("ERROR: not enough memory\n");
    exit(EXIT_FAILURE);
  }
  
  return p;
}

#if 0

/* This is example for how to allocate memory using           *
 * the individual typedef struct... not currently implemented */

int** malloc_int2d(int i_size,int j_size)
  
  individual *id;
  int i = 0;

  if (( id = (individual *) malloc( sizeof(individual) ) ) == NULL) {
     printf("ERROR ALLOCATING individual\n");
     exit;
  }

  id->dmutn = (int **) malloc( sizeof(int) * i_size);

  for(i=0; i < i_size; i++) 
     id->dmutn[i] = (int *) malloc (sizeof(int) * j_size);

  id->dmutn[0][0] = 5;
  printf("%d\n", id->dmutn[0][0]);
  return 0;

}
#endif
