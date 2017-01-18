#include <stdio.h>
#define SIZE 10

void quick_sort(double []);
void quick_sort_engine(double [], int, int);
void myheapsort(double a[], int n);

/*
int main()
{
	int i=0;
	double items[SIZE] = {0};

	for(i=0; i< SIZE; i++){
		items[i] = (double) random()/1e8;
		printf("%d. %g\n", i, items[i]);
	}
        //quick_sort(items);
        quick_sort_engine(items, 0, sizeof(items) / sizeof(items[0]));
	
	for(i=0; i< SIZE; i++){
		printf("%d. %g\n", i, items[i]);
	}
	return 0;
}
*/

void quick_sort(double items[]) { 
        int size;
        size = sizeof(items) / sizeof(items[0]);
	quick_sort_engine(items, 0, size-1);
}

void quick_sort_engine(double *items, int left, int right){
	int i, j;
	double x, y;
        extern char version_qsort[80];

        sprintf(version_qsort,
          "$Id: sort.c,v 1.3 2009/07/23 19:00:09 wes Exp $");

	i = left;
	j = right;
	x = items[(left+right)/2];

	do{
		while((items[i] < x) && (i < right))
			i++;
		while((x < items[j]) && (j > left))
			j--;
		
		if(i <= j){
			y = items[i];
			items[i] = items[j];
			items[j] = y;
			i++;
			j--;
		}
	}while(i <= j);

	if(i < right)
		quick_sort_engine(items, i, right);
	if(left < j)
		quick_sort_engine(items, left, j);
}

void myheapsort(double *ra, int n) {
  /* This routine applies the heapsort algorithm to sort array a of  *
   * length n into ascending numerical order. Array a is replaced on *
   * output by its sorted rearrangement.                             *
   * Taken from Press, et al., Numerical Recipes, 1986, p. 231.      */
  int i, j;
  int ir = n;
  int l = (n >> 1) + 1;
  double rra;
  
  for (;;) {
    if (l > 1) {         // Still in hiring phase.
      rra = ra[--l];
    } else {             // In retirement-and-promotion phase. 
      rra = ra[ir];    // Clear a space at the end of the array.
      ra[ir] = ra[1];  // Retire the top of the heap into it.
      if (--ir == 1) { // Decrease the size of the corporation.
	ra[1] = rra; // Completed the last promotion.
	return;      // Identify the least competent worker.
      }
    }
    i = l;       // Set up to sift element ra to its proper level
    j = l << 1;
    while (j <= ir) {      // Compare to better underling.
      if (j < ir && ra[j] < ra[j+1]) { ++j; }
      if (rra < ra[j]) { // Demote ra.
	ra[i] = ra[j];
	j += (i = j);
      } else {           // This is ra's level. Set j to terminate
	j = ir + 1;    // sift-down.
      }
    }
    ra[i] = rra;
  }
}

