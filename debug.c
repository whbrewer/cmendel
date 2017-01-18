#include "mendel.h"

void showInt_3dArray(int d1,int d2,int d3,int *** arr){
  int i1=0;
  int i2=0;
  int i3=0;
  printf("d1 = %d, d2 = %d, d3 = %d\n",d1,d2,d3);
  for (i1=0;i1 < d1;i1++)
    for (i2=0;i2 < d2;i2++)
      for (i3=0;i3 < d3;i3++)
        printf("array[%d][%d][%d] == %d\n",i1,i2,i3,arr[i1][i2][i3]);
}

void showInt_4dArray(int d1,int d2,int d3,int d4,int **** arr){
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  int i4 = 0;
  for (;i1 < d1;i1++)
    for (;i2 < d2;i2++)
      for (;i3 < d3;i3++)
        for (;i4 < d4;i4++)
          printf("array[%d][%d][%d][%d] == %d\n",i1,i2,i3,i4,arr[i1][i2][i3][i4]);
}



