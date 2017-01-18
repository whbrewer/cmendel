#include "mendel.h"

#ifdef MPICH
#include "/usr/local/include/mpi.h"
#endif

int mpi_myinit(int argc, char **argv)
{
  int ownid;
#ifdef MPICH
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ownid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif
  return ownid;
}

void mpi_myabort()
{
#ifdef MPICH
  MPI_Abort(MPI_COMM_WORLD, ierr);
#endif
}

void mpi_mybcast(double x, int n)
{
#ifdef MPICH
  MPI_Bcast(&x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
#endif
}

void mpi_myfinalize()
{
#ifdef MPICH
  MPI_Finalize();
#endif
}

void mpi_davg(double *x, double *xavg, int n)
{
#ifdef MPICH
  MPI_Reduce(x, xavg, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
  *xavg = *xavg/num_procs;
#endif
}

void mpi_dsum(double *x, double *xsum, int n)
{
#ifdef MPICH
  MPI_Reduce(x, xsum, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

void mpi_ravg(float *x, float *xavg, int n)
{
#ifdef MPICH
  MPI_Reduce(x, xavg, n, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  *xavg = *xavg/num_procs;
#endif
}

void mpi_isum(int *x, int *xavg, int n)
{
#ifdef MPICH
  MPI_Reduce(x, xavg, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

void mpi_migration(int ***dmutn, int ***fmutn, double ***linkage_block_fitness, 
                   int ****lb_mutn_count, int gen, int ownid, int msg_num, 
                   int ierr, bool *available)
{
#ifdef MPICH
  int dbuffs[max_del_mutn_per_indiv*num_indiv_exchanged];
  int dbuffr[max_del_mutn_per_indiv*num_indiv_exchanged];
  int fbuffs[max_fav_mutn_per_indiv*num_indiv_exchanged];
  int fbuffr[max_fav_mutn_per_indiv*num_indiv_exchanged];
  int i, j, k, l, m, n, num_receiving_tribes=0, nie=0;
  int src, dest;
  MPI_Request send_requests[5], recv_requests[5];
  MPI_Status recv_status[5];
  MPI_Status send_status[5];
  int sender;
  int max_num_dmutn, max_num_fmutn;
  int max_num_dmutn_recvd = 0, max_num_fmutn_recvd = 0;
  double lbuffs[2*num_linkage_subunits*num_indiv_exchanged];
  double lbuffr[2*num_linkage_subunits*num_indiv_exchanged];
  double t0, t1, time;
  int cbuff1s[2*num_linkage_subunits*num_indiv_exchanged];
  int cbuff1r[2*num_linkage_subunits*num_indiv_exchanged];
  int cbuff2s[2*num_linkage_subunits*num_indiv_exchanged];
  int cbuff2r[2*num_linkage_subunits*num_indiv_exchanged];
  int msg_size_dbuff, msg_size_fbuff, msg_size_lbuff;
  int msg_size_rdbuff,msg_size_rfbuff;
  int swap_list[num_indiv_exchanged*num_procs];
  int id;
  int p,q;
  float x;
  bool debug, timing;
  FILE *fp; 
  char filename[80];
  debug = false;
  timing = false;
  
  sprintf(version_mpi,
          "$Id: mpi.c,v 1.13 2009/11/06 19:20:58 wes Exp $");

  if (debug) {
     sprintf(filename,"%s.%03d.mpi",case_id,ownid);
     fp = fopen(filename,"w");
  }
  
  // swap individuals between processors
  
  for(i = 0; i < pop_size; i++) 
    available[i] = true;
  
  if (migration_model == 1) {
    num_receiving_tribes = 1;
    nie = num_indiv_exchanged;
  } else if (migration_model == 2) {
    num_receiving_tribes = 2;
    nie = num_indiv_exchanged*2;
  } else if (migration_model == 3) {
    num_receiving_tribes = num_procs - 1;
    nie = num_indiv_exchanged*(num_procs - 1);
  } else {
    printf("migration_model %d not supported\n", migration_model);
  } 
  
  //randomly select individuals to be swapped
  for(i = 0 ; i < nie; i++ )
    {
      x = (int)(current_pop_size*randomnum());
      j = min(current_pop_size,x);
      while( !available[j] ) {
	j = (j+1)%current_pop_size;
      }
      available[j] = false;
      swap_list[i] = j;
    }

  if(tracking_threshold < 1){
    //Find max number of mutations in dmutn and fmutn
    max_num_dmutn = 0;
    max_num_fmutn = 0;
    
    for(i = 0 ; i < nie; i++ ){
      max_num_dmutn = max(max_num_dmutn,max(dmutn[swap_list[i]][0][0],
                                            dmutn[swap_list[i]][1][0]));
      max_num_fmutn = max(max_num_fmutn,max(fmutn[swap_list[i]][0][0],
                                            fmutn[swap_list[i]][1][0]));
    }

    max_num_dmutn += 2;
    max_num_fmutn += 2;

    if(max_num_dmutn > max_del_mutn_per_indiv/2){
      fprintf(stderr, "ERROR: need to increase max_del_mutn_per_indiv\n");
      fprintf(stderr,"max_num_dmutn, max_del_mutn_per_indiv is %d %d\n",
              max_num_dmutn,max_del_mutn_per_indiv);
      mpi_myabort();
    }
  }
  
  //populate buffer arrays for message passing
  for(m=0; m < num_receiving_tribes; m++){
    k = 0;
    n = 0; 
    l = 0;
    for(i=0; i < num_indiv_exchanged; i++){
      id = m*num_indiv_exchanged + i;
      //printf("gen, m, i, id, sid: %d, %d, %d, %d, %d\n" ,gen, m,i,id,swap_list[id]);
      if(tracking_threshold < 1){
	for(p=0; p < max_num_dmutn; p++)
	  dbuffs[k+p] = dmutn[swap_list[id]][0][p];
	k += max_num_dmutn;

	for(p=0; p < max_num_dmutn; p++)
	  dbuffs[k+p] = dmutn[swap_list[id]][1][p];
	k += max_num_dmutn;

	for(p=0; p < max_num_fmutn; p++)
	  fbuffs[n+p] = fmutn[swap_list[id]][0][p];
	n += max_num_fmutn;

	for(p=0; p < max_num_fmutn; p++)
	  fbuffs[n+p] = fmutn[swap_list[id]][1][p];
	n += max_num_fmutn;
      }
      for(p=0; p < num_linkage_subunits; p++){
	lbuffs[l+p] = linkage_block_fitness[swap_list[id]][0][p];
	cbuff1s[l+p] = lb_mutn_count[swap_list[id]][0][0][p];
	cbuff2s[l+p] = lb_mutn_count[swap_list[id]][1][0][p];
      }
      l += num_linkage_subunits;

      for(p=0; p < num_linkage_subunits; p++){
	lbuffs[l+p] = linkage_block_fitness[swap_list[id]][1][p];
	cbuff1s[l+p] = lb_mutn_count[swap_list[id]][0][1][p];
	cbuff2s[l+p] = lb_mutn_count[swap_list[id]][1][1][p];
      }
      l += num_linkage_subunits;
    } 

    msg_size_dbuff = k; 
    msg_size_fbuff = n;
    msg_size_lbuff = l;
    
    if ( debug && ownid==0 ){
      fprintf(fp,"------ GENERATION: %d ----------------------\n", gen);
      fprintf(fp,"swap_list for proc:%d\n",ownid);
      for (i = 0; i < num_indiv_exchanged; i++)
          fprintf(fp,"%d ",swap_list[i]);
      //fprintf(fp,"\ndbuff sent: %d\n",msg_size_dbuff);
      //for (i = 1; i <= msg_size_dbuff; i++)
      //    fprintf(fp,"%d ",dbuffs[i]);
      for(j=0; j < num_indiv_exchanged; j++){ 
        fprintf(fp,"\nindividual #%d\n",swap_list[j]);
        fprintf(fp,"\ndmutn sent : %d\n",max_num_dmutn);
        for (i = 0; i < max_num_dmutn; i++)
            fprintf(fp,"%d ",dmutn[swap_list[j]][0][i]);
        fprintf(fp,"\nfmutn sent : %d\n",max_num_fmutn);
        for (i = 0; i < max_num_fmutn; i++)
            fprintf(fp,"%d ",fmutn[swap_list[j]][0][i]);
        fprintf(fp,"\nlinkage_block_fitness sent - hap1:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%lf ",linkage_block_fitness[swap_list[j]][0][i]);
        fprintf(fp,"\nlinkage_block_fitness sent - hap2:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%lf ",linkage_block_fitness[swap_list[j]][1][i]);
        fprintf(fp,"\nlb_mutn_count sent - hap1:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%d ",lb_mutn_count[swap_list[j]][0][0][i]);
        fprintf(fp,"\nlb_mutn_count - hap2:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%d ",lb_mutn_count[swap_list[j]][0][1][i]);
        fprintf(fp,"\n");
      }
    }

    t0 = MPI_Wtime();
    
    //compute destination tribe based on migration model
    if(migration_model == 1 || migration_model == 3)
      dest = (ownid+m+1)%num_procs;
    else if (migration_model == 2){
      if ((m+1)%2) dest = (ownid+m+1) % num_procs;
      else dest = (ownid+num_procs-1) % num_procs;
    } else {
      fprintf(stderr,"ERROR: migration_model not supported\n");
      exit(EXIT_FAILURE);
    }
  
    // Intead of using MPI_ANY_SOURCE would be better to use src in future
    //src = myid - 1;
    //if(src < 0) src = num_procs - 1;

    if(tracking_threshold < 1){
      
      //send the maximum number of deleterious mutations
      MPI_Irecv(&max_num_dmutn_recvd,1,MPI_INT,MPI_ANY_SOURCE,msg_num,
                MPI_COMM_WORLD,recv_requests);
      MPI_Isend(&max_num_dmutn,1,MPI_INT,dest,msg_num,MPI_COMM_WORLD,send_requests);
      //complete the nonblocking receives before processing
      MPI_Waitany(1,recv_requests,&sender,recv_status);
      //complete the nonblocking sends before proceeding
      MPI_Waitall(1,send_requests,send_status);
      msg_num++;
 
      //communicate buffer of deleterious mutations
      msg_size_rdbuff = 2*max_num_dmutn_recvd*num_indiv_exchanged;
      MPI_Irecv(&dbuffr,msg_size_rdbuff,MPI_INT,MPI_ANY_SOURCE,msg_num,
                MPI_COMM_WORLD,recv_requests);
      MPI_Isend(&dbuffs,msg_size_dbuff,MPI_INT,dest,msg_num,
                MPI_COMM_WORLD,send_requests);
      MPI_Waitany(1,recv_requests,&sender,recv_status);
      MPI_Waitall(1,send_requests,send_status);
      msg_num++;

      //send the size of the fbuff message buffer
      msg_size_rfbuff = 2*max_num_fmutn_recvd*num_indiv_exchanged;
      MPI_Irecv(&msg_size_rfbuff,1,MPI_INT,MPI_ANY_SOURCE,msg_num,
                MPI_COMM_WORLD,recv_requests);
      MPI_Isend(&msg_size_fbuff,1,MPI_INT,dest,msg_num,
                MPI_COMM_WORLD,send_requests);
      MPI_Waitany(1,recv_requests,&sender,recv_status);
      MPI_Waitall(1,send_requests,send_status);
      msg_num++;
      
      //communicate buffer of favorable mutations
      MPI_Irecv(&fbuffr,msg_size_rfbuff,MPI_INT,MPI_ANY_SOURCE,msg_num,
                MPI_COMM_WORLD,recv_requests);
      MPI_Isend(&fbuffs,msg_size_fbuff,MPI_INT,dest,msg_num,
                MPI_COMM_WORLD,send_requests);
      MPI_Waitany(1,recv_requests,&sender,recv_status);
      MPI_Waitall(1,send_requests,send_status);
      msg_num++;
    }
    
    //communicate buffer of linkage_block_fitness 
    MPI_Irecv(&lbuffr,msg_size_lbuff,MPI_DOUBLE_PRECISION,
              MPI_ANY_SOURCE,msg_num,MPI_COMM_WORLD,recv_requests);
    MPI_Isend(&lbuffs,msg_size_lbuff,MPI_DOUBLE_PRECISION,dest,msg_num,
              MPI_COMM_WORLD,send_requests);
    MPI_Waitany(1,recv_requests,&sender,recv_status);
    MPI_Waitall(1,send_requests,send_status);
    msg_num++;

    //communicate buffer of lb_mutn_count
    MPI_Irecv(&cbuff1r,msg_size_lbuff,MPI_INT,MPI_ANY_SOURCE,msg_num,
              MPI_COMM_WORLD,recv_requests);
    MPI_Isend(&cbuff1s,msg_size_lbuff,MPI_INT,dest,msg_num,
              MPI_COMM_WORLD,send_requests);
    MPI_Waitany(1,recv_requests,&sender,recv_status);
    MPI_Waitall(1,send_requests,send_status);
    msg_num++;
    
    //communicate buffer of lb_mutn_count
    MPI_Irecv(&cbuff2r,msg_size_lbuff,MPI_INT,MPI_ANY_SOURCE,msg_num,
              MPI_COMM_WORLD,recv_requests);
    MPI_Isend(&cbuff2s,msg_size_lbuff,MPI_INT,dest,msg_num,
              MPI_COMM_WORLD,send_requests);
    MPI_Waitany(1,recv_requests,&sender,recv_status); 
    MPI_Waitall(1,send_requests,send_status);
    msg_num++;

    t1 = MPI_Wtime();
    time = t1 - t0; 
    if(gen%10==0 && timing)
      printf("communication time = %d seconds \n",(int)time);
    
    //Unpack transferred buffers
    k = 0;
    n = 0;
    l = 0; 
    
    if (debug) 
      if(max_num_dmutn == max_num_dmutn_recvd)
	printf("WARNING: max_num_dmutn same for 2 tribes %d\n", max_num_dmutn);
    
    for(i=0; i < num_indiv_exchanged; i++){ 
      id = m*num_indiv_exchanged + i; 
      if(tracking_threshold < 1){ 

	for(q=0; q < 2; q++){
	   for(p=0; p < max_del_mutn_per_indiv/2; p++)
	     dmutn[swap_list[id]][q][p] = num_linkage_subunits*lb_modulo+1;
	   for(p=0; p < max_fav_mutn_per_indiv/2; p++)
	     fmutn[swap_list[id]][q][p] = num_linkage_subunits*lb_modulo+1;
	}

	for(p=0; p < max_num_dmutn_recvd; p++) 
           dmutn[swap_list[id]][0][p] = dbuffr[k+p];
	k += max_num_dmutn_recvd;
	
	for(p=0; p < max_num_dmutn_recvd; p++) 
           dmutn[swap_list[id]][1][p] = dbuffr[k+p];
	k += max_num_dmutn_recvd;
	
	for(p=0; p < max_num_fmutn_recvd; p++) 
           fmutn[swap_list[id]][0][p] = fbuffr[n+p];
	n += max_num_fmutn_recvd;
	
	for(p=0; p < max_num_fmutn_recvd; p++) 
           fmutn[swap_list[id]][1][p] = fbuffr[n+p];
	n += max_num_fmutn_recvd; 
      }

      for(p=0; p < num_linkage_subunits; p++){
	linkage_block_fitness[swap_list[id]][0][p] = lbuffr[l+p];
	lb_mutn_count[swap_list[id]][0][0][p] = cbuff1r[l+p];
	lb_mutn_count[swap_list[id]][1][0][p] = cbuff2r[l+p];
      } 
      l += num_linkage_subunits;
      
      for(p=0; p < num_linkage_subunits; p++){
	linkage_block_fitness[swap_list[id]][1][p] = lbuffr[l+p];
	lb_mutn_count[swap_list[id]][0][1][p] = cbuff1r[l+p];
	lb_mutn_count[swap_list[id]][1][1][p] = cbuff2r[l+p];
      } 
      l += num_linkage_subunits;
    }

    if (debug && ownid == 1){
      fprintf(fp,"------ GENERATION: %d ----------------------\n", gen);
      fprintf(fp,"swap_list for proc:%d\n",ownid);
      for (i = 0; i < num_indiv_exchanged; i++)
          fprintf(fp,"%d ",swap_list[i]);
      //fprintf(fp,"\ndbuff received: %d\n",msg_size_dbuff);
      //for (i = 1; i <= msg_size_dbuff; i++)
      //    fprintf(fp,"%d ",dbuffr[i]);
      for(j=0; j < num_indiv_exchanged; j++){ 
        fprintf(fp,"\nindividual #%d\n",swap_list[j]);
        fprintf(fp,"\ndmutn received: %d\n",max_num_dmutn_recvd);
        for (i = 0; i < max_num_dmutn_recvd; i++)
            fprintf(fp,"%d ",dmutn[swap_list[j]][0][i]);
        fprintf(fp,"\nfmutn received: %d\n",max_num_fmutn_recvd);
        for (i = 0; i < max_num_fmutn_recvd; i++)
            fprintf(fp,"%d ",fmutn[swap_list[j]][0][i]);
        fprintf(fp,"\nlinkage_block_fitness received - hap1:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%lf ", linkage_block_fitness[swap_list[j]][0][i]);
        fprintf(fp,"\nlinkage_block_fitness received - hap2:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%lf ", linkage_block_fitness[swap_list[j]][1][i]);
        fprintf(fp,"\nlb_mutn_count received - hap1:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%d ", lb_mutn_count[swap_list[j]][0][0][i]);
        fprintf(fp,"\nlb_mutn_count received - hap2:\n");
        for (i = 0; i < 9; i++)
            fprintf(fp,"%d ", lb_mutn_count[swap_list[j]][0][1][i]);
        fprintf(fp,"\n");
      }
    }
  }

  if(debug) fclose(fp);
#endif
}

