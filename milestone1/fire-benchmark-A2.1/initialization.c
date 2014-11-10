/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "util_read_files.h"
#include "initialization.h"
#include "mmalloc.h"
#include <mpi.h>

#include "metis.h"

#define MIN(a,b) ((a)<(b)?(a):(b))



void call_metis(int nintci, int nintcf, int points_count, int* elems, int nprocs, int my_rank, idx_t** epart){
	int i;
		idx_t nelem = (idx_t) nintcf+1;
		idx_t nnodes = (idx_t) points_count;


		// how to structure them??

		idx_t *eptr = malloc( sizeof(idx_t) * ( nelem + 1 ) );
		idx_t *eind = malloc( sizeof(idx_t) * ( nelem  ) * 8); // every cell: six neighbours


		/* every elem HAS 8 POINTS !*/
		for (i=0; i< nelem + 1; i++){
			eptr[i] = i*8;
		}

		for (i=0; i< (nelem )*8; i++){
			eind[i] = (int)elems[i];
		}



		idx_t nparts = (idx_t) nprocs;
		idx_t ncommon = 4;
		idx_t objval;
		*epart = malloc( sizeof(idx_t) * nelem );
		idx_t *npart = malloc( sizeof(idx_t) * nnodes );



		if( METIS_PartMeshDual( &nelem, &nnodes,
			 eptr,  eind,
			NULL , NULL, &ncommon,
			&nparts, NULL , NULL, &objval,
			*epart, npart) != METIS_OK )
		/*if( METIS_PartMeshNodal( &nelem, &nnodes,
					 eptr,  eind,
					NULL , NULL,
					&nparts, NULL , NULL, &objval,
					*epart, npart) != METIS_OK )*/
		{
			printf("ERROR IN METIS\n");
		}
		else
		{
			printf("METIS SUCCESSFFUL\n");
		}

}


void get_points(int nprocs, int rank, int quotient, int remainder, int* points_count,
        int*** points, int** elems, int g_points_count, int* g_elems, int** g_points) {

			int i,j,k;

			int nintci = 0;
			int nintcf = quotient + ( rank < remainder ? 1 : 0 ) - 1;


			int local_ix=0, global_ix;
			int * points_to_local = (int*) malloc(g_points_count*sizeof(int));



			for(i=0;i<g_points_count;i++){
				points_to_local[i]=-1;
			}



			int base = rank*quotient+MIN(remainder, rank);
			for(i=nintci; i<nintcf+1; i++){
				global_ix = i+base;


				for(j=0; j<8; j++){
					//global_ix = i*8+j +myrank*quotient+MIN(remainder, myrank);

					if( points_to_local[  g_elems[ global_ix*8 + j ]   ] == -1 ) {
						points_to_local[  g_elems[ global_ix*8 + j ]   ] = local_ix;

						for(k=0;k<3;k++)
							(*points)[local_ix][k] = g_points[ g_elems[ global_ix*8 + j ] ][k];

						local_ix++;


					}

					(*elems)[i*8+j] = points_to_local[ g_elems[ global_ix*8 + j ] ];

					//(*elems)[i*8+j] = g_elems[ global_ix*8 + j ];

				}

			}

			*points_count=local_ix;

			free(points_to_local);
}





void compute_local_indices(int nprocs, int myrank,
        int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, int* points_count,
        int*** points, int** elems, int g_points_count, int g_nextcf,
        int** g_lcc, int* g_elems, int** g_points, int quotient, int remainder, int get_points){

			int i, j, k;
			int local_ix=0, global_ix, local_ext_ix=0;
			int * points_to_local = (int*) malloc(g_points_count*sizeof(int));
			// points_to_local: will have index from 0 to max_index of the point contained in the domain
			// the rest will be -1 --> point with global index [i] not in the local domain
			int * elems_to_local = (int*) malloc(g_nextcf*sizeof(int));
			//*local_global_index = (int*) malloc((*nintcf-*nintci+1)*sizeof(int));

			for(i=0;i<g_points_count;i++)
				points_to_local[i]=-1;

			for(i=0;i<g_nextcf;i++)
				elems_to_local[i]=-1;

			int my_base = myrank*quotient+MIN(remainder, myrank);
			for(i=*nintci; i<*nintcf+1; i++){
				//global_ix = i+my_base;
				//(*local_global_index)[i] = global_ix;


				for(j=0; j<6; j++){
					if( g_lcc[ i ][ j ] <= *nintcf+my_base &&  g_lcc[ i ][ j ] >= my_base)
						(*lcc)[ i ][ j ] = g_lcc[ i ][ j ]-my_base;
					else{
						if( elems_to_local[  g_lcc[ i ][ j ]   ] == -1 ) {
							elems_to_local[  g_lcc[ i ][ j ]   ] = local_ext_ix;
							local_ext_ix++;
						}
						(*lcc)[ i ][ j ] = *nintcf+1 +elems_to_local[  g_lcc[ i ][ j ]   ];
					}
				}

				for(j=0; j<8; j++){

					if( points_to_local[  g_elems[ i*8 + j ]   ] == -1 ) {
						points_to_local[  g_elems[ i*8 + j ]   ] = local_ix;

						if(get_points!=0)
							for(k=0;k<3;k++)
								(*points)[local_ix][k] = g_points[ g_elems[ i*8 + j ] ][k];

						local_ix++;
					}

					(*elems)[i*8+j] = points_to_local[ g_elems[ i*8 + j ] ];

					//(*elems)[i*8+j] = g_elems[ i*8 + j ];

				}

			}

			if(get_points!=0)
				*points_count=local_ix;
			*nextci = *nintcf+1;
			*nextcf = *nextci + local_ext_ix - 1;
}












int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index) {

	int g_nintci, g_nintcf, g_nextci, g_nextcf, **g_lcc, g_points_count, **g_points, *g_elems;
	double* g_bs, * g_be, * g_bn, * g_bw, * g_bl, * g_bh, * g_bp, * g_su;


    /********** START INITIALIZATION **********/
    int i = 0, j=0, k=0;


    if(strcmp(read_type, "allread") == 0){
		// read-in the input file
		int f_status = read_binary_geo(file_in, &g_nintci, &g_nintcf, &g_nextci, &g_nextcf, &g_lcc, &g_bs,
									   &g_be, &g_bn, &g_bw, &g_bl, &g_bh, &g_bp, &g_su, &g_points_count,
									   &g_points, &g_elems);

		if ( f_status != 0 ) return f_status;


		if(strcmp(part_type, "classic")==0){

			int quotient = (g_nintcf-g_nintci+1) / nprocs;
			int remainder = (g_nintcf-g_nintci+1) % nprocs;

			*nintci = g_nintci;
			*nintcf = quotient + ( myrank < remainder ? 1 : 0 ) - 1;

			int my_base = myrank*quotient+MIN(remainder, myrank);

			allocate_local_arrays(nintci, nintcf, &g_nextci, &g_nextcf, lcc, bs, be, bn, bw, bl, bh, bp, su, &g_points_count, points, elems);

			memcpy(*bw, g_bw+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*bs, g_bs+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*be, g_be+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*bl, g_bl+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*bn, g_bn+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*bh, g_bh+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*bp, g_bp+my_base, (*nintcf-*nintci+1)*sizeof(double));
			memcpy(*su, g_su+my_base, (*nintcf-*nintci+1)*sizeof(double));

			*local_global_index = (int*) malloc((*nintcf-*nintci+1)*sizeof(int));
			for(i=*nintci; i<*nintcf+1; i++)
				(*local_global_index)[i] = my_base+i;

			compute_local_indices(nprocs, myrank,
							nintci, nintcf, nextci,
							nextcf, lcc, points_count,
							 points, elems,  g_points_count, g_nextcf,
							 &(g_lcc[my_base]), &(g_elems[my_base*8]), g_points, quotient, remainder, 1);

		} else if (strcmp(part_type, "dual")==0){

			idx_t *epart;
			call_metis(g_nintci, g_nintcf, g_points_count, g_elems, nprocs, myrank, &epart);

			allocate_local_arrays(&g_nintci, &g_nintcf, &g_nextci, &g_nextcf, lcc, bs, be, bn, bw, bl, bh, bp, su, &g_points_count, points, elems);

			//modify this
			//*nextcf = g_nextcf;
			//*nextci = g_nextci;

			*local_global_index = (int*) malloc((g_nintcf-g_nintci+1)*sizeof(int));
			int local_index=0;
			for(i=g_nintci;i<=g_nintcf;i++){
				if(epart[i]==myrank){
					(*bw)[local_index] = g_bw[ i ];
					(*bs)[local_index] = g_bs[ i ];
					(*be)[local_index] = g_be[ i ];
					(*bl)[local_index] = g_bl[ i ];
					(*bn)[local_index] = g_bn[ i ];
					(*bh)[local_index] = g_bh[ i ];
					(*bp)[local_index] = g_bp[ i ];
					(*su)[local_index] = g_su[ i ];

					memcpy((*lcc)[local_index], g_lcc[ i ], 6*sizeof(int));
					memcpy(&( (*elems)[local_index*8] ), &(g_elems[ i*8 ]), 8*sizeof(int));
					//memcpy((*points)[local_index], g_points[ i ], 3*sizeof(int));
// distr punti
					(*local_global_index)[local_index] = i;

					local_index++;
				}
			}

			*nintci = g_nintci;
			*nintcf = local_index-1;

			compute_local_indices(nprocs, 0,
							        nintci, nintcf, nextci,
							        nextcf, lcc, points_count,
							         points, elems,  g_points_count, g_nextcf,
							         *lcc, *elems, *points, 0, 0, 0);

			//modify this
			//*points_count = (local_index*8);

		}

    } else if(strcmp(read_type, "oneread") == 0){
    	if(myrank==0){
    		int f_status = read_binary_geo(file_in, &g_nintci, &g_nintcf, &g_nextci, &g_nextcf, &g_lcc, &g_bs,
    											   &g_be, &g_bn, &g_bw, &g_bl, &g_bh, &g_bp, &g_su, &g_points_count,
    											   &g_points, &g_elems);

			if ( f_status != 0 ) return f_status;


//			int quotient = (g_nintcf-g_nintci+1) / nprocs;
//			int remainder = (g_nintcf-g_nintci+1) % nprocs;
//			int my_base = myrank*quotient+MIN(remainder, myrank);
			*local_global_index = (int*) malloc((g_nintcf-g_nintci+1)*sizeof(int));
			for(i=g_nintci; i<g_nintcf+1; i++)
				(*local_global_index)[i] = i;



			idx_t *epart;
			    call_metis(g_nintci, g_nintcf, g_points_count, g_elems, nprocs, myrank, &epart);

			    char fname[100];
			    sprintf(fname, "out.aa.aa.%d.0.vtk", myrank);
			    double* s = malloc((g_nintcf+1)*sizeof(double));
			    for(i=0;i<=g_nintcf;i++)
			    	s[i]=(double)(epart[i]);
			    test_distribution(file_in, fname, *local_global_index, g_nintcf, s);

//			    sprintf(fname, "out2.aa.aa.%d.0.vtk", myrank);
//			    test_distribution2(file_in, nintci, nintcf, nextci, nextcf, lcc,
//			    		bs, be,bn,bw,bh,bl,
//			            bp, su,  points_count, points,  elems,
//			            fname, *local_global_index, g_nintcf, s);

    	}

    	MPI_Bcast(&g_nintci, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&g_nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	MPI_Bcast(&g_nextci, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&g_nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&g_points_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

		int quotient = (g_nintcf-g_nintci+1) / nprocs;
		int remainder = (g_nintcf-g_nintci+1) % nprocs;

		*nintci = g_nintci;
		*nintcf = quotient + ( myrank < remainder ? 1 : 0 ) - 1;

		allocate_local_arrays(nintci, nintcf, &g_nextci, &g_nextcf, lcc, bs, be, bn, bw, bl, bh, bp, su, &g_points_count, points, elems);


		int *sendcount, *displs;
		sendcount = (int*) malloc(nprocs*sizeof(int));
		displs = (int*) malloc(nprocs*sizeof(int));

		displs[0] = 0;
		for(i=0;i<nprocs;i++)
			sendcount[i] = quotient;
		for(i=0;i<remainder;i++)
			sendcount[i] += 1;
		displs[0] = 0;
		for(i=1;i<nprocs;i++)
			displs[i] = displs[i-1] + sendcount[i-1];

		/*if(myrank == 0){
			printf("sendcount:\t");
			for(i=0;i<nprocs;i++)
				printf("%d \t", sendcount[i]);

			printf("\n displs:\t");
			for(i=0;i<nprocs;i++)
				printf("%d \t", displs[i]);
			printf("\n");
		}*/

		MPI_Scatterv(g_bw, sendcount, displs, MPI_DOUBLE, *bw, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_bs, sendcount, displs, MPI_DOUBLE, *bs, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_be, sendcount, displs, MPI_DOUBLE, *be, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_bl, sendcount, displs, MPI_DOUBLE, *bl, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_bn, sendcount, displs, MPI_DOUBLE, *bn, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_bh, sendcount, displs, MPI_DOUBLE, *bh, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_bp, sendcount, displs, MPI_DOUBLE, *bp, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(g_su, sendcount, displs, MPI_DOUBLE, *su, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

		int my_base = myrank*quotient+MIN(remainder, myrank);
		*local_global_index = (int*) malloc((*nintcf-*nintci+1)*sizeof(int));
		for(i=*nintci; i<*nintcf+1; i++)
			(*local_global_index)[i] = my_base+i;

		if(myrank==0){
			for(i=0;i<nprocs;i++)
				for(j=displs[i]; j < displs[i]+sendcount[i]; j++)
					MPI_Send(g_lcc[ j ], 6, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for(j=0; j < sendcount[myrank]; j++)
			MPI_Recv((*lcc)[ j ], 6, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for(i=0;i<nprocs;i++)
			sendcount[i] *= 8;

		displs[0] = 0;
		for(i=1;i<nprocs;i++)
			displs[i] = displs[i-1] + sendcount[i-1];

		MPI_Scatterv(g_elems, sendcount, displs, MPI_INT, *elems, sendcount[myrank], MPI_INT, 0, MPI_COMM_WORLD);

		if(myrank==0){
			for(i=nprocs-1;i>0;i--){
				get_points(nprocs, i, quotient, remainder, points_count, points, elems, g_points_count, g_elems, g_points);
				MPI_Send(points_count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

				for(j=0; j < *points_count; j++)
					MPI_Send((*points)[ j ], 3, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			get_points(nprocs, 0, quotient, remainder, points_count, points, elems, g_points_count, g_elems, g_points);
		} else {
			MPI_Recv(points_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(j=0; j < *points_count; j++)
				MPI_Recv((*points)[ j ], 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		compute_local_indices(nprocs, myrank,
				        nintci, nintcf, nextci,
				        nextcf, lcc, points_count,
				         points, elems,  g_points_count, g_nextcf,
				         *lcc, *elems, *points, quotient, remainder, 0);
    }

    realloc_local_array(nextcf, bs,  be,  bn, bw, bl, bh, bp, su);

    printf("id %d: nintci %d nintcf %d nextci %d nextcf %d \nid a: g_nintci %d g_nintcf %d g_nextci %d g_nextcf %d point_count %d \n", myrank, *nintci, *nintcf, *nextci, *nextcf, g_nintci, g_nintcf, g_nextci, g_nextcf, *points_count);




    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*cnorm)[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*var)[i] = 0.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*cgup)[i] = 1.0 / ((*bp)[i]);
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bh)[i] = 0.0;
        (*bl)[i] = 0.0;
    }

    return 0;
}


