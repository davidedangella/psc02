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
#include "test_functions.h"
#include "mmalloc.h"
#include <mpi.h>
#include <libgen.h>
#include <papi.h>

#include <metis.h>

#define MIN(a,b) ((a)<(b)?(a):(b))


typedef struct {
	int nintci;
	int nintcf;
	int nextci;
	int nextcf;
	int** lcc;
	double* bs;
	double* be;
	double* bn;
	double* bw;
	double* bl;
	double* bh;
	double* bp;
	double* su;
	int points_count;
	int** points;
	int* elems;

} MeshData;



//void call_metis(int nintci, int nintcf, int points_count, int* elems, int nprocs, int my_rank, idx_t** epart){
void call_metis(int part_key, int nprocs, int my_rank, MeshData *global, idx_t** epart){
		int i;
		idx_t nelem = (idx_t) global->nintcf+1;
		idx_t nnodes = (idx_t) global->points_count;


		// how to structure them??

		idx_t *eptr = malloc( sizeof(idx_t) * ( nelem + 1 ) );
		idx_t *eind = malloc( sizeof(idx_t) * ( nelem  ) * 8); // every cell: six neighbours


		/* every elem HAS 8 POINTS !*/
		for (i=0; i< nelem + 1; i++){
			eptr[i] = i*8;
		}

		for (i=0; i< (nelem )*8; i++){
			eind[i] = (idx_t)global->elems[i];
		}



		idx_t nparts = (idx_t) nprocs;
		idx_t ncommon = 4;
		idx_t objval;
		*epart = malloc( sizeof(idx_t) * nelem );
		idx_t *npart = malloc( sizeof(idx_t) * nnodes );


		METIS_API(int) status ;
		if(part_key == PART_DUAL)
			status = METIS_PartMeshDual( &nelem, &nnodes, eptr,  eind, NULL , NULL, &ncommon, &nparts, NULL , NULL, &objval, *epart, npart);
		 else if (part_key == PART_NODAL)
			status = METIS_PartMeshNodal( &nelem, &nnodes, eptr,  eind, NULL , NULL, &nparts, NULL , NULL, &objval, *epart, npart);


		if( status != METIS_OK ) 	printf("ERROR IN METIS\n");
		else						printf("METIS SUCCESSFFUL\n");

}



void compute_local_lcc( MeshData* local, MeshData* global, int* local_to_global, int* global_to_local){
	int i, j;
	int global_elem;
	int local_external_node_ix = local->nintcf+1;

	for(i=local->nintci; i<local->nintcf+1; i++)

		for(j=0; j<6; j++){

			global_elem = global->lcc[ local_to_global[i] ][ j ];


			if( global_to_local[ global_elem ] == -1 ){

					if(global_elem < global->nextci){
						// internal cell not in the domain of the current processor
						// Which value should we put here?
						global_to_local[ global_elem ] = -global_elem;
					}
					else{
						local_to_global[local_external_node_ix] = global_elem;
						global_to_local[ global_elem ] = local_external_node_ix;
						local_external_node_ix++;
					}
			}

			local->lcc[ i ][ j ] = global_to_local[  global_elem   ];
		}

	local->nextci = local->nintcf+1;
	local->nextcf = local_external_node_ix -1;
}

void compute_local_lcc2( MeshData* local, MeshData* global, int* local_to_global, int* global_to_local){
	int i, j;
	int global_elem;
	int local_external_node_ix = local->nintcf+1;

	for(i=local->nintci; i<local->nintcf+1; i++)

		for(j=0; j<6; j++){

			global_elem = global->lcc[ i ][ j ];


			if( global_to_local[ global_elem ] == -1 ){

					if(global_elem < global->nextci){
						// internal cell not in the domain of the current processor
						// Which value should we put here?
						global_to_local[ global_elem ] = -global_elem;
					}
					else{
						//local_to_global[local_external_node_ix] = global_elem;
						global_to_local[ global_elem ] = local_external_node_ix;
						local_external_node_ix++;
					}
			}

			local->lcc[ i ][ j ] = global_to_local[  global_elem   ];
		}

	local->nextci = local->nintcf+1;
	local->nextcf = local_external_node_ix -1;
}



int get_global_elements_data(char* file_in, int part_key, int read_key, int nprocs, int myrank, MeshData* global){
	if(read_key == READ_ALLREAD){
			// read-in the input file
			int f_status = read_binary_geo(file_in, &(global->nintci), &(global->nintcf), &(global->nextci), &(global->nextcf), &(global->lcc), &(global->bs),
											&(global->be), &(global->bn), &(global->bw), &(global->bl), &(global->bh), &(global->bp), &(global->su), &(global->points_count),
											&(global->points), &(global->elems));

			if ( f_status != 0 ) return f_status;
	} else if(read_key == READ_ONEREAD) {
		if(myrank==0){

			int f_status = read_binary_geo(file_in, &(global->nintci), &(global->nintcf), &(global->nextci), &(global->nextcf), &(global->lcc), &(global->bs),
											&(global->be), &(global->bn), &(global->bw), &(global->bl), &(global->bh), &(global->bp), &(global->su), &(global->points_count),
											&(global->points), &(global->elems));

			if ( f_status != 0 ) return f_status;

		}
	}

	return 0;
}

void get_local_elements_data( int nprocs, int myrank, int part_key, int read_key, MeshData* global, MeshData* local, int** local_to_global){
	int i,j;




	if(read_key == READ_ALLREAD) {
		if(part_key == PART_CLASSIC){

			int g_ncells = global->nintcf-global->nintci+1;

			int quotient = g_ncells / nprocs;
			int remainder = g_ncells % nprocs;

			local->nintci = global->nintci;
			local->nintcf = quotient + ( myrank < remainder ? 1 : 0 ) - 1;
			local->points_count =  global->points_count;

			int ncells = local->nintcf-local->nintci+1;

			allocate_local_arrays(&(local->nintci), &(local->nintcf), &(global->nextci), &(global->nextcf), &(local->lcc), &(local->bs), &(local->be),
					&(local->bn), &(local->bw), &(local->bl), &(local->bh), &(local->bp), &(local->su), &(global->points_count),
					&(local->points), &(local->elems), local_to_global);

			int my_base = myrank*quotient+MIN(remainder, myrank);

			memcpy(local->bw, &(global->bw[my_base]), ncells*sizeof(double));
			memcpy(local->bs, &(global->bs[my_base]), ncells*sizeof(double));
			memcpy(local->be, &(global->be[my_base]), ncells*sizeof(double));
			memcpy(local->bl, &(global->bl[my_base]), ncells*sizeof(double));
			memcpy(local->bn, &(global->bn[my_base]), ncells*sizeof(double));
			memcpy(local->bh, &(global->bh[my_base]), ncells*sizeof(double));
			memcpy(local->bp, &(global->bp[my_base]), ncells*sizeof(double));
			memcpy(local->su, &(global->su[my_base]), ncells*sizeof(double));

			// it seems that elems and points must not be distributed to the processors

			for(i=local->nintci; i<local->nintcf+1; i++)
				(*local_to_global)[i] = my_base+i;

			int *global_to_local = (int*) malloc((global->nextcf+1)*sizeof(int));
			for(i=global->nintci; i<global->nextcf+1; i++)
				global_to_local[i] = -1;
			for(i=my_base; i<my_base+local->nintcf+1; i++)
				global_to_local[i] = i-my_base;

			compute_local_lcc(local, global, *local_to_global, global_to_local);

		} else if (part_key == PART_NODAL || part_key == PART_DUAL){

			idx_t *epart;
			call_metis(part_key, nprocs, myrank, global, &epart);

			allocate_local_arrays(&(global->nintci), &(global->nintcf), &(global->nextci), &(global->nextcf), &(local->lcc), &(local->bs), &(local->be), &(local->bn), &(local->bw), &(local->bl), &(local->bh), &(local->bp), &(local->su), &(global->points_count),
					&(local->points), &(local->elems), local_to_global);

			int *global_to_local = (int*) malloc((global->nextcf+1)*sizeof(int));

			int local_index=0;
			for(i=global->nintci;i<=global->nintcf;i++){

				if(epart[i]!=myrank){
					global_to_local[i] = -1;
				} else {
					local->bw[local_index] = global->bw[ i ];
					local->bs[local_index] = global->bs[ i ];
					local->be[local_index] = global->be[ i ];
					local->bl[local_index] = global->bl[ i ];
					local->bn[local_index] = global->bn[ i ];
					local->bh[local_index] = global->bh[ i ];
					local->bp[local_index] = global->bp[ i ];
					local->su[local_index] = global->su[ i ];

					memcpy(local->lcc[local_index], global->lcc[ i ], 6*sizeof(int));

					(*local_to_global)[local_index] = i;
					global_to_local[i] = local_index;

					local_index++;
				}
			}
			local->nintci = global->nintci;
			local->nintcf = local_index-1;
			local->points_count =  global->points_count;

			compute_local_lcc(local, global, *local_to_global, global_to_local);
		}
	} else if(read_key == READ_ONEREAD){
		if(part_key == PART_CLASSIC){
			MPI_Bcast(&(global->nintci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nintcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->points_count), 1, MPI_INT, 0, MPI_COMM_WORLD);


			int g_ncells = global->nintcf-global->nintci+1;

			int quotient = g_ncells / nprocs;
			int remainder = g_ncells % nprocs;

			local->nintci = global->nintci;
			local->nintcf = quotient + ( myrank < remainder ? 1 : 0 ) - 1;
			local->points_count =  global->points_count;

			//			int ncells = local->nintcf-local->nintci+1;

			allocate_local_arrays(&(local->nintci), &(local->nintcf), &(global->nextci), &(global->nextcf), &(local->lcc), &(local->bs), &(local->be), &(local->bn), &(local->bw), &(local->bl), &(local->bh), &(local->bp), &(local->su), &(global->points_count),
					&(local->points), &(local->elems), local_to_global);


			int *sendcount, *displs;
			sendcount = (int*) malloc(nprocs*sizeof(int));
			displs = (int*) malloc(nprocs*sizeof(int));


			for(i=0;i<nprocs;i++)
				sendcount[i] = quotient;
			for(i=0;i<remainder;i++)
				sendcount[i] += 1;

			displs[0] = 0;
			for(i=1;i<nprocs;i++)
				displs[i] = displs[i-1] + sendcount[i-1];

			MPI_Scatterv(global->bw, sendcount, displs, MPI_DOUBLE, local->bw, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->bs, sendcount, displs, MPI_DOUBLE, local->bs, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->be, sendcount, displs, MPI_DOUBLE, local->be, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->bl, sendcount, displs, MPI_DOUBLE, local->bl, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->bn, sendcount, displs, MPI_DOUBLE, local->bn, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->bh, sendcount, displs, MPI_DOUBLE, local->bh, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->bp, sendcount, displs, MPI_DOUBLE, local->bp, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatterv(global->su, sendcount, displs, MPI_DOUBLE, local->su, sendcount[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

			int my_base = myrank*quotient+MIN(remainder, myrank);

			if(myrank==0){
				for(i=0;i<nprocs;i++)
					for(j=displs[i]; j < displs[i]+sendcount[i]; j++)
						MPI_Send(global->lcc[ j ], 6, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			for(j=0; j < sendcount[myrank]; j++)
				MPI_Recv(local->lcc[ j ], 6, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


			for(i=local->nintci; i<local->nintcf+1; i++)
				(*local_to_global)[i] = i+my_base; // modify


			int *global_to_local = (int*) malloc((global->nextcf+1)*sizeof(int));
			for(i=global->nintci; i<global->nextcf+1; i++)
				global_to_local[i] = -1;
			for(i=my_base; i<my_base+local->nintcf+1; i++)
				global_to_local[i] = i-my_base;

			local->nextci=global->nextci;

			compute_local_lcc2(local, local, *local_to_global, global_to_local);
		}else if (part_key==PART_DUAL || part_key==PART_NODAL){

			MPI_Bcast(&(global->nintci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nintcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->points_count), 1, MPI_INT, 0, MPI_COMM_WORLD);

			if(myrank==0){
				idx_t *epart;
				call_metis(part_key, nprocs, myrank, global, &epart);

				int *size = (int*) calloc(nprocs, sizeof(int));
				int *local_index = (int*) calloc(nprocs, sizeof(int));

				for(i=global->nintci;i<=global->nintcf;i++)
					size[epart[i]]++;

				int **global_to_local =  (int**) malloc(nprocs*sizeof(int*));
				int **local_to_global_all =  (int**) malloc(nprocs*sizeof(int*));

				int **bw = (int**) malloc(nprocs*sizeof(int*));
				int **bs = (int**) malloc(nprocs*sizeof(int*));
				int **be = (int**) malloc(nprocs*sizeof(int*));
				int **bl = (int**) malloc(nprocs*sizeof(int*));
				int **bn = (int**) malloc(nprocs*sizeof(int*));
				int **bh = (int**) malloc(nprocs*sizeof(int*));
				int **bp = (int**) malloc(nprocs*sizeof(int*));
				int **su = (int**) malloc(nprocs*sizeof(int*));
				int **lcc = (int**) malloc(nprocs*sizeof(int*));
				for(i=0; i<nprocs;i++){
					bw[i] =  (int*) malloc(size[i]*sizeof(int));
					bs[i] =  (int*) malloc(size[i]*sizeof(int));
					be[i] =  (int*) malloc(size[i]*sizeof(int));
					bl[i] =  (int*) malloc(size[i]*sizeof(int));
					bn[i] =  (int*) malloc(size[i]*sizeof(int));
					bh[i] =  (int*) malloc(size[i]*sizeof(int));
					bp[i] =  (int*) malloc(size[i]*sizeof(int));
					su[i] =  (int*) malloc(size[i]*sizeof(int));

					lcc[i] =  (int*) malloc(size[i]*6*sizeof(int));
					global_to_local[i]  =  (int*) malloc((global->nextcf+1)*sizeof(int));
					local_to_global_all[i]  =  (int*) malloc(size[i]*sizeof(int));
				}

				for(i=global->nintci;i<=global->nintcf;i++){

					for(j=0;j<nprocs;j++)
						global_to_local[j][i] = -1;

					int proc = epart[i];
					bw[proc][local_index[proc]] = global->bw[ i ];
					bs[proc][local_index[proc]] = global->bs[ i ];
					be[proc][local_index[proc]] = global->be[ i ];
					bl[proc][local_index[proc]] = global->bl[ i ];
					bn[proc][local_index[proc]] = global->bn[ i ];
					bh[proc][local_index[proc]] = global->bh[ i ];
					bp[proc][local_index[proc]] = global->bp[ i ];
					su[proc][local_index[proc]] = global->su[ i ];

					memcpy(&(lcc[proc][local_index[proc]*6]), global->lcc[ i ], 6*sizeof(int));

					local_to_global_all[proc][local_index[proc]] = i;
					global_to_local[proc][i] = local_index[proc];

					(local_index[proc])++;

				}

				for(i=1;i<nprocs;i++){
					MPI_Send(size+i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

					MPI_Send(bw[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(bs[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(be[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(bl[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(bn[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(bh[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(bp[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(su[i], size[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

					MPI_Send(lcc[i], 6*size[i], MPI_INT, i, 0, MPI_COMM_WORLD);

					MPI_Send(local_to_global_all[i], size[i], MPI_INT, i, 0, MPI_COMM_WORLD);
					MPI_Send(global_to_local[i], global->nextcf+1, MPI_INT, i, 0, MPI_COMM_WORLD);
				}

				local->nintci = global->nintci;
				local->nintcf = size[0]-1;
				local->points_count =  global->points_count;

				allocate_local_arrays(&(local->nintci), &(local->nintcf), &(global->nextci), &(global->nextcf), &(local->lcc), &(local->bs), &(local->be),
						&(local->bn), &(local->bw), &(local->bl), &(local->bh), &(local->bp), &(local->su), &(global->points_count),
						&(local->points), &(local->elems), local_to_global);

				memcpy(*local_to_global, local_to_global_all[0], size[0]*sizeof(int));

				memcpy(local->bw, bw[0], size[0]*sizeof(double));
				memcpy(local->bs,bs[0], size[0]*sizeof(double));
				memcpy(local->be,be[0], size[0]*sizeof(double));
				memcpy(local->bl,bl[0], size[0]*sizeof(double));
				memcpy(local->bn,bn[0], size[0]*sizeof(double));
				memcpy(local->bh,bh[0], size[0]*sizeof(double));
				memcpy(local->bp,bp[0], size[0]*sizeof(double));
				memcpy(local->su,su[0], size[0]*sizeof(double));

				compute_local_lcc(local, global, *local_to_global, global_to_local[0]);

			    	char fname[100];
				sprintf(fname, "out.part.%d.vtk",  0);
				double* s = malloc((global->nintcf+1)*sizeof(double));
				int* identity = malloc((global->nintcf+1)*sizeof(int));
				for(i=0;i<=global->nintcf;i++){
					s[i]=(double)(epart[i]+1);
					identity[i]=i;
				}
				test_distribution("drall.geo.bin", fname, identity, global->nintcf+1, s);
			} else {

				int size;

				MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				local->nintci = global->nintci;
				local->nintcf = size-1;
				local->points_count =  global->points_count;

				allocate_local_arrays(&(local->nintci), &(local->nintcf), &(global->nextci), &(global->nextcf), &(local->lcc), &(local->bs), &(local->be), &(local->bn), &(local->bw), &(local->bl), &(local->bh), &(local->bp), &(local->su), &(global->points_count),
						&(local->points), &(local->elems), local_to_global);

				MPI_Recv(local->bw, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->bs, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->be, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->bl, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->bn, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->bh, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->bp, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(local->su, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				int* lcc_linear = (int*) malloc(size*6*sizeof(int));
				MPI_Recv(lcc_linear, 6*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				global->lcc=(int**) malloc(size*sizeof(int*));
				for(i=0;i<size;i++)
					global->lcc[i]=(int*) malloc(6*sizeof(int));
				for(i=0;i<size;i++)
					memcpy(global->lcc[i], lcc_linear+i*6, sizeof(int)*6);



				MPI_Recv(*local_to_global, size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				int *global_to_local = (int*) malloc((global->nextcf+1)*sizeof(int));
				MPI_Recv(global_to_local, global->nextcf+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				compute_local_lcc2(local, global, *local_to_global, global_to_local);
			}
		}
	}
}




void get_keys(char* file_in, char* part_type, char* read_type, int* input_key, int* part_key, int* read_key){
	if(strcmp(basename(file_in), "tjunc.geo.bin")==0)
		*input_key=INPUT_TJUNC;
	else if(strcmp(basename(file_in), "drall.geo.bin")==0)
		*input_key=INPUT_DRALL;
	else if(strcmp(basename(file_in), "pent.geo.bin")==0)
		*input_key=INPUT_PENT;
	else if(strcmp(basename(file_in), "cojack.geo.bin")==0)
		*input_key=INPUT_COJACK;

	if(strcmp(part_type, "classic")==0)
		*part_key=PART_CLASSIC;
	else if(strcmp(part_type, "dual")==0)
		*part_key=PART_DUAL;
	else if(strcmp(part_type, "nodal")==0)
		*part_key=PART_NODAL;

	if(strcmp(read_type, "oneread")==0)
		*read_key=READ_ONEREAD;
	else if(strcmp(read_type, "allread")==0)
		*read_key=READ_ALLREAD;
}


int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index) {
	MeshData global;
	MeshData local;

    /********** START INITIALIZATION **********/
    int i = 0;//, j=0, k=0;

    int input_key, part_key, read_key;

    get_keys(file_in, part_type, read_type, &input_key, &part_key, &read_key);


    long long start = PAPI_get_virt_usec();

    get_global_elements_data(file_in, part_key, read_key, nprocs, myrank, &global );
    get_local_elements_data(nprocs, myrank, part_key, read_key, &global, &local, local_global_index);

	*nintci= local.nintci ;
	*nintcf= local.nintcf ;
	*nextci= local.nextci ;
	*nextcf= local.nextcf ;
	* lcc= local.lcc ;
	* bs= local.bs ;
	* be= local.be ;
	* bn= local.bn ;
	* bw= local.bw ;
	* bl= local.bl ;
	* bh= local.bh ;
	* bp= local.bp ;
	* su= local.su ;
	* points_count= local.points_count ;
	* points= local.points ;
	* elems= local.elems ;

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

    long long end = PAPI_get_virt_usec();
	if((nprocs == 4 || nprocs == 12) && (input_key==INPUT_DRALL || input_key==INPUT_COJACK))
	    write_pstats_exectime( input_key, part_key, read_key, myrank, end-start );
	if((nprocs == 9) && (input_key==INPUT_PENT || input_key==INPUT_COJACK))
	    write_pstats_partition( input_key, part_key, myrank, *nintcf-*nintci+1, *nextcf-*nextci+1);

    return 0;
}


