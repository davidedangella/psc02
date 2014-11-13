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
void call_metis(char* part_type, int nprocs, int my_rank, MeshData *global, idx_t** epart){
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
		if(strcmp(part_type, "dual")==0)
			status = METIS_PartMeshDual( &nelem, &nnodes, eptr,  eind, NULL , NULL, &ncommon, &nparts, NULL , NULL, &objval, *epart, npart);
		 else if (strcmp(part_type, "nodal")==0)
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



int get_global_elements_data(char* file_in, char* part_type, char* read_type, int nprocs, int myrank, MeshData* global){
	if(strcmp(read_type, "allread") == 0){
			// read-in the input file
			int f_status = read_binary_geo(file_in, &(global->nintci), &(global->nintcf), &(global->nextci), &(global->nextcf), &(global->lcc), &(global->bs),
											&(global->be), &(global->bn), &(global->bw), &(global->bl), &(global->bh), &(global->bp), &(global->su), &(global->points_count),
											&(global->points), &(global->elems));

			if ( f_status != 0 ) return f_status;
	} else if(strcmp(read_type, "oneread") == 0) {
		if(myrank==0){

			int f_status = read_binary_geo(file_in, &(global->nintci), &(global->nintcf), &(global->nextci), &(global->nextcf), &(global->lcc), &(global->bs),
											&(global->be), &(global->bn), &(global->bw), &(global->bl), &(global->bh), &(global->bp), &(global->su), &(global->points_count),
											&(global->points), &(global->elems));

			if ( f_status != 0 ) return f_status;

		}
	}

	return 0;
}

void get_local_elements_data( int nprocs, int myrank, char* part_type, char* read_type, MeshData* global, MeshData* local, int** local_to_global){
	int i,j;




	if(strcmp(read_type, "allread") == 0) {
		if(strcmp(part_type, "classic")==0){

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

		} else if (strcmp(part_type, "dual")==0 || strcmp(part_type, "nodal")==0){

			idx_t *epart;
			call_metis(part_type, nprocs, myrank, global, &epart);

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
	} else if(strcmp(read_type, "oneread") == 0){
		if(strcmp(part_type, "classic")==0){
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
		}else if (strcmp(part_type, "dual")==0 || strcmp(part_type, "nodal")==0){

			MPI_Bcast(&(global->nintci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nintcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextci), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->nextcf), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(global->points_count), 1, MPI_INT, 0, MPI_COMM_WORLD);

			if(myrank==0){
				idx_t *epart;
				call_metis(part_type, nprocs, myrank, global, &epart);

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






int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index) {

	//int [my_base])nintci, g_nintcf, g_nextci, g_nextcf, **g_lcc, g_points_count, **g_points, *g_elems;
	//double* g_bs, * g_be, * g_bn, * g_bw, * g_bl, * g_bh, * g_bp, * g_su;

	MeshData global;
	MeshData local;

    /********** START INITIALIZATION **********/
    int i = 0;//, j=0, k=0;


    get_global_elements_data(file_in, part_type, read_type, nprocs, myrank, &global );
    get_local_elements_data(nprocs, myrank, part_type, read_type, &global, &local, local_global_index);

    /*if(strcmp(read_type, "allread") == 0){
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
			// *nextcf = g_nextcf;
			// *nextci = g_nextci;

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

			int quotient = (g_nintcf-g_nintci+1) / nprocs;
			int remainder = (g_nintcf-g_nintci+1) % nprocs;
			compute_local_indices(nprocs, myrank,
							        nintci, nintcf, nextci,
							        nextcf, lcc, points_count,
							         points, elems,  g_points_count, g_nextcf,
							         *lcc, *elems, *points, quotient, remainder, 1);

			//modify this
			// *points_count = (local_index*8);

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



//			idx_t *epart;
//			    call_metis(g_nintci, g_nintcf, g_points_count, g_elems, nprocs, myrank, &epart);
//
//			    char fname[100];
//			    sprintf(fname, "out.aa.aa.%d.0.vtk", myrank);
//			    double* s = malloc((g_nintcf+1)*sizeof(double));
//			    for(i=0;i<=g_nintcf;i++)
//			    	s[i]=(double)(epart[i]);
//			    test_distribution(file_in, fname, *local_global_index, g_nintcf, s);

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

//		if(myrank == 0){
//			printf("sendcount:\t");
//			for(i=0;i<nprocs;i++)
//				printf("%d \t", sendcount[i]);
//
//			printf("\n displs:\t");
//			for(i=0;i<nprocs;i++)
//				printf("%d \t", displs[i]);
//			printf("\n");
//		}

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
*/



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

    return 0;
}


