/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

#define INPUT_TJUNC 1
#define INPUT_DRALL 2
#define INPUT_PENT 3
#define INPUT_COJACK 4

#define PART_CLASSIC 1
#define PART_DUAL 2
#define PART_NODAL 3

#define READ_ONEREAD 1
#define READ_ALLREAD 2

int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index);

void get_keys(char* file_in, char* part_type, char* read_type, int* input_key, int* part_key, int* read_key);

#endif /* INITIALIZATION_H_ */

