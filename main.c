#include<stdio.h>
#include<stdlib.h>
#include<time.h>


#include "papi.h"
#define NUM_EVENTS 2

void fillArray(double*);
void printMatrix(double* matrix);
void resetMatrix(double*);
void resetAndFillMatrices();
double calculateGFLOPS(double);

//matrix multiplication functions
void dgemmIJK();
void dgemmJIK();
void dgemmKIJ();
void dgemmIKJ();
void dgemmJKI();
void dgemmKJI();


//Global Variables

int N = 0;

int sizeArray[] = { 1,2,3,4,5,6,7,8,9,10, 11, 12 };

double* matrixA;

double * matrixB;

double * matrixC;


int main() {

	//srand initialization
	srand(time(NULL));


	//memory allocation for array pointers
	matrixA = (double *)(malloc(N*N * sizeof(double)));
	matrixB = (double *)(malloc(N*N * sizeof(double)));
	matrixC = (double *)(malloc(N*N * sizeof(double)));

	//filling the matrix
	fillArray(matrixA);
	fillArray(matrixB);

	//reset matrixC
	resetMatrix(matrixC);




	//calling dgemm function to calculate C += A*B
	//dgemmIJK();
	//resetAndFillMatrices();



	//dgemmJIK();
	//resetAndFillMatrices();

	//dgemmKIJ();
	//resetAndFillMatrices();


	dgemmIKJ();
	resetAndFillMatrices();


	//dgemmJKI();
	//resetAndFillMatrices();

	// dgemmKJI();
	// resetAndFillMatrices();






	//freeing the dynamic memory
	free(matrixA);
	free(matrixB);
	free(matrixC);

	system("pause");
	return 0;
}


void resetAndFillMatrices(){

	resetMatrix(matrixA);
	resetMatrix(matrixB);
	resetMatrix(matrixC);

	//filling the matrix
	fillArray(matrixA);
	fillArray(matrixB);

}

void printMatrix(double * matrix) {


	printf("**************************************************************************************************************\n\n");
	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			printf("%lf  |", matrix[i*N + j]);
		}

		printf("\n");
	}
	printf("**************************************************************************************************************\n\n");
}


//Function fill array

void  fillArray(double * matrix) {

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			double r = (double)rand() / RAND_MAX*2.0;      //float in range -1 to 1
			matrix[i*N + j] = r;
		}

	}

}



//Function to initialize the matrix to 0
void resetMatrix(double* matrix) {

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			matrix[i*N + j] = 0;
		}

	}
}

/***********************************************************
 * Function: calculateGFLOPS
 * Calculates GFlops
************************************************************/
double	calculateGFLOPS(double time) {

	double gflop = (2.0 * N*N*N) / (time * 1e+9);

	return gflop;
}






/*******************************************************************************************************
Matrix Multiplication using IJK algorithm

********************************************************************************************************/
void dgemmIJK() {

	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 10; i++) {

		N = sizeArray[i];

		int Events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_TOT_CYC };

		int EventSet = PAPI_NULL;

		long long values[NUM_EVENTS];

		int retval;


		/* Initialize the Library*/

		retval = PAPI_library_init(PAPI_VER_CURRENT);


		/* Allocate space for the new eventset and do setup */
		retval = PAPI_create_eventset(&EventSet);


		/* Add Flops and total cycles to the eventset*/
		retval = PAPI_add_events(EventSet, Events, NUM_EVENTS);



		/* Start the counters */
		retval = PAPI_start(EventSet);

		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication

			for (int i = 0; i < N; i++) {

				for (int j = 0; j < N; j++) {

					double cij = C[i*N + j];

					for (int k = 0; k < N; k++) {

						cij = cij + A[i*N + k] * B[k*N + j];
						C[i*N + j] = cij;
					}
				}
			}

			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}

		/*Stop counters and store results in values */
		retval = PAPI_stop(EventSet, values);

		printf("**************************************************************************************************************\n\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t|| %lf\n", N, (sum / 3.0));
		sum = 0;

		printf("___________\n\n");

		printf("PAPI Data\n");

		for (int ctr = 0; ctr < NUM_EVENTS; ctr++) {

			printf("%lld \n", values[ctr]);
		}



		printf("**************************************************************************************************************\n\n");

	}
}

/*******************************************************************************************************
Matrix Multiplication using JIK algorithm

********************************************************************************************************/
void dgemmJIK() {


	int Events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_TOT_CYC };

	int EventSet = PAPI_NULL;

	long long values[NUM_EVENTS];


	int retval;

	/* Initialize the Library*/

	retval = PAPI_library_init(PAPI_VER_CURRENT);


	/* Allocate space for the new eventset and do setup */
	retval = PAPI_create_eventset(&EventSet);


	/* Add Flops and total cycles to the eventset*/
	retval = PAPI_add_events(EventSet, Events, NUM_EVENTS);


	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 10; i++) {

		N = sizeArray[i];

		/* Start the counters */
		retval = PAPI_start(EventSet);

		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication


			for (int j = 0; j < N; j++) {

				for (int i = 0; i < N; i++) {

					double cji = C[j*N + i];

					for (int k = 0; k < N; k++) {

						cji = cji + A[j*N + k] * B[k*N + i];
						C[j*N + i] = cji;
					}
				}
			}

			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}
		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf", N, (sum / 3.0));
		sum = 0;

		/*Stop counters and store results in values */
		retval = PAPI_stop(EventSet, values);

		printf("**************************************************************************************************************\n\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t|| %lf", N, (sum / 3.0));
		sum = 0;

		printf("___________\n\n");

		printf("PAPI Data\n");

		for (int ctr = 0; ctr < NUM_EVENTS; ctr++) {

			printf("%lld\n", values[ctr]);
		}



		printf("**************************************************************************************************************\n\n");
	}
}

/*******************************************************************************************************
Matrix Multiplication using KIJ algorithm

********************************************************************************************************/

void dgemmKIJ() {

	int Events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_TOT_CYC };

	int EventSet = PAPI_NULL;

	int retval;

	long long values[NUM_EVENTS];


	/* Initialize the Library*/

	retval = PAPI_library_init(PAPI_VER_CURRENT);


	/* Allocate space for the new eventset and do setup */
	retval = PAPI_create_eventset(&EventSet);


	/* Add Flops and total cycles to the eventset*/
	retval = PAPI_add_events(EventSet, Events, NUM_EVENTS);


	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 10; i++) {

		N = sizeArray[i];

		/* Start the counters */
		retval = PAPI_start(EventSet);

		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication


			for (int k = 0; k < N; k++) {

				for (int i = 0; i < N; i++) {

					double cki = C[k*N + i];

					for (int j = 0; j < N; j++) {

						cki = cki + A[k*N + j] * B[j*N + i];
						C[k*N + i] = cki;
					}
				}
			}
			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}
		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf", N, (sum / 3.0));
		sum = 0;

		/*Stop counters and store results in values */
		retval = PAPI_stop(EventSet, values);

		printf("**************************************************************************************************************\n\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t|| %lf", N, (sum / 3.0));
		sum = 0;

		printf("____________\n\n");

		printf("PAPI Data\n");

		for (int ctr = 0; ctr < NUM_EVENTS; ctr++) {

			printf("%lld\n", values[ctr]);
		}



		printf("**************************************************************************************************************\n\n");
	}
}

/*******************************************************************************************************
Matrix Multiplication using IKJ algorithm

********************************************************************************************************/
void dgemmIKJ() {




	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 12; i++) {

		N = sizeArray[i];


		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication


			for (int i = 0; i < N; i++) {

				for (int k = 0; k < N; k++) {

					double cik = C[i*N + k];

					for (int j = 0; j < N; j++) {

						cik = cik + A[i*N + j] * B[j*N + k];
						C[i*N + k] = cik;
					}
				}
			}

			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}
		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf", N, (sum / 3.0));
		printf("GFLOPS \t\t\t %lf\n\n",calculateGFLOPS(sum / 3.0));
		sum = 0;








		printf("**************************************************************************************************************\n\n");
	}
}

/*******************************************************************************************************
Matrix Multiplication using JKI algorithm

********************************************************************************************************/
void dgemmJKI() {

	int Events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_TOT_CYC };

	int EventSet = PAPI_NULL;

	long long values[NUM_EVENTS];

	int retval;


	/* Initialize the Library*/

	retval = PAPI_library_init(PAPI_VER_CURRENT);


	/* Allocate space for the new eventset and do setup */
	retval = PAPI_create_eventset(&EventSet);


	/* Add Flops and total cycles to the eventset*/
	retval = PAPI_add_events(EventSet, Events, NUM_EVENTS);


	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 10; i++) {

		N = sizeArray[i];

		/* Start the counters */
		retval = PAPI_start(EventSet);

		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication

			for (int j = 0; j < N; j++) {

				for (int k = 0; k < N; k++) {

					double cjk = C[j*N + k];

					for (int i = 0; i < N; i++) {

						cjk = cjk + A[j*N + i] * B[i*N + k];
						C[j*N + k] = cjk;
					}
				}
			}
			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}
		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf", N, (sum / 3.0));
		sum = 0;

		/*Stop counters and store results in values */
		retval = PAPI_stop(EventSet, values);

		printf("**************************************************************************************************************\n\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t|| %lf", N, (sum / 3.0));
		sum = 0;

		printf("_____________\n\n");

		printf("PAPI Data\n");

		for (int ctr = 0; ctr < NUM_EVENTS; ctr++) {

			printf("%lld\n", values[ctr]);
		}



		printf("**************************************************************************************************************\n\n");
	}
}

/*******************************************************************************************************
Matrix Multiplication using KJI algorithm

********************************************************************************************************/

void dgemmKJI() {

	int Events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_TOT_CYC };

	int EventSet = PAPI_NULL;

	long long values[NUM_EVENTS];


	int retval;

	/* Initialize the Library*/

	retval = PAPI_library_init(PAPI_VER_CURRENT);


	/* Allocate space for the new eventset and do setup */
	retval = PAPI_create_eventset(&EventSet);


	/* Add Flops and total cycles to the eventset*/
	retval = PAPI_add_events(EventSet, Events, NUM_EVENTS);


	double *A;
	double *B;
	double *C;

	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;

	A = matrixA;
	B = matrixB;
	C = matrixC;


	for (int i = 0; i < 10; i++) {

		N = sizeArray[i];

		/* Start the counters */
		retval = PAPI_start(EventSet);

		for (int i = 0; i < 3; i++) {

			start = clock();

			//matrix multiplication

			for (int k = 0; k < N; k++) {

				for (int j = 0; j < N; j++) {

					double ckj = C[k*N + j];

					for (int i = 0; i < N; i++) {

						ckj = ckj + A[k*N + i] * B[i*N + j];
						C[k*N + j] = ckj;
					}
				}
			}
			end = clock();

			cpu_time_used = ((double)(end - start));

			sum += cpu_time_used;
		}
		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf", N, (sum / 3.0));
		sum = 0;

		/*Stop counters and store results in values */
		retval = PAPI_stop(EventSet, values);

		printf("**************************************************************************************************************\n\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t|| %lf", N, (sum / 3.0));
		sum = 0;

		printf("_____________\n\n");

		printf("PAPI Data\n");

		for (int ctr = 0; ctr < NUM_EVENTS; ctr++) {

			printf("%lld\n", values[ctr]);
		}



		printf("**************************************************************************************************************\n\n");
	}
}
