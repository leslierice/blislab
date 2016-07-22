//timming dgemm gcc -o time_dgemm time_dgemm.c /your/path/libopenblas.a ./time_dgemm <m> <n> <k> e.g. ./time_dgemm 1000 1000 1000


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);

int main(int argc, char* argv[])
{
    int i;
    if(argc<4){
        printf("Input Error\n");
        return 1;
    }

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int k = atoi(argv[3]);

    int sizeofa = m * k;
    int sizeofb = k * n;
    int sizeofc = m * n;
    char ta = 'N';
    char tb = 'N';
    double alpha = 1.0;
    //double beta = 0.001;
    double beta = 1.0;

    int r;
    int nrepeats = 3;

    struct timeval start,finish;
    double duration, min_time;

    double* A = (double*)malloc(sizeof(double) * sizeofa);
    double* B = (double*)malloc(sizeof(double) * sizeofb);
    double* C = (double*)malloc(sizeof(double) * sizeofc);

    srand((unsigned)time(NULL));

    for (i=0; i<sizeofa; i++)
        A[i] = i%3+1;//(rand()%100)/10.0;

    for (i=0; i<sizeofb; i++)
        B[i] = i%3+1;//(rand()%100)/10.0;

    for (i=0; i<sizeofc; i++)
        C[i] = i%3+1;//(rand()%100)/10.0;
    //#if 0

    //printf("m=%d,n=%d,k=%d,alpha=%lf,beta=%lf,sizeofc=%d\n",m,n,k,alpha,beta,sizeofc);


    for ( int r = 0; r < nrepeats; r++ ) {

        gettimeofday(&start, NULL);
        dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
        gettimeofday(&finish, NULL);
        duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
        if ( r == 0 ) {
            min_time = duration;
        } else {
            min_time = min_time > duration ? duration : min_time;
        }
    }

    double gflops = 2.0 * m * n * k;
    gflops = gflops/min_time*1.0e-9;

    //printf("%dx%dx%d\t%lf second, %lf GFLOPS\n", m, n, k, duration, gflops);
    //printf( "%6lu %6lu %6lu %10.3e %6.3f\n", m, k, n, min_time, gflops );
    printf( "%6lu %6lu %6lu %6.3f\n", m, k, n, gflops );

    fflush(stdout);

    free(A);
    free(B);
    free(C);
    return 0;
}
