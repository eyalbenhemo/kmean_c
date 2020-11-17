#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*Get array of pointer to observations and calc their avg*/
float *calc_centroid(float **observations, const int len, const int d) {
    float *centroid = calloc(d, sizeof(float));
    int i = 0;
    int j = 0;
    for (; j < d; j++) {
        centroid[j] = 0;
    }
    for (; i < len; i++) {
        j = 0;
        for (; j < d; j++) {
            centroid[j] += observations[i][j];
        }
    }
    for (j = 0; j < d; j++) {
        centroid[j] /= len;
    }
    printf("%f, %f", centroid[0], centroid[1]);
    return centroid;
}

/*Get 2 observations and calc their distance*/
float euclidian_distance(const float a[], const float b[], const int d) {
    float dist = 0;
    int i = 0;
    float temp;
    for (; i < d; i++) {
        temp = (a[i] - b[i]) * (a[i] - b[i]);
        dist += temp;
    }
    return dist;
}

/*Get pointer to observation and pointer to array of centroid and return the index of closest centroid*/
int find_closest_centroid(const float a[], float **centroids, const int K, const int d) {
    int min_dist = -1;
    int min_cent = 0;
    int k = 0;
    int distance;
    for (; k < K; k++) {
        distance = euclidian_distance(a, centroids[k], d);
        if (distance < min_dist || min_dist == -1) {
            min_dist = distance;
            min_cent = k;
        }
    }
    return min_cent;
}

/*Get 2 pointers to centroids and check if they equal*//*
int check_if_equals(float new_centroids[], float centroids[]){

}
*//*Get centroids, MAX_ITER and observations
 * Calc centroids while num of iter <= MAX_ITER and last(centroids) != centroids
 * return centroids*//*
float** approximation_loop(float **centroids, int MAX_ITER ,float **observations){

}*/

int main(int argc, char *argv[]) {

    /*Parse arguments*/
    const int K = atoi(argv[1]);
    const int N = atoi(argv[2]);
    const int d = atoi(argv[3]);
    const int MAX_ITER = atoi(argv[4]);
    char c;
    int i = 0;
    double **centroids, **observations;

    /*Assertions*/
    assert(argc == 5 && "Need to get 4 args");
    assert(K > 0 && N > 0 && d > 0 && MAX_ITER > 0 && "Args should be positive");
    assert(K < N && "K need to be smaller than N");

    /*Define variables*/
    centroids = calloc(K, sizeof(double *));
    observations = calloc(N, sizeof(double *));
    for (; i < N; i++) {
        observations[i] = calloc(d, sizeof(double));
    }
    for (i = 0; i < K; i++) {
        centroids[i] = observations[i];
    }

    /*Read input*/
    i = 0;
    while (scanf("%lf%c", &observations[i / d][i % d], &c) == 2) { i++; }

    return 0;
}
