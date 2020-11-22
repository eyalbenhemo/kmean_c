#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int K, N, d, MAX_ITER;

/*Get array of pointer to observations and calc their avg*/
double **
calc_centroids(double **observations, const int *clusterAllocations, double **new_centroids, int *clustersLengths) {
    int i, j; /* looping variables */

    for (i = 0; i < K; i++) { /* initialize clusters lengths and values to 0 */
        clustersLengths[i] = 0;
        for (j = 0; j < d; j++) {
            new_centroids[i][j] = 0;
        }
    }

    for (i = 0; i < N; i++) { /* add all the values of the of the observations to the relevant centroid */
        clustersLengths[clusterAllocations[i]]++; /* update number of observations in vector */
        for (j = 0; j < d; j++) {
            new_centroids[clusterAllocations[i]][j] += observations[i][j]; /* add observation to vector */
        }
    }

    for (i = 0; i < K; i++) { /* Update the means to all centroids */
        for (j = 0; j < d; j++) {
            if (clustersLengths[i]) {
                new_centroids[i][j] /= clustersLengths[i];
            } else {
                new_centroids[i][j] = 0;
            }
        }
    }

    return new_centroids;
}

/*Get 2 observations and calc their distance*/
double euclidian_distance(const double a[], const double b[]) {
    double dist = 0;
    int i = 0;
    double temp;

    for (; i < d; i++) {
        temp = (a[i] - b[i]) * (a[i] - b[i]);
        dist += temp;
    }

    return dist;
}

/*Get pointer to observation and pointer to array of centroid and return the index of closest centroid*/
int find_closest_centroid(const double a[], double **centroids) {
    double min_dist = -1;
    int min_cent = 0;
    int k = 0;
    double distance;

    for (; k < K; k++) {
        distance = euclidian_distance(a, centroids[k]);
        if (distance < min_dist || min_dist == -1) {
            min_dist = distance;
            min_cent = k;
        }
    }

    return min_cent;
}

/*Get 2 array of centroids and check if they equal*/
int check_if_equals(double **new_centroids, double **centroids) {
    int i, j = 0;

    for (i = 0; i < K; i++) {
        for (; j <= d; j++) {
            if (centroids[i][j] != new_centroids[i][j]) {
                return 0;
            }
        }
    }

    return 1;
}

/*Get centroids, MAX_ITER and observations
 * Calc centroids while num of iter <= MAX_ITER and last(centroids) != centroids
 * return centroids*/
double **approximation_loop(double **observations) {
    int i, j;
    double **centroids = observations;
    double **newCentroids = malloc(K * sizeof(double *)); /* new centroids to be returned */
    int *clusterAllocations = malloc(N * sizeof(int)); /* create an array of where every observation in mapped to*/
    int *clustersLengths = calloc(K, sizeof(int)); /*create array of how many observations go to each centroid*/
    double **temp; /*swap variable*/
    assert(newCentroids != NULL && clusterAllocations != NULL && clustersLengths != NULL && "Allocation failed");

    for (i = 0; i < K; i++) { /* initialize clusters lengths and values to 0 */
        newCentroids[i] = calloc(d, sizeof(double));
        assert(newCentroids[i] != NULL && "Allocation failed");
    }

    for (j = 0; j < MAX_ITER; j++) {
        for (i = 0; i < N; i++) {
            clusterAllocations[i] = find_closest_centroid(observations[i], centroids);
        }
        calc_centroids(observations, clusterAllocations, newCentroids, clustersLengths);
        if (check_if_equals(centroids, newCentroids)) {
            break;
        }

        if (!j) {
            centroids = newCentroids;
            newCentroids = malloc(K * sizeof(double *));
            assert(newCentroids != NULL && "Allocation failed");
            for (i = 0; i < K; i++) {
                newCentroids[i] = calloc(d, sizeof(double));
                assert(newCentroids[i] != NULL && "Allocation failed");
            }
        } else {
            temp = centroids;
            centroids = newCentroids;
            newCentroids = temp;
        }
    }

    if (newCentroids != centroids) {
        for (i = 0; i < K; i++) {
            free(newCentroids[i]);
        }
        free(newCentroids);
    }
    free(clustersLengths);
    free(clusterAllocations);
    return centroids;
}

int main(int argc, char *argv[]) {
    char c;
    int i, j;
    double **centroids, **observations;

    /*Parse arguments*/
    K = atoi(argv[1]);
    N = atoi(argv[2]);
    d = atoi(argv[3]);
    MAX_ITER = atoi(argv[4]);

    /*Assertions*/
    assert(argc == 5 && "Need to get 4 args");
    assert(K > 0 && N > 0 && d > 0 && MAX_ITER > 0 && "Args should be positive");
    assert(K < N && "K need to be smaller than N");

    /*Define variables*/
    observations = malloc(N * sizeof(double *));
    assert(observations != NULL && "Allocation failed");
    for (i = 0; i < N; i++) {
        observations[i] = malloc(d * sizeof(double));
        assert(observations[i] != NULL && "Allocation failed");
    }

    /*Read input*/
    i = 0;
    while (scanf("%lf%c", &observations[i / d][i % d], &c) == 2) { i++; }

    /*Calc centroids*/
    centroids = approximation_loop(observations);

    /*Print centroids*/
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            if (j == d - 1) {
                printf("%.2f", centroids[i][j]);
            } else {
                printf("%.2f,", centroids[i][j]);
            }
        }
        if (i < K - 1) {
            printf("\n");
        }
    }

    for (i = 0; i < N; i++) {
        free(observations[i]);
    }
    free(observations);

    return 0;
}
