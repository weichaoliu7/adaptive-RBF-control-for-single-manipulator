#ifndef COSINE_H
#define COSINE_H
#include <math.h>

#define PI 3.14159

typedef struct Data1 {
    double* y;
    int size;
} Data1;

void cosine(Data1 *output, double m, double Ts, double t0, double t1, double amp) {
    int ARRAY_SIZE = (t1 - t0) / Ts;
    output->size = ARRAY_SIZE;
    output->y = malloc(sizeof(double) * ARRAY_SIZE);
    double t = t0;
    int i = 0;
    while (t0 <= t1 && i < ARRAY_SIZE) {
        double value = amp * cos(m * t);
        output->y[i] = value;
        t += Ts;
        i++;
    }
}

#endif
