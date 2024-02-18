#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "my_math/limit.h"

static double input_data[] = 
{
    0,
    1,
    2,
    3,
    4,
    5,
};

static double output_data[] = 
{
    0, 
    2, 
    4,
    6,
    8,
    10,
};

static const unsigned int train_count = sizeof(input_data) / sizeof(input_data[0]);
static const unsigned int generation_times = 100;

static inline double rand_number(void)
{
    return (double)rand() / RAND_MAX;
}

double cost(double w)
{
    double summaiton = 0.0;
    for(size_t i = 0; i < train_count; ++i)
    {
        double x = input_data[i];
        double y = x * w;
        
        double diff = y - output_data[i];
        /* You can replace the following statement with fabs(diff),
         * but its performance is not that good. */
        summaiton += diff * diff;
        
        //printf("Actual: %lf\t,Expected: %lf\n", y, train[i][1]);
    }
    summaiton /= train_count;

    return summaiton;
}

int main(void)
{
    srand(time(NULL));
    
    // y = x * w; w: weight; x: input; y = output;
    double w = rand_number() * 10;
    double rate = 1e-2;
    double deri = 0.0;
    
    // Train single cell with input data.
    printf("Summation of difference: %lf\n", cost(w));
    for(int i = 0; i < generation_times; ++i)
    {
        // The derivative (The rate is identical to *h* in numerical derivative).
        deri = first_symmetric_derivative(w, rate, cost);

        // The differentiation (The differentiation is equal to the production of the derivative and *h*).
        // The differentiation = The derivative * h = ((f(x+h) - f(x)) / h) * h = f(x+h) - f(x).
        w -= deri * rate;
    }

    printf("Summation of difference: %lf\n", cost(w));
    
    // Test the weight with input data.
    printf("Weight: %lf\n", w);
    for(size_t i = 0; i < train_count; ++i)
    {
        double x = input_data[i];
        double y = x * w;

        printf("Actual: %lf\t,Expected: %lf\n", y, output_data[i]);
    }
}
