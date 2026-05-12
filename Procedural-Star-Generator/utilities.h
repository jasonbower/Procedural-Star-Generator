#ifndef UTILITIES_H
#define UTILITIES_H

#define SIZE(array) (sizeof(array) / sizeof((array)[0]))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef enum boolean { FALSE, TRUE } Boolean;

typedef enum status { FAILURE, SUCCESS } Status;

double big_random_unit(void);

double my_rand_double_bias(double min, double max, float bias);

double my_rand_double(double min, double max);

char* decimal_to_sci_notation(double decimal);

double my_rand_normal(double mean, double stddev);

double my_log_interpolate(double x, const double* x_vals, const double* y_vals, int size);

void build_cumulative_table(const double* weights, double* cumulative, int size);

#endif 