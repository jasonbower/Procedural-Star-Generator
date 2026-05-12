#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilities.h"

double big_random_unit(void)
{
	const double base = (double)RAND_MAX + 1.0;
	double x = (double)rand();

	x = x * base + (double)rand();
	x = x * base + (double)rand();

	return x / (base * base * base);
}

/*  Bias < 1 pushes values closer to max
	Bias = 1 linear distribution
	Bias > 1 pushes values closer to min      */
double my_rand_double_bias(double min, double max, float bias)
{
	double r = pow(big_random_unit(), bias);

	return min + r * (max - min);
}

double my_rand_double(double min, double max)
{
	return my_rand_double_bias(min, max, 1.0);
}

char* decimal_to_sci_notation(double decimal)
{
	int i;
	char* str = malloc(32 * sizeof(char));

	if (str == NULL)
	{
		fprintf(stderr, "Error 4: program failed to allocate memory for string variable in decimal_to_sci_notation\n");
		exit(1);
	}

	if (decimal < 0.01 && decimal > 0)
	{
		for (i = 0; decimal < 1; i++)	decimal *= 10;

		snprintf(str, 20, "%.2f * 10^-%d", decimal, i);
	}
	else if (decimal >= 10000)
	{
		for (i = 0; decimal > 10; i++)	decimal /= 10;

		snprintf(str, 20, "%.2f * 10^%d", decimal, i);
	}
	else
		snprintf(str, 20, "%.2f", decimal);

	return str;
}

double my_rand_normal(double mean, double stddev)
{
	const double r1 = my_rand_double(0.000000000000001, 1.0);
	const double r2 = my_rand_double(0.0, 1.0);

	return mean + stddev * sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2);	// Box-Muller Transform
}

double my_log_interpolate(double x, const double* x_vals, const double* y_vals, int size)
{
	int i;
	double x1, x2;
	double log_y1, log_y2;
	double fraction;
	double log_y;

	if (x_vals == NULL || y_vals == NULL)
	{
		fprintf(stderr, "Error 5: null parameter(s) passed through my_log_interpolate\n");
		exit(1);
	}

	if (size < 2)
	{
		fprintf(stderr, "Error 6: illegal parameter(s) passed through my_log_interpolate\n");
		exit(1);
	}

	if (x <= x_vals[0])			return y_vals[0];

	if (x >= x_vals[size - 1])	return y_vals[size - 1];

	for (i = 0; i < size - 1; i++)
	{
		if (x >= x_vals[i] && x <= x_vals[i + 1])
		{
			x1 = x_vals[i];
			x2 = x_vals[i + 1];

			if (x2 == x1)
			{
				fprintf(stderr, "Error 7: illegal parameter(s) passed through my_log_interpolate\n");
				exit(1);
			}

			log_y1 = log10(y_vals[i]);
			log_y2 = log10(y_vals[i + 1]);
			fraction = (x - x1) / (x2 - x1);
			log_y = log_y1 + fraction * (log_y2 - log_y1);

			return pow(10.0, log_y);
		}
	}

	fprintf(stderr, "Error 8: my_log_interpolate failed\n");
	exit(1);
}

void build_cumulative_table(const double* weights, double* cumulative, int size)
{
	int i;
	double total = 0.0;

	for (i = 0; i < size; i++)
	{
		total += weights[i];
		cumulative[i] = total;
	}
}