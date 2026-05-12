#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "star_properties.h"
#include "utilities.h"

double get_msq_lifetime(double mass)
{
	const double mass_table[] = { 0.079,   0.10,   0.20,  0.45,  1.00, 2.00, 5.00, 10.00, 20.00, 60.00,  150.00, 200.00 };
	const double lifetime_table[] = { 12000.0, 6000.0, 560.0, 80.0,  10.0, 1.2,  0.10, 0.020, 0.006, 0.0035, 0.0025, 0.0022 };

	return my_log_interpolate(mass, mass_table, lifetime_table, SIZE(mass_table));
}

double get_total_lifetime(double mass)
{
	const double mass_table[] = { 0.079, 0.50, 1.00, 2.00, 8.00, 20.00, 60.00, 200.0 };
	const double extra_fraction_table[] = { 1.20,  1.18, 1.15, 1.12, 1.08, 1.05,  1.03,  1.02 };

	return get_msq_lifetime(mass) * my_log_interpolate(mass, mass_table, extra_fraction_table, SIZE(mass_table));
}

double generate_metallicity(void)
{
	const double roll = my_rand_double(0.0, 100.0);
	double fe_h;

	if (roll < 0.6)	fe_h = my_rand_normal(0.30, 0.13);
	else if (roll < 80.6)	fe_h = my_rand_normal(-0.05, 0.22);
	else if (roll < 94.6)	fe_h = my_rand_normal(-0.55, 0.25);
	else if (roll < 99.9)	fe_h = my_rand_normal(-1.45, 0.35);
	else					fe_h = my_rand_normal(-3.0, 0.40);

	if (fe_h < -5.5)	fe_h = -5.5;
	if (fe_h > 1.0)		fe_h = 1.0;

	return fe_h;
}

double get_density(double mass, double radius)
{
	const double SOLAR_MASS = 1.98847e33;   // grams
	const double SOLAR_RADIUS = 6.957e10;   // cm

	return (mass * SOLAR_MASS) / ((4.0 / 3.0) * M_PI * pow(radius * SOLAR_RADIUS, 3));
}

double get_luminosity(double radius, int surface_temp)
{
	return pow(radius, 2) * pow(surface_temp / 5772.0, 4.0);
}