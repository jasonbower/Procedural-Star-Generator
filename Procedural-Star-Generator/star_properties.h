#ifndef STAR_PROPERTIES_H
#define STAR_PROPERTIES_H

double get_msq_lifetime(double mass);

double get_total_lifetime(double mass);

double generate_metallicity(void);

double get_density(double mass, double radius);

double get_luminosity(double radius, int surface_temp);

#endif