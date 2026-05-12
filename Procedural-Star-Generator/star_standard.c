#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "star_standard.h"
#include "utilities.h"
#include "star_properties.h"

typedef enum temperature_class_letter
{
	TC_UNASSIGNED, M, K, G, F, A, B, O,
} TemperatureClassLetter;

typedef enum luminosity_class
{
	LC_UNASSIGNED, SD, V, IV, III, II, IB, IAB, IA, IA_PLUS
} LuminosityClass;

typedef struct temperature_class
{
	TemperatureClassLetter letter;
	int number;
} TemperatureClass;

typedef struct spectral_class
{
	TemperatureClass temperature_class;
	LuminosityClass luminosity_class;
} SpectralClass;

typedef struct star // KNOWN TYPE
{
	double mass;
	double metallicity;
	double age;
	double radius;
	int surface_temp;
	double luminosity;
	double density;
	SpectralClass class;
} Star;

static const double msq_mass_table[] =
{
	0.079, 0.086, 0.090, 0.105, 0.155, 0.24,  0.36, 0.45,  0.51,  0.59,
	0.6,   0.61,  0.64,	 0.67,  0.70,  0.72,  0.76, 0.82,  0.85,  0.86,
	0.89,  0.93,  0.95,  0.97,  0.98,  0.985, 0.99, 1.00,  1.04,  1.07,
	1.12,  1.16,  1.21,  1.25,  1.32,  1.38,  1.41, 1.43,  1.50,  1.60,
	1.63,  1.67,  1.76,  1.83,  1.86,  1.88,  1.92, 2.06,  2.17,  2.33,
	2.75,  3.38,  3.92,  4.30,  4.70,  5.10,  5.40, 7.30,  11.80, 17.70,
	20.2,  23.6,  28,	 35,	43,	   48,	  59,   100.0, 150.0, 200.0
};

static const double msq_radius_table[] =
{
	0.102, 0.114, 0.120, 0.137,  0.196,  0.274,  0.361,  0.446, 0.501, 0.588,
	0.608, 0.615, 0.630, 0.669,  0.701,  0.710,  0.755,  0.783, 0.797, 0.813,
	0.853, 0.914, 0.927, 0.949,  0.977,  0.991,  1.002,  1.012, 1.060, 1.100,
	1.167, 1.221, 1.324, 1.359,  1.473,  1.533,  1.578,  1.622, 1.679, 1.728,
	1.747, 1.747, 1.750, 1.775,  1.785,  1.794,  1.861,  2.117, 2.136, 2.193,
	2.490, 2.860, 2.940, 3.270,  3.360,  3.460,  3.610,  4.060, 5.710, 7.160,
	7.720, 8.470, 9.420, 10.270, 11.450, 12.130, 13.430, 16.5,  19.0,  21.0
};

static const double msq_temp_table[] = {
	2380,  2570,  2680,  2810,  3060,  3210,  3430,  3560,  3660,  3850,
	3930,  3990,  4090,  4300,  4440,  4600,  4870,  5100,  5170,  5290,
	5380,  5480,  5550,  5600,  5660,  5680,  5720,  5770,  5860,  5930,
	6050,  6180,  6280,  6350,  6550,  6670,  6750,  6820,  7020,  7220,
	7400,  7590,  7760,  7910,  8100,  8250,  8600,  8840,  9300,  9700,
	10700, 12300, 14000, 14500, 15700, 16400, 17000, 20600, 26000, 31400,
	33300, 35100, 37900, 39800, 41400, 42900, 45900, 52000, 57500, 61000
};

///// STATIC HELPER FUNCTIONS /////
static double generate_mass(void);
static Boolean should_generate_subdwarf(double mass);
static double generate_subdwarf_metallicity(void);
static double generate_age(double mass, double metallicity, Boolean is_subdwarf);
static double get_radius(double mass, double metallicity, double age);
static int get_surface_temp(double mass, double metallicity, double age, double radius);
static TemperatureClass get_temperature_class(double surface_temp);
static char print_temperature_letter(TemperatureClassLetter letter);
static LuminosityClass get_luminosity_class(double mass, double age, double luminosity);
static char* print_luminosity_class(LuminosityClass luminosity_class);

STAR star_init_default(void)
{
	Star* pStar = (Star*)malloc(sizeof(Star));

	if (pStar == NULL)
	{
		printf("Error 1: Program failed to dynamically allocate memory for uninitialized STAR object\n");
		exit(1);
	}

	pStar->mass = 0.0;
	pStar->metallicity = 0.0;
	pStar->age = 0.0;
	pStar->radius = 0.0;
	pStar->surface_temp = 0;
	pStar->luminosity = 0.0;
	pStar->density = 0.0;
	pStar->class.temperature_class.letter = TC_UNASSIGNED;
	pStar->class.temperature_class.number = 0;
	pStar->class.luminosity_class = LC_UNASSIGNED;

	return pStar;
}

void star_generate_random(STAR hStar)
{
	Star* pStar = (Star*)hStar;
	Boolean is_subdwarf;

	if (pStar == NULL)
	{
		printf("Error 2: Program failed to randomly generate STAR object\n");
		exit(1);
	}

	pStar->mass = generate_mass();

	is_subdwarf = should_generate_subdwarf(pStar->mass);

	if (is_subdwarf == TRUE)	pStar->metallicity = generate_subdwarf_metallicity();
	else						pStar->metallicity = generate_metallicity();

	pStar->age = generate_age(pStar->mass, pStar->metallicity, is_subdwarf);
	pStar->radius = get_radius(pStar->mass, pStar->metallicity, pStar->age);
	pStar->surface_temp = get_surface_temp(pStar->mass, pStar->metallicity, pStar->age, pStar->radius);

	if (is_subdwarf == TRUE)
	{
		pStar->radius *= my_rand_double(0.80, 0.90);
		pStar->surface_temp = (int)(pStar->surface_temp * my_rand_double(1.02, 1.08) + 0.5);
	}

	pStar->luminosity = get_luminosity(pStar->radius, pStar->surface_temp);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class.temperature_class = get_temperature_class(pStar->surface_temp);

	if (is_subdwarf == TRUE)	pStar->class.luminosity_class = SD;
	else						pStar->class.luminosity_class = get_luminosity_class(pStar->mass, pStar->age, pStar->luminosity);
}

void star_print_details(STAR hStar)
{
	Star* pStar = (Star*)hStar;
	char* density_str;
	char* luminosity_str;

	if (pStar == NULL)
	{
		printf("Error 3: Program failed to print details of STAR object\n");
		exit(1);
	}

	if (pStar->class.luminosity_class == SD)
		printf("\tSpectral classification: %s%c%d\n",
			print_luminosity_class(pStar->class.luminosity_class),
			print_temperature_letter(pStar->class.temperature_class.letter),
			pStar->class.temperature_class.number);
	else
		printf("\tSpectral classification: %c%d%s\n",
			print_temperature_letter(pStar->class.temperature_class.letter),
			pStar->class.temperature_class.number,
			print_luminosity_class(pStar->class.luminosity_class));

	printf("\tMass of star: %.3f (Solar masses)\n", pStar->mass);
	printf("\tRadius of star: %.2f (Solar radiuses)\n", pStar->radius);
	density_str = decimal_to_sci_notation(pStar->density);
	printf("\tDensity of star: %s (g/cm^3)\n", density_str);
	printf("\tSurface Temp of star: %d (K)\n", pStar->surface_temp);
	luminosity_str = decimal_to_sci_notation(pStar->luminosity);
	printf("\tLuminosity of star: %s (Solar luminosities)\n", luminosity_str);
	printf("\tMetallicity of star: %.2f (Fe/H)\n", pStar->metallicity);

	if (pStar->age < 0.001)	printf("\tAge of star: %.0f (Years)\n", pStar->age * 1000000000.0);
	else if (pStar->age < 1.0)		printf("\tAge of star: %.2f (Myr)\n", pStar->age * 1000.0);
	else							printf("\tAge of star: %.2f (Gyr)\n", pStar->age);

	free(density_str);
	free(luminosity_str);
}

void star_destroy(STAR* phStar)
{
	if (phStar == NULL || *phStar == NULL)
	{
		return;
	}

	Star* pStar = (Star*)*phStar;
	free(pStar);
	*phStar = NULL;
}

///// STATIC HELPER FUNCTIONS /////

static double generate_mass(void)
{
	const double weights[] = { 0.0,   0.00003, 0.12, 0.61, 3.0,  7.6, 12.0, 76.0 };
	const double generation_mass_table[] = { 200.0, 16,	    2.1,  1.4,  1.04, 0.8, 0.45, 0.079 };
	double prob_table[SIZE(weights)];
	double roll;

	build_cumulative_table(weights, prob_table, SIZE(weights));

	roll = my_rand_double(prob_table[0], prob_table[SIZE(weights) - 1]);

	return my_log_interpolate(roll, prob_table, generation_mass_table, SIZE(weights));
}

static Boolean should_generate_subdwarf(double mass)
{
	double chance;

	if (mass > 2.0)	return FALSE;

	if (mass < 0.45)	chance = 1.0;
	else if (mass < 0.80)	chance = 0.75;
	else if (mass < 1.20)	chance = 0.1;
	else					chance = 0.025;

	return (my_rand_double(0.0, 100.0) < chance) ? TRUE : FALSE;
}

static double generate_subdwarf_metallicity(void)
{
	double fe_h = my_rand_normal(-2.0, 0.45);

	if (fe_h < -7.0)	fe_h = -7.0;
	if (fe_h > -0.5)	fe_h = 0.5;

	return fe_h;
}

static double generate_age(double mass, double metallicity, Boolean is_subdwarf)
{
	double max_age = get_total_lifetime(mass);
	double bias = 0.25 * metallicity + 1.05;

	if (bias < 0.60)	bias = 0.60;
	else if (bias > 1.5)	bias = 1.5;

	if (is_subdwarf == TRUE)	bias -= 0.25;

	if (max_age > 13.8)	max_age = 13.8;

	return my_rand_double_bias(0.0, max_age, bias);
}

static double get_radius(double mass, double metallicity, double age)
{
	const double radius_mass_table[] = { 0.25, 0.3,  0.8,   1.0,   2.0,   5.0,   8.0,   20.0,  40.0, 60.0, 100.0, 200.0 };
	const double sg_multiplier_table[] = { 1.0,  1.5,  2.0,   2.5,   4.0,   5.0,   6.0,   6.0,   4.5,  3.0,  1.5,   1.25 };
	const double giant_multiplier_table[] = { 1.0,  30.0, 100.0, 110.0, 95.0,  80.0,  75.0,  65.0,  25.0, 7.5,  3.0,   2.5 };
	const double late_giant_multiplier_table[] = { 1.0,  80.0, 200.0, 220.0, 190.0, 160.0, 150.0, 130.0, 50.0, 15.0, 6.0,   5.0 };
	const double msq_lifetime = get_msq_lifetime(mass);
	const double post_msq_life_progress = (age - msq_lifetime) / (get_total_lifetime(mass) - msq_lifetime);
	const double sg_multiplier = my_log_interpolate(mass, radius_mass_table, sg_multiplier_table, SIZE(radius_mass_table));
	const double giant_multiplier = my_log_interpolate(mass, radius_mass_table, giant_multiplier_table, SIZE(radius_mass_table));
	const double late_giant_multiplier = my_log_interpolate(mass, radius_mass_table, late_giant_multiplier_table, SIZE(radius_mass_table));
	double msq_life_progress = age / msq_lifetime;
	double radius = my_log_interpolate(mass, msq_mass_table, msq_radius_table, SIZE(msq_mass_table));

	if (age > msq_lifetime)	msq_life_progress = 1.0;

	radius *= (1.0 + 0.15 * msq_life_progress) * (1.0 + 0.015 * metallicity);

	if (age > msq_lifetime)
	{
		if (post_msq_life_progress < 0.40)
			radius *= 1.0 + (post_msq_life_progress / 0.40) * (sg_multiplier - 1.0);
		else if (post_msq_life_progress < 0.90)
			radius *= sg_multiplier * (1.0 + ((post_msq_life_progress - 0.40) / 0.50) * ((giant_multiplier / sg_multiplier) - 1.0));
		else
			radius *= giant_multiplier * (1.0 + ((post_msq_life_progress - 0.90) / 0.10) * ((late_giant_multiplier / giant_multiplier) - 1.0));
	}

	return radius;
}

static int get_surface_temp(double mass, double metallicity, double age, double radius)
{
	const double temp_mass_table[] = { 0.25, 0.3,  0.8,  1.0,  2.0,  5.0,  8.0,  20.0, 40.0, 60.0, 100.0, 200.0 };
	const double min_temp_table[] = { 2400, 2500, 2700, 2800, 2900, 2900, 2900, 2900, 5000, 7000, 16000, 30000 };
	const double cooling_mass_table[] = { 0.25, 0.8,  2.0,  8.0,   20.0, 40.0, 60.0, 200.0 };
	const double cooling_exponent_table[] = { 0.50, 0.55, 0.63, 0.685, 0.72, 0.72, 0.66, 0.64 };
	const double cooling_exponent = my_log_interpolate(mass, cooling_mass_table, cooling_exponent_table, SIZE(cooling_mass_table));
	const double min_temp = my_log_interpolate(mass, temp_mass_table, min_temp_table, SIZE(temp_mass_table));
	const double msq_lifetime = get_msq_lifetime(mass);
	const double expansion = radius / (my_log_interpolate(mass, msq_mass_table, msq_radius_table, SIZE(msq_mass_table)) * 1.20 * (1.0 + 0.015 * metallicity));
	double temp = my_log_interpolate(mass, msq_mass_table, msq_temp_table, SIZE(msq_mass_table)) * (1.0 - 0.015 * metallicity);

	if (age <= msq_lifetime)	temp *= 1.0 + 0.05 * age / msq_lifetime;
	else						temp = min_temp + (temp - min_temp) / pow(expansion, cooling_exponent);

	return (int)(temp + 0.5);
}

static TemperatureClass get_temperature_class(double surface_temp)
{
	TemperatureClass temperature_class;
	double temp = surface_temp;
	int i;

	if (surface_temp < 2380)	temp = 2380;
	else if (surface_temp > 61000)	temp = 60999;

	for (i = 0; i < SIZE(msq_temp_table) - 1; i++)
	{
		if (temp >= msq_temp_table[i] && temp < msq_temp_table[i + 1])
		{
			if (i <= 9)     temperature_class.letter = M;
			else if (i <= 19)    temperature_class.letter = K;
			else if (i <= 29)    temperature_class.letter = G;
			else if (i <= 39)    temperature_class.letter = F;
			else if (i <= 49)    temperature_class.letter = A;
			else if (i <= 59)    temperature_class.letter = B;
			else                 temperature_class.letter = O;

			temperature_class.number = 9 - (i % 10);

			return temperature_class;
		}
	}

	temperature_class.letter = TC_UNASSIGNED;
	temperature_class.number = -67;

	return temperature_class;
}

static char print_temperature_letter(TemperatureClassLetter letter)
{
	switch (letter)
	{
	case M: return 'M';
	case K: return 'K';
	case G: return 'G';
	case F: return 'F';
	case A: return 'A';
	case B:	return 'B';
	case O: return 'O';
	}
	return '?';
}

static LuminosityClass get_luminosity_class(double mass, double age, double luminosity)
{
	const double msq_lifetime = get_msq_lifetime(mass);
	const double msq_life_progress = age / msq_lifetime;
	const double post_msq_life_progress = (age - msq_lifetime) / (get_total_lifetime(mass) - msq_lifetime);

	if (msq_life_progress < 1.00)							return V;

	if (mass < 8.0 && post_msq_life_progress < 0.4)	return IV;
	else if (mass < 2.0)									return III;
	else if (mass < 8.0)									return II;
	else
	{
		if (luminosity < 10000.0)						return II;
		else if (luminosity < 70000.0)						return IB;
		else if (luminosity < 300000.0)						return IAB;
		else if (luminosity < 1000000.0)					return IA;
		else												return IA_PLUS;
	}
}

static char* print_luminosity_class(LuminosityClass luminosity_class)
{
	switch (luminosity_class)
	{
	case SD:	  return "sd";
	case V:		  return "V";
	case IV:	  return "IV";
	case III:	  return "III";
	case II:	  return "II";
	case IB:	  return "Ib";
	case IAB:	  return "Iab";
	case IA:	  return "Ia";
	case IA_PLUS: return "Ia+";
	}
	return "?";
}