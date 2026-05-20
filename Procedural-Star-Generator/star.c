#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "star.h"
#include "utilities.h"
#include "star_properties.h"

typedef enum temperature_class_letter
{
	TC_UNASSIGNED, TC_M, TC_K, TC_G, TC_F, TC_A, TC_B, TC_O, 
} TemperatureClassLetter;

typedef enum luminosity_class
{
	LC_UNASSIGNED, LC_SD, LC_V, LC_IV, LC_III, LC_II, LC_IB, LC_IAB, LC_IA, LC_IA_PLUS
} LuminosityClass;

typedef enum wolf_rayet_emissions_class
{
	WRE_UNASSIGNED, WR_WN, WR_WC, WR_WO
} WolfRayetEmissionClass;

typedef enum star_type
{
	ST_UNASSIGNED, ST_STANDARD, ST_COOL_SUBDWARF, ST_WOLF_RAYET, ST_HOT_SUBDWARF, ST_WHITE_DWARF, ST_NEUTRON_STAR, ST_BLACK_HOLE
} StarType;

typedef enum white_dwarf_atmosphere
{
	AC_UNASSIGNED, AC_DA, AC_DB, AC_DO, AC_DQ, AC_DZ, AC_DC
} WhiteDwarfAtmosphere;

typedef enum neutron_star_type
{
	NST_UNASSIGNED, NST_ISOLATED, NST_PULSAR, NST_MAGNETAR
} NeutronStarType;

typedef struct spectral_class
{
	TemperatureClassLetter temperature_class;
	WolfRayetEmissionClass emission_class;
	WhiteDwarfAtmosphere atmosphere;
	int grade;
	LuminosityClass luminosity_class;
	NeutronStarType neutron_star_type;
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
	StarType type; 
	Boolean black_hole_is_active; 
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
static StarType should_generate_cool_subdwarf(double mass);
static Boolean should_generate_wolf_rayet(double mass, double metallicity, double age); 
static Boolean should_generate_hot_subdwarf(double age); 
static double generate_subdwarf_metallicity(void);
static double generate_age(double mass, double metallicity, StarType type);
static double get_radius(double mass, double metallicity, double age);
static int get_surface_temp(double mass, double metallicity, double age, double radius);
static SpectralClass get_standard_or_subdwarf_spectral_class(double mass, double age, int surface_temp, double luminosity, StarType type);
static SpectralClass get_wolf_rayet_spectral_class(double mass, double age, int surface_temp);
static SpectralClass get_white_dwarf_spectral_class(int surface_temp, double age);
static char print_temperature_letter(TemperatureClassLetter letter);
static char* print_luminosity_class(LuminosityClass luminosity_class);
static char* print_wolf_rayet_emission_class(WolfRayetEmissionClass emission_class);
static char* print_white_dwarf_atmospheric_composition(WhiteDwarfAtmosphere atmosphere); 
static void generate_standard_star_information(Star* pStar); 
static void generate_cool_subdwarf_information(Star* pStar); 
static void generate_wolf_rayet_information(Star* pStar); 
static void generate_hot_subdwarf_information(Star* pStar); 
static void generate_white_dwarf_information(Star* pStar); 

///// INTERFACE FUNCTIONS /////

STAR star_init_default(void)
{
	Star* pStar = (Star*)malloc(sizeof(*pStar));

	if (pStar == NULL)
	{
		fprintf(stderr, "Error: program failed to allocate memory for Star object in star_init_default\n");
		exit(1);
	}

	pStar->mass	= 0.0;
	pStar->metallicity = 0.0;
	pStar->age = 0.0;
	pStar->radius = 0.0;
	pStar->surface_temp = 0;
	pStar->luminosity = 0.0;
	pStar->density = 0.0;
	pStar->class.temperature_class = TC_UNASSIGNED;
	pStar->class.emission_class = WRE_UNASSIGNED;
	pStar->class.grade = -1;
	pStar->class.luminosity_class = LC_UNASSIGNED;
	pStar->type = ST_UNASSIGNED; 

	return pStar;
}

void star_generate_random(STAR hStar)
{
	Star* pStar = (Star*)hStar;
	Boolean is_remnant = FALSE; 
	Boolean is_hot_subdwarf = FALSE; 

	if (pStar == NULL)
	{
		fprintf(stderr, "Error: null parameter(s) passed through star_generate_random\n");
		exit(1);
	}

	if (my_rand_double(0.0, 100.0) < 8.0)	is_remnant = TRUE; 

	pStar->mass = generate_mass();
	pStar->type = should_generate_cool_subdwarf(pStar->mass);

	if (pStar->type == ST_UNASSIGNED)	pStar->type = ST_STANDARD;

	switch (pStar->type)
	{
		case ST_COOL_SUBDWARF:
			generate_cool_subdwarf_information(pStar);
			break;
		case ST_STANDARD:	// Also deals with Wolf-Rayet stars
			generate_standard_star_information(pStar);
			break;
		default:
			fprintf(stderr, "Error: star_generate_random failed\n");
			exit(1);
	}

	if (is_remnant == TRUE)
	{
		pStar->age = my_rand_double(0.0, clamp(13.8 - get_total_lifetime(pStar->mass), 0.0, 13.8));	// Random remnant age based on how long ago the progenitor star could have died.

		if		(pStar->mass < 8.0 && should_generate_hot_subdwarf(pStar->age) == TRUE)	generate_hot_subdwarf_information(pStar); 
		else if (pStar->mass < 8.0)														generate_white_dwarf_information(pStar);
		else if (pStar->mass >= 8.0 && pStar->mass < 20.0)								generate_neutron_star_information(pStar);
		else if (pStar->mass >= 20.0)													generate_black_hole_information(pStar); 
	}
}

void star_print_details(STAR hStar)
{
	Star* pStar = (Star*)hStar;
	char* density_str;
	char* luminosity_str;

	if (pStar == NULL)
	{
		fprintf(stderr, "Error: null parameter(s) passed through star_print_details\n");
		exit(1);
	}

	if (pStar->class.luminosity_class == LC_SD)
		printf("\tSpectral classification: %s%c%d\n",
			print_luminosity_class(pStar->class.luminosity_class),
			print_temperature_letter(pStar->class.temperature_class),
			pStar->class.grade);
	else if (pStar->type == ST_WOLF_RAYET)
		printf("\tSpectral classification: %s%d\n",
			print_wolf_rayet_emission_class(pStar->class.emission_class),
			pStar->class.grade);
	else if (pStar->type == ST_WHITE_DWARF)
		printf("\tSpectral classification: %s%d\n",
			print_white_dwarf_atmospheric_composition(pStar->class.atmosphere),
			pStar->class.grade);
	else
		printf("\tSpectral classification: %c%d%s\n",
			print_temperature_letter(pStar->class.temperature_class),
			pStar->class.grade,
			print_luminosity_class(pStar->class.luminosity_class));

	printf("\tMass of star: %.3f (Solar masses)\n", pStar->mass);
	printf("\tRadius of star: %.2f (Solar radiuses)\n", pStar->radius);
	density_str = decimal_to_sci_notation(pStar->density);

	if (density_str == NULL)
	{
		fprintf(stderr, "Error: null parameter(s) passed through star_print_details\n");
		exit(1);
	}

	printf("\tDensity of star: %s (g/cm^3)\n", density_str);
	printf("\tSurface Temp of star: %d (K)\n", pStar->surface_temp);
	luminosity_str = decimal_to_sci_notation(pStar->luminosity);

	if (luminosity_str == NULL)
	{
		free(density_str);
		fprintf(stderr, "Error: null parameter(s) passed through star_print_details\n");
		exit(1);
	}

	printf("\tLuminosity of star: %s (Solar luminosities)\n", luminosity_str);
	printf("\tMetallicity of star: %.2f (Fe/H)\n", pStar->metallicity);

	if		(pStar->age < 0.001)	printf("\tAge of star: %.0f (Years)\n", pStar->age * 1000000000.0);
	else if (pStar->age < 1.0)		printf("\tAge of star: %.2f (Myr)\n",	pStar->age * 1000.0);
	else							printf("\tAge of star: %.2f (Gyr)\n",	pStar->age);

	free(density_str);
	free(luminosity_str);
}

void star_destroy(STAR* phStar)
{
	Star* pStar;

	if (phStar == NULL || *phStar == NULL)	return;

	pStar = (Star*)*phStar;
	free(pStar);
	*phStar = NULL;
}

///// STATIC HELPER FUNCTIONS /////

// Weighted mass distribution: low-mass stars are far more common than high-mass stars.
static double generate_mass(void)
{
	const double weights[]				 = { 0.0,   0.00003, 0.12, 0.61, 3.0,  7.6, 12.0, 76.0 };
	const double generation_mass_table[] = { 200.0, 16.0,	 2.1,  1.4,  1.04, 0.8, 0.45, 0.079 };
	double prob_table[SIZE(weights)];
	double roll;

	build_cumulative_table(weights, prob_table, SIZE(weights));

	roll = my_rand_double(prob_table[0], prob_table[SIZE(weights) - 1]);

	return my_log_interpolate(roll, prob_table, generation_mass_table, SIZE(weights));
}

// Determines if a star should be a cool subdwarf, probability increases with mass.
static StarType should_generate_cool_subdwarf(double mass)
{
	double chance;

	if (mass > 2.0)	return ST_STANDARD;

	if		(mass < 0.45)	chance = 1.0;
	else if (mass < 0.80)	chance = 0.75;
	else if (mass < 1.20)	chance = 0.1;
	else					chance = 0.025;

	return (my_rand_double(0.0, 100.0) < chance) ? ST_COOL_SUBDWARF : ST_STANDARD;
}

// Detemines if a star should undergo a Wolf-Rayet phase at the end of its life, probability increases with progenitor mass, life progress, and metallicity.
static Boolean should_generate_wolf_rayet(double mass, double metallicity, double age)
{
	const double life_progress = age / get_total_lifetime(mass);
	double chance;

	if (mass < 20.0)	return FALSE;

	// WR stars are modeled as a late-life phase, not an early main-sequence state.
	if (life_progress < 0.65)	return FALSE;

	if		(mass < 30.0)	chance = 2.0;
	else if (mass < 60.0)	chance = 10.0;
	else if (mass < 100.0)	chance = 20.0;
	else					chance = 40.0;

	// Higher metallicity strengthens stellar winds, making WR formation more likely.
	chance *= clamp(1.0 + 0.25 * metallicity, 0.5, 1.5);

	if (life_progress > 0.85)	chance *= 1.5;

	return (my_rand_double(0.0, 100.0) < chance) ? TRUE : FALSE;
}

// Determines if a star should undergo a hot subdwarf phase, this phase is rare and unlikely
static Boolean should_generate_hot_subdwarf(double age)
{
	if (age > 0.10)	return FALSE;	// Hot subdwarf phase only allowed within 100 Myr after progenitor death.

	return (my_rand_double(0.0, 100.0) < 2.0) ? TRUE : FALSE;
}

// Generate low metallicity Fe/H values for metal-poor subdwarf stars.
static double generate_subdwarf_metallicity(void)
{
	const double fe_h = my_rand_normal(-2.0, 0.45);

	return clamp(fe_h, -7.0, -0.5);
}

// Metal-poor and subdwarf stars are biased toward older ages, younger stars are biased oppositely.
static double generate_age(double mass, double metallicity, StarType type)
{
	double max_age = get_total_lifetime(mass);
	double bias = clamp((0.25 * metallicity + 1.05), 0.6, 1.5);

	if (type == ST_COOL_SUBDWARF)	bias -= 0.25;

	if (max_age > 13.8)	max_age = 13.8;

	return my_rand_double_bias(0.0, max_age, bias);
}

// Calculates a reasonable radius for the star depending on its mass, metallicity, and age. 
static double get_radius(double mass, double metallicity, double age)
{
	// Post-main-sequence radius multipliers by mass, These are tuned approximation tables, not strict stellar evolution tracks, and are heavily simplified from reality.
	const double radius_mass_table[]		   = { 0.25,  0.3,  0.8,   1.0,   2.0,   5.0,   8.0,   20.0,  40.0, 60.0, 100.0, 200.0 };
	const double sg_multiplier_table[]		   = { 1.01,  1.5,  2.0,   2.5,   4.0,   5.0,   6.0,   6.0,   4.5,  3.0,  1.5,   1.25 };
	const double giant_multiplier_table[]	   = { 1.01,  40.0, 100.0, 110.0, 95.0,  80.0,  75.0,  65.0,  25.0, 7.5,  3.0,   2.5 };
	const double late_giant_multiplier_table[] = { 1.01,  80.0, 200.0, 220.0, 190.0, 160.0, 150.0, 130.0, 50.0, 15.0, 6.0,   5.0 };
	const double msq_lifetime = get_msq_lifetime(mass);
	const double post_msq_life_progress = clamp((age - msq_lifetime) / (get_total_lifetime(mass) - msq_lifetime), 0.0, 1.0);
	const double sg_multiplier = my_log_interpolate(mass, radius_mass_table, sg_multiplier_table, SIZE(radius_mass_table));
	const double giant_multiplier = my_log_interpolate(mass, radius_mass_table, giant_multiplier_table, SIZE(radius_mass_table));
	const double late_giant_multiplier = my_log_interpolate(mass, radius_mass_table, late_giant_multiplier_table, SIZE(radius_mass_table));
	const double msq_life_progress = clamp((age / msq_lifetime), 0.0, 1.0);
	double radius = my_log_interpolate(mass, msq_mass_table, msq_radius_table, SIZE(msq_mass_table));

	// Main-sequence stars expand slightly as they age, metallicity gives a small radius adjustment.
	radius *= (1.0 + 0.15 * msq_life_progress) * (1.0 + 0.015 * metallicity);

	if (age > msq_lifetime)
	{
		// Post-main-sequence radius growth is split into subgiant, giant, and late-giant phases.
		if (post_msq_life_progress < 0.40)
			radius *= 1.0 + (post_msq_life_progress / 0.40) * (sg_multiplier - 1.0);
		else if (post_msq_life_progress < 0.90)
			radius *= sg_multiplier * (1.0 + ((post_msq_life_progress - 0.40) / 0.50) * ((giant_multiplier / sg_multiplier) - 1.0));
		else
			radius *= giant_multiplier * (1.0 + ((post_msq_life_progress - 0.90) / 0.10) * ((late_giant_multiplier / giant_multiplier) - 1.0));
	}

	return radius;
}

// Minimum evolved-star temperatures and late-giant radius multipliers are tuned approximation tables.
static int get_surface_temp(double mass, double metallicity, double age, double radius)
{
	// Minimum evolved-star temperatures are tuned approximation tables.
	const double temp_mass_table[]			   = { 0.25,  0.3,  0.8,   1.0,   2.0,   5.0,   8.0,   20.0,  40.0, 60.0, 100.0, 200.0 };
	const double min_temp_table[]			   = { 2400,  2500, 2700,  2800,  2900,  2900,  3000,  3200,  8000, 11000, 16000, 30000 };
	const double late_giant_multiplier_table[] = { 1.01,  80.0, 200.0, 220.0, 190.0, 160.0, 150.0, 130.0, 50.0, 15.0, 6.0,   5.0 };
	const double min_temp = my_log_interpolate(mass, temp_mass_table, min_temp_table, SIZE(temp_mass_table));
	const double msq_lifetime = get_msq_lifetime(mass);
	const double msq_radius = my_log_interpolate(mass, msq_mass_table, msq_radius_table, SIZE(msq_mass_table));
	const double start_radius = msq_radius * 1.15 * (1.0 + 0.015 * metallicity);
	const double max_radius = start_radius * my_log_interpolate(mass, temp_mass_table, late_giant_multiplier_table, SIZE(temp_mass_table));
	const double radius_progress = clamp(((radius - start_radius) / (max_radius - start_radius)), 0.0, 1.0);	// Treat post-main-sequence cooling as progress between the starting subgiant radius and the maximum late-giant/supergiant radius.
	double temp = my_log_interpolate(mass, msq_mass_table, msq_temp_table, SIZE(msq_mass_table)) * (1.0 - 0.015 * metallicity);
	double start_temp = temp * 1.05;	// Slight temperature bump at the start of post-main-sequence evolution.
	double cooling_mass = mass; 
	double temp_progress;
	double base;

	if (age <= msq_lifetime)	temp *= 1.0 + 0.05 * age / msq_lifetime;
	else
	{
		if (mass > 20.0)	cooling_mass = 20.0 / pow(mass / 20.0, 3.0);	// Reduces high-mass cooling strength so very massive supergiants do not cool too aggressively.

		base = pow((4.0 / 3.0) * M_PI * cooling_mass, 3.0);								// Controls how quickly temperature approaches min_temp as radius expands.
		temp_progress = (1.0 - pow(base, -radius_progress)) / (1.0 - pow(base, -1.0));	// Nonlinear 0-1 cooling curve.
		temp = start_temp + (min_temp - start_temp) * temp_progress;
	}

	return (int)(temp + 0.5);
}

// Creates the spectral classification for standard stars or subdwarfs. 
static SpectralClass get_standard_or_subdwarf_spectral_class(double mass, double age, int surface_temp, double luminosity, StarType type)
{
	SpectralClass spectral_class;
	const double temp = clamp(surface_temp, 2380, 60999);
	const double msq_lifetime = get_msq_lifetime(mass);
	const double msq_life_progress = age / msq_lifetime;
	const double post_msq_life_progress = (age - msq_lifetime) / (get_total_lifetime(mass) - msq_lifetime);
	int i;

	spectral_class.emission_class = WRE_UNASSIGNED;
	spectral_class.atmosphere = AC_UNASSIGNED; 
	spectral_class.neutron_star_type = NST_UNASSIGNED;
	
	// Creates the luminosity class for standard, and cool subdwarf stars and in the future hot subdwarfs. 
	if		(type == ST_COOL_SUBDWARF && msq_life_progress < 1.00)	
															spectral_class.luminosity_class = LC_SD;
	// Main-sequence stars not in the post main-sequence are luminosity class V.
	else if (msq_life_progress < 1.00)						spectral_class.luminosity_class = LC_V;
	// Post-main-sequence class is approximated from mass, evolution stage, and luminosity.
	else if (mass < 8.0 && post_msq_life_progress < 0.4)	spectral_class.luminosity_class = LC_IV;
	else if (mass < 2.0)									spectral_class.luminosity_class = LC_III;
	else if (mass < 8.0)									spectral_class.luminosity_class = LC_II;
	else if (luminosity < 10000.0)							spectral_class.luminosity_class = LC_II;
	else if (luminosity < 70000.0)							spectral_class.luminosity_class = LC_IB;
	else if (luminosity < 300000.0)							spectral_class.luminosity_class = LC_IAB;
	else if (luminosity < 1000000.0)						spectral_class.luminosity_class = LC_IA;
	else													spectral_class.luminosity_class = LC_IA_PLUS;

	for (i = 0; i < SIZE(msq_temp_table) - 1; i++)
	{
		// Gets temperature class from the msq_temperature table by seeing what row the temp falls in, indirectly lol
		if (temp >= msq_temp_table[i] && temp < msq_temp_table[i + 1])
		{
			if		(i <= 9)	spectral_class.temperature_class = TC_M;
			else if (i <= 19)	spectral_class.temperature_class = TC_K;
			else if (i <= 29)	spectral_class.temperature_class = TC_G;
			else if (i <= 39)	spectral_class.temperature_class = TC_F;
			else if (i <= 49)	spectral_class.temperature_class = TC_A;
			else if (i <= 59)	spectral_class.temperature_class = TC_B;
			else				spectral_class.temperature_class = TC_O;

			spectral_class.grade = 9 - (i % 10);

			return spectral_class;
		}
	}

	fprintf(stderr, "Error: get_standard_or_subdwarf_spectral_class failed\n");
	exit(1); 
}

// Creates the spectral classification for Wolf-Rayet stars. 
static SpectralClass get_wolf_rayet_spectral_class(double mass, double age, int surface_temp)
{
	SpectralClass spectral_class;
	double life_progress = age / get_total_lifetime(mass);

	spectral_class.temperature_class = TC_UNASSIGNED;
	spectral_class.atmosphere = AC_UNASSIGNED;
	spectral_class.luminosity_class = LC_UNASSIGNED;
	spectral_class.neutron_star_type = NST_UNASSIGNED;

	if		(life_progress < 0.80 && mass >= 50.0 && surface_temp >= 70000)	spectral_class.emission_class = WR_WN;
	else if (life_progress < 0.95 || surface_temp < 100000)					spectral_class.emission_class = WR_WC;
	else																	spectral_class.emission_class = WR_WO;

	// Uses equal intervals to give a temperature, 30,000 or lower being 11, and 200,000 or higher being a 1
	spectral_class.grade = clamp_int((11 - (int)(((surface_temp - 30000) / (200000.0 - 30000.0)) * 10.0 + 0.5)), 1, 11);

	if (spectral_class.emission_class == WR_WO)	spectral_class.grade = clamp_int(spectral_class.grade, 1, 6);

	return spectral_class;
}

// Creates the spectral classification for white dwarfs. 
static SpectralClass get_white_dwarf_spectral_class(int surface_temp, double age)
{
	SpectralClass spectral_class;
	double roll = my_rand_double(0.0, 100.0);

	spectral_class.temperature_class = TC_UNASSIGNED;
	spectral_class.emission_class = WRE_UNASSIGNED;
	spectral_class.luminosity_class = LC_UNASSIGNED;
	spectral_class.neutron_star_type = NST_UNASSIGNED;

	// Not hard science here, reasonable approximations
	if		(surface_temp >= 45000 && roll < 60.0)	spectral_class.atmosphere = AC_DO;
	else if (roll < 75.0)							spectral_class.atmosphere = AC_DA;
	else if (roll < 87.0)							spectral_class.atmosphere = AC_DB;
	else if (surface_temp < 12000 && roll < 93.0)	spectral_class.atmosphere = AC_DQ;
	else if (surface_temp < 10000 && roll < 97.0)	spectral_class.atmosphere = AC_DZ;
	else											spectral_class.atmosphere = AC_DC;

	spectral_class.grade = clamp_int((int)(50400.0 / surface_temp + 0.5), 1, 9);

	return spectral_class;
}

// Takes the enum and creates an equivalent printable format
static char print_temperature_letter(TemperatureClassLetter letter)
{
	switch (letter)
	{
		case TC_M:	return 'M';
		case TC_K:	return 'K';
		case TC_G:	return 'G';
		case TC_F:	return 'F';
		case TC_A:	return 'A';
		case TC_B:	return 'B';
		case TC_O:	return 'O';
	}

	return '?';
}

// Takes the enum and creates an equivalent printable format
static char* print_luminosity_class(LuminosityClass luminosity_class)
{
	switch (luminosity_class)
	{
		case LC_SD:			return "sd";
		case LC_V:			return "V";
		case LC_IV:			return "IV";
		case LC_III:		return "III";
		case LC_II:			return "II";
		case LC_IB:			return "Ib";
		case LC_IAB:		return "Iab";
		case LC_IA:			return "Ia";
		case LC_IA_PLUS:	return "Ia+";
	}

	return "?";
}

// Takes the enum and creates an equivalent printable format
static char* print_wolf_rayet_emission_class(WolfRayetEmissionClass emission_class)
{
	switch (emission_class)
	{
		case WR_WN: return "WN";
		case WR_WC: return "WC";
		case WR_WO: return "WO";
	}

	return "?"; 
}

static char* print_white_dwarf_atmospheric_composition(WhiteDwarfAtmosphere atmosphere)
{
	switch (atmosphere)
	{
		case AC_DA: return "DA"; 
		case AC_DB: return "DB";
		case AC_DO: return "DO";
		case AC_DQ: return "DQ";
		case AC_DZ: return "DZ";
		case AC_DC: return "DC";
	}

	return "?"; 
}

// Generates the star information for a standard star.
static void generate_standard_star_information(Star* pStar)
{
	pStar->metallicity = generate_metallicity();
	pStar->age = generate_age(pStar->mass, pStar->metallicity, pStar->type);
	pStar->radius = get_radius(pStar->mass, pStar->metallicity, pStar->age);
	pStar->surface_temp = get_surface_temp(pStar->mass, pStar->metallicity, pStar->age, pStar->radius);
	
	// Checks to see if where the star currently is, if its a candidate to become a Wolf-Rayet, if so pivot the generation towards Wolf-Rayet.
	if (should_generate_wolf_rayet(pStar->mass, pStar->metallicity, pStar->age) == TRUE)
	{
		pStar->type = ST_WOLF_RAYET;
		generate_wolf_rayet_information(pStar);
		return;
	}

	pStar->luminosity = get_luminosity(pStar->radius, pStar->surface_temp);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class = get_standard_or_subdwarf_spectral_class(pStar->mass, pStar->age, pStar->surface_temp, pStar->luminosity, pStar->type);
}

// Generates the star information for a cool subdwarf
static void generate_cool_subdwarf_information(Star* pStar)
{
	pStar->metallicity = generate_subdwarf_metallicity();
	pStar->age = generate_age(pStar->mass, pStar->metallicity, pStar->type);
	pStar->radius = get_radius(pStar->mass, pStar->metallicity, pStar->age);
	pStar->surface_temp = get_surface_temp(pStar->mass, pStar->metallicity, pStar->age, pStar->radius);
	pStar->radius *= my_rand_double(0.80, 0.90);											// Applies cool_subdwarf radius modifier
	pStar->surface_temp = (int)(pStar->surface_temp * my_rand_double(1.02, 1.08) + 0.5);	// Applies cool_subdwarf surface temp modifier
	pStar->luminosity = get_luminosity(pStar->radius, pStar->surface_temp);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class = get_standard_or_subdwarf_spectral_class(pStar->mass, pStar->age, pStar->surface_temp, pStar->luminosity, pStar->type); 
}

// Generates the star information for a Wolf-Rayet star
static void generate_wolf_rayet_information(Star* pStar)
{
	const double progenitor_mass = pStar->mass;	// Store the original mass before WR mass loss so classification/life progress still use the progenitor mass, not the stripped Wolf-Rayet mass.
	const double wr_progress = ((clamp(pStar->age / get_total_lifetime(progenitor_mass), 0.65, 1.0)) - 0.65) / 0.35;

	pStar->mass *= (1.0 - (1.0/2.0 * wr_progress));																		// Gradually strips mass as the star evolves through the WR phase.
	// Rebuilds a compact WR radius while allowing rare massive WN stars to remain large.
	pStar->radius = clamp(((4.0 + 0.18 * pStar->mass + my_rand_normal(0.0, 3.0)) * (1.0 - 0.35 * wr_progress)), 1.5, (25.0 + 0.10 * progenitor_mass));
	pStar->surface_temp = (int)(clamp(pStar->surface_temp * (1.0 + (1.75 * wr_progress)), 30000.0, 200000.0) + 0.5);	// Increases by up to 175%, reaching 2.75x its starting temperature.
	pStar->luminosity = get_luminosity(pStar->radius, pStar->surface_temp);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class = get_wolf_rayet_spectral_class(progenitor_mass, pStar->age, pStar->surface_temp);
}

// Generates the star information for a hot subdwarf
static void generate_hot_subdwarf_information(Star* pStar)
{
	const double progenitor_mass = pStar->mass;	// Store the original mass of the star before it gets modified
	double start_temp;

	pStar->type = ST_HOT_SUBDWARF;
	// Scales hot subdwarf mass from progenitor mass, with slight random variation.
	pStar->mass = clamp(((0.079 + ((progenitor_mass - 0.079) / (8.0 - 0.079)) * (1.2 - 0.079)) * my_rand_double(0.90, 1.10)), 0.079, 1.2);
	pStar->radius = 0.10 + ((pStar->mass - 0.079) / (1.2 - 0.079)) * (0.35 - 0.10);	// Scales hot subdwarf radius from its final mass, with slight random variation.
	start_temp = 20000.0 + clamp((pStar->mass - 0.079) / (1.2 - 0.079), 0.0, 1.0) * (100000.0 - 20000.0);	// Sets the starting temperature for a hot subdwarf. 
	// Hot subdwarfs start extremely hot at young remnant age, then cool sharply as age increases.
	pStar->surface_temp = (int)(20000.0 + (start_temp - 20000.0) * (1.0 - clamp(pStar->age / 0.10, 0.0, 1.0)) + 0.5);
	pStar->luminosity = get_luminosity(pStar->radius, pStar->surface_temp);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class = get_standard_or_subdwarf_spectral_class(pStar->mass, pStar->age, pStar->surface_temp, pStar->luminosity, pStar->type);
}

// Generates the star informatiom for a white dwarf
static void generate_white_dwarf_information(Star* pStar)
{
	const double progenitor_mass = pStar->mass;

	pStar->type = ST_WHITE_DWARF;
	pStar->mass = clamp(((0.125 * progenitor_mass + 0.40) * my_rand_double(0.96, 1.04)), 0.05, 1.44);
	// White dwarf radius shrinks as mass increases, approaching zero near the Chandrasekhar limit (1.44 Solar Masses).
	pStar->radius = clamp((0.012 * sqrt(pow(1.44 / pStar->mass, 2.0 / 3.0) - pow(pStar->mass / 1.44, 2.0 / 3.0))), 0.003, 0.035);
	pStar->luminosity = clamp(0.1 / pow(pStar->age + 0.01, 5.0 / 7.0), 0.00001, 100.0);	// White dwarf luminosity decreases over time using a simplified Mestel cooling law approximation.
	// Uses the Stefan-Boltzmann law to derive surface temperature from luminosity and radius.
	pStar->surface_temp = (int)(clamp((5778.0 * pow(pStar->luminosity / (pStar->radius * pStar->radius), 0.25)), 3000.0, 150000.0) + 0.5);
	pStar->density = get_density(pStar->mass, pStar->radius);
	pStar->class = get_white_dwarf_spectral_class(pStar->surface_temp, pStar->age);
}