# Procedural-Star-Generator

A procedural stellar generation engine written in C that creates scientifically-inspired stars spanning from tiny red dwarfs to massive blue hypergiants.

The project simulates stellar mass, radius, luminosity, surface temperature, metallicity, density, age, spectral classification, luminosity classes, stellar evolution, and metal-poor subdwarf populations using interpolation tables and probabilistic models.

This project is intended to become part of a larger procedural star system generator capable of generating realistic stars, stellar remnants, exotic stellar objects, planets, and full solar systems.

--- 

## Features 

- Procedural generation of scientifically-inspired stars from 0.079–200 solar masses
- Support for:
    - Main-sequence stars
    - Subgiants
    - Giants
    - Bright giants
    - Supergiants
    - Hypergiants
    - Cool metal-poor subdwarfs
- Realistic stellar properties:
    - Mass
    - Radius
    - Surface temperature
    - Luminosity
    - Density
    - Metallicity
    - Age
    - Stellar lifetime
- Spectral classification system:
    - O, B, A, F, G, K, M temperature classes
    - sd, V, IV, III, II, Ib, Iab, Ia, Ia+ luminosity classes
    - 0–9 temperature grading
- Metallicity-based population variation
- Stellar evolution from the main sequence through late giant phases
- Logarithmic interpolation system for smooth scaling across extreme stellar ranges
- Probabilistic generation models for realistic stellar distributions
- Randomized seed generation for unique outputs every run

---

## Future Updates / Known Issues

### Planned

- Remnant star generation
    - White dwarfs
    - Neutron stars
    - Hot "Blue" subdwarfs
- Exotic star/object generation
    - Wolf Rayet Stars
    - Blue stragglers
    - Black holes
- Brown Dwarfs

### Limitations

- Generated stars may occasionally produce unrealistic cases
- Stellar evolution is simplified and does not fully model real stellar physics
- Hypergiants are represented through very high luminosity classes rather than a dedicated hypergiant model
- No user input or configuration options are available yet
- Some values are tuned for believable procedural output rather than strict scientific accuracy

_This is part of a planned project for a planetary system generator_

---

## Demo

### Example Output

```
SEED: 1828599120

Star 1:
        Spectral classification: M6V
        Mass of star: 0.141 (Solar masses)
        Radius of star: 0.18 (Solar radiuses)
        Density of star: 35.42 (g/cm^3)
        Surface Temp of star: 2990 (K)
        Luminosity of star: 2.28 * 10^-3 (Solar luminosities)
        Metallicity of star: 0.01 (Fe/H)
        Age of star: 11.73 (Gyr)
Star 2:
        Spectral classification: M2V
        Mass of star: 0.474 (Solar masses)
        Radius of star: 0.48 (Solar radiuses)
        Density of star: 6.15 (g/cm^3)
        Surface Temp of star: 3652 (K)
        Luminosity of star: 0.04 (Solar luminosities)
        Metallicity of star: -0.37 (Fe/H)
        Age of star: 12.79 (Gyr)
Star 3:
        Spectral classification: F7V
        Mass of star: 1.176 (Solar masses)
        Radius of star: 1.26 (Solar radiuses)
        Density of star: 0.83 (g/cm^3)
        Surface Temp of star: 6328 (K)
        Luminosity of star: 2.29 (Solar luminosities)
        Metallicity of star: -0.83 (Fe/H)
        Age of star: 880.01 (Myr)
```

--- 

## How to Build / Run

### Requirements

- C compiler with C89/C90 support
- Microsoft Visual Studio 2022/2026 or GCC / Clang

No external libraries are required

### Build Instructions (Visual Studio) 

1. Clone the repository: `git clone https://github.com/jasonbower/Procedural-Star-Generator.git`
2. Open the .sln file in Visual Studio.
3. Build the project: `Build -> Build Solution`
4. Run the executable.

### Build Instructions (GCC)

```
gcc *.c -o procedural_star_generator -lm
./procedural_star_generator
```

_Important note: because the program uses `<math.h>` functions such as `pow`, `sqrt`, `log`, and `cos`, GCC may require `-lm` at the end. Visual Studio handles this automatically, and is preferred for this application_

--- 

## Usage

When the program starts, it automatically generates a randomized set of stars.

For each generated star, the program displays:

- Spectral classification
- Mass
- Radius
- Density
- Surface temperature
- Luminosity
- Metallicity
- Age

The generator currently:

- Creates a randomized seed for each run
- Generates stars across a wide stellar mass range
- Applies metallicity and age variation
- Calculates stellar properties using interpolation tables
- Classifies each star by temperature class and luminosity class
- Prints the generated star data to the console

This prototype currently runs as a console-based generator and does not require user input.

---

## Project Structure

```
Procedural-Star-Generator/
├── main.c
├── star_standard.h
├── star_standard.c
├── star_properties.h
├── star_properties.c
├── utilities.h
├── utilities.c
└── README.md
```

## Important Components

- `main.c`
    - Program entry point and test driver
- `star_standard.*`
    - "Standard" star generation logic
- `star_properties.*`
    - Helper functions for all types of stars
- `utilities.*`
    - General helper functions across all sorts of matters

--- 

## Concepts Used / What I Learned

- Opaque object design
- Object-oriented design concepts
- Applying math formulas to code
- Random procedural generation logic
- Generic programming in C
- Modular software design

---

## Technical Highlights

- Implemented scientifically plausible stars
- Simulated stellar evolution
- Includes custom random generation helpers for uniform, biased, and normal distributions
- Uses probabilistic weighted tables for stellar mass generation
- Applies metallicity-based variation to stellar age, radius, and temperature

--- 

## Acknowledgments / Credits

This project was inspired by real-world stellar astronomy, stellar classification systems, and publicly available astrophysical reference data.

Special thanks to the astronomy and astrophysics communities whose research and educational resources helped inform the procedural generation models used in this project.

---

## License

A formal license has not yet been selected for this project.

Until a license is chosen, all rights are reserved by the author.  
Please contact me if you are interested in using or contributing to this project.

--- 

Author: Jason Bower

Project Started: 3/17/26
