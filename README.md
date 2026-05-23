# Procedural Star Generator [BETA]

A procedural star generation engine written in C that creates scientifically inspired stars spanning from tiny red dwarfs to massive blue hypergiants.

The project simulates stellar mass, radius, luminosity, surface temperature, metallicity, density, age, spectral classification, luminosity classes, stellar evolution, and metal-poor subdwarf populations using astrophysical approximations, stellar evolution modeling, and randomized generation techniques.

The project aims to produce scientifically believable stars while still being suitable for procedural worldbuilding, games, and simulation projects.

---

## Installation & How To Compile

**NOTE:** It is highly recommended that you download the program from the **Releases** tab, preferably the latest release, rather than downloading directly from the repository code view. Official releases are stable, usable versions of the program, while files in the repository code view may be experimental, unfinished, or buggy.

### Visual Studio: 

1. Open `Procedural-Star-Generator.slnx` in Visual Studio.
2. Select the desired build configuration (`Debug` or `Release`).
3. Build the project with `Ctrl + Shift + B`.
4. Run the program with `Ctrl + F5`.

### GCC: 

Clone the repository: 

```bash
git clone https://github.com/yourname/procedural-star-generator.git
cd procedural-star-generator
```

Compile with GCC:

```bash
gcc main.c -o star_generator -lm
```

Run the program:

```bash
./star_generator
```

---

## Usage Example

Example generated output:

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
```

--- 

## License

This project is source-available under the Bower Source License (BSL) v1.0

You may use, study, modify, and distribute this software with attribution, but you may not directly monetize or commercially exploit the software or derivative works without explicit permission from the author. 

See the LICENSE file for full details.

---

Author: Jason Bower

Contact Details: `jasonbower.mail@gmail.com`

Project Started: 3/17/2026
