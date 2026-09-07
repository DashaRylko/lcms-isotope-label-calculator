# LC-MS/MS Isotope Label Calculator

An interactive **R Shiny application** for calculating molecular isotopic distributions and evaluating potential isotopic interference when selecting isotope-labelled internal standards for **LC-MS/MS analysis**.

## Overview

Stable isotope-labelled internal standards are widely used in quantitative LC-MS/MS analysis. When selecting an isotope-labelled analogue, the mass difference between the analyte and the internal standard should be sufficiently large to minimize interference from naturally occurring isotopes of the analyte.

This application allows users to:

* Enter the elemental composition of a molecule
* Calculate the average molecular mass
* Determine the nominal mass of the most abundant isotopic composition
* Calculate the natural isotopic distribution
* Evaluate the abundance of an isotopic peak at a selected mass difference (Δm)
* Assess potential isotopic interference at the proposed mass shift

The tool can support the selection of an appropriate number of isotope labels for isotope-labelled internal standards used in LC-MS/MS workflows.

## Supported Elements

The application currently supports:

* Carbon (C)
* Hydrogen (H)
* Nitrogen (N)
* Oxygen (O)
* Sulfur (S)
* Chlorine (Cl)
* Bromine (Br)
* Fluorine (F)
* Phosphorus (P)

## How It Works

The user enters the molecular composition of the analyte and specifies a mass difference:

```text
Δm = Mass of labelled standard − Mass of analyte
```

The application calculates the probability of naturally occurring isotopic species of the analyte appearing at the corresponding nominal mass.

For example:

```text
Analyte: M
Labelled internal standard: M + 3
```

The application estimates the abundance of the analyte's natural isotopic contribution at **M + 3**, which can help evaluate the potential for isotopic interference.

## Technologies

* R
* Shiny
* bslib
* tidyverse

## Installation

Clone the repository:

```r
git clone https://github.com/YOUR-USERNAME/lcms-isotope-label-calculator.git
```

Open the project in RStudio and install the required packages:

```r
install.packages(c("shiny", "tidyverse", "bslib"))
```

Run the application:

```r
shiny::runApp()
```


## Author

**Darya Rylko**

Chemist | Biostatistics and Medical Data Analysis | R Programming


