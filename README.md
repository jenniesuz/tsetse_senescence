# <span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">tsetse_senescence code repository</span>

This repository contains all files needed to run the analyses in:

> <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Lord JS, Leyland R, Haines LR, Barreaux A, Bonsall M, Torr SJ, English S. Title to be decided. *Manuscript in preparation*.

The code is made available under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>. You are free to reuse this code provided that you give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use. Giving appropriate credit includes citation of the above publication and providing a link to this repository:

<a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/jenniesuz/tsetse_senescence" rel="dct:source">https://github.com/jenniesuz/tsetse_senescence</a>

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />

We suggest you open [this project file](tsetse_senescence.Rproj) file in [R Studio](rstudio.org). The files that form part of this project file are detailed below. We suggest that you view and run the files in the following order:
- [**r1_db_funcs.R**](r1_db_funcs.R)
- [**r1_queries.R**](r1_queries.R) 
- [**r2_motherDeaths.R**](r2_motherDeaths.R)
- [**r3_abortions.R**](r3_abortions.R) 
- [**r4_pupalWetWeight.R**](r4_pupalWetWeight.R)
- [**r5_offspringEmergence.R**](r5_offspringEmergence.R)
- [**r6_offspringSurvival.R**](r5_offspringSurvival.R)

These files plus the database are explained briefly below:

## Data
- **tsetse2018.db** - this sqlite database file contains the raw data

## Data preparation
- **r1_queries.R** - queries tsetse2018.db and retrieves the required data to produce the figures

## Mother deaths
- **r2_motherDeaths.R** - analysis of mother deaths and Figure 2 of the manuscript

## Larval abortion
- **r3_abortions.R** - analysis of the probability of abortion and Figure 3 of the manuscript

## Pupal wet weight
- **r4_pupalWetWeight.R** - analysis of pupal wet weight and Figure 4 of the manuscript

## Offspring probability of emergence and survival
- **r5_offspringEmergence.R** - analysis of the probability of offspring emergence and Figure 5a
- **r6_offspringSurvival.R** - analysis of the survival of offspring and Figure 5b
