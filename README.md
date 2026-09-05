# Monte Carlo Option Pricing & Variance Reduction

A C++ implementation and numerical study of Monte Carlo methods for option pricing, with a particular focus on variance reduction techniques and their application to European and Asian options.

## Overview

This project investigates the use of Monte Carlo simulation for pricing derivatives under the Black–Scholes framework. The primary objective is to study how classical variance reduction techniques can improve the efficiency and accuracy of Monte Carlo option pricing.

The project implements and compares:

- Standard Monte Carlo simulation
- Antithetic variates
- Control variates
- Combined antithetic + control variates
- Black–Scholes analytical pricing for validation
- Numerical estimation of option Greeks

The project also extends the framework to Asian options, where Monte Carlo simulation is particularly useful because arithmetic Asian options do not generally admit a simple closed-form pricing formula.

---

## Reproducibility

All numerical experiments use a **fixed pseudorandom seed** to ensure reproducibility.

Using a fixed seed ensures that the same sequence of pseudorandom numbers is generated each time the experiments are run, allowing numerical results and plots to be reproduced consistently.

The experiments therefore provide a controlled basis for comparing the different Monte Carlo estimators.

---

## Project Structure

```text
Monte-Carlo-Option-Pricing/
│
├── src/
│   ├── asian_main.cpp
│   └── asian_main.cpp
│
├── data/
│   ├── Asian_mc_results.csv
│   ├── Asian_mc_results_antithetic.csv
│   ├── Asian_mc_results_antithetic_CV.csv
│   ├── mc_results.csv
│   ├── mc_results_antithetic.csv
│   └── mc_results_antithetic_CV.csv
│
├── plots/
│   ├── European.png
│   └── Asian.png
│
├── report/
│   └── Monte_Carlo_Option_Pricing.pdf
│
└── README.md
