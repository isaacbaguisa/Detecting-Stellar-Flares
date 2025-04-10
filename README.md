# Detecting Stellar Flares

This project focuses on the detection and analysis of stellar flares using time-series data from NASA's Transiting Exoplanet Survey Satellite (TESS) mission. Stellar flares are intense bursts of energy resulting from magnetic reconnection, indicated by a sudden increase in brightness followed by a gradual decay. Understanding these flares is essential for insights into stellar behavior and their potential impacts on surrounding environments.

## Project Overview

The primary objective of this project is to develop and evaluate models to accurately detect stellar flares in photometric flux data.

## Data

The dataset consists of time-series measurements of stellar brightness and associated errors for three stars:

- **TIC 031381302 (Star 1)**
- **TIC 129646813 (Star 2)**
- **TIC 0131799991 (Star 3)**

Each dataset includes the following key variables:

- `time`: Observation timestamp in days (Barycentric Julian Date).
- `PDCSAP_FLUX`: Pre-search Data Conditioning Simple Aperture Photometry flux values, representing the star's brightness.
- `PDCSAP_FLUX_ERR`: Associated errors with the `PDCSAP_FLUX` measurements.

## Methods

The approach to detecting stellar flares involves these key steps:

1. **Exploratory Data Analysis (EDA)**: Understanding the characteristics of the light curves and identifying patterns indicative of flares.
2. **Simulations**: Generating synthetic flares to mimic the brightness and improve/score the models.
3. **Model Development**: Implementing detection algorithms and machine learning methods such as Gaussian Process Regression, Poisson Process Regression, and One-Class Support Vector Machines.
4. **Model Evaluation**: Assessing model performance using metrics such as sensitivity and false positive rates to score flare detection.

## Repository Structure

The repository is organized as follows:

- `Data/`: Contains the datasets used for analysis, and GPR model 
- `EDA/`: Script and report related to exploratory data analysis.
- `Methods/`: Implementation of detection algorithms and models.
- `Model_Evaluation/`: Scripts and results for the evaluation of model performance.
- `Reports/`: All reports from beginning of project outlining steps, and a reproducible example.
- `Simulations/`: Code for simulating brightness and injecting synthetic flares.
- `STA2453_Detecting_Stellar_Flares.Rproj`: R project file.
