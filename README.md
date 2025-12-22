# THESIS codebase

This repository contains the full analysis pipeline for my thesis on  
**value-based decision-making using eye-tracking and drift-diffusion modelling**.

## Project overview

The project consists of two studies:

### Study 1 – Garcia quasi-replication  
All files and folders prefixed with `Garcia`

### Study 2 – Garcia quasi-replication with different OV levels during learning  
All files and folders prefixed with `OV`

Both studies follow a similar processing and analysis structure.

---

## Repository structure

### `Mat_to_Pandas`
Code for converting raw MATLAB files into CSV files

### `Eye_tracking_code`
Eye-tracking preprocessing code, including:
- parsing EyeLink files
- extracting fixations
- assigning ROIs
- computing gaze metrics

### `Garcia_Analysis` (Study 1)
Analysis pipeline for the Garcia replication:
- behavioural dataset construction
- eye-tracking integration
- computation of aDDM regressors and required columns
- behavioural and gaze analyses

### `OV_Analysis` (Study 2)
Analysis pipeline for the OV experiment:
- behavioural and eye-tracking preprocessing
- attentional regressor construction
- analyses and visualisations
