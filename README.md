# A-SRUKF-SOC-Estimation-for-Lithium-Ion-Batteries
This project presents a comprehensive methodology for accurate State of Charge (SOC) estimation of lithium-ion batteries using an Adaptive Square Root Unscented Kalman Filter (A-SRUKF), validated within a Software in the loop(SIL) configuration.

# Adaptive Square Root Unscented Kalman Filter Based State of Charge Estimation for Lithium-Ion Batteries
## Project Overview

This project presents a complete framework for lithium-ion battery State of Charge (SOC) estimation using an Adaptive Square Root Unscented Kalman Filter (ASRUKF) integrated with a second-order Equivalent Circuit Model (2-RC ECM).

The project focuses heavily on accurate estimation of battery internal RC parameters using an Extended Kalman Filter (EKF) from Hybrid Pulse Power Characterization (HPPC) test data.

The extracted parameters:
- R0
- R1
- R2
- C1
- C2

are dynamically estimated region-wise across different SOC intervals and integrated into lookup tables for real-time SOC estimation.

The complete framework was implemented in MATLAB/Simulink and validated in a Software-in-the-Loop (SIL) environment.

# EKF-Based RC Parameter Estimation
## Why RC Parameter Estimation is Important

Lithium-ion batteries exhibit highly nonlinear electrical behavior during charging and discharging operations.

A simple ideal voltage source model cannot accurately represent:
- transient voltage response
- polarization effects
- diffusion dynamics
- internal resistance variations
- SOC-dependent voltage behavior

To accurately model these dynamics, a second-order Equivalent Circuit Model (2-RC ECM) was used.

The ECM consists of:
- Open Circuit Voltage (OCV)
- Ohmic resistance (R0)
- Fast dynamic RC branch (R1-C1)
- Slow diffusion RC branch (R2-C2)

These parameters are not constant.

They vary continuously with:
- State of Charge (SOC)
- Temperature
- Load current
- Battery aging

Accurate estimation of these RC parameters is critical because they directly influence:
- terminal voltage prediction
- SOC estimation accuracy
- filter convergence
- battery state tracking
- real-time Battery Management System (BMS) performance

  ## Second Order Equivalent Circuit Model (2-RC ECM)

The battery is modeled using:
Vt = OCV(SOC) - R0·I - V1 - V2
Where:
- Vt = terminal voltage
- OCV(SOC) = open circuit voltage
- R0 = ohmic resistance
- V1 = voltage across first RC branch
- V2 = voltage across second RC branch
- I = battery current

  ## Physical Meaning of RC Parameters

### R0 — Ohmic Resistance
Represents immediate voltage drop caused by:
- electrolyte resistance
- electrode resistance
- connector/contact resistance

R0 increases significantly at lower SOC levels due to reduced ion mobility.

---

### R1-C1 Network — Fast Transient Dynamics
Models:
- surface charge effects
- electrochemical reaction dynamics
- short-term voltage relaxation

This branch captures rapid voltage transients immediately after load changes.

---

### R2-C2 Network — Slow Diffusion Dynamics
Models:
- lithium-ion diffusion behavior
- deeper electrochemical polarization
- long-term relaxation effects

This branch captures slow stabilization behavior after current pulses.

## HPPC Dataset and Experimental Setup

Hybrid Pulse Power Characterization (HPPC) testing was used to collect real battery response data.

The battery was discharged using:
- 100 A pulse current
- 20 second discharge pulses
- 12 minute relaxation periods

The entire SOC range was divided into 9 regions:
- 100–90%
- 90–80%
- ...
- 10–0%

Each region contains:
- time
- current
- measured voltage
- simulated voltage

This region-wise segmentation improves parameter estimation accuracy by accounting for nonlinear battery behavior across SOC ranges.

## My Contribution — Region-Wise EKF Parameter Estimation

My primary contribution in this project focused on:

- preprocessing HPPC experimental datasets
- SOC region segmentation
- Extended Kalman Filter implementation
- estimation of:
  - R0
  - R1
  - R2
  - C1
  - C2
- tuning EKF process and measurement covariance matrices
- generation of region-wise lookup tables
- validation of estimated parameters against measured voltage data

The RC parameter estimation framework was designed to dynamically capture battery behavior under varying operating conditions, enabling improved SOC estimation accuracy in the ASRUKF stage.

# Extended Kalman Filter (EKF) Implementation
## EKF State Vector

The EKF operates using a 7-dimensional state vector:

x = [VRC1, VRC2, R0, R1, C1, R2, C2]
## State Prediction Equations

VRC1[k] = exp(-Δt/(R1*C1)) * VRC1[k−1]
          + (1 - exp(-Δt/(R1*C1))) * R1 * I[k−1]

VRC2[k] = exp(-Δt/(R2*C2)) * VRC2[k−1]
          + (1 - exp(-Δt/(R2*C2))) * R2 * I[k−1]

## Terminal Voltage Estimation

Vpred = OCV(SOC) - VRC1 - VRC2 - R0 * I

## Why EKF Was Used

The battery system is nonlinear because:
- OCV depends nonlinearly on SOC
- RC dynamics are exponential
- battery parameters vary with operating conditions

The Extended Kalman Filter linearizes the nonlinear system around the current operating point using Jacobian matrices.

Advantages:
- recursive online estimation
- low computational cost
- real-time feasibility
- robust parameter convergence

  # Region-Wise Parameter Estimation

Instead of estimating one global parameter set for the entire battery range, the SOC range was divided into multiple operating regions.

Each region used:
- independent EKF tuning
- separate Q and R matrices
- dedicated RC parameter extraction

Benefits:
- improved modeling accuracy
- better transient response
- lower voltage estimation error
- reduced filter divergence
- adaptive representation of battery nonlinearities

# Key Observations

- R0 increased as SOC decreased
- R1 and R2 remained relatively stable
- C1 and C2 showed minimal variation
- EKF converged successfully in all SOC regions
- terminal voltage estimation closely matched measured voltage
- low RMSE achieved during validation

  # Tools and Technologies

- MATLAB R2018a
- Simulink
- Extended Kalman Filter (EKF)
- Adaptive Square Root Unscented Kalman Filter (ASRUKF)
- HPPC Experimental Testing
- SIL Validation
- OPAL-RT
