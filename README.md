# System-Identification
## Objective
This project aims to identify a mathematical model for a given physical system based on the parameters identified at resonance.
## Context
The input is a voltage controlled oscilator, whose purpose is to find the resonance frequency and determine the mathematical model of the system through experimental data collection 
## Implementation steps
  * Collecting the input/output data after applying the signal
  * Finding the mathematical model of the system using non-parametrical methods (based on manual measurements on the collected data)
  * Finding the mathematical model of the system using parametrical system identification methods (ARX, ARMAX)
  * Comparing and deciding which model best describes the behaviour of the system.
## Validation of the model
  * The validation procces consists of comparing the collected data with the frequency response of the identified model.
  * For the non-parametrical method, validation sould be obtained through NMSE (Normalised Mean Square Error). The error should be less than 10% for the identified system to be considered validated.
  * For the parametrical approach, 2 methods were used for autocorrelation validation and 2 methods were used for intercorrelation validation.
