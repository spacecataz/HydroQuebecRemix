# HydroQuebec Simulation Files

These are the input files for the SWMF simulation of the March 1989
"HydroQuebec" storm.  Solar wind input is reconstructed via a combination
of sparse observation, empirical relationships, and data reconstruction
algorithms.

## Notes on inputs and configuration:

- Start time was selected to start on low-magnitude IMF |B| with 10 hours of preconditioning time before storm arrival.
- The conductance model is the Conductance Model for Extreme Events (CMEE, Mukhopadhyay et al., 2020, Space Weather).
- F10.7 for this event was obtained from LISRD (http://lasp.colorado.edu/lisird/data/noaa_radio_flux/).  The values for March 12, 13, and 14 are 237.6, 253.0, and 263.8 sfu, respectively.  253 is selected for the simulation.
- The new "Boris Region" functionality is employed so that the Boris factor is only applied close to the inner boundary and not everywhere.



## Future considerations:

- Non magnotometer output is very basic right now.  What other files would support this study?

