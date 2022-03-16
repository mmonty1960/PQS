# PQS
Particle Qualification Software

The software was developed in the framework of the SolarPACES project "Characterization of optical properties of particles for CSP"

The software offers several methods (#0, #1,...,#5) for correcting hemispherical reflectance spectra of particles for the presence of the transparent window used to confine particles in a box.

The different methods are obtained by the general optical model of the window-particles system under specific hypothesis

QUICK USER GUIDE:

1) Choice the standard to be used for computing mean reflectance. In UV-VIS-NIR: IEC 60904-3, ASTM G173, ISO 9050, ISO 9845-1, E 891; in IR range only IRnoWeight

2) For each item, set the filename of the spectrophotometer measurument NOT YET MULTIPLIED for the abolute reflectance of the used reference
      Rws  = reflectance of window-particles
      Tw   = window transmittance
      Rw   = window reflectance
      Rw#1 = reference Ref#1 beyond the window
      ZeroLine = spectrum acquired without anything on the port of the integrating sphere
      Ref#1 = absolute reflectance of the reference Ref#1
      Meas Ref#1 = only for IR - spectrum acquired with Ref#1 on the port
      Rw#2 = optional - reference Ref#2 beyond the window
      Rref#2 = absolute reflectance of the reference Ref#2
      
3) Select the appropriate method for the used window 
4) Push "Calculate & Plot R_particle" to plot the spetrum corrected for the set method; the light-grey curve is the window absorption.
5) If necessary, modify the limits of the Y-range
6) If the method is change, by pushing "Calculate & Plot R_particle" the old curves are mantained. Push "Refresh plot" to show only the new.

Please Note:
A) The ZeroLine is differently used in UV-VIS-NIR and IR cases: in the first it is just substracted to Rw because the signal is due to the cap covering the port; differently in IR the ZeroLine signal is due to the beam tails, so it is subtracted at all the spectra
B) In UV-VIS-NIR the baseline is recorded with the reference on the port. Conversely for IR, for each kind of measurement a new baseline is set by sending the beam onto the internal of the sphere to correct for the substitution error 
