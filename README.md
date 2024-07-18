# Introduction
REPOSITORY UNDER CONSTRUCTION
<br>
<br>This MATLAB program computes Earth-Moon low-energy transfers and optimizes the total delta v, given the orbital elements of a departure orbit, a departure date and a desired lunar circular orbit radius. The current version requires a good initial guess that has to be manually tweaked before starting the optimization and also has some convergence problems. Future versions will automate the preliminary initial guess optimization process, add the possibility of selecting a desired inclination for the lunar orbit, and improve convergence.

<br>The following kernels (can be found here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/) are required to run this program:
<br>de440.bsp
<br>gm_de440.tpc
<br>naif0012.tls
<br>pck00011.tpc

<br>
Sample plot:

![lowEnergyTransf](https://github.com/user-attachments/assets/782db5de-412b-4f7a-8fae-00da14ed1bd2)

# References
This program is mostly based on E. Belbruno's work on weak stability boundary transfers. More detailed references will be added later.
