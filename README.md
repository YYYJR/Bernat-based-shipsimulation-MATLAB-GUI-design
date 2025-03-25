# Bernat-based-shipsimulation-MATLAB-GUI-design
The design of this GUI interface was modified from the SIC forum report by Matheus V. Bernat and Ludvig H. Granström, available at: https://github.com/LINK-SIC-2021-Bernat-Granstrom/ship-simulator.

I replaced one of the Bretschneider spectral models with a **JONSWAP** spectral model to incorporate the effects of wind on waves and ships.

Of course, you can select other wave spectrum models in the wavespec function to meet your simulation needs.

**NOTE：**
  You need to create a folder called **‘wave-files’** in the runtime directory to store the wave data.
  In addition, I encountered the integration interface after the wave does not fluctuate the problem, if you have any ideas welcome to leave a message!
