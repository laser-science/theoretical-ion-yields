# theoretical-ion-yields

For the mathematica primer for the ADK calulations [click here.](https://github.com/laser-science/theoretical-ion-yields/blob/main/Primer%20for%20ADK%20code.nb)

[Carbon ADK rate calculation .cpp file](https://github.com/laser-science/theoretical-ion-yields/blob/main/adk_calcC.cpp)

[Oxygen ADK rate calculation .cpp file](https://github.com/laser-science/theoretical-ion-yields/blob/main/adk_calcO.cpp)

[Carbon Monoxide ADK rate calculation .cpp file](https://github.com/laser-science/theoretical-ion-yields/blob/main/adk_calcCO.cpp)

[Carbon ion yield calculation](https://github.com/laser-science/theoretical-ion-yields/blob/main/GaussianBeamIntegrationC.cpp)

[Oxygen ion yield calculation](https://github.com/laser-science/theoretical-ion-yields/blob/main/GaussianBeamIntegrationO.cpp)

[Carbon Monoxide ion yield calculation](https://github.com/laser-science/theoretical-ion-yields/blob/main/GaussianBeamIntegrationCO.cpp)


# How to use these codes:

For either C, O, or CO, begin with the respective adk_calc.cpp file. In this file, there are variables to edit the characteristics of the Gaussian beam, including the maximum 
intensity amplitude, the wavelength, and the FWHM. The code starts by calculating the ADK rate as a function of intensity and integrating the rates to produce the ion population curves. This code outputs the ADK rates and populations as functions of intensity in two separate .dat files.

Once you run this code, to get the total ion yield the GaussianBeamIntegration.cpp code is next. This code takes the populations and integrates based on the intensity-volume dependence of the beam. Volumes are calculated based on intensity, and the yields are updated based on the populations of the ions at those intensities. To run this code,
take the population .dat file from the previous code and edit the input of the beam integration code based on the location and name of this file. This code also outputs two files.
One is the volume-intensity dependence as a reference, and the second outputs the ion yields for each ion as a single number each. 

Running this two-code process for different maximum beam intensities will allow one to recreate the ion yield curves central to this study and reproduced in the laboratory. 

Suggested C++ compiler: mingw64 [download here](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe/download)

Suggested IDE: Eclipse [download here](https://www.eclipse.org/downloads/)

For help setting up the compiler in Eclipse, [click here](https://www.eclipse.org/4diac/documentation/html/installation/minGW.html)

Download this entire repository as a zip folder to run locally
