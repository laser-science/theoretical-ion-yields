/*ADK rate and ion population curve calculations for Carbon

Written by Liam Kelley, last updated 6/8/2021

First produces a time-dependent intensity profile of a Gaussian beam based on the FWHM and max intensity,
which a user can edit. Then uses the intensity arrays to calculate the ADK ionization rates for every species
up to carbon 4+. Then, following the rate equations for ion populations, integrates the rates to create population curves.
Also sets an upper limit on the ADK rate (1e-17 /s) to account for electron orbital period.

The code outputs two files. The first includes the ion population curves with their associated intensities.
The second is the ADK rates (in atomic units for the ion species) before the integration occurs.
*/


/****************************************************************************/
/*Headers*/
/****************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>

using namespace std;

/****************************************************************************/
/*variables and constants*/
/****************************************************************************/
double c = 2.998e8; //speed of light in SI
double e0 = 8.854e-12; //epsilon naught in SI

const int size = 58000; //size of the arrays

double dt = 0.005e-15; //time interval for the integration
double tstart = -100e-15; //starting time for the integration
double tend = 100e-15; //ending time for the integration
double tnow2 = 0;

double ipC0 = 11.2602880 / 27.2; // Factors needed to calculate the ADK rate for C0 (atomic units) (see Fittinghoff thesis)
double nstarC0 = 1.0 / sqrt(2.0 * ipC0);
double npowerC0 = 2.0 * nstarC0 - 1;
double lstarC0 = nstarC0 - 1.0;
double flmC0 = 3.0;
double epsilonC0 = pow(2.0 * ipC0, 1.5);
double c2nlC0 =
		pow(2.0, 2.0 * nstarC0)
				/ (nstarC0 * tgamma(nstarC0 + lstarC0 + 1.0)
						* tgamma(nstarC0 - lstarC0));

double ipC1 = 24.383154 / 27.2; // Factors needed to calculate the ADK rate for C1 (atomic units)
double nstarC1 = 2.0 / sqrt(2.0 * ipC1);
double lstarC1 = nstarC1 - 1.0;
double flmC1 = 3.0;
double epsilonC1 = pow(2.0 * ipC1, 1.5);
double c2nlC1 =
		pow(2.0, 2.0 * nstarC1)
				/ (nstarC1 * tgamma(nstarC1 + lstarC1 + 1.0)
						* tgamma(nstarC1 - lstarC1));

double ipC2 = 47.88778 / 27.2; // Factors needed to calculate the ADK rate for C2 (atomic units)
double nstarC2 = 3.0 / sqrt(2.0 * ipC2);
double lstarC2 = nstarC2 - 1.0;
double flmC2 = 1.0;
double epsilonC2 = pow(2.0 * ipC2, 1.5);
double c2nlC2 =
		pow(2.0, 2.0 * nstarC2)
				/ (nstarC2 * tgamma(nstarC2 + lstarC2 + 1.0)
						* tgamma(nstarC2 - lstarC2));

double ipC3 = 64.49352 / 27.2; // Factors needed to calculate the ADK rate for C3 (atomic units)
double nstarC3 = 4.0 / sqrt(2.0 * ipC3);
double lstarC3 = nstarC3 - 1.0;
double flmC3 = 1.0;
double epsilonC3 = pow(2.0 * ipC3, 1.5);
double c2nlC3 =
		pow(2.0, 2.0 * nstarC3)
				/ (nstarC3 * tgamma(nstarC3 + lstarC3 + 1.0)
						* tgamma(nstarC3 - lstarC3));

double ipC4 = 392.090515 / 27.2; // Factors needed to calculate the ADK rate for C4 (atomic units)
double nstarC4 = 5.0 / sqrt(2.0 * ipC4);
double lstarC4 = nstarC4 - 1.0;
double flmC4 = 1.0;
double epsilonC4 = pow(2.0 * ipC4, 1.5);
double c2nlC4 =
		pow(2.0, 2.0 * nstarC4)
				/ (nstarC4 * tgamma(nstarC4 + lstarC4 + 1.0)
						* tgamma(nstarC4 - lstarC4));

double IntensitySI[size]; //arrays for storing the intensity and field strength information, in SI units
double EFieldSI[size];
double EFieldau[size];
double timeArray[size];

//Beam parameters
double IMaxAmpSI = 2.4e17; //max intensity of the pulse, in W/m^2 EDIT THIS
double fwhm = 40e-15; //full width half max of the pulse, in seconds
double sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
double lambda = 800e-9;

int fieldCounter = 0;

//CONTAINERS FOR ION POPS

double ADKC0[size];
double PhiC0[size];
double PC0[size];

double ADKC1[size];
double PhiC1[size];
double PC1[size];

double ADKC2[size];
double PhiC2[size];
double PC2[size];

double ADKC3[size];
double PhiC3[size];
double PC3[size];

double ADKC4[size];
double PhiC4[size];
double PC4[size];

int popCounter = 0;

double stepC0test[size];
double stepC1test[size];
double stepC2test[size];
double stepC3test[size];
double stepC4test[size];

ofstream outFile;

/****************************************************************************/
/*global*/
/****************************************************************************/

double integrate_for_phi(double ADK[], double t) { //function for calculating the integral involved in phi[t]

	int tcounter = 0;
	double tnowInt = tstart;
	double sum = 0.0;
	double deltaSum;

	while (tnowInt <= t) {
		if (ADK[tcounter] < 1e-140) {
			deltaSum = 0.0;
		} else {
			deltaSum = ADK[tcounter] / (2.4188e-17) * dt;
		}
		tnowInt += dt;
		tcounter++;

		sum += deltaSum;
	}

	return sum;

}

double integrate_for_pop(double ADKLess[], double PopLess[], double PhiNow[],
		double t) {
	//function for calculating the integral in the expression for p[t]

	int tcounterPop = 0;
	double tnowPop = tstart;
	double sum = 0.0;
	double deltaSum;

	while (tnowPop <= t) {
		if (ADKLess[tcounterPop] < 1e-140 || PopLess[tcounterPop] < 1e-140) {
			deltaSum = 0.0;
		} else {
			deltaSum = (exp(PhiNow[tcounterPop])
					* (ADKLess[tcounterPop] / (2.4188e-17))
					* PopLess[tcounterPop]) * dt;
		}
		tnowPop += dt;
		tcounterPop++;
		sum += deltaSum;
	}

	return sum;

}

void ADKAdjust(double ADK[]) { //adjusts ADK rates to account for maximum rate (corresponding to electron orbital time)

	int hits[size];
	int hitscounter = 0;

	for (int counter = 0; counter < size; counter++) { //loop for taking values in between first and last over 1e17
		if (ADK[counter] / (2.4188e-17) >= 1e17) {
			hits[hitscounter] = counter;
			hitscounter++;
		}
	}

	int ii = 0;

	while (hits[ii] != 0) { //int for getting last over 1e17
		ii++;
	}

	int lastcounter = hits[0];

	while ((lastcounter >= hits[0]) && (lastcounter <= hits[ii - 1])) {
		ADK[lastcounter] = 1e17 * (2.4188e-17);
		lastcounter++;
	}

}

int main() {

	//POPULATING THE PULSE ARRAYS
	for (double tnow = tstart; tnow < tend; tnow += dt) {

		IntensitySI[fieldCounter] = IMaxAmpSI
				* exp(-tnow * tnow / (2.0 * sigma * sigma)); //Gaussian Distribution
		EFieldSI[fieldCounter] = sqrt(
				2.0 * IntensitySI[fieldCounter] / (c * e0))
				* abs(sin(2.0 * M_PI * c / lambda * tnow)); //Field from the intensity

		EFieldau[fieldCounter] = EFieldSI[fieldCounter] / 5.142e11;

		//ADK RATE ARRAYS

		double argC0 = -2.0 * epsilonC0 / (3.0 * EFieldau[fieldCounter]);
		ADKC0[fieldCounter] = c2nlC0 * ipC0 * flmC0
				* pow(2.0 * epsilonC0 / EFieldau[fieldCounter], npowerC0)
				* exp(argC0); //ADK Rate C1 (used atomic units for field)

		double argC1 = -2.0 * epsilonC1 / (3.0 * EFieldau[fieldCounter]);
		ADKC1[fieldCounter] = c2nlC1 * ipC1 * flmC1
				* pow(2.0 * epsilonC1 / EFieldau[fieldCounter],
						2.0 * nstarC1 - 1.0) * exp(argC1); //ADK Rate C1 (used atomic units for field)

		double argC2 = -2.0 * epsilonC2 / (3.0 * EFieldau[fieldCounter]);
		ADKC2[fieldCounter] = c2nlC2 * ipC2 * flmC2
				* pow(2.0 * epsilonC2 / EFieldau[fieldCounter],
						2.0 * nstarC2 - 1.0) * exp(argC2); //ADK Rate C2(used atomic units for field)

		double argC3 = -2.0 * epsilonC3 / (3.0 * EFieldau[fieldCounter]);
		ADKC3[fieldCounter] = c2nlC3 * ipC3 * flmC3
				* pow(2.0 * epsilonC3 / EFieldau[fieldCounter],
						2.0 * nstarC3 - 1.0) * exp(argC3); //ADK Rate C3(used atomic units for field)

		double argC4 = -2.0 * epsilonC4 / (3.0 * EFieldau[fieldCounter]);
		ADKC4[fieldCounter] = c2nlC4 * ipC4 * flmC4
				* pow(2.0 * epsilonC4 / EFieldau[fieldCounter],
						2.0 * nstarC4 - 1.0) * exp(argC4); //ADK Rate C4(used atomic units for field)

		timeArray[fieldCounter] = tnow;
		fieldCounter++;

	}

	ADKAdjust(ADKC0);
	ADKAdjust(ADKC1);
	ADKAdjust(ADKC2);
	ADKAdjust(ADKC3);
	ADKAdjust(ADKC4);

	//Main loop for ion population curves

	for (double tnow = tstart; tnow < tend; tnow += dt) {

		PhiC0[popCounter] = integrate_for_phi(ADKC0, tnow);
		PC0[popCounter] = exp(-PhiC0[popCounter]);

		PhiC1[popCounter] = integrate_for_phi(ADKC1, tnow);
		PC1[popCounter] = exp(-PhiC1[popCounter])
				* integrate_for_pop(ADKC0, PC0, PhiC1, tnow);

		PhiC2[popCounter] = integrate_for_phi(ADKC2, tnow);
		PC2[popCounter] = exp(-PhiC2[popCounter])
				* integrate_for_pop(ADKC1, PC1, PhiC2, tnow);

		PhiC3[popCounter] = integrate_for_phi(ADKC3, tnow);
		PC3[popCounter] = exp(-PhiC3[popCounter])
				* integrate_for_pop(ADKC2, PC2, PhiC3, tnow);

		PhiC4[popCounter] = integrate_for_phi(ADKC4, tnow);
		PC4[popCounter] = exp(-PhiC4[popCounter])
				* integrate_for_pop(ADKC3, PC3, PhiC4, tnow);

		popCounter++;

	}

	//OUTPUTTING TO FILES

	outFile.open("outputPopulations.dat");

	outFile << "Intensity (SI)";
	outFile << "\t";
	outFile << "Intensity (au)";
	outFile << "\t";
	outFile << "Carbon 0+ population";
	outFile << "\t";
	outFile << "Carbon 1+ population";
	outFile << "\t";
	outFile << "Carbon 2+ population";
	outFile << "\t";
	outFile << "Carbon 3+ population";
	outFile << "\t";
	outFile << "Carbon 4+ population";
	outFile << "\n";

	for (int k = 0; k < size; k++) { //writing to file the intensity - population data (W/cm^2)
		outFile << IntensitySI[k] / 10000.0;
		outFile << "\t";
		outFile << IntensitySI[k] / 10000.0 /(3.50945e16);
		outFile << "\t";
		outFile << PC0[k];
		outFile << "\t";
		outFile << PC1[k];
		outFile << "\t";
		outFile << PC2[k];
		outFile << "\t";
		outFile << PC3[k];
		outFile << "\t";
		outFile << PC4[k];
		outFile << "\n";

	}

	outFile.close();

	outFile.open("outputADKRates.dat");

	outFile << "Intensity (SI)";
	outFile << "\t";
	outFile << "Intensity (au)";
	outFile << "\t";
	outFile << "Carbon 0+ ADK Rate (au)";
	outFile << "\t";
	outFile << "Carbon 1+ ADK Rate (au)";
	outFile << "\t";
	outFile << "Carbon 2+ ADK Rate (au)";
	outFile << "\t";
	outFile << "Carbon 3+ ADK Rate (au)";
	outFile << "\t";
	outFile << "Carbon 4+ ADK Rate (au)";
	outFile << "\n";

	for (int k = 0; k < size; k++) { //writing to file the intensity - population data (W/cm^2)
		outFile << IntensitySI[k] / 10000.0;
		outFile << "\t";
		outFile << IntensitySI[k] / 10000.0 / (3.50945e16);
		outFile << "\t";
		outFile << ADKC0[k];
		outFile << "\t";
		outFile << ADKC1[k];
		outFile << "\t";
		outFile << ADKC2[k];
		outFile << "\t";
		outFile << ADKC3[k];
		outFile << "\t";
		outFile << ADKC4[k];
		outFile << "\n";

	}

	outFile.close();

	return 0;

}
