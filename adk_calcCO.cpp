/*ADK rate and ion population curve calculations for Carbon Monoxide

Written by Liam Kelley, last updated 6/8/2021

First produces a time-dependent intensity profile of a Gaussian beam based on the FWHM and max intensity,
which a user can edit. Then uses the intensity arrays to calculate the ADK ionization rates for every species
up to carbon monoxide 1+. Then, following the rate equations for ion populations, integrates the rates to create population curves.
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

const int size = 116000; //size of the arrays

double dt = 0.005e-15; //time interval for the integration
double tstart = -100e-15; //starting time for the integration
double tend = 100e-15; //ending time for the integration

double ipCO0 = 14.01 / 27.2; // Factors needed to calculate the ADK rate for C0 (atomic units)
double ZCO0 = 1.0;
double KCO0 = sqrt(2.0 * ipCO0);
double Ql0 = sqrt(1.0 / 2.0);
double Ql1 = sqrt(3.0 / 2.0);
double Bsum = 1.43 * Ql0 + 0.76 * Ql1;

double IntensitySI[size]; //arrays for storing the intensity and field strength information, in SI units
double EFieldSI[size];
double EFieldau[size];
double IMaxAmpSI = 1e22; //max intensity of the pulse, in W/m^2
double fwhm = 40e-15; //full width half max of the pulse, in seconds
double sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
double timeArray[size];
double lambda = 800e-9;

int fieldCounter = 0;

//CONTAINERS FOR ION POPS

double ADKCO0[size];
double PhiCO0[size];
double PCO0[size];

double PCO1[size];

int popCounter = 0;

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

double integrate_for_last(double ADKLess[], double PopLess[], double t) { //function for getting the last ion population in a series

	int tCounterLast = 0;
	double tNow = tstart;
	double sum = 0.0;
	double deltaSum;

	while (tNow <= t) {
		if (ADKLess[tCounterLast] < 1e-100 || PopLess[tCounterLast] < 1e-100) {
			deltaSum = 0.0;
		} else {
			deltaSum = ((ADKLess[tCounterLast] / (2.4188e-17))
					* PopLess[tCounterLast]) * dt;
		}

		tNow += dt;
		tCounterLast++;
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

		ADKCO0[fieldCounter] = pow(Bsum, 2.0) * 1.0
				/ (pow(KCO0, 2.0 * ZCO0 / KCO0 - 1.0))
				* pow(2.0 * pow(KCO0, 3.0) / EFieldau[fieldCounter],
						2.0 * ZCO0 / KCO0 - 1.0)
				* exp(-2.0 * pow(KCO0, 3.0) / (3.0 * EFieldau[fieldCounter])); //ADK Rate C1 (used atomic units for field)

		timeArray[fieldCounter] = tnow;
		fieldCounter++;

	}

	ADKAdjust(ADKCO0);

	//The main loop for population curve generation

	for (double tnow = tstart; tnow < tend; tnow += dt) {

		PhiCO0[popCounter] = integrate_for_phi(ADKCO0, tnow);
		PCO0[popCounter] = exp(-PhiCO0[popCounter]);

		PCO1[popCounter] = integrate_for_last(ADKCO0, PCO0, tnow);

		popCounter++;

	}

	//OUTPUTTING TO FILES

	outFile.open("outputPopulations.dat");

	outFile << "Intensity (SI)";
	outFile << "\t";
	outFile << "Intensity (au)";
	outFile << "\t";
	outFile << "Carbon Monoxide 0+ population";
	outFile << "\t";
	outFile << "Carbon Monoxide 1+ population";
	outFile << "\n";

	for (int k = 0; k < size; k++) { //writing to file the intensity - population data (W/cm^2)
		outFile << IntensitySI[k] / 10000.0;
		outFile << "\t";
		outFile << IntensitySI[k] / 10000.0 /(3.50945e16);
		outFile << "\t";
		outFile << PCO0[k];
		outFile << "\t";
		outFile << PCO1[k];
		outFile << "\n";

	}

	outFile.close();

	outFile.open("outputADKRates.dat");

	outFile << "Intensity (SI)";
	outFile << "\t";
	outFile << "Intensity (au)";
	outFile << "\t";
	outFile << "Carbon Monoxide 0+ ADK Rate (au)";
	outFile << "\t";
	outFile << "Carbon Monoxide 1+ ADK Rate (au)";
	outFile << "\n";

	for (int k = 0; k < size; k++) { //writing to file the intensity - population data (W/cm^2)
		outFile << IntensitySI[k] / 10000.0;
		outFile << "\t";
		outFile << IntensitySI[k] / 10000.0 /(3.50945e16);
		outFile << "\t";
		outFile << ADKCO0[k];
		outFile << "\t";
		outFile << ADKCO1[k];
		outFile << "\n";

	}

	outFile.close();

	return 0;

}
