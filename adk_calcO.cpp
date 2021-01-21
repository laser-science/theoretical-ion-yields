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
double e0 = 8.854e-12; //epsilon nauht in SI

const int size = 58000; //size of the arrays

double dt = 0.005e-15; //time interval for the integration
double tstart = -100e-15; //starting time for the integration
double tend = 100e-15; //ending time for the integration

double ipO0 = 13.618055 / 27.2; // Factors needed to calculate the ADK rate for O1 (atomic units)
double nstarO0 = 1.0 / sqrt(2.0 * ipO0);
double lstarO0 = nstarO0 - 1.0;
double flmO0 = 3.0;
double epsilonO0 = pow(2.0 * ipO0, 1.5);
double c2nlO0 = pow(2.0, 2.0 * nstarO0)
		/ (nstarO0 * tgamma(nstarO0 + lstarO0 + 1.0) * tgamma(nstarO0 - lstarO0));

double ipO1 = 35.12112 / 27.2; // Factors needed to calculate the ADK rate for O1 (atomic units)
double nstarO1 = 2.0 / sqrt(2.0 * ipO1);
double lstarO1 = nstarO1 - 1.0;
double flmO1 = 3.0;
double epsilonO1 = pow(2.0 * ipO1, 1.5);
double c2nlO1 = pow(2.0, 2.0 * nstarO1)
		/ (nstarO1 * tgamma(nstarO1 + lstarO1 + 1.0) * tgamma(nstarO1 - lstarO1));

double ipO2 = 54.93554 / 27.2; // Factors needed to calculate the ADK rate for O2 (atomic units)
double nstarO2 = 3.0 / sqrt(2.0 * ipO2);
double lstarO2 = nstarO2 - 1.0;
double flmO2 = 3.0;
double epsilonO2 = pow(2.0 * ipO2, 1.5);
double c2nlO2 = pow(2.0, 2.0 * nstarO2)
		/ (nstarO2 * tgamma(nstarO2 + lstarO2 + 1.0) * tgamma(nstarO2 - lstarO2));

double ipO3 = 77.41350 / 27.2; // Factors needed to calculate the ADK rate for O3 (atomic units)
double nstarO3 = 4.0 / sqrt(2.0 * ipO3);
double lstarO3 = nstarO3 - 1.0;
double flmO3 = 3.0;
double epsilonO3 = pow(2.0 * ipO3, 1.5);
double c2nlO3 = pow(2.0, 2.0 * nstarO3)
		/ (nstarO3 * tgamma(nstarO3 + lstarO3 + 1.0) * tgamma(nstarO3 - lstarO3));

double ipO4 = 113.8990 / 27.2; // Factors needed to calculate the ADK rate for O4 (atomic units)
double nstarO4 = 5.0 / sqrt(2.0 * ipO4);
double lstarO4 = nstarO4 - 1.0;
double flmO4 = 1.0;
double epsilonO4 = pow(2.0 * ipO4, 1.5);
double c2nlO4 = pow(2.0, 2.0 * nstarO4)
		/ (nstarO4 * tgamma(nstarO4 + lstarO4 + 1.0) * tgamma(nstarO4 - lstarO4));

double ipO5 = 138.1189 / 27.2; // Factors needed to calculate the ADK rate for O5 (atomic units)
double nstarO5 = 6.0 / sqrt(2.0 * ipO5);
double lstarO5 = nstarO5 - 1.0;
double flmO5 = 1.0;
double epsilonO5 = pow(2.0 * ipO5, 1.5);
double c2nlO5 = pow(2.0, 2.0 * nstarO5)
		/ (nstarO5 * tgamma(nstarO5 + lstarO5 + 1.0) * tgamma(nstarO5 - lstarO5));

double ipO6 = 739.32682 / 27.2; // Factors needed to calculate the ADK rate for O5 (atomic units)
double nstarO6 = 7.0 / sqrt(2.0 * ipO6);
double lstarO6 = nstarO6 - 1.0;
double flmO6 = 1.0;
double epsilonO6 = pow(2.0 * ipO6, 1.5);
double c2nlO6 = pow(2.0, 2.0 * nstarO6)
		/ (nstarO6 * tgamma(nstarO6 + lstarO6 + 1.0) * tgamma(nstarO6 - lstarO6));

double IntensitySI[size]; //arrays for storing the intensity and field strength information, in SI units
double EFieldSI[size];
double EFieldau[size];
double IMaxAmpSI = 6.8e18; //max intensity of the pulse, in W/m^2
double fwhm = 40e-15; //full width half max of the pulse, in seconds
double sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
double timeArray[size];
double lambda = 800e-9;

int fieldCounter = 0;

//CONTAINERS FOR ION POPS

double ADKO0[size];
double PhiO0[size];
double PO0[size];

double ADKO1[size];
double PhiO1[size];
double PO1[size];

double ADKO2[size];
double PhiO2[size];
double PO2[size];

double ADKO3[size];
double PhiO3[size];
double PO3[size];

double ADKO4[size];
double PhiO4[size];
double PO4[size];

double ADKO5[size];
double PhiO5[size];
double PO5[size];

double ADKO6[size];
double PhiO6[size];
double PO6[size];

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
		if (ADK[tcounter] < 1e-100) {
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
		if (ADKLess[tcounterPop < 1e-100]
				|| PopLess[tcounterPop] < 1e-100) {
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

void ADKAdjust(double ADK[]) {

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

		EFieldSI[fieldCounter] = sqrt(2.0 * IntensitySI[fieldCounter] / (c * e0))
				* abs(sin(2.0 * M_PI * c / lambda * tnow)); //Field from the intensity

		EFieldau[fieldCounter] = EFieldSI[fieldCounter] / 5.142e11;

		//ADK RATE ARRAYS

		ADKO0[fieldCounter] = c2nlO0 * ipO0 * flmO0
				* pow(2.0 * epsilonO0 / EFieldau[fieldCounter],
						2.0 * nstarO0 - 1.0)
				* exp(
						-2.0 * epsilonO0
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O1 (used atomic units for field)

		ADKO1[fieldCounter] = c2nlO1 * ipO1 * flmO1
				* pow(2.0 * epsilonO1 / EFieldau[fieldCounter],
						2.0 * nstarO1 - 1.0)
				* exp(
						-2.0 * epsilonO1
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O1 (used atomic units for field)

		ADKO2[fieldCounter] = c2nlO2 * ipO2 * flmO2
				* pow(2.0 * epsilonO2 / EFieldau[fieldCounter],
						2.0 * nstarO2 - 1.0)
				* exp(
						-2.0 * epsilonO2
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O2(used atomic units for field)

		ADKO3[fieldCounter] = c2nlO3 * ipO3 * flmO3
				* pow(2.0 * epsilonO3 / EFieldau[fieldCounter],
						2.0 * nstarO3 - 1.0)
				* exp(
						-2.0 * epsilonO3
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O3(used atomic units for field)

		ADKO4[fieldCounter] = c2nlO4 * ipO4 * flmO4
				* pow(2.0 * epsilonO4 / EFieldau[fieldCounter],
						2.0 * nstarO4 - 1.0)
				* exp(
						-2.0 * epsilonO4
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O4(used atomic units for field)

		ADKO5[fieldCounter] = c2nlO5 * ipO5 * flmO5
				* pow(2.0 * epsilonO5 / EFieldau[fieldCounter],
						2.0 * nstarO5 - 1.0)
				* exp(
						-2.0 * epsilonO5
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O5(used atomic units for field)

		ADKO6[fieldCounter] = c2nlO6 * ipO6 * flmO6
				* pow(2.0 * epsilonO6 / EFieldau[fieldCounter],
						2.0 * nstarO6 - 1.0)
				* exp(
						-2.0 * epsilonO6
								/ (3.0 * EFieldau[fieldCounter])); //ADK Rate O6(used atomic units for field)

		timeArray[fieldCounter] = tnow;
		fieldCounter++;
	}

	ADKAdjust(ADKO0);
	ADKAdjust(ADKO1);
	ADKAdjust(ADKO2);
	ADKAdjust(ADKO3);
	ADKAdjust(ADKO4);
	ADKAdjust(ADKO5);
	ADKAdjust(ADKO6);

	//INTEGRATING FOR PHI AND POPULATIONS
	for (double tnow = tstart; tnow < tend; tnow += dt) {

		PhiO0[popCounter] = integrate_for_phi(ADKO0, tnow);
		PO0[popCounter] = exp(-PhiO0[popCounter]);

		PhiO1[popCounter] = integrate_for_phi(ADKO1, tnow);
		PO1[popCounter] = exp(-PhiO1[popCounter])
				* integrate_for_pop(ADKO0, PO0, PhiO1, tnow);

		PhiO2[popCounter] = integrate_for_phi(ADKO2, tnow);
		PO2[popCounter] = exp(-PhiO2[popCounter])
				* integrate_for_pop(ADKO1, PO1, PhiO2, tnow);

		PhiO3[popCounter] = integrate_for_phi(ADKO3, tnow);
		PO3[popCounter] = exp(-PhiO3[popCounter])
				* integrate_for_pop(ADKO2, PO2, PhiO3, tnow);

		PhiO4[popCounter] = integrate_for_phi(ADKO4, tnow);
		PO4[popCounter] = exp(-PhiO4[popCounter])
				* integrate_for_pop(ADKO3, PO3, PhiO4, tnow);

		PhiO5[popCounter] = integrate_for_phi(ADKO5, tnow);
		PO5[popCounter] = exp(-PhiO5[popCounter])
				* integrate_for_pop(ADKO4, PO4, PhiO5, tnow);

		PhiO6[popCounter] = integrate_for_phi(ADKO6, tnow);
		PO6[popCounter] = exp(-PhiO6[popCounter])
				* integrate_for_pop(ADKO5, PO5, PhiO6, tnow);

		popCounter++;

	}

	//OUTPUTTING TO FILES

	outFile.open("output.dat");

	for (int k = 0; k < size; k++) { //writing to file the intensity - population data (W/cm^2)
		outFile << IntensitySI[k] / 10000;
		outFile << "\t";
		outFile << PO0[k];
		outFile << "\t";
		outFile << PO1[k];
		outFile << "\t";
		outFile << PO2[k];
		outFile << "\t";
		outFile << PO3[k];
		outFile << "\t";
		outFile << PO4[k];
		outFile << "\t";
		outFile << PO5[k];
		outFile << "\t";
		outFile << PO6[k];
		outFile << "\n";

	}

	outFile.close();

	return 0;

}
