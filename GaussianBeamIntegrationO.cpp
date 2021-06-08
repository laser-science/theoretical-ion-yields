/*
 * Volume integration of the temporally integrate populations for O
 *
 * last updated 6/8/21 by Liam Kelley
 *
 * Takes the ion population curves for O and does the volume integration for the ion
 * yields based on the intensities present in the input file.
 *
 * Outputs two files. The first is the volume correlation to the intensities just as a reference.
 * The second file is the total yield for each ion.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <array>

using namespace std;

double lambda = 800e-9; //wavelength in meters
double w0 = 2e-6; //beam waist at the center in meters
double I0 = 7e20; // peaks intensity at the center of the pulse, in W/m^2
double step = I0 / 10000; //step size of the integration (set to 10000 steps)

const int calcsize = 10000;

//arrays for the volumes and intensities of this code
array<double, calcsize> possibleInts; //has units of W/m^2
array<double, calcsize> volumes; //has units of m^3
array<double, calcsize> IntsFromCalc; // has units of W/cm^2

const int sizefile = 28570; //size of the files from the population calculation

array<double, sizefile> IntsFromFile; //array of intensities from the pop calculation
string tempintsi; //storage variable for intensities
string tempintau;

array<double, sizefile> O0pop; //arrays for the populations of each ion of carbon
array<double, sizefile> O1pop;
array<double, sizefile> O2pop;
array<double, sizefile> O3pop;
array<double, sizefile> O4pop;
array<double, sizefile> O5pop;
array<double, sizefile> O6pop;

string tempo0; //storage variables for the populations of carbon
string tempo1;
string tempo2;
string tempo3;
string tempo4;
string tempo5;
string tempo6;

double o0yield = 0; //ion yields for the carbon ions
double o1yield = 0;
double o2yield = 0;
double o3yield = 0;
double o4yield = 0;
double o5yield = 0;
double o6yield = 0;

int popindex = 0;

double Xi(double intensity) { //function for calculating Xi (see Sui's paper)
	return sqrt(I0 / intensity - 1.0);
}

double VolumeCalc(double xi) { //volume as a function of Xi, which is a funtion of intensity
	return M_PI * pow(w0, 4.0) / lambda
			* (2.0 / 9.0 * pow(xi, 3.0) + 4.0 / 3.0 * xi - 4.0 / 3.0 * atan(xi));
}

int main() {

	//generating the intensities to be calculated in the volume function
	for (int i = 0; i < calcsize; i++) {
		possibleInts[i] = step + i * step;
	}

	ifstream inFile;
	ofstream outFile;

	//READING IN TO THE POPULATION ARRAYS

	inFile.open(
			"C:\\Summer 2020 Lab Work\\C++\\VolumeIntegrationO\\Debug\\input.dat");

	while (popindex < sizefile) {

		inFile >> tempintsi;
		inFile >> tempintau;
		inFile >> tempo0;
		inFile >> tempo1;
		inFile >> tempo2;
		inFile >> tempo3;
		inFile >> tempo4;
		inFile >> tempo5;
		inFile >> tempo6;

		IntsFromFile[popindex] = atof(tempint.c_str());
		O0pop[popindex] = atof(tempo0.c_str());
		O1pop[popindex] = atof(tempo1.c_str());
		O2pop[popindex] = atof(tempo2.c_str());
		O3pop[popindex] = atof(tempo3.c_str());
		O4pop[popindex] = atof(tempo4.c_str());
		O5pop[popindex] = atof(tempo5.c_str());
		O6pop[popindex] = atof(tempo6.c_str());

		popindex++;
	}

	//first value in the intensity and volume arrays

	volumes[0] = VolumeCalc(Xi(possibleInts[0]));
	IntsFromCalc[0] = possibleInts[0] / 10000;

	//POPULATING THE VOLUME AND INTENSITY ARRAYS WITH THE VOLUME FUNCTION
	for (int i = 1; i < calcsize; i++) {

		double IntNow = possibleInts[i];
		double IntPrev = possibleInts[i - 1];

		double volumeUpTo = VolumeCalc(Xi(IntPrev));

		volumes[i] = volumeUpTo - VolumeCalc(Xi(IntNow));
		//each volume is the difference of two volume 'shells' based on intensity

		IntsFromCalc[i] = (IntNow) / 10000;

	}

	// * For this next loop, the program runs through every intensity in the calculated, evenly spaced
	// * intensities and matches it up with a close intensity from the input file.
	// *
	//* From there, the additions to the yields are calculated and added to the totals

	for (int searcher = 0; searcher < calcsize; searcher++) {

		double calcval = IntsFromCalc[searcher];

		for (int j = 0; j < sizefile; j++) {

			double fileval = IntsFromFile[j];

			if (abs(calcval - fileval) <= .01 * calcval && calcval < fileval) {

				o1yield += (volumes[searcher] * O1pop[j] * calcsize);
				o2yield += (volumes[searcher] * O2pop[j] * calcsize);
				o3yield += (volumes[searcher] * O3pop[j] * calcsize);
				o4yield += (volumes[searcher] * O4pop[j] * calcsize);
				o5yield += (volumes[searcher] * O5pop[j] * calcsize);
				o6yield += (volumes[searcher] * O6pop[j] * calcsize);

				break;

			}
		}

	}

	//printing out the yields

	cout << o1yield << "\t";
	cout << o2yield << "\t";
	cout << o3yield << "\t";
	cout << o4yield << "\t";
	cout << o5yield << "\t";
	cout << o6yield << endl;


	//Writing out intensities and volumes to file

outFile.open("outputVolumes.dat");

	 outFile << "Intensity (SI)";
	 outFile << "\t";
	 outFile << "Volume (SI)";
	 outFile << "\n";

	 for (int k = 0; k < calcsize; k++) {
	 outFile << IntsFromCalc[k];
	 outFile << "\t";
	 outFile << volumes[k];
	 outFile << "\n";
	 }

	 outFile.close();

	 outFile.open("outputYields.dat");

	 outFile << "Oxygen 1+ Yield";
	 outFile << "\t";
	 outFile << "Oxygen 2+ Yield";
	 outFile << "\t";
	 outFile << "Oxygen 3+ Yield";
	 outFile << "\t";
	 outFile << "Oxygen 4+ Yield";
	 outFile << "\t";
	 outFile << "Oxygen 5+ Yield";
	 outFile << "\t";
	 outFile << "Oxygen 6+ Yield";
	 outFile << "\n";

     outFile << o1yield;
	 outFile << "\t";
	 outFile << o2yield;
	 outFile << "\t";
	 outFile << o3yield;
	 outFile << "\t";
	 outFile << o4yield;
	 outFile << "\t";
	 outFile << o5yield;
	 outFile << "\t";
	 outFile << o6yield;
	 outFile << "\n";

	 outFile.close();


	return 0;

}
