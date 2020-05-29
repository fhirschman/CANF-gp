using namespace std;

#include <thread>

#define GNUPLOT_ENABLE_PTY
#define GNUPLOT_ENABLE_FEEDBACK
#include "./../headers/gnuplot-iostream.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "./../headers/data_structures.h"
#include "./../headers/general_functions.h"
#include <thread>

void freeAll(){
	for( int i = 0 ; i < N_Y ; i++ ){
	    delete[] cnumber_matrix[i];
	}
	delete[] cnumber_matrix;
	// gp.close();     
	return;
}

void saveIterMatrix(){
	cout << endl << " ///////////////////// " << endl << " START WRITING TO FILE  " << endl;
	ofstream res;
	res.open("iterMatrix.dat");
	for(int i = 0; i < N_Y; i++){
		for(int j = 0; j < N_X; j++){
			#ifdef GPU_MODE
			res << cuCreal(cnumber_matrix_GPU[i*N_X+j]) << " " << cuCimag(cnumber_matrix_GPU[i*N_X+j]) << " " << log(iter_matrix_GPU[i*N_X+j]) << endl;
			#else
			res << real(cnumber_matrix[i][j]) << " " << imag(cnumber_matrix[i][j]) << " " << log(iter_matrix[i][j]) << endl;

			#endif
		}
		res << endl;
	}
	res.close();
	cout << " FINISH WRITING TO FILE \n" << " ///////////////////// \n " << endl;
	return;
}

void setupGP(){
	stream << "set palette rgbformulae " << colorMatrix[COLOR_CODE][0]  << "," <<  colorMatrix[COLOR_CODE][1]   << "," << colorMatrix[COLOR_CODE][2] << ";";
	gp << stream.str();
	stream.str("");
	gp << "unset key; unset xtics; unset ytics; unset colorbox; unset label; \n";
	gp << "set lmargin screen 0; set rmargin screen 1; set tmargin screen 1; set bmargin screen 0; \n";
	gp << "set autoscale fix \n";
	if(PROD_SERIES){
		stream << "set term png size "  << N_X << "," << N_Y  << ";\n" ;
	}
	else{
		stream << "set term wxt size "  << N_X << "," << N_Y  << ";\n" ;
	}
	gp << stream.str();
	stream.str("");
	// if(SAVE_IMAGE){
		// gp << "set output \"fractal.png\"; \n";
	// }
	return;	
}

int limitsChanged(){
	if(abs(x_max-x_max_SAVE)>epsilon  || abs(x_min-x_min_SAVE)>epsilon  ||abs(y_max-y_max_SAVE)>epsilon  ||abs(y_min-y_min_SAVE)>epsilon){
		return 1;
	}
	else{
		return 0;
	}
}

int deactivatedZoomMode(){
	if(abs(abs(x_max-x_min) - abs(x_max_SAVE-x_min_SAVE)) < 0.05*N_X && abs(abs(y_max-y_min) - abs(y_max_SAVE-y_min_SAVE)) < 0.05*N_Y){
		return 1;
	}
	else{
		return 0;
	}
}


void calcRandomPolynomial(){
	int rand_n_terms = (rand() % (MAX_RANDOM_N_TERMS + 1 - 3)) + 3;
	for (int i = 0; i < MAX_RANDOM_N_TERMS; i++ ){
		double randExponent = (-1+2*((double)rand())/RAND_MAX)*RAND_EXPONENT_RANGE;
		double randCoeff = (-1+2*((double)rand())/RAND_MAX)*RAND_COEFF_RANGE;
		if(i>= rand_n_terms){
			randCoeff=0;
			randExponent=0;
		}
		#ifdef GPU_MODE
		POLY_COEFFS_GPU[i] = randCoeff;
		POLY_EXPONS_GPU[i] = randExponent;
		#endif

		POLY_COEFFS[i] = randCoeff;
		POLY_EXPONS[i] = randExponent;
	}

	cout << endl << " RAND_N_TERMS " << rand_n_terms << endl;
	cout << " RAND COEFFS ";

	if(FRAC_DER_MODE==0){
		for (int i = 0; i < MAX_RANDOM_N_TERMS; i++ ){
			cout << setprecision(16) << POLY_COEFFS[i]<< " " ;
		}
		cout << endl;
		cout << " RAND EXPONS ";
		for (int i = 0; i < MAX_RANDOM_N_TERMS; i++ ){
			cout << setprecision(16) << POLY_EXPONS[i] << " ";
		}
	}
	else{
		for (int i = 0; i < MAX_RANDOM_N_TERMS; i++ ){
			cout << setprecision(16) << POLY_COEFFS[i]<< " " ;
		}
		cout << endl;
		cout << " RAND EXPONS ";
		for (int i = 0; i < MAX_RANDOM_N_TERMS; i++ ){
			cout << setprecision(16) << POLY_EXPONS[i] << " ";
		}
		cout << endl << " FRAC_DER_MODE " << FRAC_DER_MODE << "\n";
		cout << " ALPHA_FRAC_DER " << ALPHA_FRAC_DER ;
	}

	int rand_color_ind = (rand() % (N_RAND_COLORS + 1 - 0)) + 0;
	cout << endl << " RAND_COLOR " << rand_color_ind << "\n";
	stream << "set palette rgbformulae " << colorMatrix[rand_color_ind][0]  << "," <<  colorMatrix[rand_color_ind][1]   << "," << colorMatrix[rand_color_ind][2] << ";";
	gp << stream.str();
	stream.str("");

	return;
}

void load_params(){
	ifstream inFile;
	inFile.open("./PARAMS.dat");
	if (!inFile) {
	cerr << "Unable to open file PARAMS.dat";
   		exit(1);   // call system to stop
	}	
	string line;
	int line_index = 0;
	while (getline(inFile, line)) {
		if(line[0] == '#' || line[0] == ' '|| line.length() == 0 ){
			continue;
		}
		// Get number of word in line
		istringstream test_ss(line);
		string test_word;
		int N_words=-1;
		while(test_ss){
			test_ss >> test_word;
			N_words+=1;
		}

		istringstream ss(line); 
		string word;
		int word_counter=1;
		while(ss){  
			ss >> word;
			// Account for last two identical bits in sting stream
			if(word_counter==N_words){
				ss >> word;
			}
			switch (line_index) {
				case 0: {
					MODE=stoi(word);
					break;
				}
				case 1: {
					REP_MODE=stoi(word);
					break;
				}
				case 2:	{	
					POLY_N_TERMS= stoi(word);
					break; }
				case 3: {
					POLY_COEFFS.push_back(stof(word));
					break; }
				case 4: {
					POLY_EXPONS.push_back(stof(word));
					break; }
				case 5: {
					FRAC_DER_MODE = stoi(word);
					break; }	
				case 6: {
					ALPHA_FRAC_DER = stof(word);
					break; }
				case 7: {
					COLOR_CODE = stoi(word);
					break; }
				case 8: {
					ACCURACY = stof(word);
					break; }
				case 9: {
					N_ITER = stoi(word);
					break; }
				case 10: {
					N_X = stoi(word);
					break; }
				case 11: {
					N_Y = stoi(word);
					break; }
				case 12: {
					if(word_counter==1){
						START_X = stof(word);
					}
					if(word_counter==2){
						START_Y = stof(word);
					}
					break; }
				case 13: {
					if(word_counter==1){
						INTERVAL_LENGTH_X = stof(word);
					}
					if(word_counter==2){
						INTERVAL_LENGTH_Y = stof(word);
					}
					break; }
				case 14: {
					SAVE_DATA_FILE = stoi(word);
					break; }
				case 15: {
					BLOCK_SIZE = stoi(word);
					break; }
				case 16: {
					MAX_RANDOM_N_TERMS = stoi(word);
					break; }
				case 17: {
					RAND_EXPONENT_RANGE = stof(word);
					break; }
				case 18: {
					RAND_COEFF_RANGE = stof(word);
					break; }
				case 19: {
					if(word_counter==1){
						PROD_SERIES = stoi(word);
					}
					if(word_counter==2){
						AUTO_ZOOM = stoi(word);
					}
					if(word_counter==3){
						ALPHA_VAR = stoi(word);
					}
					break; }
				case 20: {
					N_FRAMES_SERIES = stoi(word);
					break; }
				case 21: {
					if(word_counter==1){
						LOWER_ALPHA = stof(word);
					}
					if(word_counter==2){
						UPPER_ALPHA = stof(word);
					}
					break; }
				case 22: {
					ZOOM_FAC = stof(word);
					break; }
			}
			word_counter+=1;
		}
		line_index+=1;
	}
	if(MODE==1){
		POLY_N_TERMS=MAX_RANDOM_N_TERMS;
	}
	return;
}


void saveLimits(){
	gp.getLimits(x_min,x_max,y_min,y_max);
	x_max_SAVE=x_max;
	x_min_SAVE=x_min;
	y_max_SAVE=y_max;
	y_min_SAVE=y_min;
	return;
}

void setup( ){
	// Allocate Unified Memory – accessible from CPU or GPU
	cnumber_matrix = new complex<double>*[N_Y];
	for(int i=0; i<N_Y; i++){
		cnumber_matrix[i] = new complex<double>[N_X];
	};
	
	iter_matrix.resize(N_Y);
	for (int i=0; i<N_Y;i++){
		iter_matrix[i].resize(N_X);
	}

	roots_matrix.resize(N_Y);
	for (int i=0; i<N_Y;i++){
		roots_matrix[i].resize(N_X);
	}

	N_RAND_COLORS =  sizeof colorMatrix / sizeof colorMatrix[0];
	return;
}

void fillCnumberMatrix(complex<double>** cnumber_matrix,complex<double>lower_left,complex<double> upper_right)
{	
	dx_help = (real(upper_right)-real(lower_left))/(N_X);
	dx = {dx_help,0.0};
	dy_help = (imag(upper_right)-imag(lower_left))/(N_Y);
	dy = {0.0,dy_help};
	cnumber_matrix[0][0] = lower_left;
	for(int i = 0; i < N_Y; i++){
		for(int j = 0; j < N_X; j++){
			cnumber_matrix[i][j] = lower_left + double(j)*dx +double(i)*dy;
			// cout << dy << endl;
		}
	}
	return;
}

void sendToGP(vector<vector<float>>& image, bool saveImage, int frameID){
		string appdx;
		if(frameID >=100){
		}
		else if(frameID <100 and frameID >=10){
			appdx = "0";
		}
		else if(frameID <10){
			appdx = "00";
		}

		if(saveImage && PROD_SERIES){
			stream.str("");
			stream << "set output \"frames/series_" << appdx  << frameID << ".png\"; \n";
			gp << stream.str();
			stream.str("");
		}
		else if(saveImage && PROD_SERIES==0){
			stream.str("");
			stream << "set output \"frames/fractal.png\"; \n";
			gp << stream.str();
			stream.str("");
		}

		gp << "plot '-'  binary" << gp.binFmt2d(image, "array") << " with image \n";
		gp.sendBinary2d(image);

		stream << "set xrange [" << 0 << ":" << N_X  << "]; \n";
		gp << stream.str();
		stream.str("");
		
		stream << "set yrange [" << 0 << ":" << N_Y  << "]; \n";
		gp << stream.str();
		stream.str("");

		gp << "replot; \n";

}