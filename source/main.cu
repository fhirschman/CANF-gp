/*

A program to generate fractal images by using the Newton method for root 
finding of complex polynomials. The generated data is being piped to gnuplot 
for visualization. CPU and GPU support (for CUDA enabled graphic cards newer 
than Pascal generation due to unified memory) is provided.

#######

A slightly modified version of gnuplot-iostream is used in headers/gnuplot-iostream.h:

https://github.com/dstahlke/gnuplot-iostream
Copyright (c) 2013 Daniel Stahlke (dan@stahlke.org)
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

#######

*/

using namespace std;
#include <thread>

#define GNUPLOT_ENABLE_PTY
#include "./../headers/gnuplot-iostream.h"

#ifdef GPU_MODE
#include <cuda.h>
#include "./../headers/CUDA_functions.cuh"
#else
#include "./../headers/CPU_functions.h"
#endif
#include "./../headers/data_structures.h"
#include "./../headers/general_functions.h"
#include "./../headers/colors.h"

int main(void)
{	
	// Randomize seed
 	srand(time(NULL));	
 	// Load PARAMS.dat
	load_params();
	// Allocate data structures
	setup();
	// Set up gnuplot
	setupGP();

	#ifdef GPU_MODE
	// Allocate GPU data structures
	setup_GPU();
	#endif

	// Set up complex plane limits
	lower_left = { (-INTERVAL_LENGTH)*0.5 + START_X, (-INTERVAL_LENGTH)*0.5 +START_Y};
	upper_right ={ ( INTERVAL_LENGTH)*0.5 + START_X, ( INTERVAL_LENGTH)*0.5 +START_Y};


	// Interactive mode, no series of frames get produced
	if(PROD_SERIES==0){

		ZOOM_MODE_ACTIVE=0;
		RECENTLY_DEACTIVATED_ZOOM_MODE=0;

		// Big endless loop
		while(true){

			// If MODE==1 ("random mode"): generate a random complex polynomial
			if(MODE==1 && CONTINUE_TO_NEXT_POLY){
				calcRandomPolynomial();

				// If coming out of zoom mode, reset boundaries of complex plane
				if(RECENTLY_DEACTIVATED_ZOOM_MODE){
					lower_left = { (-INTERVAL_LENGTH)*0.5 + START_X, (-INTERVAL_LENGTH)*0.5 +START_Y};
					upper_right ={ ( INTERVAL_LENGTH)*0.5 + START_X, ( INTERVAL_LENGTH)*0.5 +START_Y};
					RECENTLY_DEACTIVATED_ZOOM_MODE=0;
				}
			}

			// Set up complex plane 
			fillCnumberMatrix(cnumber_matrix,lower_left,upper_right);
			#ifdef GPU_MODE
			fillCnumberMatrixGPU();
			#endif

			// Calculate image
			chrono::steady_clock::time_point begin = chrono::steady_clock::now();
			vector<vector<float> > image;
			#ifdef GPU_MODE
			// REP_MODE=0: calculate convergence matrix, REP_MODE==1 : calculate basins of attractions
			if(REP_MODE==1){
				getRootsMatrixGPU(image);
			}
			else if(REP_MODE==0){
				getIterMatrixGPU(image);
			}
			#else
			if(REP_MODE==1){
				getRootsMatrix(image);
			}
			else if(REP_MODE==0){
				getIterMatrix(image);
			}
			#endif 
			chrono::steady_clock::time_point end = chrono::steady_clock::now(); 


			if(MODE==0){

					cout << endl << " POLY_N_TERMS " << POLY_N_TERMS << endl;
					cout << " COEFFS ";

				if(FRAC_DER_MODE==0){
					for (int i = 0; i < POLY_N_TERMS; i++ ){
						cout << setprecision(16) << POLY_COEFFS[i]<< " " ;
					}
					cout << endl;
					cout << " RAND EXPONS ";
					for (int i = 0; i < POLY_N_TERMS; i++ ){
						cout << setprecision(16) << POLY_EXPONS[i] << " ";
					}
				}
				else{
					for (int i = 0; i < POLY_N_TERMS; i++ ){
						cout << setprecision(16) << POLY_COEFFS[i]<< " " ;
					}
					cout << endl;
					cout << " RAND EXPONS ";
					for (int i = 0; i < POLY_N_TERMS; i++ ){
						cout << setprecision(16) << POLY_EXPONS[i] << " ";
					}
					cout << endl << " FRAC_DER_MODE " << FRAC_DER_MODE << "\n";
					cout << " ALPHA_FRAC_DER " << ALPHA_FRAC_DER ;
				}
				cout << endl;
			}

			cout << " CALC TIME = " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " [ms]" << endl;

			// Sent iteration matrix to gnuplot
			sendToGP(image,SAVE_DATA_FILE,0);

			// Save current gnuplot limits
			saveLimits();
			
			// Check for user inputs
			CONTINUE_CHECK_LOOP=1;
			RECENTLY_DEACTIVATED_ZOOM_MODE=0;
			while(CONTINUE_CHECK_LOOP){

				JUMP_POINT:

				// Get current gnuplot limits
				gp.getLimits(x_min,x_max,y_min,y_max);

				// Act, if limits have changed
				if(limitsChanged()){
					cout << endl << " NEW LIMITS: " << x_min << " " << x_max << " " << y_min <<" " << y_max << endl;

					// For translation in x-direction
					if(x_max>x_max_SAVE || x_min<x_min_SAVE ){

						lower_left =  { real(lower_left) +  (x_max-x_max_SAVE)*dx_help , imag(lower_left)  };
						upper_right = { real(upper_right) + (x_max-x_max_SAVE)*dx_help , imag(upper_right) };

						CONTINUE_CHECK_LOOP=0;
						CONTINUE_TO_NEXT_POLY=0;
					}
					// For translation in y-direction
					else if(y_max>y_max_SAVE || y_min<y_min_SAVE ){

						lower_left =  { real(lower_left)  , imag(lower_left) +  (y_max-y_max_SAVE)*dy_help };
						upper_right = { real(upper_right) , imag(upper_right) + (y_max-y_max_SAVE)*dy_help };

						CONTINUE_CHECK_LOOP=0;
						CONTINUE_TO_NEXT_POLY=0;
					}
					// User deactivated zoom mode
					else if(deactivatedZoomMode()){
						CONTINUE_CHECK_LOOP=1;
						CONTINUE_TO_NEXT_POLY=1;
						ZOOM_MODE_ACTIVE=0;
						RECENTLY_DEACTIVATED_ZOOM_MODE=1; // Important, to reset viewing range for next polynomial!
						cout << "\n ZOOM MODE DEACTIVATED: LEFT MOUSE = Next fractal, RIGHT MOUSE = Reactivate ZOOM MODE \n " << endl;
					}
					// User wants to zoom in
					else{
						lower_left  = { real(cnumber_matrix[int(y_min)][int(x_min)]) , imag(cnumber_matrix[int(y_min)][int(x_min)])};
						upper_right = { real(cnumber_matrix[int(y_max)][int(x_max)]) , imag(cnumber_matrix[int(y_max)][int(x_max)])};
						CONTINUE_CHECK_LOOP=0;
						CONTINUE_TO_NEXT_POLY=0;
					}
					cout << " INTERVAL LENGTH " <<  real(upper_right) - real(lower_left) << " " << imag(upper_right) - imag(lower_left) << endl;
					saveLimits();
				}

				// Pass some time
				this_thread::sleep_for(200ms);
				
				// Save iter matrix if wanted
				if(SAVE_DATA_FILE){
					saveIterMatrix();
				}

				if(not(ZOOM_MODE_ACTIVE)){
					// Wait for user input
					gp.getMouse(mx, my, mb, " LEFT MOUSE = Next fractal, RIGHT MOUSE = Stay, activate ZOOM MODE");
					printf("Pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
					// If user activates zoom mode
					if(mb==3){
						ZOOM_MODE_ACTIVE=1;
						cout << "\n ZOOM MODE ACTIVATED: RIGHT MOUSE = Zoom to selection. Deactivate with very big selection! \n " << endl;
						goto JUMP_POINT;
					}
					break;
				}
			}
		}
	}

	// Production mode, series of frames get produced. Create gif with convert -delay 1 -loop 0 ./frames/*.png animated.gif, or ffmpeg -i ./frames/series_%03d.png -vf scale=512x512 output.gif 


	else if(PROD_SERIES==1){

		cout << " Coefficients: ";
		for (int i = 0; i < POLY_N_TERMS; i++ ){
			cout << setprecision(8) << POLY_COEFFS[i]<< " " ;
		}
		cout << endl;
		cout << " Exponents: ";
		for (int i = 0; i < POLY_N_TERMS; i++ ){
			cout << setprecision(8) << POLY_EXPONS[i] << " ";
		}
		cout << endl << endl;
		
		for (int frameInd=0; frameInd<N_FRAMES_SERIES; frameInd++) {
			cout << " Frame " << frameInd+1 << " of " << N_FRAMES_SERIES << endl;

			cout << lower_left << " " << upper_right << endl;

			if(AUTO_ZOOM==1){
				lower_left = { (-INTERVAL_LENGTH*pow(ZOOM_FAC,frameInd))*0.5  + START_X ,(-INTERVAL_LENGTH*pow(ZOOM_FAC,frameInd))*0.5 + START_Y};
				upper_right ={ ( INTERVAL_LENGTH*pow(ZOOM_FAC,frameInd))*0.5  + START_X ,( INTERVAL_LENGTH*pow(ZOOM_FAC,frameInd))*0.5 + START_Y};
				cout << " Zoom Fac " << pow(ZOOM_FAC,frameInd) << endl;
			} 

			if(ALPHA_VAR && FRAC_DER_MODE){
				ALPHA_FRAC_DER = LOWER_ALPHA + frameInd*(UPPER_ALPHA-LOWER_ALPHA)/(N_FRAMES_SERIES-1.0);
				cout << " Alpha  " << ALPHA_FRAC_DER << endl;
			}

			// Set up complex plane 
			fillCnumberMatrix(cnumber_matrix,lower_left,upper_right);
			#ifdef GPU_MODE
			fillCnumberMatrixGPU();
			#endif

			// Calculate image
			chrono::steady_clock::time_point begin = chrono::steady_clock::now();
			vector<vector<float> > image;
			#ifdef GPU_MODE
			if(REP_MODE==true){
				getRootsMatrixGPU(image);
			}
			else{
				getIterMatrixGPU(image);
			}
			#else
			if(REP_MODE==true){
				getRootsMatrix(image);
			}
			else{
				getIterMatrix(image);
			}
			#endif 
			chrono::steady_clock::time_point end = chrono::steady_clock::now(); 
			cout << " CALC TIME = " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " [ms]" << endl << endl;

			// Sent iteration matrix to gnuplot
			sendToGP(image,1,frameInd);

			// Save current gnuplot limits
			saveLimits();
		}
	}
	else{}

	// Free, close everything
	freeAll();
	#ifdef GPU_MODE
	freeAll_GPU();
	#endif

	return 0;
}
