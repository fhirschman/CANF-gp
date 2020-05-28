#pragma once
using namespace std;
#include<vector>
#include "./../headers/gnuplot-iostream.h"

extern int colorMatrix[53][3];
extern int POLY_N_TERMS;
extern int N_ITER;
extern int N_X;
extern int N_Y;
extern int BLOCK_SIZE;
extern int MAX_RANDOM_N_TERMS;
extern int color_1;
extern int color_2;
extern int color_3;
extern int N_RAND_COLORS;
extern int COLOR_CODE;
extern int FRAC_DER_MODE;
extern int N_FRAMES_SERIES;
extern int mb;
extern int MODE;
extern int REP_MODE;
extern double START_X;
extern double START_Y;
extern double LOWER_ALPHA;
extern double UPPER_ALPHA;
extern double x_min,x_max,y_min,y_max,x_min_SAVE,x_max_SAVE,y_min_SAVE,y_max_SAVE;
extern double ACCURACY;
extern bool USE_ZOOM;
extern bool SAVE_DATA_FILE;
extern bool SAVE_IMAGE;
extern bool PRESSED_RIGHT_MOUSE;
extern bool PROD_SERIES;
extern bool AUTO_ZOOM;
extern bool ALPHA_VAR;
extern bool ZOOM_MODE_ACTIVE;
extern bool RECENTLY_DEACTIVATED_ZOOM_MODE;
extern bool CONTINUE_TO_NEXT_POLY;
extern bool CONTINUE_CHECK_LOOP;
extern double ZOOM_FAC;
extern double INTERVAL_LENGTH;
extern double epsilon;
extern double RAND_EXPONENT_RANGE;
extern double RAND_COEFF_RANGE;
extern double passed_seconds;
extern double mx;
extern double my;
extern double dx_help;
extern double dy_help;
extern double ALPHA_FRAC_DER;
extern vector<double> POLY_COEFFS;
extern vector<double> POLY_EXPONS;
extern vector<vector<int>> iter_matrix;
extern vector<vector<int>> roots_matrix;
extern complex<double> lower_left;
extern complex<double> upper_right;
extern complex<double> dx;
extern complex<double> dy;
extern complex<double>** cnumber_matrix;
extern Gnuplot gp;
extern ostringstream stream;

#ifdef GPU_MODE
#include <cuComplex.h>
extern cuDoubleComplex* cnumber_matrix_GPU;
extern cuDoubleComplex* roots_matrix_GPU;
extern double* POLY_COEFFS_GPU;
extern double* POLY_EXPONS_GPU;
extern int* iter_matrix_GPU;
extern int numBlocks;
#endif