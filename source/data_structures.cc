using namespace std;

#include"./../headers/data_structures.h"
#include<vector>

int POLY_N_TERMS;
int N_ITER;
int N_X;
int N_Y;
int BLOCK_SIZE;
int MAX_RANDOM_N_TERMS;
int color_1;
int color_2;
int color_3;
int N_RAND_COLORS;
int COLOR_CODE;
int FRAC_DER_MODE;
int N_FRAMES_SERIES;
int mb=1;
int MODE;
int REP_MODE;
double START_X;
double START_Y;
double LOWER_ALPHA;
double UPPER_ALPHA;
double x_min,x_max,y_min,y_max,x_min_SAVE,x_max_SAVE,y_min_SAVE,y_max_SAVE;
double ACCURACY;
bool USE_ZOOM;
bool SAVE_DATA_FILE;
bool SAVE_IMAGE;
bool PRESSED_RIGHT_MOUSE;
bool PROD_SERIES;
bool AUTO_ZOOM;
bool ALPHA_VAR;
bool ZOOM_MODE_ACTIVE;
bool RECENTLY_DEACTIVATED_ZOOM_MODE=0;
bool CONTINUE_TO_NEXT_POLY=1;
bool CONTINUE_CHECK_LOOP;
double ZOOM_FAC;
double INTERVAL_LENGTH;
double epsilon=0.0000001;
double RAND_EXPONENT_RANGE;
double RAND_COEFF_RANGE;
double passed_seconds;
double mx=0.5;
double my=0.5;
double dx_help;
double dy_help;
double ALPHA_FRAC_DER;
vector<double> POLY_COEFFS;
vector<double> POLY_EXPONS;
vector<vector<int>> iter_matrix;
vector<vector<int>> roots_matrix;
complex<double> lower_left;
complex<double> upper_right;
complex<double> dx;
complex<double> dy;
complex<double>** cnumber_matrix;
Gnuplot gp;
ostringstream stream;

#ifdef GPU_MODE
cuDoubleComplex* cnumber_matrix_GPU;
cuDoubleComplex* roots_matrix_GPU;
double* POLY_COEFFS_GPU;
double* POLY_EXPONS_GPU;
int* iter_matrix_GPU;
int numBlocks;
#endif
