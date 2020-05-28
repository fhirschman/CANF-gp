#pragma once

complex<double> cgammln(complex<double> xx);

complex<double> get_fdf(complex<double>& z, vector<double>& params, int n_terms, vector<double>& exponents);

int solve_cnewton_iter(complex<double>& start, double accuracy, int max_iterations, vector<double>& params, int n_terms, vector<double>& rationalExponents);

complex<double> solve_cnewton_roots(complex<double>& start, double accuracy, int max_iterations, vector<double>& params, int n_terms, vector<double>& rationalExponents);

void getIterMatrix(vector<vector<float>>& iterImage);

void getRootsMatrix(vector<vector<float>>& rootsImage);