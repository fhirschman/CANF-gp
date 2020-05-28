#pragma once

void freeAll();

void saveIterMatrix();

void load_params();

void setupGP();

void setup();

int limitsChanged();

int deactivatedZoomMode();

void saveLimits();

void calcRandomPolynomial();

void fillCnumberMatrix(complex<double>** cnumber_matrix,complex<double> lower_left,complex<double> upper_right);

void sendToGP(vector<vector<float>>& image, bool saveImage, int frameID);
