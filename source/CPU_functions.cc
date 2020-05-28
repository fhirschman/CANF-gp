using namespace std;

#include<complex>
#include<vector>
#include"./../headers/CPU_functions.h"
#include"./../headers/data_structures.h"


// Complex gamma function
complex<double> cgammln(complex<double> xx)
{
	complex<double> x,y,tmp,ser,res;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;

	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++){
		y+=1;
		ser += cof[j]/y;
	}
	res=-tmp+log(2.5066282746310005*ser/x);
	res={real(res),0.0};
	return res;
}

complex<double> get_fdf(complex<double>& z, vector<double>& params, int n_terms, vector<double>& exponents){
	complex<double> res = 0;  
	complex<double> f = 0;
	double rationalExponent = 0;
	for(int i = n_terms-1; i>-1;i--){
		rationalExponent = exponents[i];
		complex<double> temp1 = (pow(z,rationalExponent));
		f += double(params[i])*temp1;
	}
	complex<double> df = 0;
	for(int i = n_terms-1; i>-1;i--){
		rationalExponent = exponents[i];
		// Normal derivative
		if(FRAC_DER_MODE==0){
			if(params[i] != 0){
				complex<double> temp2 = (pow(z,rationalExponent-1));
				df += double(params[i]*(rationalExponent))*temp2;
			}
		}
		// Riemann-Liouville derivative
		else{
			rationalExponent = abs(int(exponents[i]));
			if((abs(params[i]) > 0.00000000001)){
				if( (abs(rationalExponent)  > 0.00000000001)) {
				 	complex<double> temp2 = (pow(z,int(rationalExponent)-ALPHA_FRAC_DER));
					df +=  ( exp(cgammln(rationalExponent)+1.0)/exp(cgammln(rationalExponent-ALPHA_FRAC_DER+1.0)))*temp2;
				}
				else{
				 	complex<double> temp2 = params[i]*(pow(z,-ALPHA_FRAC_DER));
					df +=  temp2/exp(cgammln(1.0-ALPHA_FRAC_DER));
				}
			}
		}
	}
	res = -(f/df);
	return res;
}

int solve_cnewton_iter(complex<double>& start, double accuracy, int max_iterations, vector<double>& params, int n_terms, vector<double>& rationalExponents){
	int res=0;
	complex<double> z, help;
	z = start;
	for(int i = 0; i < max_iterations; i++){
		help = get_fdf(z,params,n_terms,rationalExponents);
		z = z + help;
		if(abs(help) < accuracy){
			res= i;
			break;
		}
		else{
			res= max_iterations;
		}
	}
	return res;
}

void getIterMatrix(vector<vector<float>>& iterImage){
	for(int i=0; i<N_Y; i++){
		vector<float> row;
		for(int j=0; j<N_X; j++){
			iter_matrix[i][j]=solve_cnewton_iter(cnumber_matrix[i][j], ACCURACY, N_ITER, POLY_COEFFS, POLY_N_TERMS, POLY_EXPONS);
			row.push_back(log(iter_matrix[i][j]));
		}
		iterImage.push_back(row);
	}
}


complex<double> solve_cnewton_roots(complex<double>& start, double accuracy, int max_iterations, vector<double>& params, int n_terms, vector<double>& rationalExponents){
	int res=0;
	complex<double> z, help;
	z = start;
	for(int i = 0; i < max_iterations; i++){
		help = get_fdf(z,params,n_terms,rationalExponents);
		z = z + help;
		if(abs(help) < accuracy){
			res= i;
			break;
		}
		else{
			res= max_iterations;
		}
	}
	return z;
}


void getRootsMatrix(vector<vector<float>>& rootsImage){
	for(int i=0; i<N_Y; i++){
		for(int j=0; j<N_X; j++){
			iter_matrix[i][j]=solve_cnewton_iter(cnumber_matrix[i][j], ACCURACY, N_ITER, POLY_COEFFS, POLY_N_TERMS, POLY_EXPONS);
		}
	}
	vector<complex<double>> roots;
	for(int i=0; i<N_Y; i++){
		vector<float> row;
		for(int j=0; j<N_X; j++){
			complex<double> root = solve_cnewton_roots(cnumber_matrix[i][j], ACCURACY, N_ITER, POLY_COEFFS, POLY_N_TERMS, POLY_EXPONS);
			if(abs(real(root)) > 10*INTERVAL_LENGTH || abs(imag(root)) > 10*INTERVAL_LENGTH ){
				root = {0,0};
			} 
			if(i==0 && j==0){
				roots.push_back(root);
			}
			bool new_root = true;
			for (int k=0; k< roots.size(); k++){
				if( abs(real(root)-real(roots[k]))<0.01 && abs(imag(root)-imag(roots[k]))<0.01){
					roots_matrix[i][j]=k;
					row.push_back(roots_matrix[i][j]);
					new_root=false;
					break;
				}
			}
			if(new_root==true){
				roots_matrix[i][j]=roots.size()+1;
				row.push_back(roots_matrix[i][j]);
				roots.push_back(root);
			}
		}
		rootsImage.push_back(row);
	}
	cout << "Found " << roots.size() << " roots:" << endl;
	for (int k=0; k< roots.size(); k++){
		cout << " " << real(roots[k]) << " " << imag(roots[k]) << endl;
	}
	return;
}