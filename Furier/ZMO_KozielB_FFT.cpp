////////////////////////////////////
//Fast Fourier Tranform
//
//Bart³omiej Kozie³
//8.01.2020
////////////////////////////////////

#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>

double M_PI = 3.141592653589793238460;
using namespace std;

void fft(complex<double>* x, int N);
void ifft(complex<double>* x, int N);

int main() {

		int N;
		double im, re;
		
		//Reading input data
		cin >> N;
		complex<double>* input = new complex<double>[N];
		complex<double>* output = new complex<double>[N];

		for (int k = 0; k < N; ++k) {
			cin >> re;
			cin >> im;
			input[k] = complex<double>(re, im);
			output[k] = complex<double>(re, im);
		}
		
		fft(output, N);
		//Displaying results
		for (int n = 0; n < N; ++n) {
			cout << n << " " << fixed << setprecision(4) << output[n] << endl;
		}
		cout << "\n";
		
		ifft(output, N);
		//Displaying results
		for (int n = 0; n < N; ++n) {
			cout << n << " " << fixed << setprecision(4) << output[n] << endl;
		}
	return 0;
}
//Fast Fourier Tranform
void fft(complex<double> * x, int N) {
	if (N <= 1) return;

	complex<double>* A0 = new complex<double>[N / 2];
	complex<double>* A1 = new complex<double>[N / 2];
	
	//Divide collection in half
	for (int i = 0; i < N / 2; ++i) {
		A0[i] = x[i * 2];
		A1[i] = x[i * 2 + 1];
	}

	//Recursive function call
	fft(A0, N / 2);
	fft(A1, N / 2);

	//Discrete Fourier Transform
	complex<double> f;

	for (int k = 0; k < N / 2; ++k) {
		f = exp(complex<double>(0, -2.0 * M_PI * k / N));
		x[k] = A0[k] + A1[k] * f;
		x[N / 2 + k] = A0[k] - A1[k] * f;
	}
}
//Inverse Fast Fourier Tranform
void ifft(complex<double>* x, int N) {
	//conjugate the complex numbers
	for (int i = 0; i < N; ++i) {
		x[i] *= complex<double>(0, -1);
	}
	
	fft(x, N);

	//conjugate the complex numbers
	for (int i = 0; i < N; ++i) {
		x[i] *= complex<double>(0, -1);
		x[i] /= N;
	}
}