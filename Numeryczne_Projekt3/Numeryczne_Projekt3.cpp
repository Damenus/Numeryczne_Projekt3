//Damian Darczuk 154824
#include <iostream>
#include <math.h>
#include <fstream>

#define K 20
#define J (K - 1)
#define N 101
#define IND 4
#define START_F 2.16
#define END_F 2.66

using namespace std;

double freq[K] = { 2.160913, 2.184642, 2.208656, 2.232956, 2.257543, 2.282417, 2.307579, 2.333029,
2.358767, 2.384794, 2.411110, 2.437714, 2.464608, 2.491789, 2.519259, 2.547017, 2.575062,
2.603393, 2.632010, 2.660913 };

double s21[5][K] = {
	//indeksy 0 lub 5
	{ 0.0154473160, 0.0182357086, 0.0194462133, 0.0118513804, 0.0414972492, 0.4124997372,
	0.9972000658, 0.9942401537, 0.9938543943, 0.9845163563, 0.9975239226, 0.9947690590, 0.9844258408,
	0.9994965668, 0.8073331021, 0.4250811874, 0.2196173998, 0.1257758825, 0.0785875430, 0.0524021994 },

	//indeksy 1 lub 6
	{ 0.0536239235, 0.0802905173, 0.1282381167, 0.2232705886, 0.4300353568, 0.8098952423, 0.9994597768,
	0.9844939866, 0.9944789352, 0.9977955581, 0.9846761685, 0.9933531543, 0.9947537723, 0.9968944459,
	0.3900713171, 0.0294907390, 0.0181270647, 0.0231763396, 0.0206478723, 0.0171024633 },

	//indeksy 2 lub 7
	{ 0.1233957541, 0.1699148622, 0.2457313981, 0.3760970395, 0.5982299829, 0.8762715691, 0.9978801826,
	0.9883261673, 0.9870309473, 0.9993667568, 0.9925572200, 0.9849207469, 0.9999956457, 0.9915374521,
	0.1937142444, 0.0034410015, 0.0071060741, 0.0005756823, 0.0040972271, 0.0067222053 },

	//indeksy 3 lub 8
	{ 0.0385736053, 0.0418602337, 0.0406238069, 0.0222182757, 0.0681165327, 0.5274988276,
	0.9992090072, 0.9877714749, 0.9997958915, 0.9872738013, 0.9872738013, 0.9997958915, 0.9877714749,
	0.9992090072, 0.5274988276, 0.0681165327, 0.0222182757, 0.0406238069, 0.0418602337, 0.0385736053 },

	//indeksy 4 lub 9
	{ 0.0385736053, 0.0418602337, 0.0406238069, 0.0222182757, 0.0681165327, 0.5274988276,
	0.9992090072, 0.9877714749, 0.9997958915, 0.9872738013, 0.9872738013, 0.9997958915, 0.9877714749,
	0.9992090072, 0.5274988276, 0.0681165327, 0.0222182757, 0.0406238069, 0.0418602337, 0.0385736053 }
};

double S_dB(double S) {
	return 20 * log10(S);
}

//wyliczyæ h(z definicji), mi[1] do mi[N - 1], tak samo lambda i delta, zerowe i Nte elementy z war.brzegowych, 
//potem rozwi¹zaæ uk³ad równañ liniowych, czyli wyznaczyæ M, podstawiæ do wzoru na szukane a, b, c i d.

double h[K - 1];
double mi[K - 1];
double lambda[K - 1];
double delta[K - 1];

void wylicz_h() {
	for (int i = 0; i < K; i++) {
		h[i] = freq[i + 1] - freq[i];
	}
}

void wylicz_mi() {
	mi[0] = 0;
	mi[K - 1] = 0;
	for (int i = 1; i < K; i++) {
		mi[i] = h[i] / (h[i] + h[i + 1]);
	}
}

void wylicz_lambda() {
	lambda[0] = 0;
	for (int i = 1; i < K; i++) {
		lambda[i] = h[i + 1] / (h[i] + h[i + 1]);
	}
}

void wylicz_delta() {
	delta[0] = 0;
	delta[K - 1] = 0;
	for (int i = 1; i < K; i++) {
		delta[i] = (6.0 / (h[i] + h[i + 1])) * (((s21[IND][i + 1] - s21[IND][i]) / h[i + 1]) - ((s21[IND][i] - s21[IND][i - 1]) / h[i]));
	}
}

void warunki_brzegowe_a() {
	lambda[0] = 0.0;
	delta[0] = 0.0;
	mi[K - 1] = 0.0;
	delta[K - 1] = 0.0;
}

void warunki_brzegowe_b() {
	lambda[0] = 1.0;
	mi[K - 1] = 1.0;

	//delta[0] = 6.0 / h[1] * (((s21[IND][1] - s21[IND][0]) / h[1]) - 1); //pochodna
	//delta[K - 1] = 6.0 / h[K - 1] * (1 - ((s21[IND][K-1] - s21[IND][K-2]) / h[K-1])); //pochodna
}

double a[K - 1];
double b[K - 1];
double c[K - 1];
double d[K - 1];

//niewiadoma!!!
struct M{
	double wartosc = 0.0;
	double wspolczynik = 1.0;
	// b + 2 * c * ( x1 -x0) + 3 * d * (x1 - x0) ^ 2
}Mn[K];

double M[K];
double macierz_wspolczynikow[K][K+1];

void stworz_uklad_rownan() {

	//wyzerowanie pól, wpisanie wartoœci do kolumny wyrazów wolnych
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < K; j++) {
			macierz_wspolczynikow[i][j] = 0.0;
		}
		macierz_wspolczynikow[i][K] = delta[i];
	}

	//wype³nienie wspó³czynikami
	for (int i = 0; i < K; i++) {
		if (i != K-1)
			macierz_wspolczynikow[i][i + 1] = lambda[i]; //dobrze?

		macierz_wspolczynikow[i][i] = 2.0;

		if (i != 0)
			macierz_wspolczynikow[i][i - 1] = mi[i]; //dobrze?
	}

}

void gauss() {

	//eliminacja
	double m;
	for (int i = 0; i < K; i++) { // ka¿da kolumna pokolei

		for (int j = i + 1; j < K; j++) {
			m = macierz_wspolczynikow[j][i] / macierz_wspolczynikow[i][i]; //wspo³czynnik mno¿enia pierwszego weirsza, by odj¹æ od kolejnych kolumn, by wyzerowaæ dolny trójk¹t
			for (int k = i; k < K + 1; k++) { // zerowanie kolumny
				macierz_wspolczynikow[j][k] = macierz_wspolczynikow[j][k] - m * macierz_wspolczynikow[i][k];
			}
		}
	}

	double s;
	//obliczanie
	//int n = K - 1;
	M[K - 1] = macierz_wspolczynikow[K - 1][K] / macierz_wspolczynikow[K - 1][K - 1];

	for (int i = K - 1 - 1; i >= 0; i--){
		s = 0;
		for (int k = i + 1; k < K; k++) {
			s = s + macierz_wspolczynikow[i][k] * M[k];
			M[i] = (macierz_wspolczynikow[i][K] - s) / macierz_wspolczynikow[i][i];
		}
	}

}

void wyznacz_wspolczyniki_abcd() {
	for (int i = 0; i < K; i++) {
		a[i] = s21[IND][i];
		b[i] = ((s21[IND][i + 1] - s21[IND][i]) / h[i + 1]) - (((2.0 * M[i] - M[i + 1]) / 6.0) * h[i + 1]);
		c[i] = M[i] / 2.0;
		d[i] = (M[i + 1] - M[i]) / (6.0 * h[i + 1]);
	}
}

double X[N];
double Y[N];

void wzynacz_X() {
	//liczye co ile bedzie X
	double krok = (END_F - START_F) / N;
	//wyznaczam wszystkie X
	for (int i = 0; i < N; i++) { //i mo¿e zaczynaæ siê od 1
		X[i] = START_F + (i * krok);
	}
}

void oblicz_Y() {
	for (int i = 0; i < N; i++) { 
		for (int j = 0; j < K; j++) {
			if (X[i] > freq[j] && X[i] < freq[j + i]) {
				Y[i] = a[j] * pow(X[i], 3) + b[j] * pow(X[i], 2) + c[j] * X[i] + d[j];
				break;
			}
		}
	}

}

int main()
{
	wylicz_h();
	wylicz_mi();
	wylicz_lambda();
	wylicz_delta();
	warunki_brzegowe_a();

	stworz_uklad_rownan();
	gauss();

	wyznacz_wspolczyniki_abcd();

	wzynacz_X();
	oblicz_Y();

	for (int i = 0; i < K; i++) {
		cout << i << ": " << M[i] << endl;
	}

	return 0;
}

