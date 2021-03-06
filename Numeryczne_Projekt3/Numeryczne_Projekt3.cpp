//Damian Darczuk 154824
#include <iostream>
#include <math.h>
#include <fstream>

#define K 20 //liczba w�z��w
#define N 101 //liczba warto�ci charakterystyki przybli�onej
#define IND 4  //ostatnia cyfra indeksu do wybrania dobrego zestawu danych
#define START_F 2.16 //pocz�tek zakresu
#define END_F 2.66 //koniec zakresu

using namespace std;
//X, cz�stotliwo�ci w kt�rych zmierzono warto�ci, inaczej w�z�y
double freq[K] = { 2.160913, 2.184642, 2.208656, 2.232956, 2.257543, 2.282417, 2.307579, 2.333029,
2.358767, 2.384794, 2.411110, 2.437714, 2.464608, 2.491789, 2.519259, 2.547017, 2.575062,
2.603393, 2.632010, 2.660913 };
//Y, mzierzone warto�ci
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
// funkcja logarytmiczna, pozwalaj�ca zobrazowa� wyniki na wykresie
double S_dB(double S) {
	return 20 * log10(S);
}

//wyliczy� h(z definicji), mi[1] do mi[N - 1], tak samo lambda i delta, zerowe i Nte elementy z war.brzegowych, 
//potem rozwi�za� uk�ad r�wna� liniowych, czyli wyznaczy� M, podstawi� do wzoru na szukane a, b, c i d.

//wsp�czynniki do stworzenia uk�adu r�wna�, K-1(19) r�na� mi�dzy w�z�ami
double h[K]; //odleg�o�ci mi�dzy w�z�ami, indeksy od 1
double mi[K - 1]; 
double lambda[K - 1];
double delta[K - 1];

void wylicz_h() {
	for (int i = 0; i < K; i++) {
		h[i+1] = freq[i + 1] - freq[i];
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

double a[K - 1];
double b[K - 1];
double c[K - 1];
double d[K - 1];

double M[K]; //niewiadoma
double macierz_wspolczynikow[K][K+1];

void stworz_uklad_rownan() {

	//wyzerowanie p�l oraz wpisanie warto�ci do kolumny wyraz�w wolnych
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < K; j++) {
			macierz_wspolczynikow[i][j] = 0.0;
		}
		macierz_wspolczynikow[i][K] = delta[i];
	}

	//wype�nienie wsp�czynikami
	for (int i = 0; i < K; i++) {
		if (i != K-1)
			macierz_wspolczynikow[i][i + 1] = lambda[i];

		macierz_wspolczynikow[i][i] = 2.0;

		if (i != 0)
			macierz_wspolczynikow[i][i - 1] = mi[i]; 
	}

}

void gauss() {

	//eliminacja
	double m;
	for (int i = 0; i < K; i++) { 
		for (int j = i + 1; j < K; j++) {
			m = macierz_wspolczynikow[j][i] / macierz_wspolczynikow[i][i]; //wspo�czynnik mno�enia pierwszego weirsza, by odj�� od kolejnych kolumn, by wyzerowa� dolny tr�jk�t
			for (int k = i; k < K + 1; k++) { // zerowanie kolumny
				macierz_wspolczynikow[j][k] = macierz_wspolczynikow[j][k] - m * macierz_wspolczynikow[i][k];
			}
		}
	}

	//obliczanie
	double s;	
	M[K - 1] = macierz_wspolczynikow[K - 1][K] / macierz_wspolczynikow[K - 1][K - 1]; //w macierzy z wyzerowanej dolnym tr�jk�tem, w ostatnim wierszu na przk�tnej(ostatnia kolumna) zostaje jeden wsp�czynnik, co pozwala wyliczenie ostatniej niewiadomej  

	for (int i = K - 1 - 1; i >= 0; i--){ //liczenie od ko�ca kolejnych niewiadmoych
		s = 0;
		for (int k = i + 1; k < K; k++) {
			s = s + macierz_wspolczynikow[i][k] * M[k];
			M[i] = (macierz_wspolczynikow[i][K] - s) / macierz_wspolczynikow[i][i];
		}
	}

}

void wyznacz_wspolczyniki_abcd() {
	for (int i = 0; i < K - 1; i++) {
		a[i] = s21[IND][i];
		b[i] = ((s21[IND][i + 1] - s21[IND][i]) / h[i + 1]) - ((((2.0 * M[i]) + M[i + 1]) / 6.0) * h[i + 1]);
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
	for (int i = 0; i < N; i++) {  
		X[i] = START_F + ((1+i) * krok);
	}
}

void oblicz_Y() {
	for (int i = 0; i < N; i++) { 
		for (int j = 0; j < K; j++) {
			if (X[i] >= freq[j] && X[i] < freq[j + 1]) { //znajd� kt�rej funckji przybli�aj�cej u�y� i wylicz warto��
				Y[i] = (d[j] * pow(X[i] - freq[j], 3)) + (c[j] * pow(X[i] - freq[j], 2)) + (b[j] * (X[i] - freq[j])) + a[j];				
				break;
			}
		}
	}

}

void zapis_do_pliku() {

	fstream plik("plik.txt", ios::out);
	
	if (plik.good())
	{
		for (int i = 0; i < N; i++)
		{
			plik << S_dB(Y[i]) << endl; //u�yj logarytmu by uzyska� wynik do pokazania na wykresie
			plik.flush();
		}

		plik << endl;

		for (int i = 0; i < N; i++)
		{
			plik << X[i] << endl;
			plik.flush();
		}

		plik.close();
	}
}

void wypisz_M() { //wypisz M (niewiadome) do por�wnania wynik�w
	for (int i = 0; i < K; i++) {
		cout << i << ": " << M[i] << endl;
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
	zapis_do_pliku();
	
	return 0;
}