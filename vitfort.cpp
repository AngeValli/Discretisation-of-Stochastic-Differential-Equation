#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>

using namespace std;

const double T = 1; /* horizon */
const double sig = 0.2; /* volatilité */
const double r = 0.05; /* taux d'intérêt */
const double S0 = 100; /* valeur initial du sous-jacent */
const int M = 100000; /* nombre de simulations indépendantes */
const int P = 6; /* nombre d'itérations sur la valeur d pas */
double N = 1; /* nombre initial de pas de discrétisation */
double current_sample;
double someul;
double careul;
double sommil;
double carmil;
double dt;
double sigdt;
double rdt;
double der;
double erreul[P];
double liceul[P];
double errmil[P];
double licmil[P];
double Npas[P];
double S[M];
double Se[M];
double Sm[M];
double maxerreul[M];
double maxerrmil[M];
double p;
int l;
int n;
int k;
int i;
int j;

double sum_pow(double * arr, int arrSize, int power)
{
    double sum = 0;
    for(int i = 0;i < arrSize; i++)
    {
        sum += std::pow(arr[i],power);
    }

    return sum;
};

void create(double * erreul, double * liceul, double * errmil, double * licmil, double * Npas)
{
    // file pointer
    std::fstream fout;

    // opens an existing csv file or creates a new file.
    fout.open("vitfort.csv", ios::out | ios::app);

    // Read the input
    fout << "erreul" << " ,"
         << "liceul" << " ,"
         << "errmil" << " ,"
         << "licmil" << " ,"
         << "Npas" << " ,"
         << "\n";
    
    for (i = 0; i < P; i++) {
        // Insert the data to file
        fout << erreul[i] << ", "
             << liceul[i] << ", "
             << errmil[i] << ", "
             << licmil[i] << ", "
             << Npas[i]
             << "\n";
        };
};

int main()
{
    std::random_device rd;
    std::mt19937 gen(rd()); 

    std::fill(std::begin(erreul), std::end(erreul), 0);
    std::fill(std::begin(liceul), std::end(liceul), 0);
    std::fill(std::begin(errmil), std::end(errmil), 0);
    std::fill(std::begin(licmil), std::end(licmil), 0);
    std::fill(std::begin(Npas), std::end(Npas), 0);

    for (i=0;i<P;++i)
    {
        dt = T/N;
        sigdt = sig*std::sqrt(dt);
        rdt = r*dt;
        der = rdt - std::pow(sigdt,2)/2;

        std::fill(std::begin(S), std::end(S), S0);
        std::fill(std::begin(Se), std::end(Se), S0);
        std::fill(std::begin(Sm), std::end(Sm), S0);
        std::fill(std::begin(maxerreul), std::end(maxerreul), 0);
        std::fill(std::begin(maxerrmil), std::end(maxerrmil), 0);

        for (k=0; k<N; k++) {
            std::normal_distribution<double> d(0, 1); 
            for(l = 0; l < M; l++) {
                current_sample = d(gen);
                S[l] = S[l]*std::exp(sigdt*current_sample+der); /* filling S */
                Se[l] = Se[l]*(1+sigdt*current_sample+rdt); /* filling Se */
                Sm[l] = Sm[l]*(1+sigdt*current_sample+ der + std::pow(sigdt*current_sample,2)/2); /* filling Sm */
                maxerreul[l] = std::max(maxerreul[l],std::abs(S[l]-Se[l])); /* filling maxerreul */
                maxerrmil[l] = std::max(maxerrmil[l],std::abs(S[l]-Sm[l])); /* filling maxerrmil */
            };
        };

        someul = sum_pow(maxerreul, M, 2);
        careul = sum_pow(maxerreul, M, 4);
        sommil = sum_pow(maxerrmil, M, 2);
        carmil = sum_pow(maxerrmil, M, 4);
        erreul[i] = someul/M;
        liceul[i] = 1.96*std::sqrt((careul/M-std::pow(someul/M,2))/M);
        errmil[i] = sommil/M;
        licmil[i] = 1.96*std::sqrt((carmil/M-std::pow(sommil/M,2))/M);
        Npas[i] = N;
        N = 2*N;

    };

    create(erreul, liceul, errmil, licmil, Npas);

    return 0;
};