#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>

using namespace std;

const double T = 1; /* Time horizon */
const double sig = 0.2; /* Volatility */
const double r = 0.05; /* Interest rates */
const double K = 110;
const double S0 = 100; /* Initial value of the underlying */
const double R = 10;
const int M = 1000000; /* Number of independent simulations */
const double RM = R*M;
const int P = 4; /* Number of iterations of d steps */
double N = 1; /* Number of initial discretisation steps */
double current_sample;
double current_sample_1;
double current_sample_2;
double payeul;
double careul;
double paymil;
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
double S2e[M];
double Sm[M];
double S2m[M];
double eul[M];
double mil[M];
double p;
int l;
int n;
int k;
int i;
int j;

double sum(double *arr, int arrSize)
{
    double sum = 0;
    for(int i = 0;i < arrSize; i++)
    {
        sum += arr[i];
    }

    return sum;
};


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
    fout.open("romberg.csv", ios::out | ios::app);

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

        payeul=0;
        careul=0;
        paymil=0;
        carmil=0;
        for (j=0;j<R;++j) {
            std::fill(std::begin(S), std::end(S), S0);
            std::fill(std::begin(Se), std::end(Se), S0);
            std::fill(std::begin(S2e), std::end(S2e), S0);
            std::fill(std::begin(Sm), std::end(Sm), S0);
            std::fill(std::begin(S2m), std::end(S2m), S0);
            for (k=0; k<N; k++) {
                std::normal_distribution<double> d1(0, 1); 
                std::normal_distribution<double> d2(0, 1); 
                for(l = 0; l < M; l++) {
                    current_sample_1 = d1(gen)/std::sqrt(2);
                    current_sample_2 = d2(gen)/std::sqrt(2);
                    current_sample = current_sample_1 + current_sample_2;
                    S[l] = S[l]*std::exp(sigdt*current_sample+der); /* filling S */
                    Se[l] = Se[l]*(1+sigdt*current_sample+rdt); /* filling Se */
                    S2e[l] = S2e[l]*(1+sigdt*current_sample_1+rdt/2)*(1+sigdt*current_sample_2+rdt/2); /* filling S2e */
                    Sm[l] = Sm[l]*(1+sigdt*current_sample+ der + std::pow(sigdt*current_sample,2)/2); /* filling Sm */
                    S2m[l] = S2m[l]*(1+sigdt*current_sample_1 + std::pow(sigdt*current_sample_1,2)/2 + der/2)*(1+sigdt*current_sample_2 + std::pow(sigdt*current_sample_2,2)/2 + der/2); /* filling S2m */
                };
            };
            
            for(l = 0; l < M; l++) {
                eul[l] = 2*std::max(K-S2e[l],double(0)) - std::max(K-Se[l],double(0)) - std::max(K-S[l],double(0));
                mil[l] = 2*std::max(K-S2m[l],double(0)) - std::max(K-Sm[l],double(0)) - std::max(K-S[l],double(0));
            };
            payeul = payeul + sum(eul,M)/M;
            careul = careul + sum_pow(eul, M, 2)/M;
            paymil = paymil + sum(mil,M)/M;
            carmil = carmil + sum_pow(mil, M, 2)/M;
        };

        erreul[i] = std::exp(-r*T)*payeul/R;
        liceul[i] = 1.96*std::exp(-r*T)*std::sqrt((careul/R-std::pow(payeul/R,2))/RM);
        errmil[i] = std::exp(-r*T)*paymil/R;
        licmil[i] = 1.96*std::exp(-r*T)*std::sqrt((carmil/R -std::pow(paymil/R,2))/RM);
        Npas[i] = N;
        N = 2*N;
    };

    create(erreul, liceul, errmil, licmil, Npas);

    return 0;
};