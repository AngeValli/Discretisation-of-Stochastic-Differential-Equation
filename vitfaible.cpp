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
const int M = 100000; /* Number of independent simulations */
const int P = 6; /* Number of iterations of d steps */
double N = 1; /* Number of initial discretisation steps */
const double d1 = (std::log(S0/K)+r*T)/(sig*std::sqrt(T)) + sig*std::sqrt(T)/2;
const double d2 = d1 - sig*std::sqrt(T);
double BS;
double current_sample;
double contcareul;
double contputeul;
double puteul;
double careul;
double contcarmil;
double contputmil;
double putmil;
double carmil;
double dt;
double sigdt;
double rdt;
double der;
double erreul[P];
double liceul[P];
double conterreul[P];
double contliceul[P];
double errmil[P];
double licmil[P];
double conterrmil[P];
double contlicmil[P];
double Npas[P];
double S[M];
double Se[M];
double Sm[M];
double paymc[M];
double payeul[M];
double paymil[M];
double diffpaymilpaymc[M];
double diffpayeulpaymc[M];
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

double normalCDF(double value)
{
   return 0.5 * std::erfc(-value * M_SQRT1_2);
}

double sum_pow(double * arr, int arrSize, int power)
{
    double sum = 0;
    for(int i = 0;i < arrSize; i++)
    {
        sum += std::pow(arr[i],power);
    }

    return sum;
};

double * diff_array(double * new_arr, double * arr1, double * arr2, int arrSize) {

    for(int i = 0;i < arrSize; i++)
    {
        new_arr[i] = arr1[i] - arr2[i];
    };

    return new_arr;
}

void create(double * erreul, double * liceul, double * conterreul, double * contliceul, double * errmil, double * licmil, double * conterrmil, double * contlicmil, double * Npas)
{
    // file pointer
    std::fstream fout;

    // opens an existing csv file or creates a new file.
    fout.open("vitfaible.csv", ios::out | ios::app);

    // Read the input
    fout << "erreul" << " ,"
         << "liceul" << " ,"
         << "conterreul" << " ,"
         << "contliceul" << " ,"
         << "errmil" << " ,"
         << "licmil" << " ,"
         << "conterrmil" << " ,"
         << "contlicmil" << " ,"
         << "Npas" << " ,"
         
         << "\n";
    
    for (i = 0; i < P; i++) {
        // Insert the data to file
        fout << erreul[i] << ", "
             << liceul[i] << ", "
             << conterreul[i] << " ,"
             << contliceul[i] << " ,"
             << errmil[i] << ", "
             << licmil[i] << ", "
             << conterrmil[i] << ", "
             << contlicmil[i] << ", "
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
    std::fill(std::begin(conterreul), std::end(conterreul), 0);
    std::fill(std::begin(contliceul), std::end(contliceul), 0);
    std::fill(std::begin(errmil), std::end(errmil), 0);
    std::fill(std::begin(licmil), std::end(licmil), 0);
    std::fill(std::begin(conterrmil), std::end(conterrmil), 0);
    std::fill(std::begin(contlicmil), std::end(contlicmil), 0);
    std::fill(std::begin(Npas), std::end(Npas), 0);

    BS = K*std::exp(-r*T)*normalCDF(-d2)-S0*normalCDF(-d1);

    for (i=0;i<P;++i)
    {
        dt = T/N;
        sigdt = sig*std::sqrt(dt);
        rdt = r*dt;
        der = rdt - std::pow(sigdt,2)/2;

        std::fill(std::begin(S), std::end(S), S0);
        std::fill(std::begin(Se), std::end(Se), S0);
        std::fill(std::begin(Sm), std::end(Sm), S0);

        for (k=0; k<N; k++) {
            std::normal_distribution<double> d(0, 1);
            for(l = 0; l < M; l++) {
                current_sample = d(gen);
                S[l] = S[l]*std::exp(sigdt*current_sample+der); /* filling S */
                Se[l] = Se[l]*(1+sigdt*current_sample+rdt); /* filling Se */
                Sm[l] = Sm[l]*(1+sigdt*current_sample+der+std::pow(sigdt*current_sample,2)/2); /* filling Sm */
            };
        };

        for(l = 0; l < M; l++) {
            paymc[l] = std::max(K-S[l],double(0));
            payeul[l] = std::max(K-Se[l],double(0));
            paymil[l] = std::max(K-Sm[l],double(0));
        };

        puteul = sum(payeul, M);
        careul = sum_pow(payeul, M, 2);

        diff_array(diffpayeulpaymc, payeul, paymc, M);
        contputeul = sum(diffpayeulpaymc, M);
        contcareul = sum_pow(diffpayeulpaymc, M, 2);

        putmil = sum(paymil,M);
        carmil = sum_pow(paymil,M,2);

        erreul[i] = std::exp(-r*T) * puteul/M - BS;
        liceul[i] = 1.96*std::exp(-r*T)*std::sqrt((careul/M-std::pow(puteul/M,2))/M);
        conterreul[i] = std::exp(-r*T)*contputeul/M;
        contliceul[i] = 1.96*std::exp(-r*T)*std::sqrt((contcareul/M -std::pow(contputeul/M,2))/M);

        diff_array(diffpaymilpaymc, paymil, paymc, M);
        contputmil = sum(diffpaymilpaymc,M);
        contcarmil = sum_pow(diffpaymilpaymc,M,2);

        errmil[i] = std::exp(-r*T) * putmil/M - BS;
        licmil[i] = 1.96*std::exp(-r*T)*std::sqrt((carmil/M-std::pow(putmil/M,2))/M);
        conterrmil[i] = std::exp(-r*T)*contputmil/M;
        contlicmil[i] = 1.96*std::exp(-r*T)*std::sqrt((contcarmil/M -std::pow(contputmil/M,2))/M);

        Npas[i] = N;
        N = 2*N;
    };

    create(erreul, liceul, conterreul, contliceul, errmil, licmil, conterrmil, contlicmil, Npas);

    return 0;
};