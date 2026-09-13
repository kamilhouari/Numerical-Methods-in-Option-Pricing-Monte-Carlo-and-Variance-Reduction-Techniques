#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;
double N(double z) {
    return 0.5* (1+std::erf(z/sqrt(2.0)));
}
double St(double aS, double ar, double as, double aT, double aZ) {
    return aS*exp((ar-0.5*as*as)*aT+as*sqrt(aT)*aZ);
}
double bsC(double aK, double aS, double ar, double as, double aT) {
    double d1=(log(aS/aK)+(ar+0.5*as*as)*aT)/(as*sqrt(aT));
    double d2=d1-as*sqrt(aT);
    return aS*N(d1)-aK*exp(-1*ar*aT)*N(d2);

}

double se(double sumsqr, double sum,int N, double r, double T) {
    double mean =sum/N;

    double variance =
        (sumsqr/N -mean*mean)/N;

    variance=std::max(variance, 0.0);

    return std::exp(-r*T) * std::sqrt(variance);
}
double betaOptimal(
    double S, double K, double r, double sigma, double T,
    const std::vector<double>& Zs)
{
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXY = 0.0;
    double sumY2 = 0.0;

    int N = Zs.size();

    for (double Z : Zs)
    {
        double ST1 = St(S, r, sigma, T, Z);
        double ST2 = St(S, r, sigma, T, -Z);

        double X = 0.5 * (std::max(ST1 - K, 0.0) +
                          std::max(ST2 - K, 0.0));

        double Y = 0.5 * (ST1 + ST2);

        sumX  += X;
        sumY  += Y;
        sumXY += X * Y;
        sumY2 += Y * Y;
    }

    double EX = sumX / N;
    double EY = sumY / N;

    double cov = (sumXY / N) - EX * EY;
    double var = (sumY2 / N) - EY * EY;

    if (var < 1e-14) return 0.0;

    return cov / var;
}
double payoff_C(double pK, double pS, double pr, double ps, double pT,double pZ) {
    double ST = St(pS,pr,ps,pT,pZ);
    return std::max(ST - pK, 0.0);

}
double payoff_C_ANT(double pK, double pS, double pr, double ps, double pT,double pZ) {
    double ST1 = St(pS,pr,ps,pT,pZ);
    double ST2 = St(pS,pr,ps,pT,-pZ);
    double payoff1 = std::max(ST1 - pK, 0.0);
    double payoff2 = std::max(ST2 - pK, 0.0);
    return 0.5 * (payoff1 + payoff2);
}
double payoff_C_CV(double pK, double pS, double pr, double ps, double pT,double pZ,double bta) {
    double ST1 = St(pS,pr,ps,pT,pZ);
    double ST2 = St(pS,pr,ps,pT,-pZ);
    double payoff1 = std::max(ST1 - pK, 0.0);
    double payoff2 = std::max(ST2 - pK, 0.0);
    double payoff = 0.5*(payoff1+payoff2);
    double STavg = 0.5*(ST1+ST2);
    double expectedST = pS*exp(pr*pT);
    return payoff - bta*(STavg-expectedST);
}


double CallPriceCv(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs,int pN,double bta) {
    double payoff=0;
    double sumPayoff = 0;

    for (double Z : Zs) {
        payoff =payoff_C_CV( (pK),  (pS),  (pr),  (ps),  (pT), (Z),bta);
        sumPayoff += payoff;
    }
    return std::exp(-pr*pT)*(sumPayoff/pN);
}
double CallPriceSt(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs,int pN) {
    double payoff=0;
    double sumPayoff = 0;

    for (double Z : Zs) {
        payoff =payoff_C( (pK),  (pS),  (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff;
    }
    return std::exp(-pr*pT)*(sumPayoff/pN);
}
double CallPriceAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs,int pN) {
    double payoff=0;
    double sumPayoff = 0;

    for (double Z : Zs) {
        payoff =payoff_C_ANT( (pK),  (pS),  (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff;
    }
    return std::exp(-pr*pT)*(sumPayoff/pN);
}
double pyfSumsSqr(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : Zs) {
        payoff =payoff_C((pK), (pS), (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSumsSqrAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : Zs) {
        payoff =payoff_C_ANT((pK), (pS), (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSumsSqrCv(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs,double bta) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : Zs) {
        payoff =payoff_C_CV((pK), (pS), (pr),  (ps),  (pT), (Z),(bta));
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSum(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& pZs) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : pZs) {
        payoff =payoff_C((pK), (pS), (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff;
    }
    return sumPayoff;
}
double pyfSumAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : Zs) {
        payoff =payoff_C_ANT((pK), (pS), (pr),  (ps),  (pT), (Z));
        sumPayoff += payoff;
    }
    return sumPayoff;
}
double pyfSumCv(double pS,double pK,double pr,double ps,double pT,const std::vector<double>& Zs,double bta) {
    double payoff=0;
    double sumPayoff = 0;
    for (double Z : Zs) {
        payoff =payoff_C_CV((pK), (pS), (pr),  (ps),  (pT), (Z),(bta));
        sumPayoff += payoff;
    }
    return sumPayoff;
}

double dlta(double pS,double pK,double pr,double ps,double pT,std::vector<double>&Zs,int pN,double bta) {
    return (CallPriceCv(pS+pS*0.01,pK,pr,ps,pT,Zs,pN,bta)-CallPriceCv(pS,pK,pr,ps,pT,Zs,pN,bta))/(pS*0.01);
}
double gmm(double pS,double pK,double pr,double ps,double pT,std::vector<double>&Zs,int pN,double bta) {
    return (CallPriceCv(pS+pS*0.01,pK,pr,ps,pT,Zs,pN,bta)-2*CallPriceCv(pS,pK,pr,ps,pT,Zs,pN,bta)+CallPriceCv(pS-pS*0.01,pK,pr,ps,pT,Zs,pN,bta))/(pS*pS*0.01*0.01);
}
double tht(double pS,double pK,double pr,double ps,double pT,std::vector<double>&Zs,int pN,double bta) {
return ((CallPriceCv(pS,pK,pr,ps,pT-pT/252,Zs,pN,bta)-CallPriceCv(pS,pK,pr,ps,pT,Zs,pN,bta))/(pT/252))/252  ;
}
double vga(double pS,double pK,double pr,double ps,double pT,std::vector<double>&Zs,int pN,double bta){
    return (CallPriceCv(pS,pK,pr,ps*0.01+ps,pT,Zs,pN,bta)-CallPriceCv(pS,pK,pr,ps,pT,Zs,pN,bta))/(ps*0.01);
}
double rho (double pS,double pK,double pr,double ps,double pT,std::vector<double>&Zs,int pN,double bta) {
    return (CallPriceCv(pS,pK,pr*0.0001+pr,ps,pT,Zs,pN,bta)-CallPriceCv(pS,pK,pr,ps,pT,Zs,pN,bta))/(pr*0.0001);
}
    int main() {
    std::random_device rd;
    std::mt19937 gen(805);
    std::normal_distribution<double> dist(0.0, 1.0);
    std::ofstream file1("mc_results.csv");
    file1 << "N,Call price,serr,err\n";
    std::ofstream file2("mc_results_antithetic.csv");
    file2 << "N,Call price,serr,err\n";
    std::ofstream file3("mc_results_antithetic_CV.csv");
    file3 << "N,Call price,serr,err\n";
    int Ns[7]= {1000,5000,10000,50000,100000,500000,1000000} ;
    double S =100;
    double r= 0.02;
    double K=100;
    double s=0.2;
    double T=1;
    double standard;
    double antithetic;
    double CV;
    for (int N: Ns ) {
        
        std::vector<double> Zs;


        for (int i = 0; i < N; i++) {
            Zs.push_back(dist(gen));
            if (i%(N/10)==0) {
                std::cout <<i<< "\n";
            }
        }
        double beta= betaOptimal(S,K,r,s,T,Zs);
        standard   = CallPriceSt( S, K, r, s, T, Zs, N);
        antithetic = CallPriceAnt( S, K, r, s, T, Zs, N);
        CV         = CallPriceCv( S, K, r, s, T, Zs, N,beta);

        std:: cout <<"N="<< N <<"\n"<<"\n";
        std::cout <<"Standard MC:                  " <<  standard << " +/- "<<se(pyfSumsSqr(S,K,r,s,T,Zs ),pyfSum(S,K,r,s,T, Zs ),N,r,T)<<std::endl;
        std::cout <<"Antithetic:                   " <<  antithetic << " +/- "<<se(pyfSumsSqrAnt(S,K,r,s,T,Zs),pyfSumAnt(S,K,r,s,T,Zs),N,r,T)<<std::endl;
        std::cout <<"Control Variate + Antithetic: " <<  CV << " +/- "<<se(pyfSumsSqrCv(S,K,r,s,T,Zs,beta),pyfSumCv(S,K,r,s,T,Zs,beta),N,r,T)<<std::endl;
        std::cout <<"Black Scholes: " << bsC(K,S,r,s,T) << std::endl;
        std::cout << "beta = " << beta << std::endl;
        file1 << N << "," <<  standard <<","<<se(pyfSumsSqr(S,K,r,s,T,Zs),pyfSum(S,K,r,s,T,Zs),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-standard,standard-bsC(K,S,r,s,T))<< "\n" ;
        file2 << N << "," <<  antithetic <<","<<se(pyfSumsSqrAnt(S,K,r,s,T,Zs),pyfSumAnt(S,K,r,s,T,Zs),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-antithetic,antithetic-bsC(K,S,r,s,T))<< "\n" ;
        file3 << N << "," <<  CV <<","<<se(pyfSumsSqrCv(S,K,r,s,T,Zs,beta),pyfSumCv(S,K,r,s,T,Zs,beta),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-CV,CV-bsC(K,S,r,s,T))<< "\n"<<"\n" ;
            // Header
            cout << left
                 << setw(15) << "Delta"
                 << setw(15) << "Gamma"
                 << setw(15) << "Theta"
                 << setw(15) << "Vega"
                 << setw(15) << "Rho"
                 << endl;

            cout << string(75, '-') << endl;

            // Single row of data
            cout << left
                 << setw(15) << dlta(S,K,r,s,T,Zs,N,beta)
                 << setw(15) << gmm(S,K,r,s,T,Zs,N,beta)
                 << setw(15) << tht(S,K,r,s,T,Zs,N,beta)
                 << setw(15) << vga(S,K,r,s,T,Zs,N,beta)
                 << setw(15) << rho(S,K,r,s,T,Zs,N,beta)
                 << endl;
            std::cout << "---------------------------------------------------------------------------\n";
    }

    file1.close();
    file2.close();
    file3.close();
    return 0;
}
