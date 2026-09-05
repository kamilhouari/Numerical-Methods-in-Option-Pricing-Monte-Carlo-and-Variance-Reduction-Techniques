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
    double sampleVar =
        (sumsqr - sum*sum/N)/(N-1);

    double se =
        exp(-r*T)*sqrt(sampleVar/N);

    return se;
}
double betaOptimal(double S,double K,double r,double sigma,double T,const std::vector<std::vector<double> >& Zs,double M)
{
    double sumX  = 0.0;
    double sumY  = 0.0;
    double sumXY = 0.0;
    double sumY2 = 0.0;
    int N = Zs.size();
    double dt = 1.0/M;
    for(const auto& pathZ : Zs)
    {
        double St = S;
        double arithmeticSum = 0.0;
        double logSum = 0.0;
        for(int k=0;k<M*T;k++)
        {
            St *= exp((r-0.5*sigma*sigma)*dt + sigma*sqrt(dt)*pathZ[k]);
            arithmeticSum += St;
            logSum += log(St);
        }
        double A = arithmeticSum/(M*T);
        double G = exp(logSum/(M*T));
        double X = std::max(A-K,0.0);
        double Y = std::max(G-K,0.0);
        sumX  += X;
        sumY  += Y;
        sumXY += X*Y;
        sumY2 += Y*Y;
    }
    double EX = sumX/N;
    double EY = sumY/N;
    double cov =sumXY/N - EX*EY;
    double var =sumY2/N - EY*EY;
    if (var < 1e-14) return 0.0;
    return cov/var;
}
double geo_pyff(double pS,double K,double pr,double ps,std::vector<double> Zs,double pM) {
    double sum=0;
    double dt=1.0/pM;
    double S=pS;
    for (double Z:Zs) {
        S=St(S,pr,ps,dt,Z);
        sum+=log(S);
    }

    return max(exp(1.0/Zs.size()*sum)-K,0.0);
}
double geoAsianClosedFormPrice(double S0,double K,double r,double sigma,double T ,int n){
    double muG =std::log(S0)+ (r - 0.5*sigma*sigma)* T * (n + 1.0) / (2.0*n);
    double sigmaG2 =sigma*sigma*T* (n + 1.0)* (2.0*n + 1.0)/ (6.0*n*n);
    double sigmaG = std::sqrt(sigmaG2);
    double d1 =(muG - std::log(K) + sigmaG2)/ sigmaG;
    double d2 = d1 - sigmaG;
    return std::exp(-r*T)* (std::exp(muG + 0.5*sigmaG2)* N(d1)- K * N(d2));
}
double payoff_C(double pK, double pS, double pr, double ps, double pT,const std::vector<double>& pZs,double pM) {

        double dt = 1.0 / pM;
        double St_path = pS;
        double average = 0.0;
        for(int k=0;k<pM*pT;k++){
            St_path *= exp((pr - 0.5*ps*ps)*dt+ ps*sqrt(dt)*pZs[k]);
            average += St_path;
        }
        average /= pM*pT;
        return max(average - pK, 0.0);
}
double payoff_C_ANT(double pK, double pS, double pr, double ps, double pT,const std::vector<double>& pZs,double pM) {
    double St1 = pS;
    double St2 = pS;
    double dt = 1.0 / pM;
    double avg1 = 0.0;
    double avg2 = 0.0;

    for(int k=0;k<pM*pT;k++)
    {
        St1 *= exp((pr-0.5*ps*ps)*dt+ ps*sqrt(dt)*pZs[k]);

        St2 *= exp((pr-0.5*ps*ps)*dt- ps*sqrt(dt)*pZs[k]);
        avg1 += St1;
        avg2 += St2;
    }
    avg1 /= pM*pT;
    avg2 /= pM*pT;
    return 0.5 * (max(avg1-pK,0.0) +max(avg2-pK,0.0));
}
double payoff_C_CV(double K,double S,double r,double sigma,double T,const std::vector<double>& pathZ,double M,double beta,double geoPrice){
    double dt = 1.0 / M;
    double S1 = S;
    double S2 = S;
    double arithmeticSum = 0.0;
    double logSum = 0.0;
    double arithmeticSumAnti = 0.0;
    double logSumAnti = 0.0;
    for(int k=0;k<M*T;k++)
    {
        S1 *= exp((r-0.5*sigma*sigma)*dt + sigma*sqrt(dt)*pathZ[k]);
        S2 *= exp((r-0.5*sigma*sigma)*dt - sigma*sqrt(dt)*pathZ[k]);
        arithmeticSum += S1;
        arithmeticSumAnti += S2;
        logSum += log(S1);
        logSumAnti += log(S2);
    }
    double A1 = arithmeticSum / (M*T);
    double A2 = arithmeticSumAnti / (M*T);
    double G1 = exp(logSum / (M*T));
    double G2 = exp(logSumAnti / (M*T));
    double X =0.5*(std::max(A1-K,0.0) + std::max(A2-K,0.0));
    double Y =0.5*(std::max(G1-K,0.0) + std::max(G2-K,0.0));
    return X - beta*(Y-exp(r*T)*geoPrice);
}

double CallPriceCV(double S,double K,double r,double sigma,double T,const std::vector<std::vector<double> >& Zs,double M,double beta,double geoPrice){
    double sumPayoff = 0.0;
    for(const auto& pathZ : Zs){
        sumPayoff += payoff_C_CV(K,S,r,sigma,T,pathZ,M,beta,geoPrice);
    }
    return exp(-r*T)*sumPayoff/Zs.size();
}
double CallPriceSt(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs,int pN,double M) {
    double payoff=0;
    double sumPayoff = 0;

    for (std::vector<double> Z : pmZs) {
        payoff =payoff_C( (pK),  (pS),  (pr),  (ps),  (pT), (Z),(M));
        sumPayoff += payoff;
    }
    return std::exp(-pr*pT)*(sumPayoff/pN);
}
double CallPriceAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs,int pN, double pM ) {
    double payoff=0;
    double sumPayoff = 0;

    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C_ANT( (pK),  (pS),  (pr),  (ps),  (pT), (Zs),(pM));
        sumPayoff += payoff;
    }
    return std::exp(-pr*pT)*(sumPayoff/pN);
}
double pyfSumsSqr(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs, double M) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C((pK), (pS), (pr),  (ps),  (pT), (Zs),(M) );
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSumsSqrAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs, double M) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C_ANT((pK), (pS), (pr),  (ps),  (pT), (Zs),(M));
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSumsSqrCv(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs,double pM ,double bta,double geoPrice) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C_CV((pK), (pS), (pr),  (ps),  (pT), (Zs),pM,(bta),geoPrice);
        sumPayoff += payoff*payoff;
    }
    return sumPayoff;
}
double pyfSum(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs, double M) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C((pK), (pS), (pr),  (ps),  (pT), (Zs),(M));
        sumPayoff += payoff;
    }
    return sumPayoff;
}
double pyfSumAnt(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs, double M) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        sumPayoff += payoff_C_ANT((pK), (pS), (pr),  (ps),  (pT), (Zs),(M));
    }
    return sumPayoff;
}
double pyfSumCv(double pS,double pK,double pr,double ps,double pT,const std::vector<std::vector<double> >& pmZs,double pM ,double bta,double geoPrice) {
    double payoff=0;
    double sumPayoff = 0;
    for (std::vector<double> Zs : pmZs) {
        payoff =payoff_C_CV((pK), (pS), (pr),  (ps),  (pT), (Zs),pM,(bta),geoPrice);
        sumPayoff += payoff;
    }
    return sumPayoff;
}

double dlta(double pS,double pK,double pr,double ps,double pT,std::vector<std::vector<double> >&pmZs,int pM,double bta,double steps) {
    return (CallPriceCV(pS+pS*0.01,pK,pr,ps,pT,pmZs,pM,betaOptimal( pS+pS*0.01, pK, pr, ps, pT,  pmZs, pM),geoAsianClosedFormPrice(pS+pS*0.01,pK,pr,ps,pT,steps))-CallPriceCV(pS,pK,pr,ps,pT,pmZs,pM,bta,geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps)))/(pS*0.01);
}
double gmm(double pS,double pK,double pr,double ps,double pT,std::vector<std::vector<double> >&pmZs,int pM,double bta,double steps) {
    return (CallPriceCV(pS+pS*0.01,pK,pr,ps,pT,pmZs,pM,betaOptimal( pS+pS*0.01, pK, pr, ps, pT,  pmZs,pM), geoAsianClosedFormPrice(pS+pS*0.01,pK,pr,ps,pT,steps))-2*CallPriceCV(pS,pK,pr,ps,pT,pmZs,pM,bta,geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps))+CallPriceCV(pS-pS*0.01,pK,pr,ps,pT,pmZs,pM,betaOptimal( pS-pS*0.01, pK, pr, ps, pT,  pmZs, pM),geoAsianClosedFormPrice(pS-pS*0.01,pK,pr,ps,pT,steps)))/(pS*pS*0.01*0.01);
}
double tht(double pS,double pK,double pr,double ps,double pT,std::vector<std::vector<double> >&pmZs,int pM,double bta,double steps) {
return ((CallPriceCV(pS,pK,pr,ps,pT,pmZs,pM,bta,geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps))-CallPriceCV(pS,pK,pr,ps,pT-pT/252,pmZs,pM,betaOptimal( pS, pK, pr, ps, pT-pT/252,  pmZs, pM),geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps)))/(pT/252))/252  ;
}
double vga(double pS,double pK,double pr,double ps,double pT,std::vector<std::vector<double> >&pmZs,int pM,double bta,double steps) {
    return (CallPriceCV(pS,pK,pr,ps*0.01+ps,pT,pmZs,pM,betaOptimal( pS, pK, pr, ps*0.01+ps, pT,  pmZs, pM), geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps))-CallPriceCV(pS,pK,pr,ps,pT,pmZs,pM,bta,geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps)))/(ps*0.01);
}
double rho (double pS,double pK,double pr,double ps,double pT,std::vector<std::vector<double> >&pmZs,int pM,double bta,double steps) {
    return (CallPriceCV(pS,pK,pr*0.0001+pr,ps,pT,pmZs,pM,betaOptimal( pS, pK, pr*0.01+pr, ps, pT,  pmZs, pM),geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps))-CallPriceCV(pS,pK,pr,ps,pT,pmZs,pM,bta,geoAsianClosedFormPrice(pS,pK,pr,ps,pT,steps)))/(pr*0.0001);
}
    int main() {
    std::random_device rd;
    std::mt19937 gen(805);
    std::normal_distribution<double> dist(0.0, 1.0);
    std::ofstream file1("Asian_mc_results.csv");
    file1 << "N,Call price,serr,err\n";
    std::ofstream file2("Asian_mc_results_antithetic.csv");
    file2 << "N,Call price,serr,err\n";
    std::ofstream file3("Asian_mc_results_antithetic_CV.csv");
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
    double M=12;
    double beta;
    double geoPrice;
    double steps=M*T;
    for (int N: Ns ) {

        std::vector<std::vector<double> > Zs;


        for (int i = 0; i < N; i++) {
            std::vector<double> z;
            for (int m = 0; m < M*T; m++) {
                z.push_back(dist(gen));
            }
            Zs.push_back(z);
            if (i%(N/10)==0) {
                std::cout <<i<< "\n";
            }
        }
        beta= betaOptimal(S,K,r,s,T,Zs,M);
        geoPrice =  geoAsianClosedFormPrice(S,K,r,s,T,steps);
        standard   = CallPriceSt( S, K, r, s, T, Zs, N , M);
        antithetic = CallPriceAnt( S, K, r, s, T, Zs, N , M);
        CV         = CallPriceCV( S, K, r, s, T, Zs, M,beta,geoPrice);

        std:: cout <<"N="<< N <<"\n"<<"\n";
        std::cout <<"Standard MC:                  " <<  standard << " +/- "<<se(pyfSumsSqr(S,K,r,s,T,Zs,M ),pyfSum(S,K,r,s,T, Zs,M ),N,r,T)<<std::endl;
        std::cout <<"Antithetic:                   " <<  antithetic << " +/- "<<se(pyfSumsSqrAnt(S,K,r,s,T,Zs,M),pyfSumAnt(S,K,r,s,T,Zs,M),N,r,T)<<std::endl;
        std::cout <<"Control Variate + Antithetic: " <<  CV << " +/- "<<se(pyfSumsSqrCv(S,K,r,s,T,Zs,M,beta,geoPrice),pyfSumCv(S,K,r,s,T,Zs,M,beta,geoPrice),N,r,T)<<std::endl;
        std::cout << "beta = " << beta << std::endl;
        double geoMC = 0.0;

        for(const auto& path : Zs)
        {
            geoMC += geo_pyff(S,K,r,s,path,M);
        }

        geoMC = exp(-r*T)*geoMC/Zs.size();

        std::cout << "Geo MC = " << geoMC << std::endl;
        std::cout << "Geo CF = " << geoPrice << std::endl;
        file1 << N << "," <<  standard <<","<<se(pyfSumsSqr(S,K,r,s,T,Zs,M),pyfSum(S,K,r,s,T,Zs,M),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-standard,standard-bsC(K,S,r,s,T))<< "\n" ;
        file2 << N << "," <<  antithetic <<","<<se(pyfSumsSqrAnt(S,K,r,s,T,Zs,M),pyfSumAnt(S,K,r,s,T,Zs,M),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-antithetic,antithetic-bsC(K,S,r,s,T))<< "\n" ;
        file3 << N << "," <<  CV <<","<<se(pyfSumsSqrCv(S,K,r,s,T,Zs,M,beta,geoPrice),pyfSumCv(S,K,r,s,T,Zs,M,beta,geoPrice),N,r,T) <<","<<fmax(bsC(K,S,r,s,T)-CV,CV-bsC(K,S,r,s,T))<< "\n"<<"\n" ;
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
                 << setw(15) << dlta(S,K,r,s,T,Zs,M,beta,steps)
                 << setw(15) << gmm(S,K,r,s,T,Zs,M,beta,steps)
                 << setw(15) << tht(S,K,r,s,T,Zs,M,beta,steps)
                 << setw(15) << vga(S,K,r,s,T,Zs,M,beta,steps)
                 << setw(15) << rho(S,K,r,s,T,Zs,M,beta,steps)
                 << endl;
            std::cout << "---------------------------------------------------------------------------\n";
    }

    file1.close();
    file2.close();
    file3.close();
    return 0;
}