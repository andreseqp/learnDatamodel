/*=============================================================================================================
                                                   random.h
===============================================================================================================
 
 Routines for the generation of pseudo-random numbers
 
 Written originally by:
    G. Sander van Doorn
    Centre for Ecological and Evolutionary Studies - Theoretical Biology Group
    University of Groningen
    the Netherlands
 
 Program version
        xx/xx/xxxx	:
 
 =============================================================================================================*/

#ifndef random_h
#define random_h
#include <random>
#include <vector>
#include "utils.h"

namespace rnd
{
    class discrete_distribution {
    public:
        discrete_distribution(const int&);
        double& operator[](const int &i) { is_accumulated = false; return pdf[i]; }
        double operator()(const int &i) const { return pdf[i]; }
        int sample();//int sample(std::mt19937 rng);
        std::vector<int> sample(int);
        int size() const { return n;}
    private:
        void accumulate();
        const int n;
        bool is_accumulated;
        std::vector<double> pdf, cdf;
    };

	//extern std::mt19937 rng_d;    
    void set_seed();
	void set_seed(const unsigned int&);
    int integer(const int&);
	bool bernoulli(const double& = 0.5);
	int binomial(const int&, const double& = 0.5);
	double pdf_binomial(const int& k, const int& n, const double& p= 0.5);
	int poisson(const double& = 1.0);
	double uniform(const double  min=0, const double  max=1);
    double normal(const double& = 0.0, const double& = 1.0);
    double exponential(const double& = 1.0);
 
	extern std::mt19937 rng;
	double gamma(const double& shape, const double& scale);
	double beta(const double&, const double&);
	double pdf_gamma(const double&, const double&, const double&);
	//double pdf_beta(const double& x, const double& shapeL, const double& shapeR);
	double r8_beta(double x, double y);
	double beta_pdf(double x, double a, double b);
}

#endif //#ifndef random_h
