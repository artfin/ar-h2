#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <gsl/gsl_sf_legendre.h>

// precomputed values of factorials
// add more if needed
std::vector<double> factorials = {1, 1, 2, 6, 24, 120, 720, 5040, 40320}; 

double tang_toennies(const int n, const double beta, const double R)
// The Tang-Toennies damping function
//
// Parameters:
// n: the order of the damping function
// r: the damping target
// beta: the damping coefficient
{
    double betaR = beta * R;
    double exp_ = std::exp(-betaR);

    double sum = 0.0;
    for ( int m = 0; m <= n; ++m ) { 
        sum += std::pow(betaR, m) / factorials[m]; 
    }

    return 1.0 - exp_ * sum;
}

double tang_toennies_derivative(const int n, const double beta, const double R)
// The derivative of Tang-Toennies dampling function with respect to R
//
// Parameters: the same as for `tang_toennies` function  
{
    double betaR = beta * R;
    double exp_ = std::exp(-betaR);

    return beta * exp_ * std::pow(betaR, n) / factorials[n]; 
}

double tang_toennies_star(const int n, const double beta, const double R)
// The ratio of the derivative of Tang-Toennies function w.r.t. to R to its value  
{
    double derivative_value = tang_toennies_derivative(n, beta, R);
    double function_value = tang_toennies(n, beta, R);
    return derivative_value / function_value;
}

double compute_A(const double epsilon, const double Re, const double C6, const double beta)
// LeRoy 1987, formula (5)
{
    double Re6 = std::pow(Re, 6);
    double Re8 = Re6 * Re * Re;

    double betaRe = beta * Re;

    double D6 = tang_toennies(6, beta, Re);
    double D8star = tang_toennies_star(8, beta, Re);
    double D6star = tang_toennies_star(6, beta, Re);

    return (epsilon * (D8star * Re - 8.0) + (D6star * Re - D8star * Re + 2) * C6 * D6 / Re6)  * std::exp(betaRe) /
        // ------------------------------------------------------------------------------------ 
                                (8.0 - betaRe - D8star * Re);
}   

double compute_C8(const double epsilon, const double Re, const double C6, const double beta, const double A)
// LeRoy 1987, formula (6)
{
    double Re6 = std::pow(Re, 6);
    double Re8 = Re6 * Re * Re;
    
    double D6 = tang_toennies(6, beta, Re);
    double D8 = tang_toennies(8, beta, Re);

    return (A * std::exp(-beta * Re) + epsilon - C6 * D6 / Re6) * Re8 / D8;
}

double Vklambda(const double R, const double beta, const double epsilon, const double C6, const double C8, const double A)
{
    double R6 = std::pow(R, 6);
    double R8 = R6 * R * R;
    double D6 = tang_toennies(6, beta, R);
    double D8 = tang_toennies(8, beta, R);

    return A * std::exp(-beta * R) - D6 * C6 / R6 - D8 * C8 / R8; 
}

double compute_V2(const double R)
{
    double beta = 3.5815; // A^(-1)
    return Vklambda(R, beta, 6.0812, 13500.0, 319951.0, 4212302.0);
    //double V21 = Vklambda(R, beta, 11.71, 29600.0, 948287.0, 15411411.0);
    //double V22 = Vklambda(R, beta, 4.035136, 5705.0, 628336.0, 11199109.0);
    //double V23 = Vklambda(R, beta, 0.0, -10395.0, 0.0, 0.0);
}

double compute_V0(const double R)
{
    double beta = 3.5815; // A^(-1)
    return Vklambda(R, beta, 51.069, 134500.0, 1202347.0, 18778229.0);
}

double pes(double R, double Theta)
{
    double V0 = compute_V0(R);
    double V2 = compute_V2(R);

    double cosTheta = std::cos(Theta);
    return V0 + V2 * gsl_sf_legendre_P2(cosTheta);
}

int main()
{
    std::cout << std::fixed << std::setprecision(16);

    int npoints = 100;
    double Rmin = 3.0;
    double Rmax = 5.0;
    double dr = (Rmax - Rmin) / npoints;

    double R = Rmin;
    for (int k = 0; k < npoints; ++k) {
        R += dr; 
    }

    return 0; 
}
