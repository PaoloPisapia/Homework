/////////////////////////////////////////////////////////////////////////
// Compile the code with:                                              //
//          g++ elliptical.cpp root_solvers.cpp quadratureRules.cpp    //
//                                                                     //
// If you wish to print the steps of the root_solvers,                 //
// modify ROOT_PRINT to 1 in root_solvers.                             //
/////////////////////////////////////////////////////////////////////////


#include "/mnt/c/Users/paolo/OneDrive/Documenti/Università/Magistrale/Algoritmi/my_header.h"

static double g_time = 0.;      
static double g_k = 0.8;        // Global variable to fix the value of k 
                                // a priori

double func(double);

double G(double);


int main()
{  
  // For more precision when printing
  cout << setiosflags(ios::scientific) << setprecision(12); 
  
  double T;               // The period (to be determined)
  double pi = M_PI;
  double phi = pi * 0.5;  // The Jacobi amplitude (given)
  
  double a, b;            // The intervals endpoints
  int N;                  // Number number of interval subdivisions
  int Ng;                 // Number of gaussian points in each interval
  
  ///////
  // 1 //
  ///////
  /* ------------Integrating using the Gauss-Legendre Rule------------*/
  
  a = 0.;
  b = phi;      // Setting the integral endpoints
  
  N = 20;
  Ng = 4;       // Setting the number of subdivisions and gaussian
                // points in each of them
  
  T = 4 * GaussQuadratureRule (func, a, b, N, Ng); // Computing the 
                                                   // period value
  
  cout << " The oscillation period T is " << T << endl;

  ///////
  // 2 //
  ///////
  /* -----Relative error compared to the small-angle approximation-----*/
  
  double T0 = 2*pi;     // Small-angle approximation
  double abs_err = fabs(T-T0);
  double rel_err = abs_err / fabs(T0);
  
  cout << "\n The oscillation period T0 is " << T0 << endl
       << " The absolute error is " << abs_err << endl
       << " The relative error is " << rel_err << endl;
       
  ///////
  // 3 //
  ///////
  /* ------------Plotting phi(t) by inverting t=F(phi, k)------------*/
  
  ofstream fdata;
    
  double tmax = 30.0;       // Using a variable to set the maximum time 
                            // at which computing phi(t)
  double h = 0.1;           // The spacing between times

  double zeroN = 0.;        // The root of the function G found with 
                            // Newton
  double zeroB = 0.;        // The root of the function G found with 
                            // Bisection

  double toll = 1.e-7;      // Tollerance for the root finders

  int checkN;               // Variable to check whether Newton worked 
                            // properly
  int checkB;               // Variable to check whether Bisection 
                            // worked properly  

  a = 0.;
  b = 21*phi;               // Setting the integral endpoints, making
                            // sure that G(a)*G(b)<=0 for every g_time
  fdata.open("phi.dat");
  
  fdata << " # 1st column: t; 2nd column: phi (using Newton as root" 
        << " finder); 3rd column: phi (using Bisection as root finder)" 
        << endl;

  while( g_time < tmax+h)   // For unkown reasons, g_time <= tmax does 
                            // not work...
  {
    checkN = Newton(G, func, a, b, toll, zeroN);
    checkB = Bisection(G, a, b, toll, zeroB);
    
    
    cout << "\n Newton check: " << checkN << endl;

    cout << " Bisection check: " << checkB << endl;

    /*cout  << G(a)   << "\t" 
            << G(b)   << "\t" 
            << g_time << "\t" 
            << zeroN  << "\t" 
            << zeroB << endl;*/ // Used simply to check that everything 
                                // run smoothly
    
    fdata << g_time << " " << zeroN << " " << zeroB << endl;
  
    g_time += h;
  }
  
  fdata.close();

  return 0;  
}

/////////////////////////////////////////////////////////////////////////
// Defining the desired functions:                                     //
//  - func: integrand of the incomplete elliptic integral of the first //
//          kind                                                       //
//  - G: function whose roots are the desired phi(t): G:=F(phi,k)-t    //
//       (required to be a function to fit into Bisection and Newton)  //
/////////////////////////////////////////////////////////////////////////

double func(double u)
{
  return 1. / sqrt(1 - g_k*g_k * sin(u)*sin(u) );
}

double G(double phi)
{
  double F = GaussQuadratureRule(func, 0., phi, 20, 4);

  return F - g_time;
}