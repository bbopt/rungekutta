#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>

bool compute_y ( int ell, long double t, long double * y );

bool compute_y_prime ( int ell, long double t, const long double * y, long double * y_prime );

bool compute_alpha_beta ( long double   beta1 , long double   beta5,
			  long double & alpha1, long double & alpha2, long double & alpha3, long double & alpha4,
			  long double & beta2 , long double & beta3 , long double & beta4 , long double & beta6   );

bool compute_error ( int ell, int step, long double t_initial, long double t_final,
		     long double alpha1, long double alpha2, long double alpha3, long double alpha4,
		     long double beta1,  long double beta2,  long double beta3,  long double beta4,
		     long double beta5,  long double beta6,
		     long double & error );
 
/*--------------------------------------*/
int main ( int argc , char ** argv ) {
/*--------------------------------------*/

  if ( argc != 2 ) {
    std::cout << 1e20 << std::endl;
    return 1;
  }

  // read the inputs beta1 and beta5:
  // --------------------------------
  
  long double beta1, beta5;

  std::ifstream in ( argv[1] );

  in >> beta1 >> beta5;

  if ( in.fail() ) {
    in.close();
    std::cout << 1e20 << std::endl;
    return 1;
  }

  in.close();

  // compute the alpha and beta values:
  // ----------------------------------
  long double alpha1, alpha2, alpha3, alpha4,
    beta2, beta3, beta4, beta6;

  if ( !compute_alpha_beta ( beta1, beta5, alpha1, alpha2, alpha3, alpha4, beta2, beta3, beta4, beta6 ) ) {
    std::cout << 1e20 << std::endl;
    return 1;
  }


  // compute all the errors:
  // -----------------------
  long double f = 0.0, error;
  int  cnt = 0;
  
  for ( int ell = 4; ell <= 7; ++ell ) {

    for ( int step = 145; step <= 150; ++step ) {

      error = 0.0;
      
      if ( !compute_error ( ell, step, 1.0, 4.0,
			    alpha1, alpha2, alpha3, alpha4, beta1, beta2, beta3, beta4, beta5, beta6,
			    error ) ) {
	std::cout << 1e20 << std::endl;
	return 1;
      }

      f = f + error*error;
    }
  }
 
  std::cout.precision(15);
  std::cout << log(f)/log(10.0) << std::endl;
  
  return 0;
}

/*-------------*/
/*  compute y  */
/*-------------*/
bool compute_y ( int ell, long double t, long double * y ) {
  for ( int i = 0 ; i < ell ; ++i ) {
    y[i] = (2.0+i)*(1.0-(i+1.0)*t*t)/(1.0+(1.0+i)*t*t);
  } 
  return true;
}

/*---------------------------*/
/*  compute y_prime given y  */
/*---------------------------*/
bool compute_y_prime ( int ell, long double t, const long double * y, long double * y_prime ) {

  long double den1, den2, sqr;
  
  for ( int i = 0 ; i < ell-1 ; ++i ) {

    den1 = t*(2.0+i);
    den2 = 3.0+i-y[i+1];
    sqr  = (2.0+i)*(3.0+i+y[i+1])/den2;  

    if ( den1 == 0.0 || den2 == 0.0 || sqr < 0.0 )
      return false;

    y_prime[i] = pow(y[i],2.0)/den1 - (2.0+i)*sqrt(sqr);    
  }

  den1 = t*(1.0+ell);
  den2 = 2.0-y[0];
  sqr  = (2.0+y[0])/den2;
  
  if ( den1 == 0.0 || den2 == 0.0 || sqr < 0.0 )
    return false;

  y_prime[ell-1] = pow(y[ell-1],2.0)/den1 - (1.0+ell)*sqrt(sqr);
  
  return true;
}

/*-------------------------------------*/
/*  compute the alpha and beta values  */
/*-------------------------------------*/
bool compute_alpha_beta ( long double   beta1 , long double   beta5,
			  long double & alpha1, long double & alpha2, long double & alpha3, long double & alpha4,
			  long double & beta2 , long double & beta3 , long double & beta4 , long double & beta6   ) {
 
  // compute the beta values:
  // ------------------------
  
  beta4 = 1.0;
  
  long double den;
  long double b12 = pow(beta1, 2.0);
  long double b13 = pow(beta1, 3.0);
  long double b22 = 1e20;
  long double b23 = 1e20;
  long double b32 = 1e20;
  
  if ( beta1 == 0.5 ) {
    beta2 = 0.5;
    if ( beta5 == 1.0 ) {
      alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
      return false;
    }
    beta3 = 1.0/(2*(1.0-beta5));
  }
  else {
    den = 8.0*(3.0*b12*beta5-2.0*beta1*beta5-beta1+1.0);
    if ( den == 0.0 ) {
      alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
      return false;
    }

    long double b14 = pow(beta1, 4.0);
    long double b52 = pow(beta5, 2.0);
    long double gamma =
      144.0 * pow(beta1,6.0)*b52 -
      384.0 * pow(beta1,5.0)*b52 +
      400.0 * b14*b52   -
      40.0  * b14*beta5 -
      192.0 * b13*b52   +
      72.0  * b13*beta5 +
      36.0  * b12*b52   +
      16.0  * b13       -
      36.0  * b12*beta5 -
      39.0  * b12       +
      4.0   * beta1*beta5 +
      30.0  * beta1       -
      7.0;

    if ( gamma < 0.0 ) {
      alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
      return false;
    }
    
    beta2 = (12.0*b13*beta5-6.0*beta1*beta5-sqrt(gamma)-5*beta1+5.0) / den;
    den   = 2.0*beta1*(2*beta1-1.0);

    if ( den == 0.0 ) {
      alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
      return false;
    }
    beta3 = (beta2*(beta1-beta2)) / den;
  }

  b32 = pow(beta3,2.0);
  b22 = pow(beta2,2.0);
  b23 = pow(beta2,3.0);
  den = 12*b13*b32-4.0*b12*beta2*beta3-8.0*b12*b32+4.0*beta1*b22*beta3+beta1*b22-b23;
  if ( den == 0.0 ) {
    alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
    return false;
  }

  beta6 = beta1*(-beta1*beta2*beta5+b22*beta5+beta1*beta3-beta3) / den;
  
  // compute the alpha values:
  // -------------------------
  
  den = 24.0*b13*beta3*(-beta1*beta2*beta5+b22*beta5+beta1*beta3-beta3);
    
  if ( den == 0.0 ) {
    alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
    return false;
  }

  alpha2 = (12.0*b12*b22*beta3*beta5-8.0*b12*beta2*beta3*beta5-4.0*b12*b32-4.0*beta1*b22*beta3+4.0*beta1*beta2*beta3+b23-b22) / den;
  den    = 24.0*b12*beta3*(beta1*beta2*beta5-b22*beta5-beta1*beta3+beta3);
  if ( den == 0.0 ) {
    alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
    return false;
  }
  alpha3 =  (12.0*b13*beta3*beta5-8.0*b12*beta3*beta5-4.0*b12*beta3+beta1*beta2+4*beta1*beta3-beta2) / den;
  alpha4 = -(12.0*b13*b32-4.0*b12*beta2*beta3-8.0*b12*b32+4.0*beta1*b22*beta3+beta1*b22-b23) / den;
  alpha1 = 1.0-alpha2-alpha3-alpha4;
 
  // check the alpha and beta values (all the "v" values must be equal to zero):
  // -------------------------------
  long double worst_violation = 0.0;
  {
    long double v = fabs ( alpha1+alpha2+alpha3+alpha4-1.0);
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha2*beta1+alpha3*beta2+alpha4*beta4-0.5 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha2*beta1*beta1+alpha3*beta2*beta2+alpha4*beta4*beta4-1.0/3.0 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha2*pow(beta1,3.0)+alpha3*pow(beta2,3.0)+alpha4*pow(beta4,3.0)-0.25 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha3*beta1*beta3+alpha4*beta1*beta5+alpha4*beta2*beta6-1.0/6.0 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha3*beta1*beta2*beta3+alpha4*beta1*beta4*beta5+alpha4*beta2*beta4*beta6-1.0/8.0 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha3*beta1*beta1*beta3+alpha4*beta1*beta1*beta5+alpha4*beta2*beta2*beta6-1.0/12.0 );
    if ( v > worst_violation ) worst_violation = v;
    v = fabs ( alpha4*beta1*beta3*beta6-1.0/24.0 );
    if ( v > worst_violation ) worst_violation = v;
  }

  if ( worst_violation > 1e-13 ) {
    alpha1 = alpha2 = alpha3 = alpha4 = beta2 = beta3 = beta4 = beta6 = 1e20;
    return false;
  }

  return true;
}

/*----------------------------------------------------------*/
/*  compute the error for one dimension (ell) and one step  */
/*----------------------------------------------------------*/
bool compute_error ( int ell, int step, long double t_initial, long double t_final,
		     long double alpha1, long double alpha2, long double alpha3, long double alpha4,
		     long double beta1,  long double beta2,  long double beta3,  long double beta4,
		     long double beta5,  long double beta6,
		     long double & error ) {

  error = 1e20;
  
  int i, j;
  long double t;
  long double h = ( t_final - t_initial ) / step;
    
  long double * y = new long double [ell];
  
  // initial solution:
  if ( !compute_y ( ell, t_initial, y ) ) {
    error = 1e20;
    delete [] y;
    return false;
  }

  long double * y_tmp = new long double [ell];
  long double * k1    = new long double [ell];
  long double * k2    = new long double [ell];
  long double * k3    = new long double [ell];
  long double * k4    = new long double [ell];

  for ( i = 0 ; i < ell ; ++i )
    k1[i] = k2[i] = k3[i] = k4[i] = y_tmp[i] = 0.0;
  
  // the RK method:
  // --------------
  for ( j = 0; j < step; ++j ) {

    t = t_initial+j*h;

    // compute k1:
    // -----------
    if ( !compute_y_prime ( ell, t, y, k1 ) ) {
      delete [] y;
      delete [] y_tmp;
      delete [] k1;
      delete [] k2;
      delete [] k3;
      delete [] k4;
      error = 1e20;
      return false;
    }
    for ( i = 0 ; i < ell ; ++i )
      k1[i] = h*k1[i];
    
      
    // compute k2:
    // -----------
    for ( i = 0 ; i < ell ; ++i )
      y_tmp[i] = y[i] + beta1*k1[i];
    if ( !compute_y_prime ( ell, t+beta1*h, y_tmp, k2 ) ) {
      delete [] y;
      delete [] y_tmp;
      delete [] k1;
      delete [] k2;
      delete [] k3;
      delete [] k4;
      error = 1e20;
      return false;
    }
    for ( i = 0 ; i < ell ; ++i )
      k2[i] = h*k2[i];  
    
    // compute k3:
    // -----------
    for ( i = 0 ; i < ell ; ++i )
      y_tmp[i] = y[i] + beta3*k2[i] + (beta2-beta3)*k1[i];
    if ( !compute_y_prime ( ell, t+beta2*h, y_tmp, k3 ) ) {
      delete [] y;
      delete [] y_tmp;
      delete [] k1;
      delete [] k2;
      delete [] k3;
      delete [] k4;
      error = 1e20;
      return false;
    }
    for ( i = 0 ; i < ell ; ++i )
      k3[i] = h*k3[i];  
    
    // compute k4:
    // -----------
    for ( i = 0 ; i < ell ; ++i ) {
      y_tmp[i] = y[i] + beta5*k2[i] + beta6*k3[i] + (beta4-beta5-beta6)*k1[i];      
    }

    if ( !compute_y_prime ( ell, t+beta4*h, y_tmp, k4 ) ) {
      delete [] y;
      delete [] y_tmp;
      delete [] k1;
      delete [] k2;
      delete [] k3;
      delete [] k4;
      error = 1e20;
      return false;
    }
    for ( i = 0 ; i < ell ; ++i )
      k4[i] = h*k4[i];  
    
    // update y:
    for ( i = 0 ; i < ell ; ++i )
      y[i] = y[i] + alpha1*k1[i] + alpha2*k2[i] + alpha3*k3[i] + alpha4*k4[i];      
  }

  delete [] k1;
  delete [] k2;
  delete [] k3;
  delete [] k4;

  // true value of y:
  compute_y ( ell, 4.0, y_tmp );
  
  // compute error:
  error = 0.0;
  for ( i = 0 ; i < ell ; ++i )
    error = error + pow ( y[i]-y_tmp[i], 2.0 );

  error = sqrt(error);
  
  delete [] y;
  delete [] y_tmp;
  
  return true;
}
