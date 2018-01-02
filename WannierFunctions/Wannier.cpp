#include <armadillo>
#include <complex>
#include <cmath>

#include "gnuplot-iostream.h"

using namespace arma;

mat Hamiltonian(int L, double V0, double q){
  vec lgrid   = linspace<vec>(-L+1,L-1,2*L-1);
  mat H       = zeros<mat>(2*L-1,2*L-1);
  H.diag()    = (q*q+V0/2.0)*ones<vec>(2*L-1) + 4.0*q*lgrid + 4.0*pow(lgrid,2);
  H.diag(-1)  = -V0/4.0*ones<vec>(2*L-2);
  H.diag(+1)  = -V0/4.0*ones<vec>(2*L-2);
  return H;
}

cx_mat calcBlochFuncs(vec& qvals, vec& xvals, double V0, int L){
  std::complex<double> cplx_i(0,1);
  mat H, eigvec, coeff(2*L-1,qvals.n_rows);
  vec eigval, lgrid = linspace<vec>(-L+1,L-1,2*L-1);
  cx_mat BlochFunc(xvals.n_rows,qvals.n_rows);

  for (size_t i = 0; i < qvals.n_rows; i++) {
    H = Hamiltonian(L,V0,qvals(i));
    eig_sym(eigval,eigvec,H);
    if (eigvec(L-1,0) > 0) {
      coeff.col(i) = eigvec.col(0);
    }
    else{
      coeff.col(i) = -eigvec.col(0);
    }

    for (size_t j = 0; j < xvals.n_rows; j++) {
      BlochFunc(j,i) = dot(exp(cplx_i*(qvals(i)+2.0*lgrid)*2.0*datum::pi*xvals(j)),coeff.col(i));
      // k_lat = 2*pi/a_lat or pi/a_lat  ??  need to compare with table
    }
  }

  return BlochFunc;
}

cx_vec calcWannierFunc(cx_mat& BlochFunc, vec& xgrid){
  cx_vec w   = sum(BlochFunc,1);

  w         -= (0.5*BlochFunc.col(0) + 0.5*BlochFunc.col(BlochFunc.n_cols-1));
  w         /= BlochFunc.n_cols;
  w         /= as_scalar(sqrt(trapz(xgrid,pow(abs(w),2))));
  return w;
}

vec getWannierAtSite_j(cx_vec& w0, vec& xgrid, int j, double a_lat){
  vec wj;
  vec site = xgrid-0.5*j*a_lat;
  vec realw = real(w0);
  interp1(xgrid,realw,site,wj,"linear",0);
  return wj;
}

vec calcPotential(vec& xgrid, double V0, double a_lat){
  return V0*pow(sin(2.0*datum::pi/a_lat*xgrid),2);
}

double calcU(vec& xgrid, cx_mat& w){
  double g   = 0.00502320473*2.0/datum::pi; //a_scat in units of a_lat
  mat normsq = real(conj(w)%w);

  normsq.each_col( [&](vec& a){
    g *= as_scalar(trapz(xgrid,pow( a ,2)));
  });
  return g;
}

double calcJ(vec& xgrid, cx_vec& w0, double V0, double a_lat){
  double dx = xgrid(1)-xgrid(0);
  cx_vec del(w0.n_rows);
  for (size_t i = 1; i < w0.n_rows-1; i++) {
    del(i) = (w0(i+1) + w0(i-1) - 2.0*w0(i))/dx/dx;
  }
  del(0)            = (-5.0*w0(1) + 4.0*w0(2) - w0(3) + 2.0*w0(0))/dx/dx;
  del(w0.n_rows-1)  = (-5.0*w0(w0.n_rows-2) + 4.0*w0(w0.n_rows-3)
                      - w0(w0.n_rows-4) + 2.0*w0(w0.n_rows-1))/dx/dx;

  cx_vec kin        = -a_lat*a_lat/4.0/datum::pi/datum::pi*del;
  cx_vec pot        = calcPotential(xgrid,V0,a_lat)%w0;
  vec w1            = getWannierAtSite_j(w0,xgrid,-1,a_lat);
  cx_vec integrand  = -(kin+pot)%w1;

  return as_scalar(trapz(xgrid,real(integrand) ));
}

mat calculateParameters(int L, double a_lat, double V0min, double V0max, size_t vals){
  double V0_T   = 20;
  auto qgrid    = linspace<vec>(-1,1,101);
  auto xgrid    = linspace<vec>(-L*a_lat,L*a_lat,1000);
  auto lgrid    = linspace<vec>(-L+1,L-1,2*L-1);
  auto Vgrid    = linspace<vec>(V0min,V0max,vals);
  mat data      = zeros<mat>(vals,3);

  for (size_t i = 0; i < vals; i++) {
    auto Bloch    = calcBlochFuncs(qgrid,xgrid,Vgrid(i),L);
    auto Wannier  = calcWannierFunc(Bloch,xgrid);
    auto Bloch_T  = calcBlochFuncs(qgrid,xgrid,V0_T,L);
    auto Wannier_T= calcWannierFunc(Bloch_T,xgrid);

    cx_mat W_tmp  = join_horiz(Wannier,Wannier_T);
    cx_mat W_full = join_horiz(W_tmp,Wannier_T);

    data(i,0)     = Vgrid(i);
    data(i,1)     = calcU(xgrid,W_full);
    data(i,2)     = calcJ(xgrid,Wannier,Vgrid(i),a_lat);
  }

  data.save("UJparams",raw_ascii);
  return data;
}

int main() {

  int L         = 36;
  double a_lat  = 1; // a_lat omitted in many formulas
  double V0     = 4.5;
  double V0_T   = 20;

  auto qgrid    = linspace<vec>(-1,1,101);
  auto xgrid    = linspace<vec>(-L*a_lat,L*a_lat,1000);
  auto lgrid    = linspace<vec>(-L+1,L-1,2*L-1);

  auto Bloch    = calcBlochFuncs(qgrid,xgrid,V0,L);
  auto Wannier  = calcWannierFunc(Bloch,xgrid);
  auto Bloch_T  = calcBlochFuncs(qgrid,xgrid,V0_T,L);
  auto Wannier_T= calcWannierFunc(Bloch_T,xgrid);

  cx_mat W_tmp  = join_horiz(Wannier,Wannier_T);
  cx_mat W_full = join_horiz(W_tmp,Wannier_T);

  auto U        = calcU(xgrid,W_full);
  auto J        = calcJ(xgrid,Wannier,V0,a_lat);

  std::cout << "U = " << U << std::endl;
  std::cout << "J = " << J << std::endl;
  std::cout << "U/J = " << U/J << std::endl;

  // calculateParameters(L,a_lat,1,20,500);

  Gnuplot gp;
  mat Potplot     = join_horiz(xgrid, calcPotential(xgrid,V0,a_lat));
  mat Blochplot   = join_horiz(xgrid, real(Bloch.col(50)));
  mat Wannierplot = join_horiz(xgrid, real(Wannier));
  mat Wannierplot2= join_horiz(xgrid, getWannierAtSite_j(Wannier,xgrid,-1,a_lat));

  // gp << "set terminal pdf\nset output 'Bloch.pdf'\n";
  gp << "set xrange [-2:2]\n";
  gp << "set ytics nomirror\nset y2tics\n";
	// Data will be sent via a temporary file.  These are erased when you call
	// gp.clearTmpfiles() or when gp goes out of scope.  If you pass a filename
	// (e.g. "gp.file1d(pts, 'mydata.dat')"), then the named file will be created
	// and won't be deleted (this is useful when creating a script).
	gp << "plot"
     << gp.file1d(Potplot) << "with lines title 'Potential' axes x1y2,"
     << gp.file1d(Blochplot) << "with lines title 'Bloch Function',"
     << gp.file1d(Wannierplot) << "with lines title 'Wannier Function',"
     << gp.file1d(Wannierplot2) << "with lines title 'Wannier Function',"
     << std::endl;


  return 0;
}
