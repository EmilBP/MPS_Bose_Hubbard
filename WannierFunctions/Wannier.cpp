#include <armadillo>
#include <complex>
#include <cmath>

#include "gnuplot-iostream.h"

using namespace arma;

mat Hamiltonian(int L, double V0, double q){
  vec lgrid   = linspace<vec>(-L,L,2*L+1);
  mat H       = zeros<mat>(2*L+1,2*L+1);
  H.diag()    = (q*q+V0/2.0)*ones<vec>(2*L+1) + 4.0*q*lgrid + 4.0*pow(lgrid,2);
  H.diag(-1)  = -V0/4.0*ones<vec>(2*L);
  H.diag(+1)  = -V0/4.0*ones<vec>(2*L);
  return H;
}

cx_mat calcBlochFuncs(vec& qvals, vec& xvals, double V0, int L){
  std::complex<double> cplx_i(0,1);
  mat H, eigvec, coeff(2*L+1,qvals.n_rows);
  vec eigval, lgrid = linspace<vec>(-L,L,2*L+1);
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
      BlochFunc(j,i) = dot(exp(cplx_i*(qvals(i)+2.0*lgrid)*datum::pi*xvals(j)),coeff.col(i));
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
  vec site = xgrid-0.5*j;
  vec realw = real(w0);
  interp1(xgrid,realw,site,wj);
  return wj;
}

vec calcPotential(vec& xgrid, double V0, double a_lat){
  return V0*pow(sin(datum::pi/a_lat*xgrid),2);
}

double calcU(vec& xgrid, cx_vec& w){
  double g = 1;
  vec normsq = cdot(w,w);
  return g*as_scalar( trapz( pow( normsq ,2)));
}

double calcJ(){
  double dx = xgrid(1)-xgrid(0);
  vec del(w.n_rows);
  for (size_t i = 1; i < w.n_rows-1; i++) {
    del(i) = (w(i+1) + w(i-1) - 2.0*w(i))/dx/dx;
  }
  del(0) = del(w.n_rows-1) = 0;
  
}

int main() {

  int L         = 10;
  double a_lat  = 1; // a_lat omitted in many formulas
  double V0     = 4;

  auto qgrid    = linspace<vec>(-1,1,101);
  auto xgrid    = linspace<vec>(-L*a_lat,L*a_lat,1000);
  auto lgrid    = linspace<vec>(-L,L,2*L+1);

  auto Bloch    = calcBlochFuncs(qgrid,xgrid,V0,L);
  auto Wannier  = calcWannierFunc(Bloch,xgrid);

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
