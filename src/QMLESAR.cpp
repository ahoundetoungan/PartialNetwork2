// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#define NDEBUG 1

using namespace Rcpp;
using namespace arma;
using namespace std;

// This function draws G for one group
//[[Rcpp::export]]
arma::mat fGm(const arma::mat& dnm,
              const int& Nm){
  arma::mat matunif(Nm, Nm, arma::fill::randu);
  arma::mat Gm = arma::normalise(arma::conv_to<mat>::from((matunif < dnm)), 1, 1);
  Gm.diag()    = arma::zeros(Nm);
  return Gm;
}


//* CASE t1 t2
//* t1 is 0 if no fixed effects and 1 otherwise
//* t2 is 0 is no contextual effects and 1 otherwise
//* Notation
//* A here is S in the paper

//*************************** Part 0: Gy and GX are observed
//********** CASES 01 and 11
//* CASE 01
void fbse0_01(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& VAy,
              arma::vec& VGy,
              arma::mat& cpV,
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const arma::mat& GX,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, dA, invdA, V;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat GXm    = GX.rows(n1, n2);
    Ay               = ym - alpha*Gym;
    V                = arma::join_rows(Xm, GXm);
    VAy             += V.t()*Ay;
    VGy             += V.t()*Gym;
    cpAy            += arma::accu(Ay%Ay);
    GyAy            += arma::accu(Gym%Ay);
    cpV             += V.t()*V;
    for(int r(0); r < R; ++r){
      dG            = fGm(dnm, Nm);
      dA            = Im - alpha*dG;
      invdA         = arma::inv(dA);
      arma::mat t2  = dG*invdA;
      Re           += arma::trace(t2);
    }
  }
  VAy    *= R;
  VGy    *= R;
  cpAy   *= R;
  GyAy   *= R;
  cpV    *= R;
}

//* CASE 11
void fbse0_11(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& VAy,
              arma::vec& VGy,
              arma::mat& cpV,
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const arma::mat& GX,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, dA, invdA, V;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat GXm    = GX.rows(n1, n2);
    Ay               = ym - alpha*Gym;
    Ay              -= mean(Ay);
    V                = arma::join_rows(Xm, GXm);
    V.each_row()    -= arma::mean(V, 0);
    VAy             += V.t()*Ay;
    VGy             += V.t()*Gym;
    cpAy            += arma::accu(Ay%Ay);
    GyAy            += arma::accu(Gym%Ay);
    cpV             += V.t()*V;
    for(int r(0); r < R; ++r){
      dG            = fGm(dnm, Nm);
      dA            = Im - alpha*dG;
      invdA         = arma::inv(dA);
      arma::mat t2  = dG*invdA;
      t2.each_row()-= arma::mean(t2, 0);
      Re           += arma::trace(t2);
    }
  }
  VAy    *= R;
  VGy    *= R;
  cpAy   *= R;
  GyAy   *= R;
  cpV    *= R;
}

//* foc alpha CASES 01 and 11
//[[Rcpp::export]]
double fcpal0_1(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& X,
                const arma::mat& GX,
                const int& Kx, 
                const int& Kx1, 
                const int& Kx2, 
                const int& M,
                const arma::vec& Ncum,
                const bool FE = false){
  int Kexo(Kx + Kx1);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec VAy(Kexo, arma::fill::zeros), VGy(Kexo, arma::fill::zeros);
  arma::mat cpV(Kexo, Kexo, arma::fill::zeros);
  if(FE){
    fbse0_11(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, GX, M, Ncum);
    Rf = R*(Ncum(M) - M);
  } else {
    fbse0_01(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, GX, M, Ncum);
    Rf = R*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpV, VAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*VAy - cpV*beta)))/Rf;
  // cout<<"se2 "<<se2<<endl;
  //if(se2 < 0) return R_PosInf;
  return (GyAy - arma::accu(beta%VGy) - se2*Re)/Rf;
}

//[[Rcpp::export]]
List efcpal0_1(const double& alpha,
         const int& R,
         List& distr, 
         List& Ilist, 
         const arma::vec& y, 
         const arma::vec& Gy,
         const arma::mat& X,
         const arma::mat& GX,
         const int& Kx, 
         const int& Kx1, 
         const int& Kx2, 
         const int& M,
         const arma::vec& Ncum,
         const bool FE = false){
  int Kexo(Kx + Kx1);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec VAy(Kexo, arma::fill::zeros), VGy(Kexo, arma::fill::zeros);
  arma::mat cpV(Kexo, Kexo, arma::fill::zeros);
  if(FE){
    fbse0_11(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, GX, M, Ncum);
    Rf = R*(Ncum(M) - M);
  } else {
    fbse0_01(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, GX, M, Ncum);
    Rf = R*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpV, VAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*VAy - cpV*beta)))/Rf;
  return List::create(Named("beta") = beta, Named("se2") = se2);
}

//********** CASES 00 and 10
//* CASE 00
void fbse0_00(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& VAy,
              arma::vec& VGy,
              arma::mat& cpV,
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, dA, invdA, V;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    Ay               = ym - alpha*Gym;
    V                = X.rows(n1, n2);
    VAy             += V.t()*Ay;
    VGy             += V.t()*Gym;
    cpAy            += arma::accu(Ay%Ay);
    GyAy            += arma::accu(Gym%Ay);
    cpV             += V.t()*V;
    for(int r(0); r < R; ++r){
      dG            = fGm(dnm, Nm);
      dA            = Im - alpha*dG;
      invdA         = arma::inv(dA);
      arma::mat t2  = dG*invdA;
      Re           += arma::trace(t2);
    }
  }
  VAy    *= R;
  VGy    *= R;
  cpAy   *= R;
  GyAy   *= R;
  cpV    *= R;
}

//* CASE 11
void fbse0_10(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& VAy,
              arma::vec& VGy,
              arma::mat& cpV,
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, dA, invdA, V;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    Ay               = ym - alpha*Gym;
    Ay              -= mean(Ay);
    V                = X.rows(n1, n2);
    V.each_row()    -= arma::mean(V, 0);
    VAy             += V.t()*Ay;
    VGy             += V.t()*Gym;
    cpAy            += arma::accu(Ay%Ay);
    GyAy            += arma::accu(Gym%Ay);
    cpV             += V.t()*V;
    for(int r(0); r < R; ++r){
      dG            = fGm(dnm, Nm);
      dA            = Im - alpha*dG;
      invdA         = arma::inv(dA);
      arma::mat t2  = dG*invdA;
      t2.each_row()-= arma::mean(t2, 0);
      Re           += arma::trace(t2);
    }
  }
  VAy    *= R;
  VGy    *= R;
  cpAy   *= R;
  GyAy   *= R;
  cpV    *= R;
}

//* foc alpha CASES 01 and 11
//[[Rcpp::export]]
double fcpal0_0(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& X,
                const int& Kx, 
                const int& M,
                const arma::vec& Ncum,
                const bool FE = false){
  int Kexo(Kx);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec VAy(Kexo, arma::fill::zeros), VGy(Kexo, arma::fill::zeros);
  arma::mat cpV(Kexo, Kexo, arma::fill::zeros);
  if(FE){
    fbse0_10(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, M, Ncum);
    Rf = R*(Ncum(M) - M);
  } else {
    fbse0_00(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, M, Ncum);
    Rf = R*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpV, VAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*VAy - cpV*beta)))/Rf;
  // cout<<"se2 "<<se2<<endl;
  //if(se2 < 0) return R_PosInf;
  return (GyAy - arma::accu(beta%VGy) - se2*Re)/Rf;
}

//[[Rcpp::export]]
List efcpal0_0(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& X,
                const int& Kx, 
                const int& M,
                const arma::vec& Ncum,
                const bool FE = false){
  int Kexo(Kx);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec VAy(Kexo, arma::fill::zeros), VGy(Kexo, arma::fill::zeros);
  arma::mat cpV(Kexo, Kexo, arma::fill::zeros);
  if(FE){
    fbse0_10(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, M, Ncum);
    Rf = R*(Ncum(M) - M);
  } else {
    fbse0_00(alpha, cpAy, GyAy, VAy, VGy, cpV, Re, R, distr, Ilist, y, Gy, X, M, Ncum);
    Rf = R*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpV, VAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*VAy - cpV*beta)))/Rf;
  return List::create(Named("beta") = beta, Named("se2") = se2);
}

//*************************** Part 1: Gy is not observed and cols in GX are observed
//********** CASES 01 and 11
//* CASE 01
void fbse1_01(const double& alpha,
              double& cpdAy, // cp is cross prod
              double& dGydAy,
              arma::vec& dVdAy,
              arma::vec& dVdGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV (dA inv(ddA) dddV - ddV)
              arma::mat& Rb,  //cp(dA inv(ddA) dddV - ddV)
              double& Rc,     //trace(cp(dA inv(ddA)))
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dA*inv(ddA)*dddV - ddV)
              double& Re,     //trace(trans(dG*inv(ddA))*dA*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X,
              const arma::mat& X1,    //GX1 is observed
              const arma::mat& GX1,
              const arma::mat& X2, 
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, dA, ddA, invddA, dGX2, ddGX1, ddGX2, dV, ddV, dddV;
  arma::vec  dGy, dAy;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat X1m    = X1.rows(n1, n2);
    arma::mat X2m    = X2.rows(n1, n2);
    arma::mat GX1m   = GX1.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dA             = Im - alpha*dG;
      dGX2           = dG*X2m;
      dGy            = dG*ym;
      dAy            = ym - alpha*dGy;
      dV             = arma::join_rows(arma::join_rows(Xm, GX1m), dGX2);
      dVdAy         += dV.t()*dAy;
      dVdGy         += dV.t()*dGy;
      cpdAy         += arma::accu(dAy%dAy);
      dGydAy        += arma::accu(dGy%dAy);
      cpdV          += dV.t()*dV;
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        ddGX1        = ddG*X1m;
        ddGX2        = ddG*X2m;
        arma::mat t1 = arma::join_rows(Xm, ddGX1);
        ddV          = arma::join_rows(t1, dGX2);        
        dddV         = arma::join_rows(t1, ddGX2);  
        arma::mat t2 = dG*invddA;
        arma::mat t3 = invddA - alpha*t2;//dA invddA
        arma::mat t4 = t3*dddV - ddV; // dA invddA dddV - ddV
        Ra          += ddV.t()*t4;
        Rb          += t4.t()*t4;
        Rc          += arma::trace(t3.t()*t3);
        Rd          += arma::trans(t2*dddV)*t4;
        Re          += arma::trace(t2.t()*t3);
      }
    }
  }
  dVdAy  *= S;
  dVdGy  *= S;
  cpdAy  *= S;
  dGydAy *= S;
  cpdV   *= S;
}

//* CASE 11
void fbse1_11(const double& alpha,
              double& cpdAy,
              double& dGydAy,
              arma::vec& dVdAy,
              arma::vec& dVdGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV dA inv(ddA) dddV
              arma::mat& Rb,  //cp(dA inv(ddA) dddV - ddV)
              double& Rc,     //trace(cp(dA inv(ddA)))
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dA*inv(ddA)*dddV - ddV)
              double& Re,     //trace(trans(dG*inv(ddA))*dA*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X,
              const arma::mat& X1, 
              const arma::mat& GX1,
              const arma::mat& X2, 
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, dA, ddA, invddA, dGX2, ddGX1, ddGX2, dV, ddV, dddV;
  arma::vec  dGy, dAy;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat X1m    = X1.rows(n1, n2);
    arma::mat X2m    = X2.rows(n1, n2);
    arma::mat GX1m   = GX1.rows(n1, n2);
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dA             = Im - alpha*dG;
      dGX2           = dG*X2m;
      dGy            = dG*ym;
      dAy            = ym - alpha*dGy;
      dV             = arma::join_rows(arma::join_rows(Xm, GX1m), dGX2); 
      dAy           -= mean(dAy);
      dV.each_row() -= arma::mean(dV, 0);
      dVdAy         += dV.t()*dAy;
      dVdGy         += dV.t()*dGy;
      cpdAy         += arma::accu(dAy%dAy);
      dGydAy        += arma::accu(dGy%dAy);
      cpdV          += dV.t()*dV;
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        ddGX1        = ddG*X1m;
        ddGX2        = ddG*X2m;
        arma::mat t1 = arma::join_rows(Xm, ddGX1);
        ddV          = arma::join_rows(t1, dGX2);        
        dddV         = arma::join_rows(t1, ddGX2);  
        arma::mat t2 = dG*invddA;
        arma::mat t3 = invddA - alpha*t2;//dA invddA
        arma::mat t4 = t3*dddV - ddV; // dA invddA dddV - ddV
        ddV.each_row() -= arma::mean(ddV, 0);
        t3.each_row()  -= arma::mean(t3, 0);
        t4.each_row()  -= arma::mean(t4, 0);
        Ra          += ddV.t()*t4;
        Rb          += t4.t()*t4;
        Rc          += arma::trace(t3.t()*t3);
        Rd          += arma::trans(t2*dddV)*t4;
        Re          += arma::trace(t2.t()*t3);
      }
    }
  }
  dVdAy  *= S;
  dVdGy  *= S;
  cpdAy  *= S;
  dGydAy *= S;
  cpdV   *= S;
}

//* foc alpha CASES 01 and 11
//[[Rcpp::export]]
double fcpal1_1(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X,
                const arma::mat& X1, 
                const arma::mat& GX1,
                const arma::mat& X2, 
                const int& Kx, 
                const int& Kx1, 
                const int& Kx2, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  int Kexo(Kx + Kx1 + Kx2);
  double cpdAy(0), dGydAy(0), Rc(0), Re(0), Rf(0);
  arma::vec dVdAy(Kexo, arma::fill::zeros), dVdGy(Kexo, arma::fill::zeros);
  arma::mat cpdV(Kexo, Kexo, arma::fill::zeros),
  Ra(Kexo, Kexo, arma::fill::zeros), Rb(Kexo, Kexo, arma::fill::zeros),
  Rd(Kexo, Kexo, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse1_11(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*(Ncum(M) - M);
  } else {
    fbse1_01(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVdAy);
  double se2      = (cpdAy - arma::accu(beta%(2*dVdAy - (cpdV - Rb)*beta)))/Rc;
  //if(se2 < 0) return R_PosInf;
  return (dGydAy - arma::accu(beta%(dVdGy + Rd*beta)) - se2*Re)/Rf;
}


//[[Rcpp::export]]
List efcpal1_1(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X,
                const arma::mat& X1, 
                const arma::mat& GX1,
                const arma::mat& X2, 
                const int& Kx, 
                const int& Kx1, 
                const int& Kx2, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  int Kexo(Kx + Kx1 + Kx2);
  double cpdAy(0), dGydAy(0), Rc(0), Re(0);
  arma::vec dVdAy(Kexo, arma::fill::zeros), dVdGy(Kexo, arma::fill::zeros);
  arma::mat cpdV(Kexo, Kexo, arma::fill::zeros),
  Ra(Kexo, Kexo, arma::fill::zeros), Rb(Kexo, Kexo, arma::fill::zeros),
  Rd(Kexo, Kexo, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse1_11(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, X1, GX1, X2, M, Ncum);
    // Rf = R*S*(Ncum(M) - M);
  } else {
    fbse1_01(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, X1, GX1, X2, M, Ncum);
    // Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVdAy);
  double se2      = (cpdAy - arma::accu(beta%(2*dVdAy - (cpdV - Rb)*beta)))/Rc;
  return List::create(Named("beta") = beta, Named("se2") = se2);
}

//********** CASES 01 and 11
//* CASE 00
void fbse1_00(const double& alpha,
              double& cpdAy, // cp is cross prod
              double& dGydAy,
              arma::vec& dVdAy,
              arma::vec& dVdGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV (dA inv(ddA) dddV - ddV)
              arma::mat& Rb,  //cp(dA inv(ddA) dddV - ddV)
              double& Rc,     //trace(cp(dA inv(ddA)))
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dA*inv(ddA)*dddV - ddV)
              double& Re,     //trace(trans(dG*inv(ddA))*dA*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, dA, ddA, invddA;
  arma::vec  dGy, dAy;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    cpdV            += Xm.t()*Xm;
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dA             = Im - alpha*dG;
      dGy            = dG*ym;
      dAy            = ym - alpha*dGy;
      dVdAy         += Xm.t()*dAy;
      dVdGy         += Xm.t()*dGy;
      cpdAy         += arma::accu(dAy%dAy);
      dGydAy        += arma::accu(dGy%dAy);
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        arma::mat t2 = dG*invddA;
        arma::mat t3 = invddA - alpha*t2;//dA invddA
        arma::mat t4 = t3*Xm - Xm; // dA invddA dddV - ddV
        Ra          += Xm.t()*t4;
        Rb          += t4.t()*t4;
        Rc          += arma::trace(t3.t()*t3);
        Rd          += arma::trans(t2*Xm)*t4;
        Re          += arma::trace(t2.t()*t3);
      }
    }
  }
  dVdAy  *= S;
  dVdGy  *= S;
  cpdAy  *= S;
  dGydAy *= S;
  cpdV   *= (R*S);
}

//* CASE 10
void fbse1_10(const double& alpha,
              double& cpdAy, // cp is cross prod
              double& dGydAy,
              arma::vec& dVdAy,
              arma::vec& dVdGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV (dA inv(ddA) dddV - ddV)
              arma::mat& Rb,  //cp(dA inv(ddA) dddV - ddV)
              double& Rc,     //trace(cp(dA inv(ddA)))
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dA*inv(ddA)*dddV - ddV)
              double& Re,     //trace(trans(dG*inv(ddA))*dA*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::mat& X,
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, dA, ddA, invddA;
  arma::vec  dGy, dAy;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    Xm.each_row()   -= arma::mean(Xm, 0);
    cpdV            += Xm.t()*Xm;
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dA             = Im - alpha*dG;
      dGy            = dG*ym;
      dAy            = ym - alpha*dGy;
      dAy           -= mean(dAy);
      dVdAy         += Xm.t()*dAy;
      dVdGy         += Xm.t()*dGy;
      cpdAy         += arma::accu(dAy%dAy);
      dGydAy        += arma::accu(dGy%dAy);
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        arma::mat t2 = dG*invddA;
        arma::mat t3 = invddA - alpha*t2;//dA invddA
        arma::mat t4 = t3*Xm - Xm; // dA invddA dddV - ddV
        t4.each_row()  -= arma::mean(t4, 0);
        t3.each_row()  -= arma::mean(t3, 0);
        Ra          += Xm.t()*t4;
        Rb          += t4.t()*t4;
        Rc          += arma::trace(t3.t()*t3);
        Rd          += arma::trans(t2*Xm)*t4;
        Re          += arma::trace(t2.t()*t3);
      }
    }
  }
  dVdAy  *= S;
  dVdGy  *= S;
  cpdAy  *= S;
  dGydAy *= S;
  cpdV   *= (R*S);
}

//* foc alpha CASES 00 and 10
//[[Rcpp::export]]
double fcpal1_0(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X,
                const int& Kx, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  double cpdAy(0), dGydAy(0), Rc(0), Re(0), Rf(0);
  arma::vec dVdAy(Kx, arma::fill::zeros), dVdGy(Kx, arma::fill::zeros);
  arma::mat cpdV(Kx, Kx, arma::fill::zeros),
  Ra(Kx, Kx, arma::fill::zeros), Rb(Kx, Kx, arma::fill::zeros),
  Rd(Kx, Kx, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse1_10(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, M, Ncum);
    Rf = R*S*(Ncum(M) - M);
  } else {
    fbse1_00(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, M, Ncum);
    Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVdAy);
  double se2      = (cpdAy - arma::accu(beta%(2*dVdAy - (cpdV - Rb)*beta)))/Rc;
  //if(se2 < 0) return R_PosInf;
  return (dGydAy - arma::accu(beta%(dVdGy + Rd*beta)) - se2*Re)/Rf;
}

//[[Rcpp::export]]
List efcpal1_0(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::mat& X,
                const int& Kx, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  double cpdAy(0), dGydAy(0), Rc(0), Re(0);
  arma::vec dVdAy(Kx, arma::fill::zeros), dVdGy(Kx, arma::fill::zeros);
  arma::mat cpdV(Kx, Kx, arma::fill::zeros),
  Ra(Kx, Kx, arma::fill::zeros), Rb(Kx, Kx, arma::fill::zeros),
  Rd(Kx, Kx, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse1_10(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, M, Ncum);
    // Rf = R*S*(Ncum(M) - M);
  } else {
    fbse1_00(alpha, cpdAy, dGydAy, dVdAy, dVdGy, cpdV, Ra, Rb, Rc, Rd, Re,
             R, S, distr, Ilist, y, X, M, Ncum);
    // Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVdAy);
  double se2      = (cpdAy - arma::accu(beta%(2*dVdAy - (cpdV - Rb)*beta)))/Rc;
  return List::create(Named("beta") = beta, Named("se2") = se2);
}

//********** Part 2: Gy is observed and cols in GX are observed
// t2 is not important
//********** CASES 01 and 11
//* CASE 01
void fbse2_01(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& dVAy,
              arma::vec& dVGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV (dddV - ddV)
              arma::mat& Rb,  //cp(dddV - ddV)
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dddV - ddV)
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const arma::mat& X1,    //GX1 is observed
              const arma::mat& GX1,
              const arma::mat& X2, 
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, ddA, invddA, dGX2, ddGX1, ddGX2, dV, ddV, dddV;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat X1m    = X1.rows(n1, n2);
    arma::mat X2m    = X2.rows(n1, n2);
    arma::mat GX1m   = GX1.rows(n1, n2);
    Ay               = ym - alpha*Gym;
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dGX2           = dG*X2m;
      dV             = arma::join_rows(arma::join_rows(Xm, GX1m), dGX2);
      dVAy          += dV.t()*Ay;
      dVGy          += dV.t()*Gym;
      cpAy          += arma::accu(Ay%Ay);
      GyAy          += arma::accu(Gym%Ay);
      cpdV          += dV.t()*dV;
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        ddGX1        = ddG*X1m;
        ddGX2        = ddG*X2m;
        arma::mat t1 = arma::join_rows(Xm, ddGX1);
        ddV          = arma::join_rows(t1, dGX2);        
        dddV         = arma::join_rows(t1, ddGX2);  
        arma::mat t2 = ddG*invddA;
        arma::mat t3 = dddV - ddV;
        Ra          += ddV.t()*t3;
        Rb          += t3.t()*t3;
        Rd          += arma::trans(t2*dddV)*t3;
        Re          += arma::trace(t2);
      }
    }
  }
  dVAy   *= S;
  dVGy   *= S;
  cpAy   *= S;
  GyAy   *= S;
  cpdV   *= S;
}

//* CASE 11
void fbse2_11(const double& alpha,
              double& cpAy, // cp is cross prod
              double& GyAy,
              arma::vec& dVAy,
              arma::vec& dVGy,
              arma::mat& cpdV,
              arma::mat& Ra,  //ddV (dddV - ddV)
              arma::mat& Rb,  //cp(dddV - ddV)
              arma::mat& Rd,  //trans(dG*inv(ddA)*dddV)*(dddV - ddV)
              double& Re,     //trace(dG*inv(ddA))
              const int& R,
              const int& S,
              List& distr, 
              List& Ilist, 
              const arma::vec& y, 
              const arma::vec& Gy,
              const arma::mat& X,
              const arma::mat& X1,    //GX1 is observed
              const arma::mat& GX1,
              const arma::mat& X2, 
              const int& M, 
              const arma::vec& Ncum){
  arma::mat dG, ddG, dA, ddA, invddA, dGX2, ddGX1, ddGX2, dV, ddV, dddV;
  arma::vec Ay;
  for(int m(0); m < M; ++m){
    int n1           = Ncum(m);
    int n2           = Ncum(m + 1) - 1;
    int Nm           = n2 - n1 + 1;
    arma::mat dnm    = distr(m);
    arma::mat Im     = Ilist(m);
    arma::vec ym     = y.subvec(n1, n2);
    arma::vec Gym    = Gy.subvec(n1, n2);
    arma::mat Xm     = X.rows(n1, n2);
    arma::mat X1m    = X1.rows(n1, n2);
    arma::mat X2m    = X2.rows(n1, n2);
    arma::mat GX1m   = GX1.rows(n1, n2);
    Ay               = ym - alpha*Gym;
    Ay              -= mean(Ay);
    for(int r(0); r < R; ++r){
      dG             = fGm(dnm, Nm);
      dA             = Im - alpha*dG;
      dGX2           = dG*X2m;
      dV             = arma::join_rows(arma::join_rows(Xm, GX1m), dGX2);
      dV.each_row() -= arma::mean(dV, 0);
      dVAy          += dV.t()*Ay;
      dVGy          += dV.t()*Gym;
      cpAy          += arma::accu(Ay%Ay);
      GyAy          += arma::accu(Gym%Ay);
      cpdV          += dV.t()*dV;
      for(int s(0); s < S; ++s){
        ddG          = fGm(dnm, Nm);
        ddA          = Im - alpha*ddG;
        invddA       = arma::inv(ddA);
        ddGX1        = ddG*X1m;
        ddGX2        = ddG*X2m;
        arma::mat t1 = arma::join_rows(Xm, ddGX1);
        ddV          = arma::join_rows(t1, dGX2);        
        dddV         = arma::join_rows(t1, ddGX2);  
        arma::mat t2 = ddG*invddA;
        arma::mat t3 = dddV - ddV;
        ddV.each_row() -= arma::mean(ddV, 0);
        t2.each_row()  -= arma::mean(t2, 0);
        t3.each_row()  -= arma::mean(t3, 0);
        Ra          += ddV.t()*t3;
        Rb          += t3.t()*t3;
        Rd          += arma::trans(t2*dddV)*t3;
        Re          += arma::trace(t2);
      }
    }
  }
  dVAy   *= S;
  dVGy   *= S;
  cpAy   *= S;
  GyAy   *= S;
  cpdV   *= S;
}

//* foc alpha CASES 01 and 11
//[[Rcpp::export]]
double fcpal2_1(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& X,
                const arma::mat& X1, 
                const arma::mat& GX1,
                const arma::mat& X2, 
                const int& Kx, 
                const int& Kx1, 
                const int& Kx2, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  int Kexo(Kx + Kx1 + Kx2);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec dVAy(Kexo, arma::fill::zeros), dVGy(Kexo, arma::fill::zeros);
  arma::mat cpdV(Kexo, Kexo, arma::fill::zeros),
  Ra(Kexo, Kexo, arma::fill::zeros), Rb(Kexo, Kexo, arma::fill::zeros),
  Rd(Kexo, Kexo, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse2_11(alpha, cpAy, GyAy, dVAy, dVGy, cpdV, Ra, Rb,  Rd, Re, R, S, distr, 
             Ilist, y, Gy, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*(Ncum(M) - M);
  } else {
    fbse2_01(alpha, cpAy, GyAy, dVAy, dVGy, cpdV, Ra, Rb,  Rd, Re, R, S, distr, 
             Ilist, y, Gy, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*dVAy - (cpdV - Rb)*beta)))/Rf;
  // cout<<"se2 "<<se2<<endl;
  //if(se2 < 0) return R_PosInf;
  return (GyAy - arma::accu(beta%(dVGy + Rd*beta)) - se2*Re)/Rf;
}


//[[Rcpp::export]]
List efcpal2_1(const double& alpha,
                const int& R,
                List& distr, 
                List& Ilist, 
                const arma::vec& y, 
                const arma::vec& Gy,
                const arma::mat& X,
                const arma::mat& X1, 
                const arma::mat& GX1,
                const arma::mat& X2, 
                const int& Kx, 
                const int& Kx1, 
                const int& Kx2, 
                const int& M, 
                const arma::vec& Ncum,
                const bool FE = false,
                const int& S = 1L){
  int Kexo(Kx + Kx1 + Kx2);
  double cpAy(0), GyAy(0), Re(0), Rf(0);
  arma::vec dVAy(Kexo, arma::fill::zeros), dVGy(Kexo, arma::fill::zeros);
  arma::mat cpdV(Kexo, Kexo, arma::fill::zeros),
  Ra(Kexo, Kexo, arma::fill::zeros), Rb(Kexo, Kexo, arma::fill::zeros),
  Rd(Kexo, Kexo, arma::fill::zeros);
  //Ra = ddV(dSinv(ddS)dddV - ddV); Rb = cp(dSinv(ddS)dddV - ddV); Rc = trace(cpdS invddS)
  //Rd = trans(dG*inv(ddS)*dddV)*(dSinv(ddS)dddV - ddV); 
  //Rf = trace(trans(dG*inv(ddS))*dS*inv(ddS))
  if(FE){
    fbse2_11(alpha, cpAy, GyAy, dVAy, dVGy, cpdV, Ra, Rb,  Rd, Re, R, S, distr, 
             Ilist, y, Gy, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*(Ncum(M) - M);
  } else {
    fbse2_01(alpha, cpAy, GyAy, dVAy, dVGy, cpdV, Ra, Rb,  Rd, Re, R, S, distr, 
             Ilist, y, Gy, X, X1, GX1, X2, M, Ncum);
    Rf = R*S*Ncum(M);
  }
  
  arma::vec beta  = arma::solve(cpdV + Ra, dVAy);
  // cout<<"beta "<<beta.t()<<endl;
  double se2      = (cpAy - arma::accu(beta%(2*dVAy - (cpdV - Rb)*beta)))/Rf;
  return List::create(Named("beta") = beta, Named("se2") = se2);
}

