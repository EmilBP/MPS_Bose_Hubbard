#ifndef OCBOSEHUBBARD_NLP_HPP
#define OCBOSEHUBBARD_NLP_HPP

#include "IpTNLP.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"

#include "OptimalControl.hpp"
#include "ControlBasis.hpp"
#include "TimeStepperTEBDfast.hpp"
#include "HamiltonianBH.hpp"
#include <assert.h>
#include <string>
#include <fstream>

using namespace Ipopt;
using OC_BH = OptimalControl<TimeStepperTEBDfast,HamiltonianBH>;

class OCBoseHubbard_nlp : public TNLP
{
private:
  OC_BH optControlProb;
  ControlBasis bControl;
  std::vector< std::vector<double> > controlCache;;

public:
  /** default constructor */
  OCBoseHubbard_nlp(OC_BH& optControlProb, ControlBasis& bControl);

  /** default destructor */
  virtual ~OCBoseHubbard_nlp();

  virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

  virtual bool get_bounds_info(Ipopt::Index n, Number* x_l, Number* x_u,
                               Ipopt::Index m, Number* g_l, Number* g_u);

  virtual bool get_starting_point(Ipopt::Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Ipopt::Index m, bool init_lambda, Number* lambda);

  virtual bool eval_f(Ipopt::Index n, const Number* x,
                      bool new_x, Number& obj_value);

  virtual bool eval_grad_f(Ipopt::Index n, const Number* x, bool new_x, Number* grad_f);

  virtual bool eval_g(Ipopt::Index n, const Number* x,
                      bool new_x, Ipopt::Index m, Number* g);

  virtual bool eval_jac_g(Ipopt::Index n, const Number* x, bool new_x,
                          Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                          Ipopt::Index *jCol, Number* values);

  virtual void finalize_solution(SolverReturn status, Ipopt::Index n,
                                 const Number* x, const Number* z_L,
                                 const Number* z_U, Ipopt::Index m, const Number* g,
                                 const Number* lambda, Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq);

  virtual bool intermediate_callback(AlgorithmMode mode,
                               Ipopt::Index iter, Number obj_value,
                               Number inf_pr, Number inf_du,
                               Number mu, Number d_norm,
                               Number regularization_size,
                               Number alpha_du, Number alpha_pr,
                               Ipopt::Index ls_trials,
                               const IpoptData* ip_data,
                               IpoptCalculatedQuantities* ip_cq);

};


#endif
