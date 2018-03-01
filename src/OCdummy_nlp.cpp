#include "OCdummy_nlp.hpp"

// constructor
OCdummy_nlp::OCdummy_nlp(OptimalControlDummy& optControlProb, ControlBasis& bControl)
 : optControlProb(optControlProb), bControl(bControl) {}

//destructor
OCdummy_nlp::~OCdummy_nlp()
{}

bool OCdummy_nlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in hs071_NLP.hpp has 4 variables, x[0] through x[3]
  n = bControl.getM();

  // one equality constraint and one inequality constraint
  m = 0;

  // in this example the Jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 0;

  // the Hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = 0;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

bool OCdummy_nlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.

  // the variables have no lower bounds
  for (Index i=0; i<n; i++)
    x_l[i] = -2e19;

  // the variables have no upper bounds
  for (Index i=0; i<n; i++)
    x_u[i] = 2e19;

  return true;
}

bool OCdummy_nlp::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish to use a warmstart option
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  auto c = bControl.getCArray();
  std::copy(c.begin(), c.end(), x);

  return true;
}

bool OCdummy_nlp::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{

  std::vector<double> input;
  input.assign(x,x+n);
  bControl.setCArray(input);

  obj_value = optControlProb.getCost(bControl);

  return true;
}

bool OCdummy_nlp::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  std::vector<double> input;
  input.assign(x,x+n);
  bControl.setCArray(input);

  auto grad = optControlProb.getAnalyticGradient(bControl);

  std::copy(grad.second.begin(), grad.second.end(), grad_f);

  return true;
}

bool OCdummy_nlp::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{

  return true;
}

bool OCdummy_nlp::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{


  return true;
}

void OCdummy_nlp::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L,
                                  const Number* z_U, Index m, const Number* g,
                                  const Number* lambda, Number obj_value,
                                  const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  printf("\n\nSolution of the primal variables, x\n");
  for (Index i=0; i<n; i++) {
    printf("x[%d] = %e\n", i, x[i]);
  }

  printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
  for (Index i=0; i<n; i++) {
    printf("z_L[%d] = %e\n", i, z_L[i]);
  }
  for (Index i=0; i<n; i++) {
    printf("z_U[%d] = %e\n", i, z_U[i]);
  }

  printf("\n\nObjective value\n");
  printf("f(x*) = %e\n", obj_value);
}
