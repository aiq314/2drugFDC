#include "dcomb_t.hpp"
#include "commons.hpp"
#include <nlopt.h>
#include <math.h>

dcomb_t::dcomb_t(){
}

dcomb_t::dcomb_t(const dcomb_t & source){
  copy(source);
}

dcomb_t& dcomb_t::operator=(const dcomb_t & source){
  if( this != &source ) copy(source);
  return *this;
}

void dcomb_t::copy(const dcomb_t &source){
  c_a = source.c_a;
  c_b = source.c_b;
  emax_a = source.emax_a;
  emax_b = source.emax_b;
  ec50_a = source.ec50_a;
  ec50_b = source.ec50_b;
  n_a = source.n_a;
  n_b = source.n_b;
  int_ab = source.int_ab;
  int_ba = source.int_ba;
  ec50_int_ab = source.ec50_int_ab;
  ec50_int_ba = source.ec50_int_ba;
  n_int_ab = source.n_int_ab;
  n_int_ba = source.n_int_ba;
  is_la = source.is_la;
  is_ec50_bi = source.is_ec50_bi;
}

void dcomb_t::init(const double c_a_val, const double c_b_val){
  c_a = c_a_val;
  c_b = c_b_val;
  emax_a = 1.0;
  emax_b = 1.0;
  ec50_a = 1.0;
  ec50_b = 1.0;
  n_a = 1.0;
  n_b = 1.0;
  int_ab = 1.0;
  int_ba = 1.0;
  ec50_int_ab = 1.0;
  ec50_int_ba = 1.0;
  n_int_ab = 1.0;
  n_int_ba = 1.0;
  is_la = true;
  is_ec50_bi = false;
}


void dcomb_t::init(const double c_a_val, const double c_b_val,
		   const double emax_a_val, const double emax_b_val,
		   const double ec50_a_val, const double ec50_b_val,
		   const double n_a_val, const double n_b_val,
		   const double int_ab_val, const double int_ba_val,
		   const double ec50_int_ab_val, const double ec50_int_ba_val,
		   const double n_int_ab_val, const double n_int_ba_val,
		   const bool is_la_val, const bool is_ec50_bi_val){
  c_a = c_a_val;
  c_b = c_b_val;
  emax_a = emax_a_val;
  emax_b = emax_b_val;
  ec50_a = ec50_a_val;
  ec50_b = ec50_b_val;
  n_a = n_a_val;
  n_b = n_b_val;
  int_ab = int_ab_val;
  int_ba = int_ba_val;
  ec50_int_ab = ec50_int_ab_val;
  ec50_int_ba = ec50_int_ba_val;
  n_int_ab = n_int_ab_val;
  n_int_ba = n_int_ba_val;
  is_la = is_la_val;
  is_ec50_bi = is_ec50_bi_val;
}
double ec50_bi(void* param_data){
  dcomb_t* d = (dcomb_t*)param_data;
  double ea, eb;
  double func;

  ea = d->emax_a / (pow(d->ec50_a * (1 + d->int_ab * pow(d->c_b, d->n_int_ab) / (pow(d->ec50_int_ab, d->n_int_ab) + pow(d->c_b, d->n_int_ab))) / d->c_a, d->n_a) + 1.0);
  eb = d->emax_b / (pow(d->ec50_b * (1 + d->int_ba * pow(d->c_a, d->n_int_ba) / (pow(d->ec50_int_ba, d->n_int_ba) + pow(d->c_a, d->n_int_ba))) / d->c_b, d->n_b) + 1.0);
  return func = ea + eb - ea * eb;
}

double emax_bi(void* param_data){
  dcomb_t* d = (dcomb_t*)param_data;
  double ea, eb;
  double func;

  ea = d->emax_a * (1 + d->int_ab * pow(d->c_b, d->n_int_ab) / (pow(d->ec50_int_ab, d->n_int_ab) + pow(d->c_b, d->n_int_ab))) / (pow(d->ec50_a / d->c_a, d->n_a) + 1.0);
    eb = d->emax_b * (1 + d->int_ba * pow(d->c_a, d->n_int_ba) / (pow(d->ec50_int_ba, d->n_int_ba) + pow(d->c_a, d->n_int_ba))) / (pow(d->ec50_b / d->c_b, d->n_b) + 1.0);
  return func = ea + eb - ea * eb;
}

double loewe_additivity(unsigned n, const double* x, double* grad, void* param_data){
  //x[0] is the combined effects
  //++count;
  dcomb_t* d = (dcomb_t*)param_data;
  double func;
    
  func = fabs(1 
    - d->c_a / (d->ec50_a * (1 + d->int_ab * pow(d->c_b, d->n_int_ab) / (pow(d->ec50_int_ab, d->n_int_ab) + pow(d->c_b, d->n_int_ab))) * (x[0] / pow(d->emax_a - x[0], 1.0 / d->n_a)))
    - d->c_b / (d->ec50_b * (1 + d->int_ba * pow(d->c_a, d->n_int_ba) / (pow(d->ec50_int_ba, d->n_int_ba) + pow(d->c_a, d->n_int_ba))) * (x[0] / pow(d->emax_b - x[0], 1.0 / d->n_b))));
  return func;
}

double calc_loewe_additivity(void* param_data){
  double lb[1] = { 0 }; /* lower bounds */
  double ub[1] = { 1 }; /* upper bounds */
  nlopt_opt opt;

  opt = nlopt_create(NLOPT_LN_COBYLA, 1); /* algorithm and dimensionality */
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  dcomb_t* data = (dcomb_t*)param_data;

  nlopt_set_min_objective(opt, loewe_additivity, data);
  nlopt_set_xtol_rel(opt, 1e-7);

  double x[1];  /* `*`some` `initial` `guess`*` */
  double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */

  /* calculating the combinational effect of drug by Loewe additivity (LA) */
  /* some initial value */
  x[0] = 0.0;//Combined effect
  if (nlopt_optimize(opt, x, &minf) < 0) {
    mpi_printf(0, "nlopt failed at C_A = %g and C_B = %g !\n", data->c_a, data->c_b);
  }
    //else {
    //    //    std::cout << data->C_A << "\t" << data->C_B << "\t" << x[0] << "\t" << minf << std::endl;
    //        //}
    
  nlopt_destroy(opt);
  return x[0];
}

double ecomb(void* param_data) {
  double result;
  dcomb_t* data = (dcomb_t*)param_data;
  if (data->is_la == true) {
    result = calc_loewe_additivity(data);
  }
  else {
    if (data->is_ec50_bi == true) {
      result = ec50_bi(data);
    }
    else {
      result = emax_bi(data);
    }
  }
  return result;
}

