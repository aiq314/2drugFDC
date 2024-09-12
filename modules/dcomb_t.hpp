#ifndef DCOMB_T_HPP
#define DCOMB_T_HPP

struct dcomb_t{
  double c_a;
  double c_b;
  double emax_a;
  double emax_b;
  double ec50_a;
  double ec50_b;
  double n_a;
  double n_b;
  double int_ab;
  double int_ba;
  double ec50_int_ab;
  double ec50_int_ba;
  double n_int_ab;
  double n_int_ba;
  double is_la;
  double is_ec50_bi;

  dcomb_t();
  dcomb_t(const dcomb_t & source);
  dcomb_t& operator=(const dcomb_t & source);
  void copy(const dcomb_t &source);
  void init(const double c_a, const double c_b); 
  void init(const double c_a_val, const double c_b_val,
	    const double emax_a_val, const double emax_b_val,
	    const double ec50_a_val, const double ec50_b_val,
	    const double n_a_val, const double n_b_val,
	    const double int_ab_val, const double int_ba_val,
	    const double ec50_int_ab_val, const double ec50_int_ba_val,
	    const double n_int_ab_val, const double n_int_ba_val,
	    const bool is_la_val, const bool is_ec50_bi_val);
};

double ec50_bi(void* param_data);
double emax_bi(void* param_data);
double loewe_additivity(unsigned n, const double* x, double* grad, void* param_data);
double calc_loewe_additivity(void* param_data);
double ecomb(void* param_data);

#endif
