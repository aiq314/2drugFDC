#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include "cellmodel.hpp"
#include "../enums/enum_Ohara_Rudy_2011.hpp"
#include "../modules/dcomb_t.hpp"

class Ohara_Rudy_2011 : public Cellmodel
{
public:
  Ohara_Rudy_2011();
  ~Ohara_Rudy_2011();
  void initConsts ();
  void initConsts (double celltype);
  void initConsts (bool is_dutta);
  void initConsts (double celltype, bool is_dutta);
  void initConsts (double celltype, double conc, double *ic50, bool is_dutta);
  void initConsts (double celltype, double conc_a, double conc_b, 
	double *ic50_a, double *ic50_b, double ddis_effect, double *cvar,
	bool is_dutta, bool is_ddis, bool is_syntopic, bool is_cvar);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
private:
  // apply Dutta conductance scaling based on Dutta et al. 2017
  void ___applyDutta();
  // apply drug-induced equation based on drug-induced equation
  // epsilon is used to avoid NaN value data in ic50
  // (denominator is zero, app is broken)
  void ___applyDrugEffect(double conc, double *ic50, double epsilon);
  // actual initial condition function
  // that will be called by public functions
  void ___initConsts();
  // prompt the info of celltype
  // 0 is endo, 1 is epi, 2 is M cell
  void ___printCelltype(int celltype);
  /*==============*/
  /* Added by ALI */
  /*==============*/
  void ___applyCvar(double *cvar);
  void ___applyDrugCombination(const double conc_a, const double conc_b, double *ic50_a, double *ic50_b, double ddis_la_effect, bool is_ddis_la, bool is_syntopic, double epsilon);
};

#endif

