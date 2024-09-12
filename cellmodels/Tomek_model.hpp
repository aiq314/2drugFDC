#ifndef TOMEK_MODEL_ENDO_HPP
#define TOMEK_MODEL_ENDO_HPP

#include "cellmodel.hpp"
#include "../enums/enum_Tomek_model.hpp"

class Tomek_model : public Cellmodel
{
public:
  Tomek_model();
  ~Tomek_model();
  void initConsts ();
  void initConsts (double celltype);
  void initConsts (double celltype, double conc, double *ic50);
  void initConsts (double celltype, double conc_a, double conc_b, double *ic50_a, double *ic50_b, double ddis_la_effect, bool is_dutta, bool is_ddis_la);
  void computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC );
  void solveAnalytical( double dt );
private:
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
  void ___applyDrugCombination(const double conc_a, const double conc_b, double *ic50_a, double *ic50_b, double ddis_la_effect, bool is_ddis_la, double epsilon);
};



#endif

