#ifndef PARAM_HPP
#define PARAM_HPP

struct param_t
{
  bool is_dutta; // TRUE if using Dutta scaling
  bool is_using_output; // TRUE if using last output file
  double bcl; // basic cycle length
  unsigned short pace_max; // maximum pace
  double celltype;  // cell types
  double dt;        // time step
  double dt_write;  // writing step
  double inet_vm_threshold; // Vm threshold for calculating inet
  void init();
  void show_val();
  
  /*==============*/
  /* Added by ALI */
  /*==============*/
  char cvar_file[1024];
  char ddis_file[1024];
  char hill_file_a[1024];
  char hill_file_b[1024];
  char drug_name_a[100];
  char drug_name_b[100];
  double cmax_a;
  double cmax_b;
  char radius[100];
  char angles[100];
  bool is_cvar; // TRUE if using cvar file
  bool is_ddis_la; // TRUE if using ddis file
  bool is_syntopic; // TRUE if using syntopic model, FALSE if using allotopic model
};

#endif
