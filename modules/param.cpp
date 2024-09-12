#include "param.hpp"

#include <cstdio>
#include "commons.hpp"

void param_t::init()
{
  is_dutta = false;
  is_using_output = false;
  bcl = 2000.;
  pace_max = 1000;
  celltype = 0.;
  dt = 0.1;
  dt_write = 1.0;
  inet_vm_threshold = -88.0;
  //snprintf(hill_file, sizeof(hill_file), "%s", "./drugs/bepridil/IC50_samples10.csv");
  //snprintf(concs, sizeof(concs), "%s", "0.0");
  
  /*==============*/
  /* Added by ALI */
  /*==============*/
  snprintf(ddis_file, sizeof(ddis_file), "%s", "./DDIs/data_LA_10.csv");
  snprintf(cvar_file, sizeof(cvar_file), "%s", "./cvars/H_M.csv");
  snprintf(hill_file_a, sizeof(hill_file_a), "%s", "./drugs/bepridil/IC50_samples10.csv");
  snprintf(hill_file_b, sizeof(hill_file_b), "%s", "./drugs/bepridil/IC50_samples10.csv");
  snprintf(drug_name_a, sizeof(drug_name_a), "%s", "bepridil");
  snprintf(drug_name_b, sizeof(drug_name_b), "%s", "verapamil");
  cmax_a = 0.0; 
  cmax_b = 0.0; 
  snprintf(radius, sizeof(radius), "%s", "0.0");
  snprintf(angles, sizeof(angles), "%s", "0.0");
  is_ddis_la = false;
  is_syntopic = false;
  is_cvar = false;
}

void param_t::show_val()
{
  //mpi_printf( 0, "%s -- %s\n", "Hill File", hill_file );
  mpi_printf( 0, "%s -- %lf\n", "Celltype", celltype);
  mpi_printf( 0, "%s -- %d\n", "Is_Dutta", is_dutta );
  mpi_printf( 0, "%s -- %d\n", "Is_Using_Output", is_using_output );
  mpi_printf( 0, "%s -- %lf\n", "Basic_Cycle_Length", bcl);
  mpi_printf( 0, "%s -- %d\n", "Number_of_Pacing", pace_max);
  mpi_printf( 0, "%s -- %lf\n", "Time_Step", dt);
  mpi_printf( 0, "%s -- %lf\n", "Inet_Vm_Threshold", inet_vm_threshold);
  mpi_printf( 0, "%s -- %lf\n", "Writing_Step", dt_write);
  //mpi_printf( 0, "%s -- %s\n", "Drug_Name", drug_name);
  //mpi_printf( 0, "%s -- %s\n", "Concentrations", concs);
  /*==============*/
  /* Added by ALI */
  /*==============*/
  mpi_printf( 0, "%s -- %s\n", "DDIs File", ddis_file );
  mpi_printf( 0, "%s -- %s\n", "Cvar File", cvar_file );
  mpi_printf( 0, "%s -- %s\n", "Hill File a", hill_file_a );
  mpi_printf( 0, "%s -- %s\n", "Hill File b", hill_file_b );
  mpi_printf( 0, "%s -- %s\n", "Drug_Name_A", drug_name_a);
  mpi_printf( 0, "%s -- %s\n", "Drug_Name_A", drug_name_b);
  mpi_printf( 0, "%s -- %lf\n", "Cmax A", cmax_a);
  mpi_printf( 0, "%s -- %lf\n", "Cmax B", cmax_b);
  mpi_printf( 0, "%s -- %s\n", "Radius", radius);
  mpi_printf( 0, "%s -- %s\n", "Angle", angles);
  mpi_printf( 0, "%s -- %d\n", "Is_DDIs", is_ddis_la );
  mpi_printf( 0, "%s -- %d\n", "Is_Syntopic", is_syntopic );
  mpi_printf( 0, "%s -- %d\n", "Is_Cvar", is_cvar );
}
