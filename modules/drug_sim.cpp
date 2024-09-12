#include "drug_sim.hpp"

#ifdef TOMEK_2019
#include "../cellmodels/Tomek_model.hpp" 
#else
#include "../cellmodels/Ohara_Rudy_2011.hpp"
#endif
#include "commons.hpp"

#include <cmath>


void do_drug_sim(const double conc_a, const double conc_b, 
std::array<double, 14> ic50_a, std::array<double, 14> ic50_b, 
std::array<double,6> ddis, std::array<double,18> cvar,
const param_t *p_param, const unsigned short sample_id, 
Cellmodel *p_cell, cvode_t *p_cvode, bool is_firsttime)
{
  // buffer for writing in snprintf() function
  char buffer[255];
  // for calculating qinward
  double inal_auc_control, ical_auc_control, inal_auc_drug, ical_auc_drug;
  // CVode variables
  double tnext, tcurr;
  int cvode_retval;
  unsigned int icount, imax;

  // files for storing results
  // time-series result
  FILE *fp_vm, *fp_ca, *fp_dvmdt, *fp_ires, *fp_inet;
  // features result
  FILE *fp_ap_profile, *fp_ca_profile, *fp_qni;

  // simulation parameters
  double dt = p_param->dt;
  double dtw = p_param->dt_write;
  const char *drug_name_a = p_param->drug_name_a;
  const char *drug_name_b = p_param->drug_name_b;
  const double bcl = p_param->bcl;
  const double inet_vm_threshold = p_param->inet_vm_threshold;
  const unsigned short pace_max = p_param->pace_max;
  const unsigned short celltype = (unsigned short)p_param->celltype;
  const unsigned short last_pace_print = 3;
  const unsigned short last_drug_check_pace = 250;
  const unsigned int print_freq = (1./dt) * dtw;
  unsigned short pace_count = 0;
  unsigned short pace_steepest = 0;

  // drug features
  double inet,qnet;
  double inal_auc, ical_auc;
  double vm_repol30, vm_repol50, vm_repol90;
  double t_depol;
  double t_ca_peak, ca_amp50, ca_amp90;
  double cad50_prev, cad50_curr, cad90_prev, cad90_curr;

  // variables to store features
  // temp_result is the result of features in 1 pace,
  // will be interchanged during the simulation.
  // cipa_result is the final result of the simulation.
  cipa_t cipa_result, temp_result;

  // TRUE if vm peak more than 0 mV
  bool is_eligible_AP;
   
  /*==============*/
  /* Added by ALI */
  /*==============*/
  // We assume that DDIs effect only applied to hERG channel 
  ddis_t ddis_la;
  dcomb_t* dcomb_param = new dcomb_t(); 
  double ddis_effect;
  if(p_param->is_ddis_la == true){
    dcomb_param->init(conc_a,conc_b, // drug concentrations in nM
	  1.0,1.0, // emax
	  ic50_a[12],ic50_b[12], // ic50
	  ic50_a[13],ic50_b[13], // Hill's coeffs
	  ddis[0],ddis[1], // int_ab and int_ba
	  pow(10.0,ddis[2]+3),pow(10.0,ddis[3]+3), // ec50_int_ab and ec50_int_ba in nM
	  ddis[4],ddis[5], // n_int_ab and n_int_ba
	  true,false // is_la and is_ec50_bi
	);
    ddis_effect = ecomb(dcomb_param);
  }
  // apply some cell initialization
  p_cell->initConsts( celltype, conc_a, conc_b, 
 	ic50_a.data(), ic50_b.data(), ddis_effect, cvar.data(), 
	p_param->is_dutta, p_param->is_ddis_la, 
	p_param->is_syntopic, p_param->is_cvar);
  p_cell->CONSTANTS[BCL] = bcl;
  
  FILE *fp_states;
  if( p_param->is_using_output > 0 ){
#ifdef TOMEK_2019
   fp_states = fopen("output_tomek.dat", "r");
#else
   if( p_param->is_dutta > 0){
     fp_states = fopen("output_orudy_dutta.dat", "r");
   }
   else {
     fp_states = fopen("output_orudy.dat", "r");
   }
#endif
   if( fp_states != NULL ){
     mpi_printf(0, "Using initial condition from steady state!\n");
     int idx = 0;
     while(fgets( buffer, sizeof(buffer), fp_states) != NULL){
       p_cell->STATES[idx++] = strtod(buffer, NULL);
     }
   }
   else{
     mpi_printf(0, "No initial file found! Skipped using initial file!\n");
   }
  }
  
  tcurr = 0.;
  tnext = dt;
  
  init_cvode(p_cvode, p_cell, tcurr, is_firsttime);

  // generate file for time-series output
  snprintf(buffer, sizeof(buffer), 
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_vmcheck_smp%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, 
	    sample_id );
  fp_vm = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer), 
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_dvmdt_smp%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, 
	    sample_id );
  fp_dvmdt = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer),  
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_ca_i_smp%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, 
	    sample_id);
  fp_ca = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer),  
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_ires_smp%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, 
	    sample_id);
  fp_ires = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer),  
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_inet_smp%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, 
	    sample_id);
  fp_inet = fopen( buffer, "w" );
  snprintf(buffer, sizeof(buffer),  
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_qnet_proc%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, mympi::rank );
  fp_qni = fopen( buffer, "a" );
  snprintf(buffer, sizeof(buffer), 
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_ap_profile_proc%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, mympi::rank );
  fp_ap_profile = fopen( buffer, "a" );
  snprintf(buffer, sizeof(buffer),  
	    "result/%.2lf_%.2lf/%s_%.2lf_%s_%.2lf_ca_profile_proc%d.plt", 
            conc_a, conc_b, drug_name_a, conc_a, drug_name_b, conc_b, mympi::rank );
  fp_ca_profile = fopen( buffer, "a" );

  fprintf(fp_vm, "%s %s\n", "Time", "Vm");
  fprintf(fp_dvmdt, "%s %s\n", "Time", "dVm/dt");
  fprintf(fp_ca, "%s %s\n", "Time", "cai");
  fprintf(fp_ires, "%s %s %s %s %s %s %s %s\n", 
              "Time", "INa", "INaL", "ICaL", 
              "Ito", "IKr", "IKs", "IK1");
  fprintf( fp_ap_profile, "%s %s %s %s %s %s %s %s %s\n",
              "Sample_ID", "Dvm/Dt_Repol", "Max_Dvm/Dt", "Vm_Peak", 
              "Vm_Resting","APD90", "APD50", "APDTri", "Steepest_Pace");
  fprintf( fp_ca_profile, "%s %s %s %s %s %s %s\n",
                 "Sample_ID", "Ca_Peak", "Ca_Diastole", "CaD90", 
                 "CaD50","Catri", "Steepest_Pace");
  fprintf( fp_qni, "%s %s %s\n","Sample_ID", "Qnet", "Qinward");
  
  icount = 0;
  imax = (unsigned int)((pace_max * bcl)/dt);
  inet = 0.;
  qnet = 0.;
  is_eligible_AP = 0.;

  cipa_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
  temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );

  while(icount < imax)
  // begin simulation loop
  {
    cvode_retval = CVode(p_cvode->cvode_mem, tnext, p_cvode->states_vec, &tcurr, CV_NORMAL);
    if( cvode_retval == CV_SUCCESS ){
      icount++;
      tnext += dt;
    }
    else{
      mpi_fprintf(0, stderr, "CVode error at sample_ID %d and concentration_A %.2lf concentration_B %.2lf at rank %d\n", 
      sample_id, conc_a, conc_b, mympi::rank);
      break;
    }
    // get the RATES array
    p_cell->computeRates(tnext, p_cell->CONSTANTS, p_cell->RATES, N_VGetArrayPointer_Serial(p_cvode->states_vec), p_cell->ALGEBRAIC);
    // calculate inet and AUC
    if(p_cell->STATES[V] > inet_vm_threshold){
      // According to CiPA protocol (Li. et al 2018), the qNet is only affected by 4 ion channels (Na, NaL, CaL, and Kr)
      //inet += (p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[Ito]+p_cell->ALGEBRAIC[IKr]+p_cell->ALGEBRAIC[IKs]+p_cell->ALGEBRAIC[IK1])*dt;
      inet += (p_cell->ALGEBRAIC[INa]+p_cell->ALGEBRAIC[INaL]+p_cell->ALGEBRAIC[ICaL]+p_cell->ALGEBRAIC[IKr])*dt;
      inal_auc += p_cell->ALGEBRAIC[INaL]*dt;
      ical_auc += p_cell->ALGEBRAIC[ICaL]*dt;
    } 

    // save temporary result
    if(pace_count >= pace_max-last_drug_check_pace && icount % print_freq == 0){
      temp_result.cai_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[cai]) );
      temp_result.vm_data.insert( std::pair<double, double> (tcurr, p_cell->STATES[V]) );
      temp_result.dvmdt_data.insert( std::pair<double, double> (tcurr, p_cell->RATES[V]) );
      snprintf( buffer, sizeof(buffer), "%lf %lf %lf %lf %lf %lf %lf", 
              p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL], p_cell->ALGEBRAIC[ICaL], p_cell->ALGEBRAIC[Ito], 
              p_cell->ALGEBRAIC[IKr], p_cell->ALGEBRAIC[IKs], p_cell->ALGEBRAIC[IK1] );
      temp_result.ires_data.insert( std::pair<double, string> (tcurr, string(buffer)) );
    }

    // execute at the beginning of new pace
    if(icount % (int)(bcl/dt) == 0 && icount > 0){

      // assuming the pace is eligible,
      // we will extract result
      if( is_eligible_AP && pace_count >= pace_max-last_drug_check_pace-1) {
        for(std::multimap<double, double>::iterator itrmap = temp_result.cai_data.begin(); 
            itrmap != temp_result.cai_data.end() ; itrmap++ ){
          // before the peak calcium
          if( itrmap->first < t_ca_peak ){
            if( itrmap->second < ca_amp50 ) cad50_prev = itrmap->first;
            if( itrmap->second < ca_amp90 ) cad90_prev = itrmap->first;
          }
          // after the peak calcium
          else{
            if( itrmap->second > ca_amp50 ) cad50_curr = itrmap->first;
            if( itrmap->second > ca_amp90 ) cad90_curr = itrmap->first;
          }
        }
        temp_result.cad50 = cad50_curr - cad50_prev;
        temp_result.cad90 = cad90_curr - cad90_prev;
        temp_result.qnet = inet/1000.0;
        temp_result.inal_auc = inal_auc;
        temp_result.ical_auc = ical_auc;
        temp_result.vm_dia = p_cell->STATES[V];
        temp_result.ca_dia = p_cell->STATES[cai];

        // replace result with steeper repolarization AP or first pace from the last 250 paces
        if( temp_result.dvmdt_repol > cipa_result.dvmdt_repol ) {
          pace_steepest = pace_count;
          cipa_result = temp_result;
        }
      }// end pace eligible
      // resetting inet and AUC values
      // and increase the pace count
      pace_count++;
      inet = 0.;
      inal_auc = 0.;
      ical_auc = 0.;
      if(pace_count >= pace_max-last_drug_check_pace){
        temp_result.init( p_cell->STATES[V], p_cell->STATES[cai] );
        t_ca_peak = tcurr;
        t_depol = (p_cell->CONSTANTS[BCL]*pace_count)+p_cell->CONSTANTS[stim_start];
        is_eligible_AP = false;
      }
    }// end beginning pace process
    // entering the last 250 paces
    if(pace_count >= pace_max-last_drug_check_pace){
    
      // get maximum dvmdt
      if( temp_result.dvmdt_max < p_cell->RATES[V] ) temp_result.dvmdt_max = p_cell->RATES[V];

      // get the peak Vm 6 secs after depolarization (when Na channel just closed after bursting)
      if( icount % (int)(((p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+6.)) / dt) == 0 ){
        temp_result.vm_peak = p_cell->STATES[V];
        if( temp_result.vm_peak > 0. ){
          vm_repol30 = temp_result.vm_peak - (0.3 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol50 = temp_result.vm_peak - (0.5 * (temp_result.vm_peak - temp_result.vm_valley));
          vm_repol90 = temp_result.vm_peak - (0.9 * (temp_result.vm_peak - temp_result.vm_valley));
          is_eligible_AP = true;
        }
        else is_eligible_AP = false;
      }
      // these operations will be executed if it's eligible AP and executed at the beginning of repolarization
      else if( is_eligible_AP && tcurr > (p_cell->CONSTANTS[BCL]*pace_count)+(p_cell->CONSTANTS[stim_start]+6.) ){
        // get the steepest dvmdt_repol around vm_repol30 and vm_repol90 if AP shape is eligible
        if( vm_repol30 >= p_cell->STATES[V] && p_cell->STATES[V] >= vm_repol90 && temp_result.dvmdt_repol < p_cell->RATES[V] ){
          temp_result.dvmdt_repol = p_cell->RATES[V];
        }
        // get the APD90, APD50, peak calcium, 50% and 90% of amplitude of Calcium, and time of peak calcium
        if( vm_repol50 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol50-2 ) temp_result.apd50 = tcurr - t_depol;
        if( vm_repol90 > p_cell->STATES[V] && p_cell->STATES[V] > vm_repol90-2 ) temp_result.apd90 = tcurr - t_depol;
        if( temp_result.ca_peak < p_cell->STATES[cai] ){  
          temp_result.ca_peak = p_cell->STATES[cai];
          ca_amp50 = temp_result.ca_peak - (0.5 * (temp_result.ca_peak - temp_result.ca_valley));
          ca_amp90 = temp_result.ca_peak - (0.9 * (temp_result.ca_peak - temp_result.ca_valley));
          t_ca_peak = tcurr;
        } 
      }   
    }// end of last 250 paces operations

  }// end simulation loop

  // print graph result and features
  // if CVode success
  if( cvode_retval == CV_SUCCESS ){
    //if(p_param->is_print_graph == 1 ) {
      for(std::multimap<double, double>::iterator itrmap = cipa_result.cai_data.begin(); itrmap != cipa_result.cai_data.end() ; itrmap++ ){
        fprintf(fp_ca, "%lf %lf\n", itrmap->first, itrmap->second); 
      }
      for(std::multimap<double, double>::iterator itrmap = cipa_result.vm_data.begin(); itrmap != cipa_result.vm_data.end() ; itrmap++ ){
        fprintf(fp_vm, "%lf %lf\n", itrmap->first, itrmap->second); 
      }
      for(std::multimap<double, double>::iterator itrmap = cipa_result.dvmdt_data.begin(); itrmap != cipa_result.dvmdt_data.end() ; itrmap++ ){
        fprintf(fp_dvmdt, "%lf %lf\n", itrmap->first, itrmap->second); 
      }
      for(std::multimap<double, string>::iterator itrmap = cipa_result.ires_data.begin(); itrmap != cipa_result.ires_data.end() ; itrmap++ ){
        fprintf(fp_ires, "%lf %s\n", itrmap->first, (itrmap->second).c_str()); 
      }
    //}
    fprintf( fp_ap_profile, "%d %lf %lf %lf %lf %lf %lf %lf %d\n",
    sample_id, cipa_result.dvmdt_repol, cipa_result.dvmdt_max, cipa_result.vm_peak,
    cipa_result.vm_dia, cipa_result.apd90, cipa_result.apd50, cipa_result.apd90-cipa_result.apd50,
             pace_steepest);
    fprintf( fp_ca_profile, "%d %lf %lf %lf %lf %lf %d\n",
    sample_id, cipa_result.ca_peak, cipa_result.ca_dia,
    cipa_result.cad90, cipa_result.cad50, cipa_result.cad90-cipa_result.cad50, pace_steepest );
    
    if( (int)ceil(conc_a) == 0 && (int)ceil(conc_b) == 0) {
      inal_auc_control = cipa_result.inal_auc;
      ical_auc_control = cipa_result.ical_auc;
      fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0, 0.0 );
    }
    else{
      inal_auc_drug = cipa_result.inal_auc;
      ical_auc_drug = cipa_result.ical_auc;
      fprintf( fp_qni, "%d %lf %lf\n", sample_id, (cipa_result.qnet)/1000.0,
               ( (inal_auc_drug/inal_auc_control) + (ical_auc_drug/ical_auc_control) ) * 0.5 );
    }
  }
  else{
    fprintf( fp_ap_profile, "%d ERR ERR ERR ERR ERR ERR\n", sample_id);
    fprintf( fp_ca_profile, "%d ERR ERR ERR ERR ERR\n", sample_id);
    fprintf( fp_qni, "%d ERR ERR\n", sample_id);
  }


  // clean the memories
  fclose(fp_states);
  fclose(fp_qni);
  fclose(fp_inet);
  fclose(fp_ires);
  fclose(fp_ca);
  fclose(fp_dvmdt);
  fclose(fp_vm);
  delete dcomb_param;
}
