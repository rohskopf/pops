#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>       /* round, floor, ceil, trunc */

#include "popsinput.h"
#include "popsmemory.h"
#include "objfun.h"

using namespace LAMMPS_NS;

using namespace POPS_NS;

Objfun::Objfun(POPS *pops) : PopsPointers(pops) {}

Objfun::~Objfun() {};

vector<double> Objfun::calcfitness(double sub_xpop_r[], int nepp)
{

  int pop_size = pops->popsinput->pop_size_in;
  int unknowns_tot = pops->popsinput->unknowns_in + pops->popsinput->unknowns_extra_in;

  // Check sub_xpop_r
  /*for (int i=0; i < nepp; i++)
  {
    debug << sub_xpop_r[i] << endl;
  }*/

  // Allocate memory for params
  pops->popsmemory->allocate(params, nepp);
  
  // Prepare error quantities
  vector<double> sub_Z_vec;
  sub_Z_vec.resize(pop_size);
  sub_Zf2_vec.resize(pop_size);
  int i_indx = 0;
  for (int i=0; i < nepp/unknowns_tot; i++)
  {
    for (int u=0; u < unknowns_tot; u++)
    {
      //debug << u << endl;
      params[u] = sub_xpop_r[u + i_indx];
      //debug << "param" << u << ": " << params2[u] << endl;
    }
    /************ Calculate Z value here *************/
    sub_Z_vec[i] = calcz(params, i);
    //debug << sub_Z_vec[i] << endl;
    i_indx = i_indx + unknowns_tot;
    //debug << i << endl;
  }

  //debug << nepp << endl;
  pops->popsmemory->deallocate(params);

  return sub_Z_vec;
 
} // void calcfitness()

double Objfun::calcz(double params[], int ind)
{

  /* Get necessary variables for this function */
  string rest = pops->popsinput->commands;
  string nontab_commands = pops->popsinput->nontab_commands;
  vector<double> d_vec = pops->popsinput->d_vec_in;
  int unknowns = pops->popsinput->unknowns_in;

  size_t found;
  size_t n;
  string substring;
  
  int sensparam_counter = 0;
  string line;
  string l_whole;
  double changedparam;
  // Put all the parameters from an individual into the file ------------------------------------------------------------------------------
  int current_line = 1;
  //cout << sensparam_counter << endl;

  //  Get the individual process ID.
  int rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  char refpot[64];
  char temppot[64];
  char objfundebug[64];
  sprintf (refpot, "REFPOT%d", rank);
  sprintf (temppot, "TEMPPOT%d", rank);
  sprintf (objfundebug, "OFDEBUG%d", rank);
  ifstream fh_final_r("REFPOT");
  ofstream fh_temp_w(temppot);
  ofstream ofdebug(objfundebug);

  //ofdebug << "This is process " << rank << endl;

  int sym_vec_counter = 0; 
  int sym_coeffs_indx = 0;
  int searchindx = 0;
  bool is_sym = false;
  string sym_nums;
  //string sym_params;
  vector<double> sym_params;
  int param_counter = 0;
  string tab_symparam_nums; // This is a string "123456..." that labels the symmetric parameters that exist
  vector<string> symparams_str;

  while (!fh_final_r.eof())
  { 
    int searchindx = 0;
    getline(fh_final_r, l_whole);
    string l_temp = l_whole;
    string str1 = l_temp;
    string str2 = "[";
    //debug << searchindx << endl;
    //debug << l_temp << endl;
    // Check to see if there are any left brackets:
    //found = str1.find(str2, searchindx);
    found = l_temp.find(str2, searchindx);
    //debug << "found: " << found << endl;
    searchindx = found+1;
    /*debug << "found: " << found << endl;
    //debug << "SI: " << searchindx << endl;*/
    // Loop through the line if there are left brackets:
    while (found != string::npos)
    {
      if (found != string::npos)
      {
        // Extract inside of bracket
        size_t right_brack_found = l_temp.find("]");
        string inside = l_temp.substr(found+1, right_brack_found - found-1);
        //debug << "inside: " << inside << endl;
        // If inside is empty, proceed like normal
        if (inside.empty() )
        {
          // Find length of inside
          const char* inside_cc = inside.c_str();
          size_t inside_length = strlen(inside_cc);
          //debug << "inside_length: " << inside_length << endl;
          // Determine length of entire portion to be removed
          int tot_length = inside_length + 2;
          //debug << "tot_length: " << tot_length << endl;
          // Replace the original substring with a parameter in place of the bracketed term
          l_temp.replace(found, tot_length, " "); // Delete the bracketed term
          // Extract the parameter
          //debug << "COUNTER: " << nontab_param_counter << endl;
          double param = params[param_counter];
          //debug << "PARAM: " << param << endl;
          // Round to proper decimal accuracy
          param = round(param*pow(10,d_vec[param_counter]))/(pow(10,d_vec[param_counter]));
          // Convert to string
          string param_str = to_string(param);
          // Insert parameter into the command
          l_temp.insert(found, param_str);
          // Increment the param_counter
          param_counter++;
          // Find the next changeable parameter in the line (if any)
          found = l_temp.find("[");
        } // if (inside.empty() )
        else // we have a symmetric parameter
        {

          //debug << "inside: " << inside << endl;

          //----- Examine the contents
          string inside_examine = inside;
          // Now examine the symmetric parameters
          string sym_str = inside_examine;
          //debug << "sym_str: " << sym_str << endl;
          // Check if a member of symparam_nums
          std::size_t sym_found = tab_symparam_nums.find(sym_str);
          //debug << "sym_found: " << sym_found << endl;
          // If the parameter is not a part of symmetric params, add it!
          if (sym_found==std::string::npos)
          {
            //debug << "NEW SYMMETRIC PARAM" << endl;
            tab_symparam_nums.append(sym_str);
            tab_symparam_nums.append(",");
            //debug << "symparam_nums: " << tab_symparam_nums << endl;
            //debug << symparams[0] << endl;
            // Find length of inside
            const char* inside_cc = inside.c_str();
            size_t inside_length = strlen(inside_cc);
            //debug << "inside_length: " << inside_length << endl;
            // Determine length of entire portion to be removed
            int tot_length = inside_length + 2;
            //debug << "tot_length: " << tot_length << endl;
            // Replace the original substring with a parameter in place of the bracketed term
            l_temp.replace(found, tot_length, " "); // Delete the bracketed term
            // Extract the parameter
            //debug << "COUNTER: " << param_counter << endl;
            double param = params[param_counter];
            //debug << "PARAM: " << param << endl;
            double newparam = param;
            // Round to proper decimal accuracy
            param = round(param*pow(10,d_vec[param_counter]))/(pow(10,d_vec[param_counter]));
            // Convert to string
            string param_str = to_string(param);
            symparams_str.push_back(param_str);
            //debug << symparams_str[0] << endl;
            // Insert parameter into the command
            l_temp.insert(found, param_str);
            // Increment the param_counter
            param_counter++;
            // Find the next changeable parameter in the line (if any)
            found = l_temp.find("[");
            // if (sym_found==std::string::npos)
          } // if (sym_found == std::string::npos)
          else // we have a symmetric parameter that has already been defined 
          {
            // Find length of inside
            const char* inside_cc = inside.c_str();
            size_t inside_length = strlen(inside_cc);
            //debug << "inside_length: " << inside_length << endl;
            // Determine length of entire portion to be removed
            int tot_length = inside_length + 2;
            //debug << "tot_length: " << tot_length << endl;
            // Replace the original substring with a parameter in place of the bracketed term
            l_temp.replace(found, tot_length, " "); // Delete the bracketed term
            // Extract the parameter
            //debug << "COUNTER: " << nontab_param_counter << endl;
            // Find which symmetric parameter to use...
            int sym_indx = atoi(inside_examine.c_str() );
            //debug << "PARAM: " << param << endl;
            // Convert to string
            string param_str = symparams_str[sym_indx - 1];
            // Insert parameter into the line
            l_temp.insert(found, param_str);
            // DON'T increment the param_counter for this symmetric parameter
            //param_counter++;
            // Find the next changeable parameter in the line (if any)
            found = l_temp.find("[");
          } // else

        } // else we have a symmetric parameter

       } // if found != nopos

    } // while (found != string::npos)
    //debug << l_temp << endl;
    fh_temp_w << l_temp << endl;

    //param_counter = 0;

    current_line ++;

  } // while (fh_final_r.eof())

  current_line = 1;
  //sensparam_counter ++;
  fh_final_r.close();
  fh_temp_w.close();

  // Loop through the pairstyle commands and input them
  int nontab_param_counter = 0;
  //size_t param_found;
  n = nontab_commands.find(',');
  rest = nontab_commands;
  bool isnegative = false;
  string symparam_nums;
  //vector<double> symparams;
  symparams_str.clear();
  string setcommands;
  int num_setcommands = 0;
  while (n != std::string::npos)
  {
    substring = rest.substr(0,n);
    //debug << "MAIN SUBSTRING: " << substring << endl;
    // Loop through the substring and find the occurences of "PARAM" and replace them with params[i] 
    size_t param_found = substring.find("[");
    size_t left_brace = substring.find("{");
    //debug << "Left brace index: " << left_brace << endl;
    if (left_brace != string::npos)
    {
      //debug << "found a rank!" << endl;
      size_t right_brace = substring.find("}");
      //debug << "Right brace index: " << right_brace << endl;
      string inside = substring.substr(left_brace+1, right_brace - left_brace-1);
      //debug << "inside string: " << inside << endl;
      if (inside.empty() ) // then continue like normal
      {
        // Find length of inside
        const char* inside_cc = inside.c_str();
        size_t inside_length = strlen(inside_cc);
        //debug << "inside_length: " << inside_length << endl;
        // Determine length of entire portion to be removed
        int tot_length = inside_length + 2;
        //debug << "tot_length: " << tot_length << endl;
        // Replace the original substring with a parameter in place of the bracketed term
        substring.replace(left_brace, tot_length, " "); // Delete the bracketed term
        //debug << "Substring after brace deletion: " << substring << endl;
        // Now append a rank integer to the end
        //debug << "rank: " << rank << endl;
        string rank_str = to_string(rank);
        // Insert rank into the command
        substring.insert(left_brace, rank_str);
        //debug << "Final substring: " << substring << endl;
        // Extract the parameter
        //debug << "COUNTER: " << nontab_param_counter << endl;
        /*double param = params[unknowns + nontab_param_counter];
        //debug << "PARAM: " << param << endl;
        // Round to proper decimal accuracy
        param = round(param*pow(10,d_vec[unknowns+nontab_param_counter]))/(pow(10,d_vec[unknowns+nontab_param_counter]));
        // Convert to string
        string param_str = to_string(param);
        // Insert parameter into the command
        substring.insert(param_found, param_str);
        // Increment the nontab_param_counter
        nontab_param_counter++;
        // Find the next changeable parameter in the line (if any)
        param_found = substring.find("[");*/
      } // if (inside.empty() )
      
    }
    while (param_found != std::string::npos)
    {
      // Extract inside of bracket
      size_t right_brack_found = substring.find("]");
      string inside = substring.substr(param_found+1, right_brack_found - param_found-1);
      /*debug << "substring: " << substring << endl;
      debug << param_found << endl;
      debug << right_brack_found << endl;
      debug << right_brack_found - param_found -1<< endl;*/
      //debug << "inside: " << inside << endl;
      
      //---- Examine contents of "inside" to decide what to do
      // Check if there is no inside -> proceed like normal
      if (inside.empty() ) // then continue like normal
      {
        // Find length of inside
        const char* inside_cc = inside.c_str();
        size_t inside_length = strlen(inside_cc);
        //debug << "inside_length: " << inside_length << endl;
        // Determine length of entire portion to be removed
        int tot_length = inside_length + 2;
        //debug << "tot_length: " << tot_length << endl;
        // Replace the original substring with a parameter in place of the bracketed term
        substring.replace(param_found, tot_length, " "); // Delete the bracketed term
        // Extract the parameter
        //debug << "COUNTER: " << nontab_param_counter << endl;
        double param = params[unknowns + nontab_param_counter];
        //debug << "PARAM: " << param << endl;
        // Round to proper decimal accuracy
        param = round(param*pow(10,d_vec[unknowns+nontab_param_counter]))/(pow(10,d_vec[unknowns+nontab_param_counter]));
        // Convert to string
        string param_str = to_string(param);
        // Insert parameter into the command
        substring.insert(param_found, param_str);
        // Increment the nontab_param_counter
        nontab_param_counter++;
        // Find the next changeable parameter in the line (if any)
        param_found = substring.find("[");
      } // if (inside.empty() )
      // If not empty, we have some modifications to make...
      else
      {

        //debug << "inside: " << inside << endl;

        //----- Examine the contents
        string inside_examine = inside;
        // Check for negative signs
        std::size_t inside_found = inside.find("-");
        //debug << inside_found << endl;
        if (inside_found!=std::string::npos)
        {
          //debug << "FACE" << endl;
          // IS_NEGATIVE = TRUE
          isnegative = true;
          // Delete the negative
          inside_examine.erase(inside_found, 1);
          //debug << "inside after delete: " << inside_examine << endl;
        }
        // Now examine the symmetric parameters
        string sym_str = inside_examine;
        //debug << "sym_str: " << sym_str << endl;
        // Check if a member of symparam_nums
        std::size_t sym_found = symparam_nums.find(sym_str);
        //debug << "sym_found: " << sym_found << endl;
        // If the parameter is not a part of symmetric params, add it!
        if (sym_found==std::string::npos)
        {
          //debug << "NEW SYMMETRIC PARAM" << endl;
          symparam_nums.append(sym_str);
          //debug << "symparam_nums: " << symparam_nums << endl;
          //debug << symparams[0] << endl;
          // Find length of inside
          const char* inside_cc = inside.c_str();
          size_t inside_length = strlen(inside_cc);
          //debug << "inside_length: " << inside_length << endl;
          // Determine length of entire portion to be removed
          int tot_length = inside_length + 2;
          //debug << "tot_length: " << tot_length << endl;
          // Replace the original substring with a parameter in place of the bracketed term
          substring.replace(param_found, tot_length, " "); // Delete the bracketed term
          // Extract the parameter
          //debug << "COUNTER: " << nontab_param_counter << endl;
          double param = params[unknowns + nontab_param_counter];
          //debug << "PARAM: " << param << endl;
          double newparam = param;
          // Round to proper decimal accuracy
          param = round(param*pow(10,d_vec[unknowns+nontab_param_counter]))/(pow(10,d_vec[unknowns+nontab_param_counter]));
          // Convert to string
          string param_str = to_string(param);
          symparams_str.push_back(param_str);
          //debug << symparams_str[0] << endl;
          // Insert parameter into the command
          substring.insert(param_found, param_str);
          // Increment the nontab_param_counter
          nontab_param_counter++;
          // Find the next changeable parameter in the line (if any)
          param_found = substring.find("[");
        } // if (sym_found==std::string::npos)
        else // we have a symmetric parameter that already exists
        {
          // Find length of inside
          const char* inside_cc = inside.c_str();
          size_t inside_length = strlen(inside_cc);
          //debug << "inside_length: " << inside_length << endl;
          // Determine length of entire portion to be removed
          int tot_length = inside_length + 2;
          //debug << "tot_length: " << tot_length << endl;
          // Replace the original substring with a parameter in place of the bracketed term
          substring.replace(param_found, tot_length, " "); // Delete the bracketed term
          // Extract the parameter
          //debug << "COUNTER: " << nontab_param_counter << endl;
          // Find which symmetric parameter to use...
          int sym_indx = atoi(inside_examine.c_str() );
          /*// Check if negative
          if (isnegative)
          {
            param = (-1)*param;
            isnegative = false;
          }*/
          //debug << "PARAM: " << param << endl;
          // Convert to string
          string param_str = symparams_str[sym_indx - 1];
          // Check if negative
          if (isnegative)
          {
            param_str.insert(0,"-");
            isnegative = false;
          }
          // Insert parameter into the command
          substring.insert(param_found, param_str);
          // DON'T increment the nontab_param_counter for this symmetric parameter
          //nontab_param_counter++;
          // Find the next changeable parameter in the line (if any)
          param_found = substring.find("[");
        } // else


      }

    }
    //debug << substring << endl;
    // Convert strings to const char * types so that lammps can read them
    size_t setfound = substring.find("set");
    const char* substring_cc = substring.c_str();
    //cout << "LAMMPS INPUT: " << substring_cc << endl;
    if (setfound != string::npos)
    {
      setcommands.append(substring);
      setcommands.append(",");
      num_setcommands++;
    }
    else
    {
    // Input the command
    lmp->input->one(substring_cc);
    }
    //cout << substring_cc << endl;
    // Find location of next command and repeat
    rest = rest.substr(n+1);
    n = rest.find(',');
  }

  // Get vectors and quantities necessary for the property prediction
  double *x_all = pops->popsinput->x_all;
  int *N_vec = pops->popsinput->N_vec;
  char *atomtypes = pops->popsinput->atomtypes;
  double *fio = pops->popsinput->fio;
  double *eio = pops->popsinput->eio;
  double *eio_s = pops->popsinput->eio_s;
  double *sio = pops->popsinput->sio;
  double *pio = pops->popsinput->pio;
  double *pio_s = pops->popsinput->pio_s;

  // Calculate MD quantities and compare to reference quantities
  double sum_f_diff=0;
  double sum_e_diff=0;
  double sum_e_diff_s=0;
  double sum_s_diff=0;
  double sum_p_diff=0;
  double sum_p_diff_s=0;
  double Z_2_f = 0; // % error in forces
  double Z_2_e = 0; // % error in energies
  // Denominators for normalization
  double sum_fio = 0;
  double sum_eio = 0;
  double sum_eio_s = 0;
  double sum_sio = 0;
  double sum_pio = 0;
  double sum_pio_s = 0;
  int M_indx=0;
  int M_indx_e=0;
  int M_indx_s = 0;
  int t_indx = 0;
  int Natoms;
  double E0_md;
  double P0_md;

  pops->popsmemory->allocate(fi, pops->popsinput->M_in*pops->popsinput->N_tot*3);
  pops->popsmemory->allocate(ei, pops->popsinput->M_in);
  pops->popsmemory->allocate(ei_s, pops->popsinput->M_in);
  pops->popsmemory->allocate(pi, pops->popsinput->M_in);
  pops->popsmemory->allocate(pi_s, pops->popsinput->M_in);

  for (int i=0; i<pops->popsinput->M_in; i++)
  {  

    // Delete atoms from previous configuration
    lmp->input->one("delete_atoms group all");
    /*// Set number of atoms
    Natoms = N_vec[i+1];
    string Natoms_str = to_string(Natoms);
    //debug << Natoms << endl;
    string create_atoms_str = "create_atoms 1 random NATOMS 69 NULL";
    size_t NATOMS_found = create_atoms_str.find("NATOMS");
    create_atoms_str.replace(NATOMS_found, 6, " "); // Delete "NATOMS"
    create_atoms_str.insert(NATOMS_found, Natoms_str);
    const char* create_atoms_str_cc = create_atoms_str.c_str();
    debug << create_atoms_str << endl;
    lmp->input->one("delete_atoms region box");
    //debug << create_atoms_str_cc << endl;
    lmp->input->one(create_atoms_str_cc);*/

    // Set number of atoms
    Natoms = N_vec[i];
    string Natoms_str = to_string(Natoms);

    pops->popsmemory->allocate(x, Natoms*3);
    // Get the x array
    for (int a=0; a<Natoms*3; a++)
    {
      x[a] = x_all[a+M_indx];
      //cout << x[a] << endl;
    }

    // Extract the cube length of the box from CONFIGS
    //double ao = aio[i];
    //double caoo = coaio[i];
    double x_length = pops->popsinput->x_lengths[i];
    double y_length = pops->popsinput->y_lengths[i];
    double z_length = pops->popsinput->z_lengths[i];
    double xy_box = pops->popsinput->xy_vec[i];
    double xz_box = pops->popsinput->xz_vec[i];
    double yz_box = pops->popsinput->yz_vec[i];
    string x_str = to_string(x_length);
    string y_str = to_string(y_length);
    string z_str = to_string(z_length);
    string xy_str = to_string(xy_box);
    string xz_str = to_string(xz_box);
    string yz_str = to_string(yz_box);
    //debug << ao << endl;

    // Convert box length to string
    //string ao_str = to_string(ao);
    //string caoo_str = to_string(caoo);
    string box_size_str = "change_box all x final 0 LENGTH y final 0 LENGTH z final 0 LENGTH xy final LENGTH xz final LENGTH yz final LENGTH units box";
    // Replace "LENGTH" with the ao_str for all three coordinates...
    // Replace x coordinate
    size_t LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, x_str);
    // Replace y coordinate
    LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, y_str);
    // Replace z coordinate
    LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, z_str);
    // Replace xy coordinate
    LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, xy_str);
    // Replace xz coordinate
    LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, xz_str);
    // Replace yz coordinate
    LENGTH_found = box_size_str.find("LENGTH");
    box_size_str.replace(LENGTH_found, 6, " "); // Delete "LENGTH"
    box_size_str.insert(LENGTH_found, yz_str);
    //debug << box_size_str << endl;
    // Convert string to const char* so that lammps can read it
    const char* box_size_str_cc = box_size_str.c_str();

    // Change the box size
    lmp->input->one(box_size_str_cc);

    // Create atoms (just to declare memory - positions don't matter yet)
    string create_atoms_str = "create_atoms 1 random NATOMS 69 NULL";
    size_t NATOMS_found = create_atoms_str.find("NATOMS");
    create_atoms_str.replace(NATOMS_found, 6, " "); // Delete "NATOMS"
    create_atoms_str.insert(NATOMS_found, Natoms_str);
    const char* create_atoms_str_cc = create_atoms_str.c_str();
    //debug << create_atoms_str << endl;
    //lmp->input->one("delete_atoms region box");
    //cout << create_atoms_str_cc << endl;
    lmp->input->one(create_atoms_str_cc);

    //vector<int> type_vec = {1,1,1,1,1,1,1,1};
    //string types = "11111111";
    string set_str = "set atom ATOMID type TYPE";
    // Loop through all atoms and change the type according to type_vec
    //cout << Natoms << endl;
    for (int t=0; t < Natoms; t++)
    {
      char type = atomtypes[t+t_indx];
      const char* type_cc = &type;
      //debug << t_indx << endl;
      //debug << "type: " << type << endl;
      int atomid = t+1;
      string atomid_str = to_string(atomid);
      //debug << "atomid_str: " << atomid_str << endl;
      // Replace ATOMID with the atomid_str
      size_t ATOMID_found = set_str.find("ATOMID");
      set_str.replace(ATOMID_found, 6, " "); // Delete "ATOMID"
      set_str.insert(ATOMID_found, atomid_str);
      //cout << "set_str: " << set_str << endl;
      // Repalce TYPE with "type"
      size_t TYPE_found = set_str.find("TYPE");
      set_str.replace(TYPE_found, 4, " "); // Delete "TYPE"
      //cout << "type: " << type << endl;
      set_str.insert(TYPE_found, 1, type);
      //cout << "set_str: " << set_str << endl;
      // Now convert the set_str to a const char* so that lammps can read it
      const char* set_str_cc = set_str.c_str();
      // Change the atom type
      lmp->input->one(set_str_cc);
      set_str = "set atom ATOMID type TYPE";
    }

    t_indx = t_indx + Natoms;
    //debug << t_indx << endl;

    // Scatter the atoms with the new positions
    lammps_scatter_atoms(lmp,"x",1,3,x);

    /******** LAMMPS SETUP COMMANDS *******/
    // Loop through LAMMPS set commands and input them 
    size_t set_indx = setcommands.find(',');
    string substring;
    string rest = setcommands;
    //writefile << commands << endl;
    while (set_indx != std::string::npos)
    {  
      substring = rest.substr(0,set_indx);
      // Convert strings to const char * types so that lammps can read them
      const char* substring_cc = substring.c_str();
      //cout << "SUBSTRING: " << substring_cc << endl;
      // Input the command
      lmp->input->one(substring_cc);
      //debug << "input succeeded" << endl;
      // Find location of next command and repeat
      rest = rest.substr(set_indx+1);
      //debug << "REST: " << rest << endl;
      set_indx = rest.find(',');
      //debug << "n: " << n << endl;
    }

    // Evaluate the potential
    lmp->input->one("run 0 pre yes post no");

    // Get the forces
    pops->popsmemory->allocate(f, Natoms*3);
    lammps_gather_atoms(lmp,"f",1,3,f);
    //cout << f[0] << endl;
    //debug << x[0] << endl;

    // Extract the energy
    e = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
    double e_md = e[0];

    // Extract the stress tensor
    s = (double *) lammps_extract_compute(lmp, "S", 0, 1); // style = 0 for global data, type = 1 for vector quantity ****make type=1 to get a vector(stress tensor)

    // Extract pressure
    press = (double *) lammps_extract_compute(lmp, "S", 0, 0); // style = 0 for global data, type = 1 for vector quantity ****make type=0 to get a scalar quantity (total pressure)
    double pressure_md = press[0];
    //debug << press_md << endl;

    double si_xx = s[0];
    double si_yy = s[1];
    double si_zz = s[2];
    double si_xy = s[3];
    double si_yz = s[4];
    double si_zx = s[5];

    // Put forces into fi and calculate sum_f_diff for the forces
    //debug << "CONFIGURATION: " << i << "------------------------------------------" << endl;
    double fi_x;
    double fi_y; 
    double fi_z;
    int forceindx = 0;
    double f_diff_mag;
    double denom;
    double fio_x, fio_y, fio_z;
    //debug << "Natoms: " << Natoms << endl;
    pops->popsmemory->allocate(fi, pops->popsinput->M_in*pops->popsinput->N_tot*3);
    for (int n=0; n < Natoms; n++)
    {
      fi[n+M_indx] = f[n];
      fi_x = f[forceindx];
      fi_y = f[forceindx+1];
      fi_z = f[forceindx+2];
      fio_x = fio[forceindx + M_indx];
      fio_y = fio[forceindx+1 + M_indx];
      fio_z = fio[forceindx+2 + M_indx];

      // Calculate the magnitude of the difference of the vectors
      f_diff_mag = sqrt( pow(fi_x - fio_x, 2) + pow(fi_y - fio_y, 2) + pow(fi_z - fio_z, 2) );
      //denom = pow(0.5, pow(2, fio_x) + pow(2, fio_y) + pow(2, fio_z) );
      denom = sqrt(pow(fio_x, 2) + pow(fio_y, 2) + pow(fio_z, 2) );
      //double denom2 = pow( ( pow(fio_x, 2) + pow(fio_y, 2) + pow(fio_z, 2) ), 0.5  );

      // Calculate sums of forces for normalization
      sum_fio = sum_fio + pow(sqrt(pow(fio_x, 2) + pow(fio_y, 2) + pow(fio_z, 2) ), 2);

      /*debug << "fiox: " << fio_x << endl;
      debug << "fioy: " << fio_y << endl;
      debug << "fioz: " << fio_z << endl;	
      debug << "fix: " << fi_x << endl;
      debug << "fiy: " << fi_y << endl;
      debug << "fiz: " << fi_z << endl;*/	
      /*debug << "magdiff: " << f_diff_mag << endl;
      debug << "denom: " << denom << endl;*/
      //debug << "denom2 " << denom2 << endl;

      if (denom != 0)
      {
        Z_2_f = Z_2_f + (f_diff_mag/denom);    
        //debug << fun << endl;
      }

      /*debug << "contribution: " << contribution << endl;
      debug << "sum_f_diff: " << sum_f_diff << endl;*/
      //sum_f_diff = sum_f_diff + (abs(fi[n+M_indx] - fio[n+M_indx])/abs(fio[n+M_indx]));
      /*debug << "fio: " << fio[n+M_indx] << endl;
      debug << "fi: " << fi[n+M_indx] << endl;
      debug << "sum_f_diff: " << sum_f_diff << endl;*/


      /*if (pow(f_diff_mag,2) < 1e-20)
      {
        debug << "f_diff_mag:" << f_diff_mag << endl;
      }*/
      sum_f_diff = sum_f_diff + pow(f_diff_mag,2);
      //debug << "sum_f_diff: " << sum_f_diff << endl;

      forceindx = forceindx + 3;
    } // for (int n=0; n < Natoms; n++)
    M_indx = M_indx + Natoms*3;
    M_indx_e = M_indx_e + Natoms;

    // Calculate portion of the cohesive energy objective function
    /*pops->popsmemory->allocate(ei, pops->popsinput->M_in);
    pops->popsmemory->allocate(ei_s, pops->popsinput->M_in);*/
    ei[i] = e_md;
    if (i == 0) // Use this as the reference energy
    {
      E0_md = e_md;
      //debug << "E0_md: " << E0_md << endl;
    }
    sum_e_diff = sum_e_diff + pow(eio[i] - ei[i],2);
    sum_eio = sum_eio + pow(eio[i],2);
    // Calculate portion of the cohesive energy SHAPE objective function
    ei_s[i] = e_md - E0_md;
    sum_e_diff_s = sum_e_diff_s + pow(eio_s[i] - ei_s[i],2);
    sum_eio_s = sum_eio_s + pow(eio_s[i],2);
    /*debug << "eio_s: " << eio_s[i] << endl;
    debug << "ei_s: " << ei_s[i] << endl;
    debug << "sum_e_diff_s: " << sum_e_diff_s << endl;
    debug << "sum_eio_s: " << sum_eio_s << endl;*/
    // Calculate percent error in energy
    double e_diff = eio[i] - ei[i];
    Z_2_e = Z_2_e + abs((e_diff)/(eio[i]));
    /*debug << "e_md: " << e_md << endl;
    debug << "e_io: " << eio[i] << endl;
    debug << "sum_e_diff: " << sum_e_diff << endl;*/

    // Calculate portions of stress objective function (and convert from  kB to B)
    double sio_xx = sio[M_indx_s]*1000;
    double sio_yy = sio[M_indx_s+1]*1000;
    double sio_zz = sio[M_indx_s+2]*1000;
    double sio_xy = sio[M_indx_s+3]*1000;
    double sio_yz = sio[M_indx_s+4]*1000;
    double sio_zx = sio[M_indx_s+5]*1000;
    sum_s_diff = sum_s_diff + pow(sio_xx-si_xx,2) + pow(sio_yy-si_yy,2) + pow(sio_zz-si_zz,2) + pow(sio_xy-si_xy,2) + pow(sio_yz-si_yz,2) + pow(sio_zx-si_zx,2);
    /*debug << "dft_xx: " << sio_xx << endl;
    debug << "dft_yy: " << sio_yy << endl;
    debug << "dft_zz: " << sio_zz << endl;
    debug << "dft_xy: " << sio_xy << endl;
    debug << "dft_yz: " << sio_yz << endl;
    debug << "dft_zx: " << sio_zx << endl;
    debug << "md_xx: " << si_xx << endl;
    debug << "md_yy: " << si_yy << endl;
    debug << "md_zz: " << si_zz << endl;
    debug << "md_xy: " << si_xy << endl;
    debug << "md_yz: " << si_yz << endl;
    debug << "md_zx: " << si_zx << endl;
    debug << "sum_s_diff: " << sum_s_diff << endl;*/
    sum_sio = sum_sio + pow(sio_xx,2) + pow(sio_yy,2) + pow(sio_zz,2) + pow(sio_xy,2) + pow(sio_yz,2) + pow(sio_zx,2);
    M_indx_s = M_indx_s+6;

    // Calculate portions of the pressure objective function (and covert from kB to B)
    /*pops->popsmemory->allocate(pi, pops->popsinput->M_in);
    pops->popsmemory->allocate(pi_s, pops->popsinput->M_in);*/
    pi[i] = pressure_md;
    if (i == 0) // Use this as the reference pressure
    {
      P0_md = pressure_md;
      //debug << "E0_md: " << E0_md << endl;
    }
    sum_p_diff = sum_p_diff + pow(pio[i]*1000 - pi[i],2);
    sum_pio = sum_pio + pow(pio[i],2);
    /*debug << "pio: " << pio[i]*1000 << endl;
    debug << "pi: " << pi[i] << endl;*/
    // Calculate portion of the pressure SHAPE objective function
    pi_s[i] = pressure_md - P0_md;
    sum_p_diff_s = sum_p_diff_s + pow(pio_s[i] - pi_s[i],2);
    sum_pio_s = sum_pio_s + pow(pio_s[i],2);
    /*debug << "eio_s: " << eio_s[i] << endl;
    debug << "ei_s: " << ei_s[i] << endl;
    debug << "sum_e_diff_s: " << sum_e_diff_s << endl;
    debug << "sum_eio_s: " << sum_eio_s << endl;*/


    /* Deallocate arrays */
    pops->popsmemory->deallocate(x);
    pops->popsmemory->deallocate(f);
    /*pops->popsmemory->deallocate(fi);
    pops->popsmemory->deallocate(ei);
    pops->popsmemory->deallocate(ei_s);
    pops->popsmemory->deallocate(pi);
    pops->popsmemory->deallocate(pi_s);*/

  } // for (int i=0; i<M; i++)
  /* Deallocate group quantitiy arrays */
  pops->popsmemory->deallocate(fi);
  pops->popsmemory->deallocate(ei);
  pops->popsmemory->deallocate(ei_s);
  pops->popsmemory->deallocate(pi);
  pops->popsmemory->deallocate(pi_s);

  // Check if all of the configurations in CONFIGS have zero forces, and just set sum_fio to 1 so that Z_f doesn't blow up to infinity
  if (sum_fio == 0)
  {
    sum_fio = 1;
  }

  //debug << "N_tot: " << N_tot << endl;
  // Plug all of this into Z
  Z_2_f = ((Z_2_f)/(pops->popsinput->N_tot)) * 100;
  Z_2_e = ((Z_2_e)/(pops->popsinput->M_in)) * 100;
  double Z_f = ((sum_f_diff)/(sum_fio*pops->popsinput->N_tot));
  double Z_e = ((sum_e_diff)/(sum_eio*pops->popsinput->M_in));
  double Z_e_s ((sum_e_diff_s)/(sum_eio_s*pops->popsinput->M_in));
  double Z_s = ((sum_s_diff)/(sum_sio*pops->popsinput->M_in*6));
  double Z_p = ((sum_p_diff)/(sum_pio*pops->popsinput->M_in));
  double Z_p_s = ((sum_p_diff_s)/(sum_pio_s*pops->popsinput->M_in));
  double w_f = pops->popsinput->w_f_in;
  double w_e = pops->popsinput->w_e_in;
  double w_s = pops->popsinput->w_s_in;
  double w_e_s = pops->popsinput->w_e_s_in;
  double w_p = pops->popsinput->w_p_in;
  double w_p_s = pops->popsinput->w_p_s_in;
  double Z = w_f*Z_f + w_e*Z_e + w_s*Z_s + w_e_s*Z_e_s + w_p*Z_p + w_p_s*Z_p_s;
  
  sub_Zf2_vec[ind] = Z_2_f;
  /*debug << "sum_f_diff: " << sum_f_diff << endl;
  debug << "Z_f: " << Z_f << endl;
  debug << "sum_e_diff: " << sum_e_diff << endl;
  debug << "Z_e: " << Z_e << endl;
  debug << "sum_e_diff_s: " << sum_e_diff << endl;
  debug << "Z_e_s: " << Z_e_s << endl;
  debug << "sum_s_diff: " << sum_s_diff << endl;
  debug << "Z_s: " << Z_s << endl;
  debug << "sum_p_diff: " << sum_s_diff << endl;
  debug << "Z_p: " << Z_p << endl;
  debug << "Z_2_f: " << Z_2_f << endl;
  debug << "Z_2_e: " << Z_2_e << endl;
  debug << "Z: " << Z << endl;*/
   
  // Check Z_2_f and the corresponding parameters
  /*debug << "Z_2_f: " << Z_2_f << endl;
    debug << params[0] << endl;
    debug << params[1] << endl;
    debug << params[2] << endl;
    debug << params[3] << endl;
    debug << params[4] << endl;
    debug << params[5] << endl;
    debug << params[6] << endl;
    debug << params[7] << endl;
    debug << params[8] << endl;
    debug << params[9] << endl;
    debug << params[10] << endl;
    debug << params[11] << endl;
  debug << "-----------------------------------" << endl;*/
  

  /*// Change mutation rate if needed
  if (abs(Z_2_e) < 0.0005) //&& Z_2_f < 0.15)
  {
    debug << "GOOD FIT ACHIEVED" << endl;
    debug << Z << endl;
    debug << params[0] << endl;
    debug << params[1] << endl;
    debug << params[2] << endl;
    debug << params[3] << endl;
    debug << params[4] << endl;
    debug << params[5] << endl;
    debug << params[6] << endl;
    debug << params[7] << endl;
    debug << params[8] << endl;
    debug << params[9] << endl;
    debug << params[10] << endl;
    for (int u=0; u <= unknowns - 1; u++)
    {
      debug << "PARAM" << endl;
      debug << params[u] << endl;
      //debug << "param" << u << ": " << params[u] << endl;
    }
    mut_rate = 0.05;
  }*/

  // Get the normalized Z functions for trivial solution analysis (trivial solutions will have Z_0 = 1)
  double Z_f0 = (sum_f_diff)/(sum_fio);
  double Z_e0 = (sum_e_diff)/(sum_eio);
  double Z_s0 = (sum_s_diff)/(sum_sio);
  /*debug << "Zf0: " << Z_f0 << endl;
  debug << "Ze0: " << Z_e0 << endl;
  debug << "Zs0: " << Z_s0 << endl;
  debug << "Z: " << Z << endl;*/

  // Catch trivial solutions
  if (0.95 < Z_e0 && Z_e0 < 1.1)
  {
    Z = pow(10,100);
  }
  if (0.95 < Z_s0 && Z_s0 < 1.1)
  {
    Z = pow(10,100);
  }
  // Catch the trivial solutions
  /*if (0.9 < Z_test && Z_test < 1.1)
  {
    Z = pow(10,10);
  }*/
  //Catch the nonreal solutions
  if (Z != Z) // In cpp, if a number is not equal to itself then it is not real
  {
    Z = pow(10,10);
  }

  // Run GULP

  //int Z = 0.0; 

  return Z;
}
