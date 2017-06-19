#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <iterator>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <math.h>       /* round, floor, ceil, trunc */
#include <algorithm>
#include "mpi.h"

#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

using namespace std;
//#include "GAfunctions.h"

// LAMMPS NAMESPACE AND POINTERS
using namespace LAMMPS_NS;
char *args1[] = {
    (char *) "lmp",
    (char *) "-log",
    (char*) "none",
     0};
LAMMPS *lmp;

const int memchunk1 = 20000; // represents M*N*3 as a maximum
const int memchunk2 = 400; // represents number of unknowns as a maxium
const int memchunk3 = 20000; // represents max population size
const int memchunk4 = 10000; // represents max sub string length

// make debug file and lastgen file for data analysis 
ofstream debug("DEBUG");
ofstream lastgen("LASTGEN");
ofstream data("HISTORY");

int gen_num;
int g; // generation variable

// Initialize and xpop
int unknowns; // in potential file
int unknowns_extra; // LJ and coulomb terms
int unknowns_tot; // total unknowns

vector<int> ssl;
int extrainds;
int pop_size;
int length; // lengths of individuals
double mut_rate;
int d; // Used by many functions to extract decimal place accuracy from "d_vec"
double eliteperc; // percentage of population to be chosen as elite
vector<double> intervals;
char pop [memchunk3][memchunk4];
double xpop[memchunk3][memchunk2];
unsigned long long intval; //for conversion to z
int z;
int test;
string line; // used by many functions for file I/O
string commands; // lammps commands
string nontab_commands; // non-tabulated parameter commands
string pairstyle_commands; // pair style commands

// Fitness vector and calculation
vector<double> Z_vec;
//vector<double> sub_Z_vec;
vector<double> Z_f_vec;
vector<double> Z_e_vec;
double Z;
double x1, x2, x3;
//vector<double> sensparam_line_vec; // Size must be greater than unknowns!
//vector<double> sensparam_col_vec; // location of sensitive parameters in lines
vector<double> d_vec; // vector of decimal place accuracy values for each parameter
vector<double> sym_vec; // vector of location of lines in potential file to abide by symmetric mixing rules
string sym_coeffs; // string of coefficients for symmetric mixing
double params[memchunk2];
double params2[memchunk2];
size_t found; // for finding parameters in a line
int N;
int M;

// Force and energy values
double fio[memchunk1];
double fi[memchunk1];
int fio_indx;
double *e; //energy vector for compute extraction
double *s; // stress vector for stress tensor extraction
double *press; // pressure vector for pressure extraction
double eio[memchunk1]; // Reference (DFT) cohesive energies
double eio_s[memchunk1]; // Reference (DFT) cohesive energies for shape fitting (E-E0)
double sio[memchunk1]; // Reference (DFT) stress tensor values (xx, yy, zz, xy, yz, zx)
double pio[memchunk1]; // Reference pressures
double pio_s[memchunk1]; // Reference pressures for shape fitting (P_md - P0_md)
// Reference box lengths
double x_lengths[memchunk1];
double y_lengths[memchunk1];
double z_lengths[memchunk1];
double xy_vec[memchunk1];
double xz_vec[memchunk1];
double yz_vec[memchunk1];
double ei[memchunk1]; // MD cohesive energies
double ei_s[memchunk1]; // MD cohesive energies for shape fitting (E_md - E0_md)
double pi[memchunk1]; // MD pressures
double pi_s[memchunk1]; // MD pressures for shape fitting (P_md - P0_md)
double w_f; // force error weight
double w_e; // energy error weight
double w_e_s; // energy shape error weight
double w_s; // stress error weight
double w_p; // pressure error weight
double w_p_s; // pressure shape error weight
int natoms;
double f[memchunk1];
double x[memchunk1];
double x_all[memchunk1];
int M_indx;
double E0; // Reference energy for EV shape fitting
double E0_md;
double P0; // Reference energy for PV shape fitting
double P0_md;
int N_vec[memchunk1];
int N_tot; //  Total number of atoms across all configurations

// atom type array
char atomtypes[memchunk1];

// Selection
const int numelite = pop_size*eliteperc;
char elite[memchunk3][memchunk4];
char rest[memchunk3][memchunk4];
char parents[memchunk3][memchunk4]; // establish parent vector here
char sorted_pop[memchunk3][memchunk4];
double sorted_xpop[memchunk3][memchunk2];

// Crossover
char children[memchunk3][memchunk4];

// --------------------------------------- GET INPUT ------------------------------------------------------------------------
void getinput()
{

  int procs;
  int rank;

  //  Get the number of processes.
  procs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group
  //  Get the individual process ID.
  rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  char filename_in[64];
  sprintf (filename_in, "GETINPUT%d", rank);
  ofstream writefile(filename_in);

  // Declare scalar inputs
  double value;
  double gen_num_in;
  double extrainds_in;
  double pop_size_in;
  double mut_rate_in;
  double eliteperc_in; 
  double unknowns_in;
  double unknowns_extra_in;
  double M_in;
  double N_in;
  double w_f_in;
  double w_e_in;
  double w_e_s_in;
  double w_s_in;
  double w_p_in;
  double w_p_s_in;
  // Declare vector inputs
  vector<double> intervals_in;
  vector<int> ssl_in;
  //vector<double> sensparam_line_vec_in;
  //vector<double> sensparam_col_vec_in;
  vector<double> d_vec_in;
  vector<double> sym_vec_in;
  // Declare string inputs
  string sym_coeffs_in;

  // Open INPUT file
  ifstream INPUT("INPUT");
  // Ignore the first line
  getline(INPUT, line); 
  string characters;
  // Loop through 18 more lines
  for (int i=1; i<=18; i++)
  {
    getline(INPUT, line);
    switch (i)
    {
      case 1:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               gen_num_in = value;
               //cout << gen_num << endl;
             }
      case 2:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               extrainds_in = value;
               //cout << extrainds << endl;
             }
      case 3:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               pop_size_in = value;
             }
      case 4:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               mut_rate_in = value;
             }
      case 5:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               eliteperc_in = value;
             }
      case 6:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               unknowns_in = value;
             }
      case 7:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               unknowns_extra_in = value;
             }
      case 8:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               M_in = value;
             }
      case 9:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               N_in = value;
             }
      case 10:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_f_in = value;
             }
      case 11:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_e_in = value;
             }
      case 12:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_e_s_in = value;
             }
      case 13:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_s_in = value;
             }
      case 14:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_p_in = value;
             }
      case 15:{ //cout << line << endl;
               stringstream ss(line);
               ss >> characters >> value;
               //cout << characters << endl;
               //cout << value << endl;
               w_p_s_in = value;
             }
      case 16:{ //cout << line << endl; // VECTOR INPUTS START HERE
                string temp;
               	vector<double> vec;
	        size_t i = 0, start = 0, end;
	        do {
		  end = line.find_first_of ( ' ', start );
		  temp = line.substr( start, end );
		  if ( isdigit ( temp[0] ) || temp[0] == '-' )
		  {
		    vec.push_back ( atof ( temp.c_str ( ) ) );
		    ++i;
		  }
		  start = end + 1;
	        } while ( start );
                 // Assign the input variable to vec
                 intervals_in = vec;
             }
      case 17:{ //cout << line << endl;
                string temp;
               	vector<double> vec;
	        size_t i = 0, start = 0, end;
	        do {
		  end = line.find_first_of ( ' ', start );
		  temp = line.substr( start, end );
		  if ( isdigit ( temp[0] ) )
		  {
		    vec.push_back ( atof ( temp.c_str ( ) ) );
		    ++i;
		  }
		  start = end + 1;
	        } while ( start );
                 // Assign the input variable to vec
                 d_vec_in = vec;
             }
      /*case 15:{ //cout << line << endl;
                string temp;
               	vector<double> vec;
	        size_t i = 0, start = 0, end;
	        do {
		  end = line.find_first_of ( ' ', start );
		  temp = line.substr( start, end );
		  if ( isdigit ( temp[0] ) )
		  {
		    vec.push_back ( atof ( temp.c_str ( ) ) );
		    ++i;
		  }
		  start = end + 1;
	        } while ( start );
                 // Assign the input variable to vec
                 sym_vec_in = vec;
             }
      case 16:{ //cout << line << endl;
                string temp;
               	// Find colon
                size_t found = line.find(":");
                string substr = line.substr(found+2, line.size() );
                sym_coeffs_in = substr;
             }*/

  
    } // switch (i)

  } // for (int i=1..)

  // Now loop through the LAMMPS setup commands, and store the commands in "commands", delimited by a ","
  getline(INPUT, line);
  while (line != "----------------------------------------LAMMPS PAIR STYLE COMMANDS")
  {
    commands.append(line);
    commands.append(",");
    getline(INPUT,line);
    //debug << commands << endl;
  }
  // Continue looping through to get the pair style commands
  getline(INPUT, line);
  while (line != "----------------------------------------END OF FILE")
  {
    nontab_commands.append(line);
    nontab_commands.append(",");
    getline(INPUT,line);
    //debug << commands << endl;
  }
  /*// Continue looping through to get the pair style commands
  getline(INPUT, line);
  char character;
  while (line != "----------------------------------------END OF FILE")
  {
    pairstyle_commands.append(line);
    pairstyle_commands.append(",");
    getline(INPUT,line);
  }*/
  INPUT.close();

  // Assign values to global variables from INPUT
  gen_num = gen_num_in;
  extrainds = extrainds_in;
  pop_size = pop_size_in + extrainds;
  N = N_in;
  M = M_in;
  mut_rate = mut_rate_in;
  eliteperc = eliteperc_in;
  unknowns = unknowns_in;
  unknowns_extra = unknowns_extra_in;
  //numberoflines = numberoflines_in + 1; // +1 due to nature of indicies and the way the code works
  w_f = w_f_in;
  w_e = w_e_in;
  w_e_s = w_e_s_in;
  w_s = w_s_in;
  w_p = w_p_in;
  w_p_s = w_p_s_in;

  // Calculation of other global variables and post-processing from INPUT data
  unknowns_tot = unknowns+unknowns_extra;
  Z_vec.resize(pop_size);

  // Resize the vectors to be used from INPUT
  ssl.resize(unknowns_tot);
  //sensparam_line_vec.resize(unknowns_tot);
  //sensparam_col_vec.resize(unknowns_tot);
  d_vec.resize(unknowns_tot);
  intervals.resize(unknowns_tot*2);
  // Assign the vectors from INPUT
  //sensparam_line_vec = sensparam_line_vec_in;
  //sensparam_col_vec = sensparam_col_vec_in;
  d_vec = d_vec_in;
  intervals = intervals_in;

  // Determine the substring lengths based on the intervals
  int interval_indx = 0;
  for (int i = 0; i < unknowns_tot; i++)
  {
    double interval1 = intervals[interval_indx];
    double interval2 = intervals[interval_indx+1];
    /*debug << interval1 << endl;
    debug << interval2 << endl;*/

    double l = 1; // initial guess for length (small as possible)
    double m = (interval2 - interval1)/(pow(2,l) - 1);
    double S = log10(m); // metric to minimize
    double d = d_vec[i];
    // Find smallest d such that S > -d
    while (S > -d)
    {
      l++;
      m = (interval2 - interval1)/(pow(2,l) - 1);
      S = log10(m); // metric to minimize
    } // while (S > -d)
    //debug << "FINAL: " << l << endl;
    ssl[i] = l;
    interval_indx = interval_indx + 2;
  } // for (int i =0; i < unknowns_tot ...)

  // Calculate binary length of indivudals
  length=0;
  for(std::vector<int>::iterator j=ssl.begin();j!=ssl.end();++j)
    length += *j;


  // Loop through CONFIGS to get "x_all", "fio", "eio", "sio", and "aio" arrays
  ifstream configs("CONFIGS");
  int atom;
  double xcoor;
  double ycoor;
  double zcoor;
  double fx, fy, fz;
  double E;
  int N_indx = 0;
  int S_indx = 0;
  double xx, yy, zz, xy, yz, zx;
  double a;
  double coa;
  double x_length, y_length, z_length, xy_box ,xz_box, yz_box;
  double pressure;
  char type;
  int type_indx = 0;
  int Natoms;
  N_tot = 0;
  for (int k=1; k<=M; k++)
  {
    //debug << "CONFIG " << k << endl;
    // Get to 2nd line of data file
    for (int i=1; i <=2; i++)
    {
      getline(configs, line);
    }
    int coor_indx = 0;
    int i = 1;
    getline(configs, line);
    while (line != "/")
    {
      //debug << line << endl;
      stringstream ss(line);
      ss >> type >> xcoor >> ycoor >> zcoor >> fx >> fy >> fz;
      x_all[N_indx+coor_indx] = xcoor;
      x_all[N_indx+coor_indx+1]=ycoor;
      x_all[N_indx+coor_indx+2]=zcoor;
      fio[N_indx+coor_indx] = fx;
      fio[N_indx+coor_indx+1]=fy;
      fio[N_indx+coor_indx+2]=fz;
      //debug << xcoor << " " << ycoor << " "  << zcoor << " " << fx << " " << fy << " " << fz <<  endl;
      //debug << type << endl;
      // Store the atom type
      atomtypes[type_indx+i-1] = type;
      //debug << type << endl;
      //debug << atomtypes[type_indx + i -1] << endl;
      coor_indx=coor_indx+3;
      i++; 
      //debug << i << endl;
      getline(configs, line);
    }
    Natoms = i -1;
    N_vec[k] = Natoms; // Add number of atoms for this configuration ot the N_vec
    N_tot = N_tot + Natoms; // Add number of atoms to N_tot
    //debug << N_tot << endl;
    /*for (int i=1; i<=N; i++)
    { 
      //debug << line << endl;
      getline(configs,line);
      //debug << line << endl;
      stringstream ss(line);
      ss >> type >> xcoor >> ycoor >> zcoor >> fx >> fy >> fz;
      x_all[N_indx+coor_indx] = xcoor;
      x_all[N_indx+coor_indx+1]=ycoor;
      x_all[N_indx+coor_indx+2]=zcoor;
      fio[N_indx+coor_indx] = fx;
      fio[N_indx+coor_indx+1]=fy;
      fio[N_indx+coor_indx+2]=fz;
      //debug << xcoor << " " << ycoor << " "  << zcoor << " " << fx << " " << fy << " " << fz <<  endl;
      // Store the atom type
      atomtypes[type_indx+i-1] = type;
      //debug << type << endl;
      coor_indx=coor_indx+3;
    }*/
    // The next line is the cohesive energy
    getline(configs, line);
    stringstream ss(line);
    ss >> E;
    //debug << E << endl;
    if (k == 1) // use this as the "reference" k for EV shape fitting
    {
      E0 = E;
      //debug << E0 << endl;
    }
    eio[k-1] = E;
    eio_s[k-1] = E - E0; // E-shape curve for shape fitting
    // The next line is the stress tensor
    getline(configs, line);
    stringstream ss2(line);
    ss2 >> xx >> yy >> zz >> xy >> yz >> zx;
    sio[S_indx] = xx;
    sio[S_indx+1] = yy;
    sio[S_indx+2] = zz;
    sio[S_indx+3] = xy;
    sio[S_indx+4] = yz;
    sio[S_indx+5] = zx;
    type_indx = type_indx + Natoms;
    N_indx = N_indx+Natoms*3;
    S_indx = S_indx+6;
    // The next line is the pressure
    getline(configs, line); 
    stringstream ss4(line);
    ss4 >> pressure;
    if (k == 1) // use this as the "reference" k for PV shape fitting
    {
      P0 = pressure;
      //debug << P0 << endl;
    }
    pio[k-1] = pressure;
    pio_s[k-1] = pressure - P0; // P-shape curve for shape fitting
    // The next line is the box size (x, y, z)
    getline(configs, line);
    stringstream ss3(line);
    ss3 >> x_length >> y_length >> z_length >> xy_box >> xz_box >> yz_box;
    x_lengths[k-1] = x_length;
    y_lengths[k-1] = y_length;
    z_lengths[k-1] = z_length;
    xy_vec[k-1] = xy_box;
    xz_vec[k-1] = xz_box;
    yz_vec[k-1] = yz_box;
  }
  configs.close();

  
  //debug << "done with configs!" << endl;
  writefile.close();
} // void getinput()


// ---------------------------------------- INITIALIZE -----------------------------------------------------------------------
void initialize(int seed)
{

  /// MAKE RANDOM BINARY POPULATION -------------------------------------------
  //char pop [pop_size][unknowns*sub_str_length];
  srand(seed); // seed the random number generator
  for (int a=0; a <= pop_size - 1; a++)
  {
    for (int b=0; b <= length - 1; b++)
    {
      // Randomly choose a "0" or a "1"
      int rando = rand() % 2;
      if (rando == 1)
      {
        //pop[a][b] = 1;
        pop[a][b] = '1';
      }
      else
      {
        //pop[a][b] = 0;
        pop[a][b] = '0';
      } 
     // debug << pop[a][b];
    }
    //debug << "\n";
  } // for int a=0; a<=pop_size;...

  // Add individuals to the population
  ifstream extrastream("EXTRAINDS");
  string extraind;
  int randoinds = pop_size - extrainds; // Number of random individuals, also serves as an index to start adding the extra individuals
  for (int i=1; i<= extrainds; i++)
  {
    getline(extrastream, extraind);
    for (int j=0; j<= length-1; j++)
    {
      pop[randoinds][j]=extraind[j];
    }
    randoinds = randoinds + 1;
  }

  extrastream.close();

} //intialize()

// ---------------------------------------- DETERMINE PHENOTYPES -----------------------------------------------------------------
double* phenotype()
{

  //cout << "Population phenotypes" << endl;
  // Make XPOP ARRAY -------------------------------------------------------------
  int interval_indx = 0;
  //double xpop[pop_size][unknowns];
  double current_interval[2];
  double r_max;
  double r_min;
  double z_max;
  double m;
  //int subpop[pop_size][sub_str_length];
  double subxpop[pop_size]; // sub population of r values
  //int row[sub_str_length];
  int int_vec[pop_size];
  int dec_vec[pop_size];
  double r; // "r" value
  //double r2;
  int nextindx = 0; // For extracting the next subpop
  for (int a=0; a <= unknowns_tot - 1; a++)
  {
    //int row[ssl[a]];
    char row[ssl[a]];
    int subpop[pop_size][ssl[a]]; 
    char subpop2[pop_size][ssl[a]];
    //debug << "Parameter: " << a << endl;
    //current_interval = intervals[interval_indx:interval_indx+1];
    current_interval[0] = intervals[interval_indx];
    current_interval[1] = intervals[interval_indx+1];
    r_min = current_interval[0];
    r_max = current_interval[1];
    z_max = pow(2, ssl[a]) - 1;
    //cout << z_max << endl;
    m = (r_max - r_min)/(z_max);
    //cout << m << endl;
    // Extract appropriate sub population
    if (a == 0) 
    {
      //subpop = pop
      for (int i=0; i <= pop_size - 1; i++)
      {
        for (int j=0; j <= ssl[a] - 1; j++)
        {
          //subpop[i][j] = pop[i][j];
          subpop[i][j] = pop[i][j];
        } // for int j=0...

      } // for int i=0; i<=pop_size;...
    } // if (a == 0)
    
    else
    {
      for (int i=0; i <= pop_size - 1; i++)
      {
        for (int j=0; j <= ssl[a] - 1; j++)  
        {
          //subpop[i][j] = pop[i][j+ a*ssl[a]];
          //subpop[i][j] = pop[i][j + nextindx];
          subpop[i][j] = pop[i][j + nextindx];
        } // for int j=0...

      } // for int i=0; i<=pop_size;...
    
    } // else
  
    // DEBUG
    /*debug << "Subpop: " << a +1 << "-------------------------------" << endl;
    for (int i=0; i <= pop_size -1; i++)
    {
      for (int j=0; j <= ssl[a] -1 ; j++)
      {
        debug << subpop[i][j];
      }
      debug << endl;
    }*/

    // Take each row of subpop and make it into a single r value, thus making a subxpop array. Then add this value into the appropriate spot in the xpop array
    // Convert to whole ints row by row
    for (int i=0; i <= pop_size -1 ; i++)
    {
      // Get the row
      for (int j=0; j <= ssl[a] - 1; j++)
      {
        row[j] = subpop[i][j];
        //row2[j] = subpop2[i][j];

      } // for (int j=0; ...)

      // Convert the row to a single int 
      stringstream bin;
      for (int k : row)
      {
        
        bin << k;

        
      }

      // Convert row to decimal
      int jj = ssl[a]-1;
      double z = 0;
      for (int j = 0; j <= ssl[a]-1; j++)
      {
        if (row[j] == '1')
        {
          z += pow(2,jj);
          //cout << pow(2,j+1) << endl;
        }
        jj--;
        //debug << row[j];
      }
      //debug << endl;
      // Convert decimal to r value
      r = m*z + r_min;
      //debug << "r: " << r << endl;
      //debug << "z: " << z << endl;
 
      // Add this to subxpop array
      subxpop[i]=r;

      // Also add to the appropriate place in the xpop array
      xpop[i][a] = r;

      //cout << xpop[i][a] << endl;
    } // for (int i=0; ...)
    nextindx = nextindx + ssl[a];
    interval_indx = interval_indx + 2;
  } // for (int a=0; ...)

  // DEBUG
  /*for(int i=0; i <= pop_size-1; i++)
  {
    debug << "INDIVIDUAL: " << i+1 << endl; // Realize that the 2nd component is the parameter
    debug << "x1: " << xpop[i][0] << endl;
    debug << "x2: " << xpop[i][1] << endl;
    debug << "x3: " << xpop[i][2] << endl;
    debug << "x4: " << xpop[i][3] << endl;
    debug << "x5: " << xpop[i][4] << endl;
    debug << "x6: " << xpop[i][5] << endl;
    debug << "x7: " << xpop[i][6] << endl;
    debug << "x8: " << xpop[i][7] << endl;
  }*/


  // Make xpop_r
  double* xpop_r = 0; // This will be the returned xpop
  xpop_r = new double[pop_size*unknowns_tot];

  int h_indx = 0;
  for (int h = 0; h < pop_size; h++)
  {
    //xpop_r[h] = new double[unknowns_tot];
    for (int w = 0; w < unknowns_tot; w++)
     {
      // fill in some initial values
      // (filling in zeros would be more logic, but this is just for the example)
      h_indx = h*unknowns_tot;
      xpop_r[w+h_indx] = xpop[h][w];
     }
   }

  return xpop_r;

} //void phenotype

// --------------------------------------- ESTABLISH OBJECTIVE FUNCTION -------------------------------------------------------
double objfun(double params[])
{

  size_t n;
  string substring;
  string rest = commands;
  
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
    const char* substring_cc = substring.c_str();
    //debug << "LAMMPS INPUT: " << substring_cc << endl;
    // Input the command
    lmp->input->one(substring_cc);
    //debug << substring_cc << endl;
    // Find location of next command and repeat
    rest = rest.substr(n+1);
    n = rest.find(',');
  }

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
  M_indx=0;
  int M_indx_e=0;
  int M_indx_s = 0;
  int t_indx = 0;
  int Natoms;
    for (int i=0; i<M; i++)
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
      Natoms = N_vec[i+1];
      string Natoms_str = to_string(Natoms);
      //debug << Natoms << endl;

      // Get the x array
      for (int a=0; a<Natoms*3; a++)
      {
        x[a] = x_all[a+M_indx];
        //debug << x[a] << endl;
      }

      // Extract the cube length of the box from CONFIGS
      //double ao = aio[i];
      //double caoo = coaio[i];
      double x_length = x_lengths[i];
      double y_length = y_lengths[i];
      double z_length = z_lengths[i];
      double xy_box = xy_vec[i];
      double xz_box = xz_vec[i];
      double yz_box = yz_vec[i];
      string x_str = to_string(x_length);
      string y_str = to_string(y_length);
      string z_str = to_string(z_length);
      string xy_str = to_string(xy_box);
      string xz_str = to_string(xz_box);
      string yz_str = to_string(yz_box);
      //debug << ao << endl;

      /*// Create atoms (just to declare memory - positions don't matter yet)
      string create_atoms_str = "create_atoms 1 random NATOMS 69 NULL";
      size_t NATOMS_found = create_atoms_str.find("NATOMS");
      create_atoms_str.replace(NATOMS_found, 6, " "); // Delete "NATOMS"
      create_atoms_str.insert(NATOMS_found, Natoms_str);
      const char* create_atoms_str_cc = create_atoms_str.c_str();
      //debug << create_atoms_str << endl;
      lmp->input->one("delete_atoms region box");
      //debug << create_atoms_str_cc << endl;
      lmp->input->one(create_atoms_str_cc);*/

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
      //debug << create_atoms_str_cc << endl;
      lmp->input->one(create_atoms_str_cc);

      //vector<int> type_vec = {1,1,1,1,1,1,1,1};
      //string types = "11111111";
      string set_str = "set atom ATOMID type TYPE";
      // Loop through all atoms and change the type according to type_vec
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
        //debug << "set_str: " << set_str << endl;
        // Repalce TYPE with "type"
        size_t TYPE_found = set_str.find("TYPE");
        set_str.replace(TYPE_found, 4, " "); // Delete "TYPE"
        set_str.insert(TYPE_found, 1, type);
        //debug << "set_str: " << set_str << endl;
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

      // Evaluate the potential
      lmp->input->one("run 0 pre yes post no");

      // Get the forces
      lammps_gather_atoms(lmp,"f",1,3,f);
      //debug << f[0] << endl;
      //debug << x[0] << endl;

      // Extract the energy
      e = (double *) lammps_extract_compute(lmp, "P", 0, 0); // style = 0 for global data, type = 0 for scalar quantity
      double e_md = e[0];
      //debug << e_md << endl;

      // Extract the stress tensor
      s = (double *) lammps_extract_compute(lmp, "S", 0, 1); // style = 0 for global data, type = 1 for vector quantity ****make type=0 to get a scalar quantity (total pressure)

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
        debug << "fiz: " << fi_z << endl;	
        debug << "magdiff: " << f_diff_mag << endl;*/
	//debug << "denom: " << denom << endl;*/
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
      debug << "e_md: " << e_md << endl;
      debug << "e_io: " << eio[i] << endl;
      //debug << "sum_e_diff: " << sum_e_diff << endl;

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
      
    } // for (int i=0; i<M; i++)

  // Check if all of the configurations in CONFIGS have zero forces, and just set sum_fio to 1 so that Z_f doesn't blow up to infinity
  if (sum_fio == 0)
  {
    sum_fio = 1;
  }

  //debug << "N_tot: " << N_tot << endl;
  // Plug all of this into Z
  Z_2_f = (Z_2_f)/(N_tot);
  Z_2_e = (Z_2_e)/(M);
  double Z_f = ((sum_f_diff)/(sum_fio*N_tot));
  double Z_e = ((sum_e_diff)/(sum_eio*M));
  double Z_e_s ((sum_e_diff_s)/(sum_eio_s*M));
  double Z_s = ((sum_s_diff)/(sum_sio*M*6));
  double Z_p = ((sum_p_diff)/(sum_pio*M));
  double Z_p_s = ((sum_p_diff_s)/(sum_pio_s*M));
  double Z = w_f*Z_f + w_e*Z_e + w_s*Z_s + w_e_s*Z_e_s + w_p*Z_p + w_p_s*Z_p_s;
  debug << "sum_f_diff: " << sum_f_diff << endl;
  debug << "Z_f: " << Z_f << endl;
  debug << "Z_f*w_f: " << Z_f*w_f << endl;
  debug << "sum_e_diff: " << sum_e_diff << endl;
  debug << "Z_e: " << Z_e << endl;
  debug << "Z_e*w_e: " << Z_e*w_e << endl;
  debug << "sum_e_diff_s: " << sum_e_diff << endl;
  debug << "Z_e_s: " << Z_e_s << endl;
  debug << "Z_e_s*w_e_s: " << Z_e_s*w_e_s << endl;
  debug << "sum_s_diff: " << sum_s_diff << endl;
  debug << "Z_s: " << Z_s << endl;
  debug << "Z_s*w_s: " << Z_s*w_s << endl;
  debug << "sum_p_diff: " << sum_s_diff << endl;
  debug << "Z_p: " << Z_p << endl;
  debug << "Z_2_f: " << Z_2_f << endl;
  debug << "Z_2_e: " << Z_2_e << endl;
  debug << "Z: " << Z << endl;
   
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

  return Z; 

} // double objfun()


// ---------------------------------------- CALCULATE FITNESS VECTOR ----------------------------------------------------------
vector<double> calcfitness(double sub_xpop_r[], int nepp)
{

  // Check sub_xpop_r
  /*for (int i=0; i < nepp; i++)
  {
    debug << sub_xpop_r[i] << endl;
  }*/

  //double params[unknowns];
  
  vector<double> sub_Z_vec;
  sub_Z_vec.resize(pop_size);
  int i_indx = 0;
  for (int i=0; i < nepp/unknowns_tot; i++)
  {
    for (int u=0; u < unknowns_tot; u++)
    {
      //debug << u << endl;
      params2[u] = sub_xpop_r[u + i_indx];
      //debug << "param" << u << ": " << params2[u] << endl;
    }
    sub_Z_vec[i] = objfun(params2);
    //debug << sub_Z_vec[i] << endl;
    i_indx = i_indx + unknowns_tot;
    //debug << i << endl;
  }

  //debug << nepp << endl;


  return sub_Z_vec;
 
} // void calcfitness()


// ---------------------------------------- SELECTION -----------------------------------------------------------------
// Purpose: Make a "parents" vector for crossover
void selection()
{
  //cout << "Selection" << endl;
  // Need to sort the Z_vec and retain the indices of the original permutation

  // Make a vector of sorted indices for Z_vec
  vector<int> sorted_indices(pop_size);
  size_t n(0);
  generate(begin(sorted_indices), end(sorted_indices), [&]{ return n++; });

  // Now sort the sorted_indices vector according to the sorting in Z_vec
  sort(begin(sorted_indices),
    end(sorted_indices),
    [&](int i1, int i2) { return Z_vec[i1] < Z_vec[i2]; } );

  // Now "sorted_indices" contains the indices of the sorted Z_vec, but we still have to sort Z_vec:
  sort(Z_vec.begin(), Z_vec.end(), [](double a, double b) {
    return b > a;
  });
  // Now Z_vec is sorted.

  //DEBUG
  //debug << "Z VALUES" << endl;
  for(int i=0; i < pop_size; i++)
  {
    //debug << Z_vec[i] << endl;
  }

  // Now we must sort xpop according to the sorted indices
  //int sorted_pop[pop_size][unknowns*sub_str_length];
  //double sorted_xpop[pop_size][unknowns];
  int u = 0;
  int indx;
  for (int i=0; i <= pop_size - 1; i++)
  {
    indx = sorted_indices[i];
    for (int u=0; u <= unknowns_tot -1; u++)
    {
      sorted_xpop[i][u] = xpop[indx][u];
    }
    for (int j=0; j <= length -1; j++)
    {
      sorted_pop[i][j] = pop[indx][j];
      //debug << sorted_pop[i][j];
    }
    //debug << endl;
  }
  // At this point, "Z_vec" is sorted and so is "sorted_xpop" and "sorted_pop"
  //debug << "GENERATION: " << g << "----------------------------------------" << endl;

  // Data analysis
  data << g << " " << Z_vec[0] << endl;

  // Last generation data
  if (g == gen_num)
  {
    double value;
    string val_str;
    lastgen << "Generation: " << g << endl;
    for(int i=0; i <= pop_size-1; i++)
    {
      lastgen << "--------------------- Individual: " << i+1 << endl;
      for (int j=0; j < length; j++)
      {
        lastgen << sorted_pop[i][j];
      }
      lastgen << endl;
      for (int u=0; u <= unknowns_tot-1; u++)
      {
        d = d_vec[u];
        value = sorted_xpop[i][u];
        value = round(value*pow(10,d))/(pow(10,d));
        val_str = to_string(value);
        lastgen << "param" << u+1 << ": " << val_str << endl;
      }
      lastgen << "Z: " << Z_vec[i] << endl;
    }
  }

    double value;
    string val_str;
    lastgen << "Generation: " << g << endl;
      lastgen << "--------------------- Individual: " << "BEST" << endl;
      for (int j=0; j < length; j++)
      {
        lastgen << sorted_pop[0][j];
      }
      lastgen << endl;
      for (int u=0; u <= unknowns_tot-1; u++)
      {
        d = d_vec[u];
        value = sorted_xpop[0][u];
        value = round(value*pow(10,d))/(pow(10,d));
        val_str = to_string(value);
        lastgen << "param" << u+1 << ": " << val_str << endl;
      }
      lastgen << "Z: " << Z_vec[0] << endl;

  // "sorted_xpop" and Z_vec are now sorted from lowest to highest, as to find the minimum
  /*for (int i=0; i <= pop_size - 1; i++)
  {
    cout << sorted_xpop[i][u] << " " << sorted_xpop[i][u+1] << " " << sorted_xpop[i][u+2] << " " << Z_vec[i] << endl;
  }*/

  // The next step is to select individuals based on whatever selection algorithm


  // Now we take the top "eliteperc" to save for the next generation
  int numelite = pop_size*eliteperc;
  /*int elite[numelite][unknowns*sub_str_length];
  int rest[pop_size-numelite][unknowns*sub_str_length];
  int parents[(pop_size/2)][unknowns*sub_str_length]; // establish parent vector here*/
  // DEBUG
  //debug << "ELITE and PARENT" << endl;
  for (int i =0; i <= pop_size*eliteperc - 1; i++) // make array of elite individuals
  {
    /*elite[i][u] = sorted_xpop[i][u];
    elite[i][u+1] = sorted_xpop[i][u+1];
    elite[i][u+2] = sorted_xpop[i][u+2];*/
    if (i <= pop_size*eliteperc - 1)
    {
      for (int j=0; j <= length - 1; j++)
      {
        elite[i][j] = sorted_pop[i][j];
        parents[i][j] = sorted_pop[i][j];
        //debug << elite[i][j];
        //cout << j << endl;
      }
      //debug << endl;
    }// if
    else
    {
      for (int j=0; j <= length - 1; j++)
      {
        rest[i][j] = sorted_pop[i][j];
        //cout << j << endl;
      }
    } //else
  } //for

  // Make the parents out of the "elite" and (0.5 - eliteperc) of the "rest" (see vec.cpp for random selection from vectors)
  // Make a vector of indicies
  vector<int> indices(pop_size - numelite);
  for (int i=0; i <= (pop_size-numelite)-1; i++)
  {
    indices[i] = i;
    //cout << indices[i] << endl;
  }
  // Pick randomly from the "rest" and add them to the "parents"
  int ii;
  int rando;
  int randoindx;
  //debug << "numelite: " << numelite << endl;
  for (int i = numelite; i <= (pop_size/2)-1; i++)
  {
    ii = i - numelite;
    // Make a random number in the size of the indices vector
    rando = rand() % indices.size();
    // Extract the random index
    randoindx = indices[rando];
    // Delete this index from the indices vector 
    indices.erase (indices.begin()+(rando));
    
    //debug << "Parents: ";
    for (int j=0; j <= length -1 ; j++)
    {
      
      parents[i][j] = pop[randoindx][j];
      //debug << parents[i][j];
    }
    //debug << endl;
  }
} // selection


// ---------------------------------------- CROSSOVER -----------------------------------------------------------------
void crossover()
{
  //cout << "Crossover" << endl;
  int f = 1;
  int child_indx = 1;
  int limit = pop_size/2;
  // Preallocate array of children
  //int children[pop_size/2][unknowns*sub_str_length];
  // Generate a vector of indices for random selection from parents
  vector<int> parent_indices(pop_size/2);
  for (int i=0; i <= (pop_size/2)-1; i++)
  {
    parent_indices[i] = i;
    //cout << indices[i] << endl;
  } // for (int=0..)

  int indx;
  //int length = unknowns*sub_str_length; // lengths of individuals
  int numb_indices; // this will decrease as parents are selected, to prevent them from being selected twice
  int rand_indx1;
  int indx1;
  int rand_indx2;
  int indx2;
  int parent1[length];
  int parent2[length];
  int child1[length];
  int child2[length];
  int mutindx;
  int addindx = 0; // for adding to the children array
  //int children[pop_size/2][unknowns*sub_str_length];
  while(f <= limit/2)
  {
    // Generate a random index, 0 <= indx <= length-1, to crossover
    indx = rand() % length;
    //cout << indx << endl;
    
    // Take 2 parents at random ---------------
    // Need number of indices in list to generate a random number from 1 to the number of indices
    numb_indices = parent_indices.size();
    // Now need a random indx between 1 and length of parent_indices
    rand_indx1 = rand() % numb_indices;
    // Extract this indx from the indx list
    indx1 = parent_indices[rand_indx1];
    // Now eliminate this option from the list of acceptable indices
    parent_indices.erase (parent_indices.begin() + (rand_indx1));
    // Update the number of indices
    numb_indices = parent_indices.size();
    // Make another random index
    rand_indx2 = rand() % numb_indices;
    // Extract this indx from the indx list
    indx2 = parent_indices[rand_indx2];
    // Now eliminate this option from the list of acceptable indices
    parent_indices.erase (parent_indices.begin() + (rand_indx2));

    // Extract the two parents
    //debug << "Parent1: ";
    for (int j=0; j <= length-1; j++)
    {
      parent1[j] = parents[indx1][j];
      parent2[j] = parents[indx2][j];
      //debug << parent1[j];
    }
    //debug << endl;

    // Make the crossover switch ------------------------
    // Extract the tails
    /*int tail1[length-indx];
    int tail2[length-indx];
    int ii = 0;
    for (int j=indx; j <= length - 1; j++)
    {
      tail1[ii] = parent1[j];
      tail2[ii] = parent2[j];
      ii++;
    }*/
    // Exchange the tails while birthing the children
    for (int j=0; j <= indx-1; j++) // Make the heads
    {
      child1[j] = parent1[j];
      child2[j] = parent2[j];
      //cout << "Head indx: " << j << endl;
    }
    for (int j=indx; j <= length - 1; j++) // Make the tails by switching the last parts of parent1 and parent2
    {
      child1[j] = parent2[j];
      child2[j] = parent1[j];
      //cout << "Tail indx: " << j << endl;
    }
    
    // Mutate identical individuals if necessary ------------
    // Make children vectors out of the children arrays for comparisons
    vector<int> child1_vec(length);
    vector<int> child2_vec(length);
    for (int j=0; j <= length-1; j++)
    {
      child1_vec[j] = child1[j];
      child2_vec[j] = child2[j];
    }
    // See if child1 equals child2 or the parents
    if( equal (child1_vec.begin(), child1_vec.end(), child2) || equal (child1_vec.begin(), child1_vec.end(), parent1) || equal (child1_vec.begin(), child1_vec.end(), parent2) )
    {
      //cout << "MUTATE CHILD1" << endl;
      // Generate a random index between 0 and (length-1)
      mutindx = rand() % length;
      // Mutate the child at this index
      //cout << child1[mutindx] << endl;
      if (child1[mutindx] == '0')
        child1[mutindx] = '1';
      else
        child1[mutindx] = '0';
    }
    // See if child2 equals child1 or the parents
    if( equal (child2_vec.begin(), child2_vec.end(), child1) || equal (child2_vec.begin(), child2_vec.end(), parent1) || equal (child2_vec.begin(), child2_vec.end(), parent2) )
    {
      //cout << "MUTATE CHILD2" << endl;
      // Generate a random index between 0 and (length-1)
      mutindx = rand() % length;
      // Mutate the child at this index
      //cout << child1[mutindx] << endl;
      if (child1[mutindx] == '0')
        child1[mutindx] = '1';
      else
        child1[mutindx] = '0';
    }

    // Add the child to the children array
    for (int j=0; j <= length-1; j++)
    {
      children[addindx][j] = child1[j];
      children[addindx+1][j] = child2[j];
    }

    addindx = addindx + 2;
    f++;
  } // while(f <= limit/2)  

} // void crossover()


// ---------------------------------------- MUTATION -------------------------------------------------------------------------------------------------------------------------------------
void mutation()
{
  //cout << "Mutation and adding individuals to new population" << endl;
  // If a random number between 1 and 100 is less than the mutation rate times 100, then apply mutation
  int num_children = pop_size/2;
  int f = 0;
  int mut_indx;
  int rando;
  while (f <= num_children - 1)
  {
    // Generate a random number between 1 and 100
    rando = rand() % 100 + 1;
    // Test to see if it is less than the mutation rate times 100
    if (rando < mut_rate*100)
    {
      // Generate a random index between 0 and (length-1) and mutate this index
      //debug << "MUTATION --------------" << endl;
      mut_indx = rand() % length;
      //debug << "f: " << f << endl;
      //debug << children[f][mut_indx] << endl;
      if (children[f][mut_indx] == '0'){
        children[f][mut_indx] = '1';
      }
      else{
        children[f][mut_indx] = '0';
      }
      //debug << children[f][mut_indx] << endl;
    } // if(rando<mut_rate*100)
    //debug << children[f][mut_indx] << endl; 
    f++; 
  }//while(f<=num_children)
  // Make the new population with the parents and the children
  for (int i=0; i <= (pop_size/2)-1; i++) // Add the parents
  {
    for (int j=0; j <= (length-1); j++)
    {
      pop[i][j] = parents[i][j];
    }
  }
  for (int i=(pop_size/2); i <= pop_size-1; i++) // Add the children
  {
    for (int j=0; j <= (length-1); j++)
    {
      pop[i][j] = children[i - (pop_size/2)][j];
    }
  }
  
  /*// DEBUG
  for (int i=0; i <= pop_size-1; i++)
  {
    debug << "Mutated: ";
    for (int j=0; j <= length-1; j++)
    {
      debug << pop[i][j];
    } 
    debug << endl;
  }*/

} //void mutation()



int main(int argc, char **argv) //-------------------------------------------------------------------- MAIN -------------------------------------------------------------------------------------
{

  int color,key,global,local;
  MPI_Comm comm;

  /* Initialize MPI */
  MPI_Init(&argc,&argv);

  bool time_seed;
  /* Extract optional random number seed */
  double seed;
  if (argc == 3) {
    char *seed_char = argv[1]; // obtain a value from the input to seed the random number generator
    seed = atoi(seed_char); // convert the seed value to an int to be used for srand()
    time_seed = false;
  }
  else {
    seed = 0;
    time_seed = true;
  }

  if (time_seed)
  {
    // Calculate current number of seconds since the new year
    time_t now;
    struct tm newyear;
    double seconds;
    time(&now);  // get current time; same as: now = time(NULL)
    newyear = *localtime(&now);
    newyear.tm_hour = 0; newyear.tm_min = 0; newyear.tm_sec = 0;
    newyear.tm_mon = 0;  newyear.tm_mday = 1;
    seconds = difftime(now,mktime(&newyear));
    // Add this number of seconds to the seed
    seed = seed + seconds;
  }

  // Split the communicators so that multiple instances of LAMMPS can be run
  MPI_Comm_rank(MPI_COMM_WORLD, &global);
  color = global / 1; // Change "1" to 2 in order to use 2 procs per instance, etc..
  key = global; 
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
  MPI_Comm_rank(comm,&local);
  
  // Clarify Open MPI quantities
  int num_elements_per_proc = atoi(argv[2]);

  int procs;
  int rank;

  //  Get the number of processes.
  procs = MPI::COMM_WORLD.Get_size ( ); //Get_size gets number of processes (np) in communicator group

  //  Get the individual process ID.
  rank = MPI::COMM_WORLD.Get_rank ( ); // Get_rank gets the rank of the calling process in the communicator

  lmp = new LAMMPS(3,args1,comm);

  char filename_in[64];
  sprintf (filename_in, "PROC%d", rank);
  ofstream writefile(filename_in);
  writefile << "This is process " << rank << endl;
  writefile << "Process global rank: " << global << endl;
  writefile << "Process local rank: " << local << endl;

  // Declare global variables from INPUT and CONFIGS
  getinput();

  if (rank == 0)
  {
    debug << "Master process: " << rank << endl;
    debug << "Number of elements per proc: " << num_elements_per_proc << endl;
    debug << "Number of procs: " << procs << endl;
    debug << "Number of unknowns: " << unknowns_tot << endl;
  }


  // Initialize the algorithm (create random population, etc..)
  if (rank == 0)
  {
    initialize(seed); 
    //debug << "Random number seed: " << seed << endl;
    //debug << "Initialization succeeded" << endl;
  }

  // Loop through LAMMPS commands and input them 
  size_t n = commands.find(',');
  string substring;
  string rest = commands;
  //writefile << commands << endl;
  while (n != std::string::npos)
  {
    substring = rest.substr(0,n);
    // Convert strings to const char * types so that lammps can read them
    const char* substring_cc = substring.c_str();
    //debug << "SUBSTRING: " << substring_cc << endl;
    // Input the command
    lmp->input->one(substring_cc);
    //debug << "input succeeded" << endl;
    // Find location of next command and repeat
    rest = rest.substr(n+1);
    //debug << "REST: " << rest << endl;
    n = rest.find(',');
    //debug << "n: " << n << endl;
  }
  //writefile << "LAMMPS SETUP COMMANDS SUCCESSFULLY INPUTTED" << endl;

  // Loop through the evolutionary process
  writefile << "Beginning evolution with " << gen_num << " generations..." <<  endl;
  double* xpop_r = NULL;
  g = 1;
  while (g <= gen_num)
  {

    writefile << "Generation " << g << " started" << endl;

    if (rank == 0)
    {
      debug << "Generation " << g << " started" << endl;
      xpop_r = phenotype();
      /*for (int i = 0; i < unknowns_tot*pop_size; i++)
      {
        debug << "xpop_r[" << i+1 << "]: " << xpop_r[i] << endl;
      }*/
    }
    //writefile << "phenotype() complete" << endl;

    // Check the xpop array
    /*for (int i = 0; i < pop_size*unknowns_tot; i++)
    {
      writefile << "xpop_r[" << i+1 << "]: " << xpop_r[i] << endl;
    }*/

    // Create sub_xpop_r
    double *sub_xpop_r = (double *)malloc(sizeof(double) * num_elements_per_proc);

    // Scatter the random numbers from the root process to all processes in the MPI world
    //MPI::COMM_WORLD.Scatter(xpop_r, num_elements_per_proc, MPI::DOUBLE, sub_xpop_r, num_elements_per_proc, MPI::DOUBLE, 0);
    //MPI::comm::Scatter(xpop_r, num_elements_per_proc, MPI::DOUBLE, sub_xpop_r, num_elements_per_proc, MPI::DOUBLE, 0);
    //writefile << "MPI_Scatter begin" << endl;
    MPI_Scatter(xpop_r, num_elements_per_proc, MPI::DOUBLE, sub_xpop_r, num_elements_per_proc, MPI::DOUBLE, 0, MPI_COMM_WORLD);
    //writefile << "MPI_Scatter complete" << endl;

    // Check the sub xpop array
    /*for (int i = 0; i < num_elements_per_proc; i++)
    {
      writefile << "sub_xpop_r[" << i+1 << "]: " << sub_xpop_r[i] << endl;
    }*/

    vector<double> sub_Z_vec = calcfitness(sub_xpop_r, num_elements_per_proc);

    // Convert vector to array for MPI_Gather
    int sub_Z_vec_size = sub_Z_vec.size();
    double sub_Z_vec_arr[sub_Z_vec_size];
    copy(sub_Z_vec.begin(), sub_Z_vec.end(), sub_Z_vec_arr);
    
    // Check the contents of sub_Z_vec_arr
    /*for (int i=0; i<num_elements_per_proc/unknowns_tot; i++)
    {
      writefile << "sub_Z_vec_arr[" << i+1 << "]: " << sub_Z_vec_arr[i] << endl;
    }*/

    double *Z_vec_arr;
    Z_vec_arr=new double [pop_size*sizeof(double)];;
    //writefile << "Z_vec_size: " << pop_size << endl;

    // Gather the sub_Z_vec vectors into Z_vec
    //MPI::COMM_WORLD.Gather(sub_Z_vec_arr,num_elements_per_proc , MPI::DOUBLE, Z_vec_arr, num_elements_per_proc, MPI::DOUBLE, 0);
    /*writefile << "nepp: " << num_elements_per_proc << endl;
    writefile << "pop_size: " << pop_size << endl;
    writefile << "nepp/pop_size: " << num_elements_per_proc/pop_size << endl;*/
    /*MPI_Gather(sub_Z_vec_arr,num_elements_per_proc/pop_size , MPI::DOUBLE, Z_vec_arr, num_elements_per_proc/pop_size, MPI::DOUBLE, 0, MPI_COMM_WORLD);*/
    MPI_Gather(sub_Z_vec_arr,pop_size/procs, MPI::DOUBLE, Z_vec_arr, pop_size/procs, MPI::DOUBLE, 0, MPI_COMM_WORLD);
    //writefile << "MPI_Gather complete" << endl;

    // Check the contents of Z_vec_arr
    /*if (rank == 0)
    {
      for (int i=0; i<pop_size; i++)
      {
        writefile << "Z_vec_arr[" << i+1 << "]: "  << Z_vec_arr[i] << endl;
      }
    }*/

    // Convert Z_vec_arr to Z_vec
    //vector<double> Z_vec(Z_vec_arr, Z_vec_arr + pop_size );
    // Check the contents of Z_vec
    /*for (int i=0; i<pop_size; i++)
    {
      writefile << Z_vec[i] << endl;
    }*/

    // Convert Z_vec_arr to Z_vec
    if (rank == 0)
    {
      for (int i=0; i<pop_size; i++)
      {
        //debug << Z_vec_arr[i] << endl;
        Z_vec[i] = Z_vec_arr[i];
      }
    }



    if (rank == 0)
    {
      selection();
      crossover();
      mutation();
      g++;
    }
  }

  // Clarify that procs completed tasks
  writefile << "Proc " << global << " finished!" << endl;

  lastgen.close();
  writefile.close();
  debug.close();

  // Delete LAMMPS
  delete lmp;

  // Close MPI
  int MPI_Comm_free(MPI_Comm *comm);
  MPI_Finalize();
} // main
