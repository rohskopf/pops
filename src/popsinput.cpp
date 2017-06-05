#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include "mpi.h"

#include "popsinput.h"
#include "popsmemory.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

using namespace std;

using namespace LAMMPS_NS;

using namespace POPS_NS;

PopsInput::PopsInput(POPS *pops) : PopsPointers(pops) {}

PopsInput::~PopsInput() 
{
  /*int *test;
  popsmemory->allocate(test, 3) ; 
  popsmemory->deallocate(test) ;

  lmp->input->one("units metal");*/

  //cout << "Input destructor" << endl;

  pops->popsmemory->deallocate(N_vec);
  pops->popsmemory->deallocate(eio);
  pops->popsmemory->deallocate(eio_s);
  pops->popsmemory->deallocate(sio);
  pops->popsmemory->deallocate(pio);
  pops->popsmemory->deallocate(pio_s);
  pops->popsmemory->deallocate(x_lengths);
  pops->popsmemory->deallocate(y_lengths);
  pops->popsmemory->deallocate(z_lengths);
  pops->popsmemory->deallocate(xy_vec);
  pops->popsmemory->deallocate(xz_vec);
  pops->popsmemory->deallocate(yz_vec);

  pops->popsmemory->deallocate(x_all);
  pops->popsmemory->deallocate(fio);
  pops->popsmemory->deallocate(atomtypes);

};

void PopsInput::readinput()
{

  //std::cout << "Reading INPUT on proc " << pops->rank << std::endl;

  string line;

  // Declare scalar inputs
  double value;

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

  INPUT.close();

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

}

void PopsInput::readconfigs()
{
  //cout << "Reading CONFIGS on proc " << pops->rank << endl;

  // Loop through CONFIGS to get "x_all", "fio", "eio", "sio", and "aio" arrays
  int M = pops->popsinput->M_in;

  /*pops->popsmemory->allocate(x_all, pops->popsinput->M_in);
  x_all[0] = 1.0;
  x_all[1] = 2.0;
  x_all[2] = 3.0;
  std::cout << "test2[0]: " << x_all[0] << std::endl;
  std::cout << "test2[1]: " << x_all[1] << std::endl;
  std::cout << "test2[2]: " << x_all[2] << std::endl;
  pops->popsmemory->deallocate(x_all);*/

  // Allocate group variables //
  pops->popsmemory->allocate(N_vec, M);
  pops->popsmemory->allocate(eio, M);
  pops->popsmemory->allocate(eio_s, M);
  pops->popsmemory->allocate(sio, M*6);
  pops->popsmemory->allocate(pio, M);
  pops->popsmemory->allocate(pio_s, M);
  pops->popsmemory->allocate(x_lengths, M);
  pops->popsmemory->allocate(y_lengths, M);
  pops->popsmemory->allocate(z_lengths, M);
  pops->popsmemory->allocate(xy_vec, M);
  pops->popsmemory->allocate(xz_vec, M);
  pops->popsmemory->allocate(yz_vec, M);

  double E;
  double E0;
  double P0;
  int N_indx = 0;
  int S_indx = 0;
  double xx, yy, zz, xy, yz, zx;
  double a;
  double coa;
  double x_length, y_length, z_length, xy_box ,xz_box, yz_box;
  double pressure;
  char type;
  int type_indx = 0;
  // Scan through CONFIGS to get the total number of atoms
  ifstream configs1("CONFIGS");
  string line;
  N_tot = 0;
  for (int k=1; k<=M; k++)
  {
    //cout << "CONFIG " << k << endl;
    // Get to 2nd line of data file
    for (int i=1; i <=2; i++)
    {
      getline(configs1, line);
    }
    int coor_indx = 0;
    int i = 1;
    getline(configs1, line);
    while (line != "/")
    {
      //debug << line << endl;
      i++; 
      getline(configs1, line);
    }
    int Natoms = i -1;
    N_vec[k-1] = Natoms; // Add number of atoms for this configuration to the N_vec
    N_tot = N_tot + Natoms; // Add number of atoms to N_tot
    //debug << N_tot << endl;
    // The next line is the cohesive energy
    getline(configs1, line);
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
    getline(configs1, line);
    stringstream ss2(line);
    ss2 >> xx >> yy >> zz >> xy >> yz >> zx;
    sio[S_indx] = xx;
    sio[S_indx+1] = yy;
    sio[S_indx+2] = zz;
    sio[S_indx+3] = xy;
    sio[S_indx+4] = yz;
    sio[S_indx+5] = zx;
    S_indx = S_indx+6;
    // The next line is the pressure
    getline(configs1, line); 
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
    getline(configs1, line);
    stringstream ss3(line);
    ss3 >> x_length >> y_length >> z_length >> xy_box >> xz_box >> yz_box;
    x_lengths[k-1] = x_length;
    y_lengths[k-1] = y_length;
    z_lengths[k-1] = z_length;
    xy_vec[k-1] = xy_box;
    xz_vec[k-1] = xz_box;
    yz_vec[k-1] = yz_box;
  }
  configs1.close();

  // Allocate memory for per-atom quantities
  pops->popsmemory->allocate(x_all, M*N_tot*3);
  pops->popsmemory->allocate(fio, M*N_tot*3);
  pops->popsmemory->allocate(atomtypes, M*N_tot);

  // Now get the per-atom quantities
  int atom;
  double xcoor;
  double ycoor;
  double zcoor;
  double fx, fy, fz;
  ifstream configs2("CONFIGS");
  for (int k=1; k<=M; k++)
  {
    //debug << "CONFIG " << k << endl;
    // Get to 2nd line of data file
    for (int i=1; i <=2; i++)
    {
      getline(configs2, line);
    }
    int coor_indx = 0;
    int i = 1;
    getline(configs2, line);
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
      //cout << xcoor << " " << ycoor << " "  << zcoor << " " << fx << " " << fy << " " << fz <<  endl;
      //debug << type << endl;
      // Store the atom type
      atomtypes[type_indx+i-1] = type;
      //cout << type << endl;
      //cout << atomtypes[type_indx + i -1] << endl;
      //cout << atomtypes[1] << endl;
      coor_indx=coor_indx+3;
      i++; 
      //debug << i << endl;
      getline(configs2, line);
    }
    int Natoms = i -1;
    type_indx = type_indx + Natoms;
    N_indx = N_indx+Natoms*3;
    // The next line is the cohesive energy
    getline(configs2, line);
    // The next line is the stress tensor
    getline(configs2, line);
    // The next line is the pressure
    getline(configs2, line); 
    // The next line is the box size (x, y, z)
    getline(configs2, line);
  }
  configs2.close();

  //cout << "atomtypes[0]: " << atomtypes[0] << endl;
  //cout << x_all[0] << endl;

}


