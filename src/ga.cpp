#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>

#include "popsinput.h"
#include "popsmemory.h"

#include "ga.h"

using namespace std;

using namespace POPS_NS;

Ga::Ga(POPS *pops) : PopsPointers(pops) 
{

}

Ga::~Ga() 
{

  //cout << "Destructing" << endl;
  /* Deallocate arrays */
  // Phenotype
  pops->popsmemory->deallocate(pop);
  //cout << "pop" << endl;
  pops->popsmemory->deallocate(xpop);
  //cout << "xpop" << endl;
  pops->popsmemory->deallocate(xpop_r);
  //cout << "xpop_r" << endl;
  pops->popsmemory->deallocate(subxpop);
  //cout << "subxpop" << endl;

  // Selection
  pops->popsmemory->deallocate(sorted_xpop);
  //cout << "sorted_xpop" << endl;
  pops->popsmemory->deallocate(sorted_pop);
  //cout << "sorted_pop" << endl;
  pops->popsmemory->deallocate(parents);
  //cout << "parents" << endl;
  pops->popsmemory->deallocate(rest);
  pops->popsmemory->deallocate(elite);

  // Crossover
  pops->popsmemory->deallocate(parent1);
  pops->popsmemory->deallocate(parent2);
  pops->popsmemory->deallocate(child1);
  pops->popsmemory->deallocate(child2);
  pops->popsmemory->deallocate(children);

  fclose(history);
  fclose(individuals);

  //cout << "Destructed" << endl;
};

// ---------------------------------------- INITIALIZE -----------------------------------------------------------------------
void Ga::initialize()
{

  //std::cout << "GA Initialization on proc " << pops->rank << std:: endl;
  /************************** Organize inputs to prepare GA variables **************************/
  int unknowns = pops->popsinput->unknowns_in;
  int unknowns_extra = pops->popsinput->unknowns_extra_in;
  unknowns_tot = unknowns+unknowns_extra;
  Z_vec.resize(pops->popsinput->pop_size_in);

  // Resize the vectors to be used from INPUT
  ssl.resize(unknowns_tot);
  d_vec.resize(unknowns_tot);
  intervals.resize(unknowns_tot*2);
  // Assign the vectors from INPUT
  d_vec = pops->popsinput->d_vec_in;
  intervals = pops->popsinput->intervals_in;

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

  // Allocate the binary population
  pops->popsmemory->allocate(pop, pops->popsinput->pop_size_in, length);

  // Initialize the random binary population
  srand(pops->seed); // seed the random number generator
  for (int a=0; a <= pops->popsinput->pop_size_in - 1; a++)
  {
    for (int b=0; b <= length - 1; b++)
    {
      // Randomly choose a "0" or a "1"
      int rando = rand() % 2;
      //cout << rando << endl;
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
  int randoinds = pops->popsinput->pop_size_in - pops->popsinput->extrainds_in; // Number of random individuals, also serves as an index to start adding the extra individuals
  for (int i=1; i<= pops->popsinput->extrainds_in; i++)
  {
    getline(extrastream, extraind);
    for (int j=0; j<= length-1; j++)
    {
      pop[randoinds][j]=extraind[j];
    }
    randoinds = randoinds + 1;
  }

  extrastream.close();

  /* Allocate arrays */
  // Phenotype
  pops->popsmemory->allocate(xpop, pops->popsinput->pop_size_in, pops->ga->unknowns_tot);
  pops->popsmemory->allocate(xpop_r, pops->popsinput->pop_size_in*pops->ga->unknowns_tot);
  pops->popsmemory->allocate(subxpop, pops->popsinput->pop_size_in);

  // Selection
  pops->popsmemory->allocate(sorted_xpop, pops->popsinput->pop_size_in, pops->ga->unknowns_tot);
  pops->popsmemory->allocate(sorted_pop, pops->popsinput->pop_size_in, pops->ga->length);
  pops->popsmemory->allocate(parents, pops->popsinput->pop_size_in, pops->ga->length);
  pops->popsmemory->allocate(rest, pops->popsinput->pop_size_in, pops->ga->length);
  pops->popsmemory->allocate(elite, pops->popsinput->pop_size_in*pops->popsinput->eliteperc_in, pops->ga->length);

  // Crossover
  pops->popsmemory->allocate(parent1, pops->ga->length);
  pops->popsmemory->allocate(parent2, pops->ga->length);
  pops->popsmemory->allocate(child1, pops->ga->length);
  pops->popsmemory->allocate(child2, pops->ga->length);
  pops->popsmemory->allocate(children, pops->popsinput->pop_size_in, pops->ga->length);

  if (pops->rank == 0) fprintf(history, "Generation      Z          %% Force Error            \n");
  

} //intialize()


// ---------------------------------------- DETERMINE PHENOTYPES -----------------------------------------------------------------
double* Ga::phenotype()
{

  //cout << "Determine Phenotypes on proc " << pops->rank << endl;
  // Get necessary variables from input
  vector<double> intervals = pops->popsinput->intervals_in;
  int pop_size = pops->popsinput->pop_size_in;
  vector<int> ssl = pops->ga->ssl;
  int unknowns_tot = pops->ga->unknowns_tot;
  

  // Make XPOP ARRAY -------------------------------------------------------------
  int interval_indx = 0;
  //double xpop[pop_size][unknowns];
  double current_interval[2];
  double r_max;
  double r_min;
  double z_max;
  double m;
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


// ---------------------------------------- SELECTION -----------------------------------------------------------------
// Purpose: Make a "parents" vector for crossover
void Ga::selection()
{


  //cout << "Selection" << endl;
  // Need to sort the Z_vec and retain the indices of the original permutation

  /* Get the necessary variables required for this function */
  vector<double> Z_vec = pops->Z_vec;
  vector<double> Zf2_vec = pops->Zf2_vec;
  int pop_size = pops->popsinput->pop_size_in; 
  int unknowns_tot = pops->ga->unknowns_tot;
  double eliteperc = pops->popsinput->eliteperc_in;


  // Make a vector of sorted indices for Z_vec
  vector<int> sorted_indices(pop_size);
  size_t n(0);
  generate(begin(sorted_indices), end(sorted_indices), [&]{ return n++; });

  // Now sort the sorted_indices vector according to the sorting in Z_vec
  sort(begin(sorted_indices),
    end(sorted_indices),
    [&](int i1, int i2) { return Z_vec[i1] < Z_vec[i2]; } );

  //cout << "Z_vec[0]: " << Z_vec[0] << endl;

  // Now "sorted_indices" contains the indices of the sorted Z_vec, but we still have to sort Z_vec:
  sort(Z_vec.begin(), Z_vec.end(), [](double a, double b) {
    return b > a;
  });
  //cout << "Z_vec[0]: " << Z_vec[0] << endl;
  // Sort Zf2_vec
  sort(Zf2_vec.begin(), Zf2_vec.end(), [](double a, double b) {
    return b > a;
  });
  // Now Z_vec is sorted.

  //DEBUG
  //cout << "Z VALUES" << endl;
  /*for(int i=0; i < pop_size; i++)
  {
    cout << Z_vec[i] << endl;
  }*/

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
  cout << pops->g << " " << Z_vec[0] << endl;

  /////////////////////////////////////////////////////////
  /////// MORE DATA ANALYSIS FOR OUTPUT CAN GO HERE ///////
  /////////////////////////////////////////////////////////
  fprintf(history, "%d             %.10f         %.6f                       \n", pops->g, Z_vec[0], Zf2_vec[0]);

  // Print the best individual genes
  for (int i=0; i <= length-1; i++)
  {
    fprintf(individuals, "%c", sorted_pop[0][i]);
  }
  fprintf(individuals, "\n");
  double value;
  string val_str;
  fprintf(individuals, "Generation: %d\n", pops->g);
  for (int u=0; u <= unknowns_tot-1; u++)
  {
    int d = d_vec[u];
    value = sorted_xpop[0][u];
    value = round(value*pow(10,d))/(pow(10,d));
    val_str = to_string(value);
    fprintf(individuals, "param %d: %.6f\n", u+1, value);
  }
  fprintf(individuals, "Z: %.6f\n", Z_vec[0]);

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


void Ga::crossover()
{

  //cout << "Crossover" << endl;

  // Establish variables for this function
  int pop_size = pops->popsinput->pop_size_in;

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


} //crossover

void Ga::mutation()
{
  //cout << "Mutation" << endl;

  /*Declare variables used in this function*/
  int pop_size = pops->popsinput->pop_size_in;
  double mut_rate = pops->popsinput->mut_rate_in;

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


}
