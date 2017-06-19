% Loop through all "OUTCAR#" files and generate a "dftforces" file
fh_debug = fopen("DEBUG", 'w');
fh_configs = fopen("CONFIGS", 'w');
% Set the number of atoms
N = 8; % atoms
eatoms = N*(-0.10839383); % Summed energy of all individual atoms 

for (num = 1:5 )
  num_str = num2str(num);
  fh = fopen(["OUTCAR" num_str], 'r');

  if fh > 0

  % Loop through the OUTCAR file until the "TOTAL-FORCE" line is found
  z = 1;
  l = fgetl(fh);
  while (ischar(l))

    if strcmp(l, "  LATTYP: Found a simple cubic cell.")
      l = fgetl(fh);
      %fputs(fh_debug, [l "\n"]);
      equal_indx = find(l == "=");
      a_str = l(equal_indx+1:end);
      %fputs(fh_debug, [a_str "\n"]);
      l = fgetl(fh);
      l = fgetl(fh);
      equal_indx = find(l == "=");
      c_over_a_str = l(equal_indx+1:end);
      %fputs(fh_debug, [c_over_a_str "\n"]);
      %fputs(fh_debug, [l "\n"]);
      box_size = str2num(l);
    end

    if strcmp(l, "  FORCE on cell =-STRESS in cart. coord.  units (eV):")
      % advance 14 lines
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      % Find the "B" in "kB"
      indx = find(l == "B");
      % Get the stress values
      stress_tensor = l(indx+1:end);
      % Get the total pressure
      l = fgetl(fh);
      indx = find(l == "=");
      indx = indx(1);
      k_indx = find(l == "k");
      k_indx = k_indx(1);
      pressure = l(indx+1:k_indx-1);

      % Now get the box size
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      % The following lines are the lattice vectors and reciprocal lattice vectors
      l = fgetl(fh);
      l_num = str2num(l);
      x_length = l_num(1);
      l = fgetl(fh);
      l_num = str2num(l);
      y_length = l_num(2);
      l = fgetl(fh);
      l_num = str2num(l);
      z_length = l_num(3);
      box_size2 = [x_length y_length z_length];
            

    end %if strcmp(l, ...)

    if strcmp(l, " POSITION                                       TOTAL-FORCE (eV/Angst)")
      % Then go to the next lines and write the forces of N atoms to "dftforces"
      l = fgetl(fh); % ignore this line
      % Write two lines to "dftforces" and "dftpos"
      fputs(fh_configs, "\n");
      fputs(fh_configs, ["--------------- CONFIGURATION " num_str " ----------------\n"]);
      atom = 1;
      while (atom <= N)
        l = fgetl(fh);
        l_num = str2num(l);
        % Get the positions and forces
        x = l_num(1);
        y = l_num(2);
        z = l_num(3);
        fx = l_num(4);
        fy = l_num(5);
        fz = l_num(6);

        % Get the positions
        pos = [x y z];
        pos_str = num2str(pos);

        % Get the forces
        forces = [fx fy fz];
        forces_str = num2str(forces);

        % Write the positions and forces to CONFIGS
        inputline = [1 pos forces];
        fprintf(fh_configs, '%1.0f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n',inputline);
 
        atom++;
      end% while (atom <= N)
      fputs(fh_configs, "/\n");

      % Now go 13 more lines
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
      l = fgetl(fh);
  
      %%% Get the energy at the end of the line
      % Find where there are equal signs in the line
      indices = find(l=="=");
      indx = indices(2); % The quantity of interest is beyond the 2nd equal sign
      energy = l(indx+1:end);
      % Subtract N*eatom to get cohesive energy
      energy = str2num(energy);
      energy = energy - eatoms;
      %E = [E energy];
      %energy = num2str(energy);
      % Put into "CONFIGS"
      inputline = energy;
      fprintf(fh_configs, '%1.15f\n',inputline);
      %fputs(fh_configs, [energy "\n"]);
      % Also put the stress tensor under the energy
      fputs(fh_configs, [stress_tensor "\n"]);
      % Put the pressure under the stress tensor
      fputs(fh_configs, [pressure "\n"]);
      % Put the box size in
      inputline = [x_length y_length z_length 0.0 0.0 0.0];
      fprintf(fh_configs, '%1.15f %1.15f %1.15f %1.15f %1.15f %1.15f\n',inputline);
      %box_size_str = num2str(box_size);
      %fputs(fh_configs, [a_str c_over_a_str "\n"] );
      
    end %if  

    l = fgetl(fh);
  end % while (ischar(l))
  num++;

end

end

fclose(fh_configs);
fclose(fh_debug);
