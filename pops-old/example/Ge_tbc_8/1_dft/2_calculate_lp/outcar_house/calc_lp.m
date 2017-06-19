% Loop through all "OUTCAR#" files and generate a "dftforces" file
fh_debug = fopen("DEBUG", 'w');
fh_configs = fopen("CONFIGS", 'w');
% Set the number of configs and atoms
M_start = 1; %first config
M_end = 50; %last config
N = 8; % atoms
n=8;
eatoms = N*(-0.10839383); % Summed energy of all individual atoms 

a = [5.6 5.62 5.64 5.65 5.655 5.66 5.67 5.68];
E = [];

num = M_start;
for (num = a )
  num_str = num2str(num);
  fh = fopen(["OUTCAR" num_str], 'r');

  if fh > 0

  % Loop through the OUTCAR file until the "TOTAL-FORCE" line is found
  z = 1;
  l = fgetl(fh);
  while (ischar(l))

    if strcmp(l, "  LATTYP: Found a simple cubic cell.")
      fputs(fh_debug, l);
      l = fgetl(fh);
      fputs(fh_debug, l);
      box_size = str2num(l);
    end

    if strcmp(l, "  FORCE on cell =-STRESS in cart. coord.  units (eV):")
      % ingore the next 12 lines
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
	input_line = [pos_str "     " forces_str "\n"];
        fputs(fh_configs, ["1 " input_line]);
 
        atom++;
      end% while (atom <= N)

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
      E = [E energy];
      energy = num2str(energy);
      % Put into "CONFIGS"
      fputs(fh_configs, [energy "\n"]);
      % Also put the stress tensor under the energy
      fputs(fh_configs, [stress_tensor "\n"]);
      % Put the box size in
      box_size_str = num2str(box_size);
      fputs(fh_configs, [box_size_str "\n"] );
      
    end %if  

    l = fgetl(fh);
  end % while (ischar(l))
  num++;

end

end

V = a.^3;

%%%% The rest should remain the same!

[P1, S1] = polyfit(a, E, 2);
[P2, S2] = polyfit(V, E/n*N, 2);
%[P2, S2] = polyfit(V, E/N, 2);

a0 = -P1(2)/(2*P1(1));
V0 = a0^3;
Ecoh = polyval(P1,a0)/n;
%Ecoh = polyval(P1,a0)/N;
B = 2*V0*P2(1) * 160.2; % conversion from eV/A^3 to GPa

disp(sprintf('a0   = %.6f (A)',a0));
disp(sprintf('Ecoh = %.6f (eV)',Ecoh));
disp(sprintf('B    = %.6f (GPa)',B));

aaxis = [0:0.01:1]*(max(a)-min(a))+min(a);
Efit = polyval(P1,aaxis);
vaxis = aaxis.^3;

figure(1)
plot(a,E,'o',aaxis,Efit,'-');


fclose(fh_configs);
fclose(fh_debug);
