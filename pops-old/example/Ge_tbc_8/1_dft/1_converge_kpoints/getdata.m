%%% Extract forces and stress tensor from OUTCAR for convergence testing
%%% Gets DFT data to test for convergence of a sample
function getforces(val)

val_str = num2str(val);

fh = fopen("OUTCAR");
fh_w = fopen(["SUMMARY" val_str] ,'w');

N = 8;

l = fgetl(fh);

while( ischar (l) )

    if strcmp(l, "  LATTYP: Found a simple cubic cell.")
      l = fgetl(fh);
      box_size = str2num(l);
      fputs(fh_w, ["LP -----\n"]);
      fputs(fh_w, [num2str(box_size) "\n"]);
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
      % Write to file
      fputs(fh_w, ["Stress tensor -----\n"]);
      fputs(fh_w, [stress_tensor "\n"]);
      
      % Get the total pressure
      l = fgetl(fh);
      press = l;
      % Find first "="
      indx_press = find(press == "=");
      indx_press = indx_press(1);
      indx_press_k = find(press == "k");
      indx_press_k = indx_press_k(1);
      press_str = press(indx_press+1:indx_press_k-1);
      % Write pressure to file
      fputs(fh_w, ["Pressure -----\n"]);
      fputs(fh_w, [press_str "\n"]);
      

    end %if strcmp(l, ...)


    if strcmp(l, " POSITION                                       TOTAL-FORCE (eV/Angst)")
      % Then go to the next lines and write the forces of N atoms to "dftforces"
      fputs(fh_w, ["Forces-----\n"]);
      l = fgetl(fh); % ignore this line
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
	input_line = [forces_str "\n"];
        fputs(fh_w, [input_line]);
 
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
      % Write energy to file
      fputs(fh_w, ["Energy -----\n"]);
      fputs(fh_w, [energy "\n"]);
    end



  l = fgetl(fh);

end % while (ischar(l))

fclose(fh);
fclose(fh_w);

end% Function
