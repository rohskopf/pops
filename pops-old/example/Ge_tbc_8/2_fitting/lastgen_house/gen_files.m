% Loop through all "OUTCAR#" files and generate a "dftforces" file
fh_debug = fopen("DEBUG", 'w');
fh_bestfits = fopen("BESTFITS", 'w');
fh_extrainds = fopen("EXTRAINDS", 'w');
% Set the number of configs and atoms
M = 100; %trials
N = 8; % atoms
P = 10; % number of parameters
G = 100; % number of generations

num = 1;
while (num <= M)
  num_str = num2str(num);
  fh = fopen(["LASTGEN" num_str], 'r');

  % Find where "Generation: ${G}"
  l = fgetl(fh);
  while (~strcmp(l, ["Generation: " num2str(G)]) )
    l = fgetl(fh);

  end

  % Ignore next line of LASTGEN file
  z = 1;
  l = fgetl(fh);
  % Loop through the next P+2 lines to write the parameters + genes + Z value to BESTFITS
  while (z <= P+2)
    l = fgetl(fh);
  %  % Put the 2nd line in the EXTRAINDS file
    if (z == 1)
      fputs(fh_extrainds, [l " \n"]);
    end % if (z==2)
    fputs(fh_bestfits, [l "\n"]);
    z++;
  end % while (z < P+2)
  num++;

end

fclose(fh_bestfits);
fclose(fh_debug);
