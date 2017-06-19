%a = 5.65; % lattice parameter
%p = 0.15; % percentage of lattice parameter to displace
u = 0.15; % upper limit
l = -0.15; % lower limit
fh = fopen("POSCAR_equilib",'r');
fh_w = fopen("POSCAR", 'w');

% get to the 8th line
a = 1;
while (a<8)
  line = fgetl(fh);

  fputs(fh_w, [line "\n"]);
  a++;
end% while (a<8)

% Get the next line
line = fgetl(fh);
% Loop through the rest of the lines
while (ischar(line))
  numline = str2num(line);
  % Add a random number to each component (between l and u)
  r1 = (u-l).*rand + l;
  r2 = (u-l).*rand + l;
  r3 = (u-l).*rand + l;

  n1 = numline(1)+r1;
  n2 = numline(2)+r2;
  n3 = numline(3)+r3;
  n1 = num2str(n1);
  n2 = num2str(n2);
  n3 = num2str(n3);

  numline = [n1 " " n2 " " n3];

  %numline = [n1 n2 n3];

  % Put the changed line into a new POSCAR
  %newline = num2str(numline);
  fputs(fh_w, [numline "\n"]);
  line = fgetl(fh);
end %while fh

fclose(fh);
fclose(fh_w);
