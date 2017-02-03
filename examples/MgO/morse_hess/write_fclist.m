%%% Use the neighlist file and Hessian to write a FC list

%%% INPUTS
N = 64; % Number of atoms
h = dlmread('mgo222.hessian');


nl = dlmread('NEIGHLIST1');
fh_ii = fopen('IILIST', 'w');
fh_ij = fopen('IJLIST', 'w');

% Now need to make an array fclist that has the same order of itag jtag as neighlist
[length_h width_h] = size(h);
[length_nl width_nl] = size(nl);

%fputs(fh_ii, [num2str(length_nl) "\n"]);
fputs(fh_ii, "ii\n");
fputs(fh_ij, "ij\n");

iilist = [];
for a = 1:N
  
  % i-i interaction
  indices = find(h(:,1) == a & h(:,3) == a);
  fcs = h(indices, :);
  iilist = [iilist; fcs];

end

ijlist = [];
for a=1:length_nl

  l = nl(a,:);
  itag = l(1);
  jtag = l(2);

  % i-j interaction
  indices = find(h(:,1) == itag & h(:,3) == jtag);
  fcs = h(indices, :);
  vals = fcs(:,5);

  fcs(:,5) = vals;
  ijlist = [ijlist; fcs];

end


% Write the fclist to a file
dlmwrite(fh_ii, iilist, ' ');
dlmwrite(fh_ij, ijlist, ' ');

fclose(fh_ii);
fclose(fh_ij);

