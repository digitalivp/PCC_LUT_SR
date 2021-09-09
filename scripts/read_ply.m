function [V C] = read_ply(filename)

fp = fopen(filename,'r');
x = char(fread(fp,Inf,'char'))';
fclose(fp);

s1 = strfind(x,['end_header'])+length('end header');
s2 = find(x(s1:end)==10,1);
% delim = x(s1-1+[1:s2]);

VC = str2num(x((s1+s2):end));
V = VC(:,1:3);
if size(VC,2) == 6
	C = VC(:,4:6);
else
	C = zeros(size(V));
end