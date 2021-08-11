function neighs = get_neighbours(V, wnd_size)
% For each voxel V(k,:), get_neighbours()
% 	indicates the presence of each neighbour around
% 	-wnd_size<=V(k,m)<=wnd_size, for all dimensions 'm'.
%
% Inputs:
%	 V: N x 3 pointcloud geometry
%	 wnd_size: the maximum displacement allowed in each direction.
%
% Output:
%	 neighs: N x 1 vector, where each bit indicates the
%	 	presence or absence of neighbours around
%	 	-wnd_size<=Vtgt(k,m)<=wnd_size, for all dimensions 'm'.
%
% The bits in 'neighs' indicate the presence
% of neighbours in the following dislocations
% (from most to less significant bit):
%    1   1   1
%    1   1  -0
%    1   1  -1
%    1  -0   1
%    1  -0  -0
%    1  -0  -1
%    1  -1   1
%    1  -1  -0
%    1  -1  -1
%   -0   1   1
%   -0   1  -0
%   -0   1  -1
%   -0  -0   1
%    0   0  -1
%    0  -1   1
%    0  -1   0
%    0  -1  -1
%   -1   1   1
%   -1   1   0
%   -1   1  -1
%   -1   0   1
%   -1   0   0
%   -1   0  -1
%   -1  -1   1
%   -1  -1   0
%   -1  -1  -1
%

[L, N] = size(V);
neighs = uint32(zeros(L,1));
deltas = dec2base(0:((wnd_size*2+1)^N-1), ...
	wnd_size*2+1) - '0' - wnd_size;

%%% Faster method
deltas = deltas(1:(find(sum(abs(deltas),2)==0)-1),:);
m = ([2^(2*size(deltas,1)-1) 1]); %mudança de tipo de variavel para 
                                        %compatibilidade com matlab
for k = 1:size(deltas,1)
	[~, ps1, ps2] = intersect(V, ...
		V - ones(L,1)*deltas(k,:), 'rows');
	neighs(ps1) = neighs(ps1) + m(2);
	neighs(ps2) = neighs(ps2) + m(1);
	m = m.*[1/2 2]; % matlab nao aceita esse produto com inteiros
end

% %%% Slower method
% deltas = deltas(sum(abs(deltas),2)>0,:);
% deltas = deltas*jmp_size;

% m = 1;
% for k = 1:size(deltas,1)
% 	[~, ~, ps] = intersect(V, V - ones(L,1)*deltas(k,:), 'rows');
% 	neighs(ps) = neighs(ps) + m;
% 	m = m*2;
% end

