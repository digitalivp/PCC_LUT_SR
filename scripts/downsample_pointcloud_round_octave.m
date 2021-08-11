function [V_dwn, C_dwn, Ns] = downsample_pointcloud_round_octave(V, C, down_fact)
% V = single(V); %this is to avoid approximations problems
% down_fact = single(down_fact);
[V_dwn, ~, Vdj] = unique(round(V/down_fact), 'rows');
C_dwn = zeros(size(V_dwn));
Ns =  accumarray(Vdj, ones(size(V(:,1))));
for k = 1:3
% 	Ns =  accumarray(Vdj, ones(size(C(:,k))));
	C_dwn(:,k) = accumarray(Vdj, C(:,k))./Ns;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%            THIS IS SHORTER, BUT TAKES MUUUCH LONGER         %%%
	% C_dwn(:,k) = accumarray(Vdj, C(:,k), [size(V_dwn,1) 1], @mean); %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
V_dwn = double(V_dwn);
C_dwn = round(C_dwn);