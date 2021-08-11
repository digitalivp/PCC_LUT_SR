function kids = child_node_occupancy_frac_round_octave(Vds, V, s)
%CHILD_NODE_OCCUPANCY(Vdwn,V)
% For each point of the downsampled point cloud Vdwn, gets which of the 8
% possible child nodes is occupied.
% Inputs:
% Vdwn: downsampled point cloud (N-by-3)
% V: original point cloud (M-by-3)
%
% Output: children, an N-by-1 matrix containing information for each child
% node from Vdwn.
%
% Child nodes represent which points from parent Vdwn exist in V, if we
% consider the parent point in [0,0,0] its child nodes, by index position
% will be:
% 1) [-0.5 -0.5 -0.5]    5) [0.5 -0.5 -0.5]
% 2) [-0.5 -0.5  0.5]    6) [0.5 -0.5  0.5]
% 3) [-0.5  0.5 -0.5]    7) [0.5  0.5 -0.5]
% 4) [-0.5  0.5  0.5]    8) [0.5  0.5  0.5]
%
% to convert from dec2bin use
% children = logical(dec2bin(children,8)-'0');


% s should be in the form of p/q, 1<s<2

[p,q] = rat(s);

% first we take one 'division period' to find values that map two values in
% the downsampled version, i.e., duplicate values in xd that would become
% multiparous coordinates
x = 0:floor(s)*(p-1);
xd = round((x)/s);

dup = find(hist(xd, x)==ceil(s))-1;
e = x-round(xd*s);
i_dup = ismember(xd,dup);
x_ch = reshape(e(i_dup),ceil(s),[])';

if q==1
    dup = 0;
end

% then we set as true for every coordinate which will generate duplicate
% children (multiparous coordinate)

r = mod(round(Vds),q);
mult_child = zeros(size(r,1),size(r,2),length(dup));
for i = 1:length(dup)
    mult_child(:,:,i) = sum(x_ch(i,:))*double(r==dup(i));
end
num_child = sum(mult_child,3);
num_child_lbl = abs(num_child)*[4;2;1];

% Label  binary     multiparous  uniparous       Children
%                   coordinate   coordinate      number 
%   7     [1 1 1]   x, y, z                      8 children
%   6     [1 1 0]   x, y         z               4 children
%   5     [1 0 1]   x, z         y               4 children
%   4     [1 0 0]   x            y, z            2 children
%   3     [0 1 1]   y, z         x               4 children
%   2     [0 1 0]   y            x, z            2 children
%   1     [0 0 1]   z            x, y            2 children
%   0     [0 0 0]                x, y, z         1 child

% conditions: 
c = cell(8,1);
c{1} = num_child_lbl == 0; % [0 0 0], parents with 1 child

c{2} = num_child_lbl == 4; % [1 0 0], parents with 2 children
c{3} = num_child_lbl == 2; % [0 1 0], parents with 2 children
c{4} = num_child_lbl == 1; % [0 0 1], parents with 2 children

c{5} = num_child_lbl == 3; % [0 1 1], parents with 4 children
c{6} = num_child_lbl == 5; % [1 0 1], parents with 4 children
c{7} = num_child_lbl == 6; % [1 1 0], parents with 4 children

c{8} = num_child_lbl == 7; % [1 1 1], parents with 8 children

% children cases for each condition:
deltas = xyz_displacements(0:1);
children = cell(8,1);

children{1} = deltas(deltas(:,1)==0 & deltas(:,2)==0 & deltas(:,3)==0,:); % single child [0 0 0]

children{2} = deltas(deltas(:,2)==0 & deltas(:,3)==0,:);% 2 children yz coordinates are uniparous
children{3} = circshift(children{2},1,2); % 2 children xz coordinates are uniparous
children{4} = circshift(children{3},1,2); % 2 children xy coordinates are uniparous

children{5} = deltas(deltas(:,1)==0,:); % 4 children, only x coordinate is uniparous
children{6} = circshift(children{5},1,2);% 4 children y coordinate is uniparous
children{7} = circshift(children{6},1,2);% 4 children z coordinate is uniparous

children{8} = deltas; % 8 children

% adding children for each condition
Vup = cell(8,1);

for i=1:8
    if sum(c{i}) > 0
        epsilon = kron(num_child(c{i},:),ones(size(children{i},1),1));
        Vup{i} = upsample_pointcloud_children_round_octave(Vds(c{i},:), 0*Vds(c{i},:), s, children{i}, epsilon);
    end
end

kids = cell(8,1);
Vp = cell(8,1);
for i=1:8
    if sum(c{i})>0
        isChildren = ismember(Vup{i}, V,'rows'); %check which child really exist in V
        l = size(children{i},1);
        isChildren = reshape(isChildren,l,[]);
    %     kids{i} = bi2de(isChildren','left-msb');
        if i == 1
            kids{i} = (isChildren .* 2.^bin2dec(char(children{i} + '0')))';
        else
            kids{i} = sum(isChildren .* 2.^bin2dec(char(children{i} + '0')))';
        end
        Vp{i} = Vds(c{i},:);
    end
end
% Vp = cell2mat(Vp);
Vp = vertcat(Vp{:});
[~, idx] = ismember(Vds,Vp, 'rows');

kids = vertcat(kids{:});
kids = kids(idx,:); % na mesma ordem de Vds
end %function