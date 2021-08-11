function [Vsr, Csr] = LUT_SR_fractional_octave(Vds, Cds, s, lut)
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
r = mod(Vds,q);
mult_child = zeros(size(r,1),size(r,2),length(dup));
for i = 1:length(dup)
    mult_child(:,:,i) = sum(x_ch(i,:))*double(r==dup(i));
end

%prepare shift for color SR
gap = round(s*dup' + x_ch)/s - dup';
[w,y] = rat(gap);
gap_even = [-q/2*ones(length(dup),1), +q/2*ones(length(dup),1)];
shift = min(abs(gap_even - (p./y).*w),[],2)/p;

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

% increments in case length(dup)>1
inc = [0 0 0
       1 0 0;
       0 1 0;
       0 0 1;
       0 1 1;
       1 0 1;
       1 1 0;
       1 1 1];
   

n = get_neighbours(Vds, 1);
if nargin < 4
    lut = build_LUT_frac_round_octave(Vds, s);
end
Vsr = cell(size(c,1),1);
Csr = cell(size(c,1),1);
neighs = xyz_displacements(-1:1);

for i=1:8
    if sum(c{i}) > 0
        if i < 2
            [Vsr{i}, Csr{i}] = upsample_pointcloud_children_round_octave(Vds(c{i},:), Cds(c{i},:), s, children{i});
        else
            neigh_state = n(c{i},:);
            max_val = sum(2.^bin2dec(char(children{i} + '0')));
            child_ok = max_val*ones(size(neigh_state));
            [ok,idx] = ismember(neigh_state, lut{i-1}(:,1));
            child_ok(ok) = lut{i-1}(idx(idx~=0),2);
%             child_ok = de2bi(child_ok,'right-msb',size(deltas,1));
            child_ok = fliplr(dec2bin(child_ok,size(deltas,1)) - '0');
            child_ok = double(reshape(child_ok',[],1));
            child_ok(child_ok==0) = NaN;
            epsilon = kron(num_child(c{i},:),ones(size(deltas,1),1));
            Vsr{i} = upsample_pointcloud_children_round_octave(Vds(c{i},:), Cds(c{i},:), s, deltas, child_ok.*epsilon);
            
            %color SR
            if nargout > 1
                a = ones(1,1,length(dup));
                a(1,1,:) = shift;
                approach_parents = sum(mult_child(c{i},:,:).*a,3);
                rep = size(deltas,1);
                approach_parents = kron(inc(i,:).*approach_parents,ones(rep,1)).*child_ok;
                approach_parents(any(isnan(approach_parents), 2), :) = [];
                d2 = inc(i,:);
                l = size(children{i},1);
                switch l
                    case 2
                        num_neighs = 14;
                    case 4
                        num_neighs = 10;
                    case 8
                        num_neighs = 7;
                end
                [~, val_father_dist] = knnsearch(neighs, 1/(s*2) * d2, 'K', num_neighs, 'NSMethod', 'kdtree', 'Distance', 'euclidean'); %there are 7 valid father & uncles for each child
                val_father_dist = double((unique(val_father_dist)));
                Vs = zeros(size(Vsr{i}));
                Cs = Vs*0;
                for k = 1:3
                    if d2(k)==0
                        Vs(:,k) = round(Vsr{i}(:,k)/s);
                    else
                        Vs(:,k) = Vsr{i}(:,k)/s - approach_parents(:,k); %shift to make possible identifing father and children
                    end
                end
                [idx, d] = knnsearch(Vds, Vs, 'K', num_neighs, 'NSMethod', 'kdtree', 'Distance', 'euclidean');              
                weights = 0*ones(size(d));
                for j = 1:size(val_father_dist,2)
                    if j == 1
                        g = 8/s;
                    else
                        g = 1;
                    end
                    weights(abs(single(d)-single(val_father_dist(j)))<=1e-6) = g/val_father_dist(j); %father % de onde tirei esse *4? acho que foi empirico
                end   
                
                for ii=1:size(Cs,1)
                    Cs(ii,:) = (weights(ii,:)*Cds(idx(ii,:),:))/sum(weights(ii,:));
                end               
                Csr{i} = round(Cs);
            end
        end
    end
end
    

[Vsr,i_u,i_dj] = unique(cell2mat(Vsr),'rows'); %acho que nao precisa desse unique ja que ja tem unique no upsample
if nargout > 1
    Csr = cell2mat(Csr);
    Csr = Csr(i_u,:);
%     Csr2 = zeros(size(Vsr));
%     Ns = accumarray(i_dj, ones(size(Csr(:,1))));
%     for k = 1:3
%         Csr(:,k) = accumarray(i_dj, Csr(:,k))./Ns; % algo pra levar em conta o numero de filhos?
%     end
end
    
end