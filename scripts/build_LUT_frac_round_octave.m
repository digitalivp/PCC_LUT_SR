function lut = build_LUT_frac_round_octave(V, s)
%%
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

if s == 2
    deltas = xyz_displacements(0:1);
else
    deltas = xyz_displacements(-1:1);
end

uncles = cell(7,1);
kids = cell(7,1);
for j = 1:size(deltas,1)
    Vin = V + deltas(j,:);
    Vp  = downsample_pointcloud_round_octave(Vin,Vin*0,s);
    varphi = get_neighbours(Vp,1);
    sigma = child_node_occupancy_frac_round_octave(Vp,Vin,s);
    r = mod(round(Vp),q);
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
    c = cell(7,1);
    % c{1} = num_child_lbl == 0; % [0 0 0], parents with 1 child
    c{1} = num_child_lbl == 4; % [1 0 0], parents with 2 children
    c{2} = num_child_lbl == 2; % [0 1 0], parents with 2 children
    c{3} = num_child_lbl == 1; % [0 0 1], parents with 2 children
    c{4} = num_child_lbl == 3; % [0 1 1], parents with 4 children
    c{5} = num_child_lbl == 5; % [1 0 1], parents with 4 children
    c{6} = num_child_lbl == 6; % [1 1 0], parents with 4 children   
    c{7} = num_child_lbl == 7; % [1 1 1], parents with 8 children
    
    for i=1:length(c)
        if sum(c{i})>0
            uncles{i} = [uncles{i}; varphi(c{i},:)];
            kids{i}   = [kids{i}; sigma(c{i},:)];
        end
    end
end

% children cases for each condition:

children = cell(7,1);
deltas = xyz_displacements(0:1);
% children{1} = deltas(deltas(:,1)==0 & deltas(:,2)==0 & deltas(:,3)==0,:); % single child [0 0 0]

children{1} = deltas(deltas(:,2)==0 & deltas(:,3)==0,:);% 2 children yz coordinates are uniparous
children{2} = circshift(children{1},1,2); % 2 children xz coordinates are uniparous
children{3} = circshift(children{2},1,2); % 2 children xy coordinates are uniparous

children{4} = deltas(deltas(:,1)==0,:); % 4 children, only x coordinate is uniparous
children{5} = circshift(children{4},1,2);% 4 children y coordinate is uniparous
children{6} = circshift(children{5},1,2);% 4 children z coordinate is uniparous

children{7} = deltas; % 8 children
      

%% Step 3: Fill in LUT
% a ideia � fazer uma lut pra cada condi�ao, ja que cada condicao vai ter
% um numero e filhos especificos de possibilidades. do jeito que estou
% fazendo pode ser que eu chegue a um par neigh_ref - child_ref que nao �
% possivel para o no pai em questao.

lut = cell (7,1);
for i = 1:7 % uniparous parents do not matter for the LUT
    if ~isempty(kids{i})
%         l = size(children{i},1);
        map = [uncles{i} kids{i}]; % one map for each condition
        [~, idx_s] = sort(map);
        sortmap = [map(idx_s(:,1),1), map(idx_s(:,1),2)]; % sort by neigh state
        x = uncles{i};
        [a,~,c] = unique(x);
        neig_count = [a, accumarray(c,1)];
        
        % find neigourhood values that don't repeat. These values are added
        % directly to the LUT.
        unicos = neig_count(:,2)==1;
        non_repeat_idx = ismember(sortmap(:,1),neig_count(unicos,1));
        lut{i} = [sortmap(non_repeat_idx,1), sortmap(non_repeat_idx,2:end)];
        
        % For values that repeat, we check the statistics of their children.
        %     repeat_idx = find(ismember(sortmap(:,1),neig_count(unicos,1)) ~= 1);
        repeat_idx = ~non_repeat_idx;
        if any(repeat_idx)
            repeat_map = sortmap(repeat_idx,:);% pensar em como lidar com os repetidos
            [a,~,c] = unique(repeat_map(:,1));
            curr_lut = zeros(size(a));
            full_lut = dec2bin(repeat_map(:,2)) - '0';
%             full_lut = de2bi(repeat_map(:,2));
            n = ones(size(full_lut(:,1)));
            for k = 1:size(full_lut,2)
                curr_lut(:,k) = round(accumarray(c,full_lut(:,k))./accumarray(c,n));
            end
%             rep_pair = [a, bi2de(curr_lut,'right-msb')];
            rep_pair = [a, sum(curr_lut.*(2.^(size(curr_lut,2)-1:-1:0)),2)];
            lut{i} = [lut{i}; rep_pair];
            [~, idx_s] = sort(lut{i});
            lut{i} = [lut{i}(idx_s(:,1),1) lut{i}(idx_s(:,1),2)];
%             max_val = sum(2.^bin2dec(char(children{i} + '0')));
            min_val = min(2.^bin2dec(char(children{i} + '0')));
            lut{i}(lut{i}(:,2)==0,2) = min_val;
        end
    end
end

end %function
