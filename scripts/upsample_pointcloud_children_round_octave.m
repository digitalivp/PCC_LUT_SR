function [Vup, Cup] = upsample_pointcloud_children_round_octave(Vd, Cd, s, children, child_ok)
up_delta = ones(size(children,1),1);
up_VC = ones(size(Vd,1),1);
if nargin < 5
    child_ok = ones(length(up_VC)*length(up_delta),1);
end
e = child_ok < 0; %this is to flipud children order that have negative increments
% for s>2 this value of e needs to be corrected...(TO BE IMPLEMENTED) 
Vup = round(kron((Vd*s), up_delta) + (abs(child_ok).*kron(up_VC,children) - e));

if nargout > 1
    Cup = kron(Cd, up_delta);
    Cup(any(isnan(Vup), 2), :) = [];
end
Vup(any(isnan(Vup), 2), :) = [];
Vup = double(Vup);

%check for values greater than the maximum allowed
L = nextpow2(max(max(Vd*s)));
t_max = Vup > (2^L-1);
idx_nok_max = any(t_max,2);
if any(idx_nok_max)
    Vup(t_max) = (2^L-1);
    if nargout > 1
%         [~,~,idj] = unique(Vup,'rows','stable');
        [~,~,idj] = unique(Vup,'rows','first');
        idj = sort(idj);
        Ns =  accumarray(idj, ones(size(idx_nok_max)));
        Cup2 = zeros(size(Ns,1),3);
        for k = 1:3
            Cup2(:,k) = accumarray(idj, Cup(:,k))./Ns;
        end
        Cup = round(Cup2(idj,:));
    end
end
%check for values lower than the minimum allowed
t_min = Vup < 0;
idx_nok_min = any(t_min,2);
if any(idx_nok_min)
    Vup(t_min) = 0;
    if nargout > 1
%         [~,~,idj] = unique(Vup,'rows','stable');
        [~,~,idj] = unique(Vup,'rows','first');
        idj = sort(idj);
        Ns =  accumarray(idj, ones(size(idx_nok_min)));
        Cup2 = zeros(size(Ns,1),3);
        for k = 1:3
            Cup2(:,k) = accumarray(idj, Cup(:,k))./Ns;
        end
        Cup = round(Cup2(idj,:));
    end
end
end %function