function C_new = pointcloud_interpolation_2(V, C, V_new, xyz_max_dist)
% V and V_new must have the same resolution (depth)
% xyz_max_dist max neighbourhood to considerer
[old_pos, old_pos_org] = ismember(V_new,V,'rows');
new_pos = ~old_pos;
V_created = V_new(new_pos,:); % V_created has only the points from V_new that
                              % don't exist in V

deltas = xyz_displacements(-xyz_max_dist:xyz_max_dist);
dists = sqrt(sum(deltas.^2,2));
val_deltas = find(dists>0);
deltas = deltas(val_deltas,:);
dists = dists(val_deltas,:);
%create the neighbourhood to be searched (deltas) and its repective distances
%to the center point
%%
C_created = V_created*0;
C_weights = C_created;
for k = 1:length(dists)
	[~,delta_ok,created_ok] = intersect(bsxfun(@plus,V,deltas(k,:)),V_created,'rows');
    %check if the points in the V_created are in the max displacement
    %allowed
	C_created(created_ok,:) = C_created(created_ok,:) + C(delta_ok,:)/dists(k);
    %for created points within the max displacement, their color will be
    %the sum of the color of each neighbour point multiplied by the inverse of the displacement distance
	C_weights(created_ok,:) = C_weights(created_ok,:) + 1/dists(k);
    %the weight is the inverse of the displacement distance
end
C_created = bsxfun(@rdivide, C_created, C_weights); % take the weighted average mean color
C_new = V_new*0;
C_new(old_pos,:) = C(old_pos_org(old_pos_org>0),:);
C_new(new_pos,:) = C_created;
end