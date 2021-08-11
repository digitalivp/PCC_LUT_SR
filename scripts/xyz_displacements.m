function deltas = xyz_displacements(interval)

[x y z] = meshgrid(interval, ...
					interval, ...
					interval);

deltas = sortrows([x(:) y(:) z(:)]);
% deltas = deltas(end:-1:1,:);