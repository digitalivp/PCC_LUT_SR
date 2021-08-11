function ply_write(filename,V,C)
% Write the pointcloud to a ply file
% Usage: ply_write(filename,V,C)
% Input:
%   filename: name of ply file in current directory
%   V: Nx3 matrix of 3D coordinates of total N vertices
%   C: Nx3 matrix of R,G,B colors on the N vertices
%		If C is not provided, only the geometry
%		is written to file
%

fid = fopen(filename,'w');

fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'comment Nothing to comment\n');
fprintf(fid,'element vertex %d\n', length(V));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');

if nargin==3
	fprintf(fid,'property uchar red\n');
	fprintf(fid,'property uchar green\n');
	fprintf(fid,'property uchar blue\n');
end

fprintf(fid,'end_header\n');

if nargin==2
	fprintf(fid,'%f %f %f \n', V');
elseif nargin==3
	fprintf(fid,'%f %f %f %d %d %d \n', [double(V), double(C)]');
end

fclose(fid);
