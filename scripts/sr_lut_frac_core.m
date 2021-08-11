function sr_lut_frac_core(filename,i_r,end_path)

[Vdec,Cdec] = read_ply([end_path '/' filename '.ply.bin.decoded.ply']);

gp = nextpow2(max(max(Vdec)));
s = 1/TMC13_positionQuantizationScale(i_r, gp);
##T = min(Vdec);
T = zeros(1,3);
Vds = round((Vdec - T)/s);

if s <= 2
    Vlut = LUT_SR_fractional_octave(Vds,Cdec,s);
    Clut = round(pointcloud_interpolation_2(round(Vdec), Cdec, Vlut, 1));
    Vlut = Vlut + T;
elseif s <= 4
    Vlut2 = LUT_SR_fractional_octave(Vds,Cdec,s/2);
    Clut2 = pointcloud_interpolation_2(round(Vdec/2), Cdec, Vlut2, 1);
    Vlut  = LUT_SR_fractional_octave(Vlut2, Clut2,s/2);
    Clut  = round(pointcloud_interpolation_2(round(Vlut2*2), Clut2, Vlut, 1));
    Vlut  = Vlut + T;
elseif s <= 8
    Vlut4 = LUT_SR_fractional_octave(Vds,Cdec,s/4);
    Clut4 = pointcloud_interpolation_2(round(Vdec/4), Cdec, Vlut4, 1);
    Vlut2 = LUT_SR_fractional_octave(Vlut4, Clut4,s/4);
    Clut2 = pointcloud_interpolation_2(round(Vlut4*2), Clut4, Vlut2, 1);
    Vlut  = LUT_SR_fractional_octave(Vlut2, Clut2,s/4);
    Clut  = round(pointcloud_interpolation_2(round(Vlut2*2), Clut2, Vlut, 1));
    Vlut  = Vlut + T;
elseif s <= 16
    Vlut8 = LUT_SR_fractional_octave(Vds,Cdec,s/8);
    Clut8 = pointcloud_interpolation_2(round(Vdec/8), Cdec, Vlut8, 1);
    Vlut4 = LUT_SR_fractional_octave(Vlut8,Clut8,s/8);
    Clut4 = pointcloud_interpolation_2(round(Vlut8*2), Clut8, Vlut4, 1);
    Vlut2 = LUT_SR_fractional_octave(Vlut4, Clut4,s/8);
    Clut2 = pointcloud_interpolation_2(round(Vlut4*2), Clut4, Vlut2, 1);
    Vlut  = LUT_SR_fractional_octave(Vlut2, Clut2,s/8);
    Clut  = round(pointcloud_interpolation_2(round(Vlut2*2), Clut2, Vlut, 1));
    Vlut  = Vlut + T;
end

out_name = [end_path '/' filename '.ply.bin.decoded.lutSR.ply'];

d=dir();

for k=1:length(d)
  if strcmp(d(k).name,out_name)==1
    delete(out_name)
    disp('file rewrite')
  endif
end

ply_write(out_name,Vlut,Clut)

