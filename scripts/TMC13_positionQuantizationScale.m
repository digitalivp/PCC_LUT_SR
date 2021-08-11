function pqs = TMC13_positionQuantizationScale(rp, gp)
rp = 6-rp;
p_min = max(gp - 9, 7);
start = min(1, gp - (p_min + 6));
step = max(1, (min(gp - 1, p_min + 7) - p_min) / 5);
y = start + round(rp * step);
div = 2 ^ (abs(y) + 1);
pqs = mod((1 - 2*(y < 0)), div) / div;
end
         
         