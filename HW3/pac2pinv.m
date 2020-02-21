function pinv = pac2pinv(pac)
pac = pac / 1000;
if pac >= 10 && pac <= 100
    pinv = 0.556*pac - 145.556;
elseif pac > 100 && pac <= 500
    pinv = 97.5 - 0.075 * pac;
else
    pinv = NaN;
end
end