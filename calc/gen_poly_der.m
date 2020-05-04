% use minpoly for actual derivation, this is just for removing the grunt
% work from the typed derivation

clear
clc

syms a x

min_poly_1 = (x + a^3)*(x + a^6)*(x+a^12)*(x+a^24)*(x+a^48)*(x+a^96);
min_poly_1 = collect(expand(min_poly_1));

pretty(min_poly_1);