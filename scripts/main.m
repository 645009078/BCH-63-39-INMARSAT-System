clear
clc

% add files in src  
addpath('../src/');

gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
code = bch(63,39,9,4,gen_poly,prim_poly );

msg = zeros(1, 39);
msg(1:2:39) = 1;

code_poly = code.encode(msg);

code_poly_altered = code_poly;
code_poly_altered(38) = 1;
code_poly_altered(1:2) = 0;
code_poly_altered(12) = 1;

r_msg = code.decode(code_poly_altered);


code_poly_altered ~= code_poly
r_msg ~= msg
