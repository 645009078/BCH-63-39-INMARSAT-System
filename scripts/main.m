clear
clc

% add files in src  
addpath('../src/');

gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
code = bch(63,39,9,4,gen_poly,prim_poly );
code2 = bch2(63,39,9,4,gen_poly,prim_poly );

msg = randi([0 1],1,39); 

code_poly = code.encode(msg);
                    err = randi([0 1],1,63); 

code_poly_altered = mod(code_poly + err, 2);
% code_poly_altered = code_poly;
% 
% 
s1 = code.calculate_syndrome(code_poly_altered);
s2 = code2.syndrome_indexes(code_poly_altered);
s1 
s2
% 
% r_msg = code.decode(code_poly_altered);
% r_msg2 = code2.decode(code_poly_altered);
% 
% code_poly_altered ~= code_poly;
% r_msg ~= msg
% r_msg2 ~= msg
% r_msg ~= r_msg2
