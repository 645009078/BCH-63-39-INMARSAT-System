clear
clc

A = zeros(1, 63);
A(9) = 2170;
A(10) = 11718;
A(11) = 32382;
A(12) = 140322;
A(13) = 628866;
A(14) = 2245950;
A(15) = 7302603;
A(16) = 21907809;
A(17) = 60355638;
A(18) = 154242186;
A(19) = 365056650;
A(20) = 803124630;
A(21) = 1648195230;
A(22) = 3146554530;
A(23) = 5596735032;
A(24) = 9327891720;
A(25) = 14579965764;
A(26) = 21309180732;
A(27) = 29146649420;
A(28) = 37474263540;
A(29) = 45314900820;
A(30) = 51356887596;
A(31) = 54561631635;
A(32) = 54561631635;
A(33) = 51356887596;
A(34) = 45314900820;
A(35) = 37474263540;
A(36) = 29146649420;
A(37) = 21309180732;
A(38) = 14579965764;
A(39) = 9327891720;
A(40) = 5596735032;
A(41) = 3146554530;
A(42) =1648195230;
A(43) = 803124630;
A(44) = 365056650;
A(45) = 154242186;
A(46) = 60355638;
A(47) = 21907809;
A(48) = 7302603;
A(49) = 2245950;
A(50) = 628866;
A(51) = 140322;
A(52) = 32382;
A(53) = 11718;
A(54) = 2170;
A(63) = 1 ;


syms p
p_undetected_error = 0;
for i = 1:63
     p_undetected_error = p_undetected_error + A(i) * p^i * (1-p)^(63-i);
end

dmin = 9;
p_estimated_undetected_error = 2170 * p^dmin;


ps = 0:0.001:0.01;
ps_errs = double(subs(p_undetected_error, p, ps));
ps_est_errs = double(subs(p_estimated_undetected_error, p, ps));
table(ps', ps_errs', ps_est_errs')