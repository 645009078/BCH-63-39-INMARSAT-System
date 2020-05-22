clear
clc

% add files in src  
addpath('../src/');

A = zeros(1, 63); A(9) = 2170; A(10) = 11718; A(11) = 32382; A(12) = 140322;
A(13) = 628866; A(14) = 2245950; A(15) = 7302603; A(16) = 21907809; A(17) = 60355638;
A(18) = 154242186; A(19) = 365056650;
A(20) = 803124630; A(21) = 1648195230; A(22) = 3146554530; A(23) = 5596735032; 
A(24) = 9327891720; A(25) = 14579965764; A(26) = 21309180732; A(27) = 29146649420; 
A(28) = 37474263540; A(29) = 45314900820; A(30) = 51356887596;
A(31) = 54561631635; A(32) = 54561631635; A(33) = 51356887596;
A(34) = 45314900820; A(35) = 37474263540; A(36) = 29146649420;
A(37) = 21309180732; A(38) = 14579965764; A(39) = 9327891720;
A(40) = 5596735032; A(41) = 3146554530; A(42) =1648195230;
A(43) = 803124630; A(44) = 365056650; A(45) = 154242186; A(46) = 60355638;
A(47) = 21907809; A(48) = 7302603; A(49) = 2245950; A(50) = 628866;
A(51) = 140322; A(52) = 32382; A(53) = 11718; A(54) = 2170; A(63) = 1 ;


gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
n = 63;
k = 39;
dmin = 9;
t = 4;
code = bch(n, k, dmin, t,gen_poly,prim_poly );

simulator = bch_simulation(code);
ps = 0.0025:0.005:0.07;
num_sym = 1000;

[bit_error_rates] = simulator.simulate_bsc_BER(num_sym, ps);

% theoretical estimations

% upper bound estimation
ber_upper_bound = zeros(size(ps));
for i = (t+1):n
    ber_upper_bound = ber_upper_bound + (t + i) * nchoosek(n, i) * ps.^(i) .* (1-ps).^(n-i);
end
ber_upper_bound = ber_upper_bound / n;

% ber approximation
ber_approx = zeros(size(ps));
for i = (t+1):n
    ber_approx = ber_approx  + i * nchoosek(n, i) * ps.^i .* (1 - ps).^(n-i);
end
ber_approx = ber_approx / n;

% ber approximation for p << 1
ber_low_p_approx = 1/n * (t+1) * nchoosek(63, (t+1)) * ps.^(t+1);

% plot ber vs theory
figure();
semilogy(ps, bit_error_rates);
hold on
semilogy(ps, ber_approx);
hold on
semilogy(ps, ber_upper_bound);
hold on
semilogy(ps, ber_low_p_approx);
set(gca,'yscale','log')
ylabel('Bit Error Probability');
xlabel('BSC Transition Probability');
grid on
legend('Simulated', 'BER Approximate', 'Upper Bound', 'BER Approximate (p << 1)');
title('BCH(63, 39) BSC BER');

% get prob undetected error 
[undetected_err_prob] = simulator.simulate_bsc_prob_undetected_error(num_sym, ps);

% approximate undetected error  probability 
approx_undetected_err_prob = 2170 * ps.^dmin;

% theoretical undetected error probability
undetected_err_prob_theoretical = zeros(1, size(ps, 2));
for i = 1:n
     undetected_err_prob_theoretical = undetected_err_prob_theoretical + A(i) * ps.^i .* (1-ps).^(63-i);
end

% plot probability undetected error estimation vs simulation
figure();
hold on
plot(ps, undetected_err_prob);
hold on
plot(ps, approx_undetected_err_prob);
hold on
plot(ps, undetected_err_prob_theoretical);
set(gca,'yscale','log')
ylabel('Undetected Error Probability');
xlabel('BSC Transition Probability');
grid on
legend('Simulated', 'Theoretical (Approximate)', 'Theoretical');
title('BCH(63, 39) BSC Undetected Error Probability');