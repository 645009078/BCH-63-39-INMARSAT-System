
clear
clc

% add files in src  
addpath('../src/');
 
gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
n = 63;
k = 39;
dmin = 9;
t = 4;
code = bch(n,k,dmin,t,gen_poly,prim_poly);
simulator = bch_simulation(code);

% normalise SNR
mod_order = 4;
SNR_norm_dB = 0:10;
num_sym = 200;
[coded_bit_error_rates, uncoded_bit_error_rates] = simulator.simulate_awgn_BER(num_sym, mod_order, SNR_norm_dB)

figure();
semilogy(SNR_norm_dB, coded_bit_error_rates);
hold on
semilogy(SNR_norm_dB, uncoded_bit_error_rates);
xlabel('E_b/N_o');
ylabel('BER');
grid on 
set(gca,'yscale','log');
legend('Coded', 'Uncoded');
title('BCH(63, 39) QPSK AWGN BER');