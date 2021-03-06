clear
clc


gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
n = 63;
k = 39;
dmin = 9;
t = 4;
code = bch(n,k,dmin,t,gen_poly,prim_poly);
simulator = bch_simulation(code);

SNR = 1:0.1:2.6; % not in db
mod_order = 4;
num_sym = 200;
[coded_bit_error_rates, uncoded_bit_error_rates] = simulator.simulate_awgn(num_sym, mod_order, SNR)

% normalise SNR
SNR_norm = SNR / ((k/n) * log2(mod_order));
SNR_norm_db = 10*log(SNR_norm);

figure();
semilogy(SNR_norm_db, coded_bit_error_rates);
hold on
semilogy(SNR_norm_db, uncoded_bit_error_rates);
xlabel('E_b/N_o');
ylabel('BER');
grid on 
set(gca,'yscale','log');
legend('Coded', 'Uncoded');
