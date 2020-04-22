% Assuming that QPSK is used as the modulation scheme, simulate the coding system
% in a BSC channel. Plot the BER versus p curve. Please note that you need to
% choose the appropriate values of p , such that the BERs at least cover the range
% from 10  5 to 10  1 . 

% Plot the undetected error probability with the same range of p
% values. 

% Compare the simulation results with your estimation produced in 2) and 3)
% and comment on the accuracy of your estimations and justify.

clear
clc


gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
prim_poly = [1,1,0,0,0,0,1];
n = 63;
k = 39;
dmin = 9;
t = 4;
code = bch(n, k, dmin, t,gen_poly,prim_poly );


simulator = bch_simulation(code);
% ps = 0.01:0.0025:0.04;
ps = 0.008:0.0025:0.02;
num_sym = 5000;

[bit_error_rates, probability_of_undetected] = simulator.simulate_bsc(num_sym, ps)


% theoretical estimations
prob_undetected_theoretical = 2170 * ps.^dmin;
ber_estimation_theoretical = 1/n * (t+1) * nchoosek(63, (t+1)) * ps.^(t+1);


% plot ber vs theory
figure();
semilogy(ps, bit_error_rates);
hold on
semilogy(ps, ber_estimation_theoretical);
set(gca,'yscale','log')
ylabel('Bit Error Probability');
xlabel('BSC Transition Probability');
grid on
legend('Simulated', 'Theoretical');


% plot probability undetected error estimation vs simulation
figure();
plot(ps, prob_undetected_theoretical);
hold on
plot(ps, probability_of_undetected);
set(gca,'yscale','log')
ylabel('Undetected Error Probability');
xlabel('BSC Transition Probability');
grid on
legend('Simulated', 'Theoretical');



