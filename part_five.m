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
code = bch(63,39,9,4,gen_poly);


simulator = bch_simulation(code);
ps = 0.005:0.00005:0.037;
num_sym = 50000;
[bit_error_rates, probability_of_undetected] = simulator.simulate_bsc(num_sym, ps);

% theoretical estimations
prob_undetected_theoretical = 2170 * ps.^9;
ber_estimation_theoretical = 1/63 * 5 * nchoosek(63, 5) * ps.^5 .* (1 - ps).^5;

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


% 
% 
% % plotting approximated probability of undetectecd error
% prob_undetected_theoretical = 2170 * ps.^9;
% figure();
% plot(ps, prob_undetected_theoretical);
% hold on
% ber_estimation_theoretical = 1/63 * 5 * nchoosek(63, 5) * ps.^5;
% plot(ps, ber_estimation_theoretical);
% 


