clear
clc

gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
code = bch(63,39,9,4,gen_poly);


% simulation
% -----------------------------------------------------------------------------------
% 
% sim = bch_simulation(code);
% 
% ps = 0:0.01:0.5;
% % ps = ps * 10e-3
% num_sym = 10;
% sim.simulate_bsc(num_sym, ps);

% 
% 
% % SNR = 2:1:3; % not in db
% SNR = 1.3:0.01:3.5; % not in db
% 
% mod_order = 4;
% num_sym = 10000;
% sim.simulate_awgn(num_sym, mod_order, SNR);
% 
% 


% ---------------------------------------------------------------------------------------------------















































% % trying to use bsc 
% 
% rand_msgs = randi([0 1],64,39); % Random matrix
% rand_msgs_decoded_my_implementation = zeros(size(rand_msgs));
% rand_msgs_decoded_their_implementation = zeros(size(rand_msgs));
% 
% for n = 1:size(rand_msgs, 1)
%     cp = bchenc(gf(rand_msgs(n, :)), 63, 39);
%     rp = bsc(double(cp.x), 0.1);
%     x = bchdec(gf(rp), 63, 39);
%     rand_msgs_decoded_their_implementation(n, :) = double(x.x);
%     
%     code_poly = code.encode(rand_msgs(n, :));
%     r_poly = bsc(code_poly, 0.1);
%     rand_msgs_decoded_my_implementation(n, :) = code.decode(r_poly);
% end
% 
% [num_errs_mine, pct_errs_mine] = biterr(rand_msgs, rand_msgs_decoded_my_implementation)
% 
% [num_errs_theirs, pct_errs_theirs] = biterr(rand_msgs, rand_msgs_decoded_their_implementation)
% 
% % rank against matlabs 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 






% 
% gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
% code = bch(63,39,9,4,gen_poly);
% 
% message_sample = zeros(1,39);
% message_sample(1:2:end) = 1;
% 
% code_poly = code.encode(message_sample);
% % introduce errors 
% code_poly(end-1) = 1;
% % code_poly(end-3) = 1;
% code_poly(end-4) = 0;
% code_poly(end-9) = 1;
% code_poly(end-15) = 1; % more than 4 errors 
% 
% 
% decoded_msg = code.decode(code_poly);
% 
% message_sample == decoded_msg
% 



%   
% % code.syndrome_indexes(gen_poly);
% % 
% % 
% % 
% % res = code.gf_poly_eval([1,1,0,0,0,0,1], 3);
% % 
% % 
% % syndrome = code.syndrome_indexes(gen_poly); % should be all 0
% % 
% % % example 
% % 
% % message = zeros(1,39);
% % message(1:20) = 1;
% % message
% % 
% % encode 
% codepoly = code.encode(message)
% % 
% % check systematic 
% codepoly(25:end) == message
% % 
% % show parity
% disp("parity bits");
% codepoly(1:24)
% % 
% % 
% % syndrome 
% msg_synd = code.syndrome_indexes(codepoly);
% disp("message syndrome, should be all 1's (0)");
% msg_synd
% 
% 
% adding an error and showing syndrome 
% disp("code poly with 1 error in position 22, syndrome");
% code_poly(22) = 1;
% inv_syndrome = code.syndrome_indexes(code_poly);
% inv_syndrome
% 
% 
% 
% 
% 
% 
% 
% 


% 
% 
% % checking, example 6.1 and 6.4 from textbook
% gen_61 = [1,0,0,0,1,0,1,1,1];
% code_61 = bch(15, 7, 5, 2, gen_61);
% syndrome_example_64 = code_61.syndrome_indexes([1,0,0,0,0,0, 0,0,1,0,0,0,0,0,0]);
% code_61.decode([1,0,0,0,0,0, 0,0,1,0,0,0,0,1,0]);
% % result is [alpha2, alpha4, alpha7, alpha8], shown in example 6.4
% % 
% % 
% % % testing power set indicies 
% % ps_inds = code.gf_power_set_indexes(3, 6); % 2 --> alpha^2 
% % % result should be [2,5,8,11,14,17,20]
% % ps_inds




%  m = [1,0,1,1];
%  code_poly = code.encode(m);
%  gfpretty(flip(code_poly));
  
  

% making stuff work -----------------------------------------
% px = [1,1,0,0,1];
% [rt, rt_tuple] = gfroots(px, 4,2);
% gfpretty(px);
% 
% p = 2;
% m = 4;
% field = gftuple([-1:p^m-2]',m,p);
% 
% Show_Poly(field);
% 
% 
% 
% function Show_Poly(Field)
%     Power=arrayfun(@(x) ['s^',num2str(x)],(1:length(Field))'-2,'UniformOutput',false);
%     Tuple=mat2cell(Field,ones(1,size(Field,1)),size(Field,2));
%     Poly=cellfun(@(x) poly2str(fliplr(x),'s'),Tuple,'UniformOutput',false);
%     disp([Power,Poly,cellfun(@(x) num2str(x),Tuple,'UniformOutput',false)])
% end