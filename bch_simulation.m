classdef bch_simulation 
    properties (Access = private)
       code; % bch code
    end
    methods
        function obj = bch_simulation(bch_code)
            obj.code = bch_code;
        end
        function [enc, dec] = init_psk_mod_and_demod(obj, modulation_order)
            enc = comm.PSKModulator();
            enc.ModulationOrder = modulation_order;
            enc.BitInput=true;
            enc.PhaseOffset = pi/4; % why? ----------------------------------------------------------------

            dec=comm.PSKDemodulator();
            dec.ModulationOrder = modulation_order;;
            dec.BitOutput = true;
            dec.PhaseOffset = pi/4; % why -------------------------------------------------------------
        end
        function [bit_error_rates, prob_undetected_errors] = simulate_bsc(obj, num_sym, ps)
            msgs = randi([0 1],num_sym,39); % Random matrix
            dec_msgs = zeros(size(msgs));
            bit_error_rates = zeros(1, length(ps));
            prob_undetected_errors = zeros(1, length(ps));
            for i = 1:length(ps)
                i
                num_undetected = 0;
                for n = 1:num_sym
                    code_poly = obj.code.encode(msgs(n, :));
                    r_code = bsc(code_poly, ps(i));
                    [r_msg, num_errors_corrected] = obj.code.decode(r_code);
                    if (num_errors_corrected ~= 0 && sum(r_msg ~= msgs(n, :)) ~= 0)
                        num_undetected = num_undetected + 1;
                    end 
                    dec_msgs(n, :) = r_msg;                    
                end
                [~, ber] = biterr(msgs, dec_msgs);
                bit_error_rates(i) = ber;
                prob_undetected_errors(i) = num_undetected / num_sym;
                % restore 
                dec_msgs(:) = 0;
            end
            return;
        end
        function [coded_bit_error_rates, uncoded_bit_error_rates] = simulate_awgn(obj, num_sym, mod_order, SNR)
            [psk_enc, psk_dec] = obj.init_psk_mod_and_demod(mod_order);
            msgs = randi([0 1],num_sym,39); % Random matrix
            dec_msgs = zeros(size(msgs));
            uncoded_msgs = zeros(size(msgs));
            coded_bit_error_rates = zeros(1, length(SNR));
            uncoded_bit_error_rates = zeros(1, length(SNR));
            SNR_db = 10*log(SNR);
            for r = 1:numel(SNR_db)
                r
                for n = 1:num_sym
                    % encode 
                    code_poly = obj.code.encode(msgs(n, :));
                    % pad to make code poly integer multiple of m length
                    code_poly(end+1) = 0; 
                    % encode code poly with psk
                    tx = psk_enc(code_poly');
                    % transmit over awgn channel
                    rx_awgn = awgn(tx, SNR_db(r), 'measured');
                    % decode psk
                    r_code = psk_dec(rx_awgn);
                    % remove padding bits
                    r_code = r_code(1:end-1); 
                    % decode recieved code poly
                    r_msg = obj.code.decode(r_code');
                    % save recieved message
                    dec_msgs(n, :) = r_msg; 
                    % save result of uncoded message
                    % code is systematic, final k bits of recieved code
                    % poly
                    uncoded_msgs(n, :) = r_code(end-obj.code.k+1: end);
                end
                [~, ber_coded] = biterr(msgs, dec_msgs);
                [~, ber_uncoded] = biterr(msgs, uncoded_msgs);
                coded_bit_error_rates(r) = ber_coded;
                uncoded_bit_error_rates(r) = ber_uncoded;
                % restore 
                dec_msgs(:) = 0;
            end
            return;
        end 
    end
end























































% 
%  gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
% 
%             message_sample = zeros(1,39);
%             message_sample(1:2:end) = 1;
% 
%             code_poly = obj.code.encode(message_sample);
%             % introduce errors 
%             code_poly(end-1) = 1;
%             code_poly(end-3) = 1;
%             code_poly(end-4) = 0;
% %             code_poly(end-9) = 1;
%             code_poly(end-15) = 1; % more than 4 errors 
% 
% 
%             decoded_msg = obj.code.decode(code_poly);
% 
%             message_sample == decoded_msg