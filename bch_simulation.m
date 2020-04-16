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
        function simulate_bsc(obj, num_sym, mod_order, ps)
            [psk_enc, psk_dec] = obj.init_psk_mod_and_demod(mod_order);
            msgs = randi([0 1],num_sym,39); % Random matrix
            dec_msgs = zeros(size(msgs));
            calc_bers = zeros(1, length(ps));
            for i = 1:length(ps)
                i
                for n = 1:num_sym
                    code_poly = obj.code.encode(msgs(n, :));
%                     code_poly(end+1) = 0; % pad
%                     rx = psk_enc(code_poly');
                    r_code = bsc(code_poly, ps(i));
%                     r_code = psk_dec(rx);
%                     r_code = r_code(1:end-1)'; % remove pad bit and transpose to row vec
                    r_msg = obj.code.decode(r_code);
                    dec_msgs(n, :) = r_msg;                    
                end
                [~, ber] = biterr(msgs, dec_msgs);
                calc_bers(i) = ber;
                % restore 
                dec_msgs(:) = 0;
            end
            calc_bers
            % plot result
            figure();
            semilogy(ps, calc_bers);
            
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