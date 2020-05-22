classdef bch_simulation 
    properties (Access = private)
       code; % bch code
    end
    methods
        function obj = bch_simulation(bch_code)
            obj.code = bch_code;
        end
        function [mod, demod] = init_psk_mod_and_demod(obj, modulation_order)
            % initialise psk modulator 
            mod = comm.PSKModulator();
            mod.ModulationOrder = modulation_order;
            mod.BitInput=true;
            mod.PhaseOffset = pi/4; 
            % initialise psk demodulator 
            demod=comm.PSKDemodulator();
            demod.ModulationOrder = modulation_order;
            demod.BitOutput = true;
            demod.PhaseOffset = pi/4; 
            return;
        end
        function [bit_error_rates] = simulate_bsc_BER(obj, num_sym, ps)
            % initialise return arrays
            bit_error_rates = zeros(1, length(ps));
            % simulate BER 
            for i = 1:length(ps)
                fprintf("simulating transition probability %d / %d\n", i, length(ps));
                n_errs = 0;
                n_bits = 0;
                while n_errs <= 100 && n_bits <= 1e7
                    % generate random msg frame 
                    msgs = randi([0 1],num_sym,obj.code.k); 
                    % send msgs over bsc 
                    for n = 1:num_sym
                        % encode row
                        codeword = obj.code.encode(msgs(n, :));
                        % send over binary symmetric channel
                        r_codeword = bsc(codeword, ps(i));
                        r_msg = obj.code.decode(r_codeword);
                        err_count = sum(r_msg ~= msgs(n, :));                        
                        % record detected codeword
                        if (err_count ~= 0)
                            n_errs = n_errs + err_count;
                        end 
                        n_bits = n_bits + obj.code.k;
                    end
                end
                bit_error_rates(i) = n_errs / n_bits;
            end
            return;
        end
        function [coded_bit_error_rates, uncoded_bit_error_rates] = simulate_awgn_BER(obj, num_sym, mod_order, SNR_norm_dB)
            % initialise psk modulator and demodulator
            [psk_mod, psk_demod] = obj.init_psk_mod_and_demod(mod_order);
            % initialise return vectors
            coded_bit_error_rates = zeros(1, length(SNR_norm_dB));
            uncoded_bit_error_rates = zeros(1, length(SNR_norm_dB));
            % generate channel 
            awgn = comm.AWGNChannel();
            awgn.NoiseMethod = 'Signal to noise ratio (Eb/No)';
            awgn.BitsPerSymbol = log2(mod_order);
            % simulate each SNR
            for r = 1:numel(SNR_norm_dB)
                fprintf("simulating SNR %d / %d\n", r, length(SNR_norm_dB));
                n_coded_errs = 0;
                n_uncoded_errs = 0;
                n_bits = 0;
                awgn.EbNo = SNR_norm_dB(r);
                while n_coded_errs < 100 && n_bits < 1e7
                    % generate random message frame 
                    msgs = randi([0 1],num_sym,obj.code.k); 
                    % send each msg over awgn channel
                    for n = 1:num_sym
                        % encode msg
                        code_poly = obj.code.encode(msgs(n, :));
                        % pad with 0 for psk mod 
                        code_poly(ceil(obj.code.n / mod_order) * mod_order) = 0;
                        % psk mod
                        tx = psk_mod(code_poly');
                        % transmit over awgn channel 
                        rx = awgn(tx);
                        % psk demod
                        r_code_poly = psk_demod(rx);
                        % remove padding bits
                        r_code_poly = r_code_poly(1:obj.code.n);
                        % decode bch poly
                        r_msg = obj.code.decode(r_code_poly');
                        % count errors 
                        n_bits = n_bits + obj.code.k;
                        n_coded_errs = n_coded_errs + sum(r_msg ~= msgs(n, :));
                        n_uncoded_errs = n_uncoded_errs + sum(r_code_poly(end-obj.code.k+1:end)' ~= msgs(n, :));
                    end
                end
                coded_bit_error_rates(r) = n_coded_errs / n_bits;
                uncoded_bit_error_rates(r) = n_uncoded_errs / n_bits;
            end
            return;
        end
        function [prob_undetected_errors] = simulate_bsc_prob_undetected_error(obj, num_sym, ps)
            % initialise return arrays
            prob_undetected_errors = zeros(1, length(ps));
            % simulate BER 
            for i = 1:length(ps)
                fprintf("simulating transition probability %d / %d\n", i, length(ps));
                n_msgs = 0;
                n_undetected = 0;
                while n_undetected <= 1
                    % generate random msg frame 
                    msgs = randi([0 1],num_sym,obj.code.k); 
                    % send msgs over bsc 
                    for n = 1:num_sym
                        % encode row
                        codeword = obj.code.encode(msgs(n, :));
                        % send over binary symmetric channel
                        r_codeword = bsc(codeword, ps(i));
                        % calc syndrome
                        syndrome = obj.code.calculate_syndrome(r_codeword);
                        if (all(syndrome == -Inf)) 
                            if (msgs(n, :) ~= r_codeword(end-obj.code.k+1: end))
                                n_undetected = n_undetected + 1;
                            end
                        end
                        n_msgs = n_msgs + 1;
                    end
                end
                prob_undetected_errors(i) = n_undetected / n_msgs;
            end
            return;
        end
    end
end





