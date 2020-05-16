classdef bch_simulation 
    properties (Access = private)
       code; % bch code
    end
    methods
        function obj = bch_simulation(bch_code)
            obj.code = bch_code;
        end
        function [enc, dec] = init_psk_mod_and_demod(obj, modulation_order)
            % initialise psk modulator 
            enc = comm.PSKModulator();
            enc.ModulationOrder = modulation_order;
            enc.BitInput=true;
            enc.PhaseOffset = pi/4; 
            % initialise psk demodulator 
            dec=comm.PSKDemodulator();
            dec.ModulationOrder = modulation_order;
            dec.BitOutput = true;
            dec.PhaseOffset = pi/4; 
        end
        function [bit_error_rates, prob_undetected_errors] = simulate_bsc(obj, num_sym, ps)
            % initialise return arrays
            bit_error_rates = zeros(1, length(ps));
            prob_undetected_errors = zeros(1, length(ps));
            % simulate BER 
            for i = 1:length(ps)
                fprintf("simulating transition probability %d / %d\n", i, length(ps));
                n_errs = 0;
                n_bits = 0;
                n_undetected = 0;
                while n_errs <= 100 || n_undetected <= 1
                    % generate random msg frame 
                    msgs = randi([0 1],num_sym,obj.code.k); 
                    % send msgs over bsc 
                    for n = 1:num_sym
                        % encode row
                        codeword = obj.code.encode(msgs(n, :));
                        % send over binary symmetric channel
                        r_codeword = bsc(codeword, ps(i));
                        [r_msg, no_error_detected] = obj.code.decode(r_codeword);
                        err_count = sum(r_msg ~= msgs(n, :));
                        % check for undetected errors  
                        if (no_error_detected == 1 && err_count ~= 0)
                            n_undetected = n_undetected + 1;
                        end
                        % record detected codeword
                        if (err_count ~= 0)
                            n_errs = n_errs + err_count;
                        end 
                        n_bits = n_bits + obj.code.k;
                    end
                end
                bit_error_rates(i) = n_errs / n_bits;
                prob_undetected_errors(i) = n_undetected / (n_bits / obj.code.k);
            end
            return;
        end
        function [coded_bit_error_rates, uncoded_bit_error_rates] = simulate_awgn(obj, num_sym, mod_order, SNR)
            % initialise psk modulator and demodulator
            [psk_enc, psk_dec] = obj.init_psk_mod_and_demod(mod_order);
            % initialise return vectors
            coded_bit_error_rates = zeros(1, length(SNR));
            uncoded_bit_error_rates = zeros(1, length(SNR));
            % convert SNR to dB
            SNR_dB = 10*log(SNR);
            % generate channel 
            awgn = comm.AWGNChannel();
            awgn.NoiseMethod = 'Signal to noise ratio (Eb/No)';
            awgn.BitsPerSymbol = log2(mod_order);
            % simulate each SNR
            for r = 1:numel(SNR_dB)
                fprintf("simulating SNR %d / %d\n", r, length(SNR));
                n_coded_errs = 0;
                n_uncoded_errs = 0;
                n_bits = 0;
                awgn.EbNo = SNR_dB(r);
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
                        tx = psk_enc(code_poly');
                        % transmit over awgn channel 
                        rx = awgn(tx);
                        % psk demod
                        r_code_poly = psk_dec(rx);
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
    end
end





