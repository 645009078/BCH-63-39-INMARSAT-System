classdef bch
    properties
        n; % code length
        k; % message length
        m; % prim poly degree
        d; % minimum distance
        t; % error correcting capability
        gen_poly;  % ascending power
        gen_poly_deg; % degree of gen poly
        gf_ext; % 2^m x (m-1) matrix, rows are elements of the extended gf in polynomial form
        prim_poly;
    end
    methods
        function obj = bch(n, k,d,t,gen_poly, prim_poly)
            arguments
                n {mustBeInteger}
                k {mustBeInteger} 
                d {mustBeInteger} 
                t {mustBeInteger} 
                gen_poly(1,:) {mustBeLessThanOrEqual(gen_poly, 1)} % ascending power 
                prim_poly(1,:) {mustBeLessThanOrEqual(prim_poly,1)} % ascending power 
            end
            obj.n = n;
            obj.k = k;
            obj.m = log2(n+1);
            obj.d = d;
            obj.t = t;
            obj.gen_poly = gen_poly;
            obj.prim_poly = prim_poly;
            obj.gf_ext = gftuple([-1:2^obj.m-2]', prim_poly); 
        end 
        function code_polynomial = encode(obj, msg)
            arguments
                obj
                msg(1,:) {mustBeLessThanOrEqual(msg,1)} 
            end
            msg_shifted = [zeros(1, obj.n -obj.k), msg]; % multiply by x^(n-k), increasing power
            [~, rem] = gfdeconv(msg_shifted, obj.gen_poly); % get remainder
            code_polynomial = gfadd(msg_shifted, rem, 2); 
        end
        function syndrome = calculate_syndrome(obj, poly_coeffs)
            arguments
                obj
                poly_coeffs (1, :) {mustBeLessThanOrEqual(poly_coeffs, 1), mustBeInteger}
            end
            % elements [alpha, ...., alpha^2t] of gf(2^m)
            exponents =  1:2*obj.t;
            % convert poly to ext gf
            poly_gf = zeros(size(poly_coeffs));
            poly_gf(poly_coeffs == 0) = -Inf;
            % evaluate polynomial for all elements 
            syndrome = obj.gf_ext_poly_eval(poly_gf,exponents)';
            return;
        end
        function bia_tbl = init_simplified_bia_table(obj)
            % t + 1 rows, for range [-1, t]
            % 5 columns for mu, sigma fn, d, l, and mu - l
            bia_tbl = cell(obj.t + 2, 5);
            bia_tbl{1,1} = -1/2; % mu
            bia_tbl{1,2} = [0]; % sigma(x) coefficients, increasing order, exponential form
            bia_tbl{1,3} = 0; % d
            bia_tbl{1,4} = 0; % l
            bia_tbl{1,5} = -1; % 2mu - l
            bia_tbl{2,1} = 0; % mu
            bia_tbl{2,2} = [0]; % sigma(x)
            bia_tbl{2,4} = 0; % l
            bia_tbl{2,5} = 0; %2mu - l
            return;
        end
        function d_mu = calculate_simplified_discrepency(obj, bia_tbl, syndrome, mu)
            arguments
                obj
                bia_tbl
                syndrome (1,:) {mustBeNumeric}
                mu {mustBeNumeric}
            end 
            sigma_mu = bia_tbl{mu+3, 2};
            l_mu = bia_tbl{mu+3, 4};
            d_mu = syndrome(2*mu+3);
            for i = 1:l_mu
                term_i = gfmul(syndrome(2*mu+3-i), sigma_mu(i+1), obj.gf_ext);
                d_mu = gfadd(term_i, d_mu, obj.gf_ext);
            end
            return;
        end 
        function rho = calculate_rho(obj, bia_tbl, max_rho)
            arguments 
                obj
                bia_tbl
                max_rho {mustBeInteger}
            end
            rho = -1/2;
            val = -1; % 2rho - l
            for mu = max_rho:-1:0
                if bia_tbl{mu+2, 3} ~= -Inf && bia_tbl{mu+2, 5} > val
                    rho = mu;
                    val = bia_tbl{mu+2, 5};
                end
            end
            return;
        end
        function sigma_x = simplified_berlekamp_massey(obj, syndrome) 
            arguments
                obj
                syndrome (1,:) {mustBeNumeric} % syndrome in exp form 
            end
            % initialise 2t + 2 x 5 table for iterative algorithm'
            % rows represent interations over [-1, 2t]
            bia_tbl = obj.init_simplified_bia_table();
            % initialise default syndrome 
            bia_tbl{2,3} = syndrome(1);
            % iterate
            for mu = 0:obj.t-1
                % calc discrepency 
                bia_tbl{mu+2, 3} = obj.calculate_simplified_discrepency(bia_tbl, syndrome, mu-1);

                d_mu = bia_tbl{mu+2, 3};
                % carry sigma to iteration m + 1
                bia_tbl{mu+3, 2} = bia_tbl{mu+2, 2};                
                % non-zero discrepency
                if d_mu ~= -Inf
                    rho = obj.calculate_rho(bia_tbl, mu-1);
                    % element in extended gf
                    d_rho = bia_tbl{floor(rho)+2,3};
                    % polynomial over extended gf
                    sigma_rho = bia_tbl{floor(rho)+2, 2};
                    % amount of shift 
                    x_deg = 2*(mu-rho);
                    % multiply sigma_rho by x_deg
                    tmp = [-Inf + zeros(1, x_deg), sigma_rho];
                    % get inverse of d_rho
                    d_rho_inv = gfdiv(0, d_rho, obj.gf_ext);
                    % multiply to get coeff
                    coeff = gfmul(d_mu, d_rho_inv, obj.gf_ext);
                    % multiply tmp by coeff
                    coeff_resized = zeros(size(tmp)) + coeff;
                    tmp = gfmul(tmp, coeff_resized, obj.gf_ext);
                    % add discrepency
                    bia_tbl{mu+3, 2} = gfadd(bia_tbl{mu+3, 2}, tmp, obj.gf_ext);
                end
                % fill table for mu + 1
                bia_tbl{mu+3, 1} = mu + 1;
                bia_tbl{mu+3, 4} = length(bia_tbl{mu+3, 2}) -1;
                bia_tbl{mu+3, 5} = 2*bia_tbl{mu+3, 1} - bia_tbl{mu+3, 4};
            end
            sigma_x = bia_tbl{obj.t+2, 2};
            return;
        end
        function power_set = gf_power_set(obj, exponents, max_pow)
            arguments
                obj
                exponents (1, :) {mustBeInteger}
                max_pow {mustBeInteger}
            end
            % restrict i to gf ext exponents 
            exponents = mod(exponents, 2^obj.m - 1);
            % init power set 
            power_set =  zeros(size(exponents, 2), max_pow + 1);
            % powers
            pows = 0:max_pow;
            power_set(:, pows + 1) = mod(exponents' .* pows, 2^obj.m-1) ;
            return;  
        end
        function res = gf_ext_poly_eval(obj, poly_coeffs, exponents) 
            arguments 
                obj
                poly_coeffs (1,:) {mustBeNumeric} % indexes of coefficients in gf_ext table
                exponents (1,:) {mustBeNumeric} % exponents of alpha in gf ext, to be evaluated
            end
            % generate power set, one row for each value of i
            power_set = obj.gf_power_set(exponents, size(poly_coeffs,2) -1);
            % rep poly coeffs
            poly_coeffs_rep = repmat(poly_coeffs, length(exponents), 1);
            % eval terms 
            eval_terms = gfmul(power_set, poly_coeffs_rep, obj.gf_ext);
            % init results 
            res = zeros(length(exponents), 1) + -Inf;
            % eval each row 
            for i = 1:length(poly_coeffs)
                res = gfadd(res, eval_terms(:, i), obj.gf_ext);
            end
            return;
        end
        function [message, no_error_detected] = decode(obj, r)
            arguments
                obj
                r (1,:) {mustBeLessThanOrEqual(r, 1), mustBeInteger}
            end
            % syndrome in exponential form of recieved code polynomial
            syndrome = obj.calculate_syndrome(r);
            % set error detection flag 
            no_error_detected = all(syndrome == -Inf);
            if (no_error_detected == 1)
                % no error detected
                message = r(end-obj.k+1:end);
                return;
            end
            % simplified berlekov iterative algorithm 
            sigma_x = obj.simplified_berlekamp_massey(syndrome);
            % find error locations
            exponents = 0:2^obj.m-2;
            eval_exponents = obj.gf_ext_poly_eval(sigma_x, exponents);
            root_exps = exponents(eval_exponents == -Inf);
            % correct errors
            if numel(root_exps) ~= 0
                % invert roots to get error locations
                error_locations = gfdiv(zeros(size(root_exps)), root_exps, obj.gf_ext) + 1;
                % correct errors
                r(error_locations) = mod(r(error_locations) + 1, 2);
            end
            % extract message bits (systematic code)
            message = r(end-obj.k+1:end);
            return;
        end 
    end
end 

