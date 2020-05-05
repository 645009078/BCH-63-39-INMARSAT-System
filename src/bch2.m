
classdef bch2
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
        function obj = bch2(n, k,d,t,gen_poly, prim_poly)
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
            obj.m = size(prim_poly, 2) - 1;
            obj.d = d;
            obj.t = t;
            obj.gen_poly = gen_poly;
            obj.gen_poly_deg = find(gen_poly, 1, 'last') - 1;
            Power = [-1:2^obj.m-2]';
            Polynomial = gftuple(Power, prim_poly); 
            Decimal = bi2de(Polynomial);
            obj.gf_ext = table(Power, Decimal, Polynomial); 
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
        function ps_indexes = gf_power_set_indexes(obj, i, max_pow)
            arguments
               obj
               i {mustBeInteger}
               max_pow {mustBeInteger}
            end
            i = mod(i, 2^obj.m-1);
            ps_indexes = zeros(size(i, 2), max_pow+1);
            ps_indexes(:, 1) = 2; 

            if max_pow == 0 
                return;
            end
            
            % calculate power set exponents 
            ps_indexes(:, 2) = i' + 2;
            pows = 2:max_pow;
            ps_indexes(:, pows+1) = mod(i' .* pows, 2^obj.m-1) + 2;
            return;
        end 
        function res = gf_poly_eval(obj, poly_coeffs, i) % i is for alpha^i
            arguments
                obj
                poly_coeffs(1, :) {mustBeLessThanOrEqual(poly_coeffs, 1)}
                i {mustBeInteger}
            end
            i = mod(i, 2^obj.m-1);
            
            
            % calculate row indexes for each term of each power set 
            power_set_indexes = obj.gf_power_set_indexes(i, size(poly_coeffs,2)-1);
            % convert to row major order for grouping (reshaping) 
            power_set_indexes = power_set_indexes';
            power_set_indexes = power_set_indexes(:)';
            
            % calculate column vector [(alpha^i)^1;...; (alpha^i)^p], for
            % all i. p is the degree of the input polynomial.
            power_set_polynomials = obj.gf_ext.Polynomial(power_set_indexes, :);
            % repeat for each power set
            repeated_coeffs = repmat(poly_coeffs, 1, size(i, 2));
            % element wise multiplication with respective coefficients 
            power_set_polynomials = power_set_polynomials .* repeated_coeffs';
            
            % reshape power set, grouping polynomial form of power set of 
            % alpha^i in power_set_polynomials(:,:,i)
            power_set_polynomials = reshape(power_set_polynomials', [obj.m, size(poly_coeffs, 2), size(i, 2)]);
            
            % eval the polynomial, adding the terms of each set multiplied
            % by their coefficients. summing c0*1 + c1*alpha^i + ... +
            % cp*(alpha^i)^p for each i
            res = mod(sum(power_set_polynomials, 2), 2);
            res = reshape(res, [obj.m,size(i,2),1])';
            res = res(:,:,1);

        end 
        function s_indexs = syndrome_indexes(obj, poly) 
            arguments
                obj
                poly {mustBeLessThanOrEqual(poly, 1)}
            end
            % elements [alpha, ...., alpha^2t] of gf(2^m)
            i =  1:2*obj.t;
            % evaluate polynomial for all elements 
            si_poly = obj.gf_poly_eval(poly,i);
            % convert result for each input alpha^i to decimal format.
            [~, s_indexs] = ismember(bi2de(si_poly), obj.gf_ext.Decimal);
            s_indexs = s_indexs';
            return;
        end
        function bia_tbl = init_simplified_bia_table(obj)
            % t + 1 rows, for range [-1, t]
            % 5 columns for mu, sigma fn, d, l, and mu - l
            bia_tbl = cell(obj.t + 2, 5);
            
            bia_tbl{1,1} = -1/2; % mu
            bia_tbl{1,2} = [2]; % sigma(x) coefficients, increasing order, indexes of gf_ext
            bia_tbl{1,3} = 2; % d
            bia_tbl{1,4} = 0; % l
            bia_tbl{1,5} = -1; % 2mu - l
            bia_tbl{2,1} = 0; % mu
            bia_tbl{2,2} = [2]; % sigma(x)
            bia_tbl{2,4} = 0; % l
            bia_tbl{2,5} = 0; %2mu - l
            return;
        end
        function d_index = calculate_simplified_discrepency(obj, bia_tbl, syndrome_indexes, mu)
            arguments
                obj
                bia_tbl
                syndrome_indexes(1, :) {mustBeInteger}
                mu {mustBeInteger}
            end
            syndrome_indexes = [syndrome_indexes, ones(1, 2*mu+3 - length(syndrome_indexes))];
            sigma_mu = bia_tbl{mu+3, 2};
            l = bia_tbl{mu+3, 4};
            d_index = syndrome_indexes(2*mu+3);
            for i = 1:l
                temp = obj.multiply_gf_ext_elements(syndrome_indexes(2*mu+3-i), sigma_mu(i+1));
                d_index = obj.add_gf_ext_elements(temp, d_index);
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
                if bia_tbl{mu+2, 3} ~= 1 && bia_tbl{mu+2, 5} > val
                    rho = mu;
                    val = bia_tbl{mu+2, 5};
                end
            end
            return;
        end
        function inv_index = invert_gf_ext_element(obj, elem_index)
            arguments
                obj
                elem_index {mustBeInteger}
            end 
            elem_exponent = elem_index - 2;            
            inv_exponent = mod(2^obj.m - 1 - elem_exponent, 2^obj.m - 1);
            inv_index = inv_exponent + 2;
        end
        function res_indexs = multiply_gf_ext_elements(obj, a_indexs, b_indexs)
            arguments
                obj
                a_indexs {mustBeInteger}
                b_indexs {mustBeInteger}
            end 
            a_exponents = a_indexs - 2;
            b_exponents = b_indexs - 2;
            res_exponents = gfmul(a_exponents', b_exponents', obj.gf_ext.Polynomial);
            % account for 0 in extended field 
            res_exponents(res_exponents == -Inf) = -1;
            res_indexs = res_exponents + 2;
        end
        function res_indexs = add_gf_ext_elements(obj, a_indexs, b_indexs)
             arguments
                obj
                a_indexs {mustBeInteger}
                b_indexs {mustBeInteger}
            end 
            a_polys = obj.gf_ext.Polynomial(a_indexs, :);
            b_polys = obj.gf_ext.Polynomial(b_indexs, :);
            res_poly = gfadd(a_polys, b_polys, 2);
            [~, res_indexs] = ismember(bi2de(res_poly), obj.gf_ext.Decimal);
            return;
        end
        function res = add_gf_ext_polynomials(obj, poly_a, poly_b)
            arguments
                obj
                poly_a (1,:) {mustBeInteger} % indexes in gf_ext table
                poly_b (1,:) {mustBeInteger} % indexes in gf_ext table
            end 
            p_a = ones(1, max(numel(poly_a), numel(poly_b)));
            p_a(1:numel(poly_a)) = poly_a;
            p_b = ones(size(p_a));
            p_b(1:numel(poly_b)) = poly_b;
            res = obj.add_gf_ext_elements(p_a, p_b)';
            return;
        end
        function sigma_x = simplified_berlekamp_massey(obj, s_inds) 
            arguments
                obj
                s_inds (1,:) {mustBeInteger} % syndrome indexes 
            end
            % initialise 2t + 2 x 5 table for iterative algorithm'
            % rows represent interations over [-1, 2t]
            bia_tbl = obj.init_simplified_bia_table();
            bia_tbl{2,3} = s_inds(1);
            
            for mu = 0:obj.t-1
                d_mu_index = bia_tbl{mu+2, 3};
                bia_tbl{mu+3, 2} = bia_tbl{mu+2, 2};                
                % non-zero discrepency
                if d_mu_index ~= 1
                    rho = obj.calculate_rho(bia_tbl, mu-1);
                    d_mu_index = bia_tbl{mu+2, 3}; 
                    % element in extended gf
                    d_rho_index = bia_tbl{floor(rho)+2,3}; 
                    % polynomial over extended gf
                    sigma_rho = bia_tbl{floor(rho)+2, 2}; 
                    x_deg = 2*(mu-rho);
                    % multiply by x^x_deg
                    % using ones because coefficients are indexes of gf_ext
                    tmp = [ones(1, x_deg) ,sigma_rho];
                    d_rho_inv_index = obj.invert_gf_ext_element(d_rho_index);
                    coeff_index = obj.multiply_gf_ext_elements(d_mu_index, d_rho_inv_index);
                    for i = 1:length(tmp)
                        tmp(i) = obj.multiply_gf_ext_elements(tmp(i), coeff_index);
                    end
                    bia_tbl{mu+3, 2} = obj.add_gf_ext_polynomials(bia_tbl{mu+3, 2}, tmp);
                end
                % fill table for mu + 1
                bia_tbl{mu+3, 1} = mu + 1;
                bia_tbl{mu+3, 4} = length(bia_tbl{mu+3, 2}) -1;
                bia_tbl{mu+3, 3} = obj.calculate_simplified_discrepency(bia_tbl, s_inds, mu)
                bia_tbl{mu+3, 5} = 2*bia_tbl{mu+3, 1} - bia_tbl{mu+3, 4};
            end
            sigma_x = bia_tbl{obj.t+2, 2};
            return;
        end
        function res = gf_ext_poly_eval(obj, poly_coeff_indexes, i) 
            arguments 
                obj
                poly_coeff_indexes (1,:) {mustBeInteger} % indexes of coefficients in gf_ext table
                i (1,:) {mustBeInteger} % exponents of alpha in gf ext, to be evaluated
            end
            i = mod(i, 2^obj.m-1);
            % generate power set, one row for each value of i
            power_set_indexes = obj.gf_power_set_indexes(i, size(poly_coeff_indexes,2)-1);
            % repeat polynomial coeffs, one row for each row in power set
            % indexes rep
            poly_coeff_indexes_rep = repmat(poly_coeff_indexes, length(i), 1);
            % multiply element wise
            eval_term_indexes = arrayfun(@obj.multiply_gf_ext_elements, power_set_indexes, poly_coeff_indexes_rep);
            % sum terms for each polynomial 
            res = ones(length(i), 1); % zero default, one result per i 
            for i = 1:length(poly_coeff_indexes)
                res = obj.add_gf_ext_elements(res, eval_term_indexes(:, i));
            end   
            return;
        end
        function root_indexes = find_roots_sigma_x(obj, sigma_x)
            arguments
                obj
                sigma_x (1,:) {mustBeInteger} % cofficients, indexes in gf_ext
            end
            exponents = 0:obj.n-1;
            % evaluate sigma_x for each value of i
            values = obj.gf_ext_poly_eval(sigma_x, exponents);
            % indexes of 1 correspond to 0 in gf ext 
            root_indexes = find(values == 1) + 1 ;
        end
        function [message, no_error_detected] = decode(obj, r)
            arguments
                obj
                r (1,:) {mustBeLessThanOrEqual(r, 1)}
            end
            % syndrome element indexes of recieved code polynomial
            s_inds = obj.syndrome_indexes(r)
            % set error detection flag
            no_error_detected = all(s_inds == 1);
            if (no_error_detected == 1)
                message = r(end-obj.k+1:end);
                return;
            end
            % simplified berlekov iterative algorithm 
            sigma_x = obj.simplified_berlekamp_massey(s_inds)
            % find error locations
            root_indexes = obj.find_roots_sigma_x(sigma_x);
            % invert roots
            inverse_root_indexes = obj.invert_gf_ext_element(root_indexes);
            % find exponents (error locations)
            error_locations = obj.gf_ext.Power(inverse_root_indexes) + 1;
            % correct errors
            r(error_locations) = mod(r(error_locations) + 1, 2);
            % extract message bits (systematic code)
            message = r(end-obj.k+1:end);
            return;
        end 
    end
end 

% 
% 
