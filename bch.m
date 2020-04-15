% potentially encode using generator matrix, make it systemaic via imrot,
% rref, then imrot

classdef bch
    properties (Access = private)
        n; % code length
        k; % message length
        m; % prim poly degree
        d; % minimum distance
        t; % error correcting capability
        gen_poly;  % decreasing power
        gen_poly_deg; % degree of gen poly
        gf_ext; % 2^m x (m-1) matrix, rows are elements of the extended gf in polynomial form
        
    end
    methods
        function obj = bch(n, k,d,t,gen_poly)
            arguments
                n {mustBeNumeric}
                k {mustBeNumeric} 
                d {mustBeNumeric} 
                t {mustBeNumeric} 
                gen_poly(1,:) {mustBeNumeric} % decreasing power
            end
            obj.n = n;
            obj.k = k;
            obj.m = log2(n + 1);
            obj.d = d;
            obj.t = t;
            obj.gen_poly = gen_poly;
            obj.gen_poly_deg = find(gen_poly, 1, 'last') - 1;
            Power = [-1:2^obj.m-2]';
            Polynomial = gftuple(Power, obj.m); 
            Decimal = bi2de(Polynomial);
            obj.gf_ext = table(Power, Decimal, Polynomial); % doesnt include 0
            
        end 
        function code_polynomial = encode(obj, msg)
             msg_shifted = [zeros(1, obj.n -obj.k), msg]; % multiply by x^(n-k), increasing power
             [~, rem] = gfdeconv(msg_shifted, obj.gen_poly); % get remainder
             code_polynomial = gfadd(msg_shifted, rem, 2); 
        end
        function ps_exponents = gf_power_set_indexes(obj, i, max_pow)
            i = mod(i, 2^obj.m-1);
            ps_exponents = zeros(size(i, 2), max_pow+1)
            ps_exponents(:, 1) = 2; 

            if max_pow == 0 
                return;
            end
            % calculate power set exponents 
            ps_exponents(:, 2) = i' + 2;
            pows = 2:max_pow;
            ps_exponents(:, pows+1) = mod(i' .* pows, 2^obj.m-1) + 2;
        end 
        function res = gf_poly_eval(obj, poly_coeffs, i) % i is for alpha^i
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
        function s_indexs = syndrome_indexes(obj, poly) % fix this to use matrix mult with parity check matrix in the future, for now its fine 
            % elements [alpha, ...., alpha^2t] of gf(2^m)
            s_indexs = ones(1, 2*obj.t); 
            i =  1:2*obj.t;
            % evaluate polynomial for all elements 
            si_poly = obj.gf_poly_eval(poly,i);
            % convert result for each input alpha^i to decimal format.
            si_dec = bi2de(si_poly);
            for n = i
                % find row of result in gf_ext table
                [s_indexs(n), ~] = find(obj.gf_ext.Decimal == si_dec(n), 1);
            end 
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
        function d_index = calculate_simplified_discrepency(obj, bia_tbl, syndrome_indexes, mu);
            syndrome_indexes = [syndrome_indexes, ones(1, 2*mu+3 - length(syndrome_indexes))];
            sigma_mu = bia_tbl{mu+3, 2};
            l = bia_tbl{mu+3, 4};
            d_index = syndrome_indexes(2*mu+3);
            
            for i = 1:l
                % reverse order 
                temp = obj.multiply_gf_ext_elements(syndrome_indexes(2*mu+3-i), sigma_mu(i+1));
                d_index = obj.add_gf_ext_elements(temp, d_index);
            end
            return;
        end 
        function rho = calculate_rho(obj, bia_tbl, max_rho)
           
            
            rho = -1/2;
            val = -1; % 2rho - l
            for mu = max_rho:-1:0
               
                if bia_tbl{mu+2, 3} ~= 1 && bia_tbl{mu+2, 5} - mu  > val
                    rho = mu;
                    val = bia_tbl{mu+2, 5};
                end
            end
            rho
            
            
            
            
            return;
        end
        function inv_index = invert_gf_ext_element(obj, elem_index)
            elem_exponent = elem_index - 2;            
            inv_exponent = mod(2^obj.m - 1 - elem_exponent, 2^obj.m - 1);
            inv_index = inv_exponent + 2;
        end
        function res_index = multiply_gf_ext_elements(obj, a_index, b_index)
            if a_index == 1 || b_index == 1
                res_index = 1;
                return;
            end
            
            if a_index == 2 
                res_index = b_index;
                return;
            end
            
            if b_index == 2
                res_index = a_index;
                return;
            end
            
            a_exponent = a_index - 2;
            b_exponent = b_index - 2;
            res_exponent = mod(a_exponent + b_exponent, 2^obj.m - 1);
            res_index = res_exponent + 2;

        end
        function res_index = add_gf_ext_elements(obj, a_index, b_index)
            a_poly = obj.gf_ext.Polynomial(a_index, :);
            b_poly = obj.gf_ext.Polynomial(b_index, :);
            res_poly = gfadd(a_poly, b_poly, 2);
            res_index = find(obj.gf_ext.Decimal == bi2de(res_poly), 1);
        end
        function res = add_gf_ext_polynomials(obj, poly_a, poly_b)
            
            % pad with ones 
            if length(poly_b) > length(poly_a)
                poly_a = [poly_a, ones(1, length(poly_b) - length(poly_a))];
            end 
            if length(poly_a) > length(poly_b)
                poly_b = [poly_b, ones(1, length(poly_a) - length(poly_b))];
            end 
            
            
            res = poly_a;
            for i = 1:length(res)
                 res(i) = obj.add_gf_ext_elements(res(i), poly_b(i));
            end
            return;
        end
        function sigma_x = simplified_berlekamp_massey(obj, r) 
            % initialise 2t + 2 x 5 table for iterative algorithm'
            % rows represent interations over [-1, 2t]
            bia_tbl = obj.init_simplified_bia_table();
            % syndrome element indexes of recieved code polynomial
            s_inds = obj.syndrome_indexes(r);
            bia_tbl{2,3} = s_inds(1);
            
            for mu = 0:obj.t-1
                d_mu_index = bia_tbl{mu+2, 3};
                bia_tbl{mu+3, 2} = bia_tbl{mu+2, 2};
                % non-zero discrepency
                if d_mu_index ~= 1
                    rho = obj.calculate_rho(bia_tbl, mu-1);  
                    d_mu_index = bia_tbl{mu+2, 3}; % element in extended gf
                    d_rho_index = bia_tbl{floor(rho)+2,3}; % element in extended gf
                    sigma_rho = bia_tbl{floor(rho)+2, 2}; % polynomial over extended gf 
                    x_deg = 2*(mu-rho);
                    
                    % multiply by x^x_deg
                    % using ones because coefficients are indexes of gf_ext
                    % field, where 0 is in index 1. 
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
                bia_tbl{mu+3, 3} = obj.calculate_simplified_discrepency(bia_tbl, s_inds, mu);
                bia_tbl{mu+3, 5} = 2*bia_tbl{mu+3, 1} - bia_tbl{mu+3, 4};
            end
            sigma_x = bia_tbl{obj.t+2, 2};
            return;
        end
        function res = gf_ext_poly_eval(obj, poly_coeff_indexes, i) % i are exponents, not indexes!
            i = mod(i, 2^obj.m-1);
            power_set_indexes = obj.gf_power_set_indexes(i, size(poly_coeff_indexes,2)-1);
            
            res = 1; % zero default
            prod_index = ones(1, length(power_set_indexes));
            for k = 1:length(power_set_indexes);
                prod_index = obj.multiply_gf_ext_elements(poly_coeff_indexes(k), power_set_indexes(k));
                res = obj.add_gf_ext_elements(res, prod_index);
            end   
        end
        function root_indexes = find_roots_sigma_x(obj, sigma_x)
            values_poly = [];
            exponents = 0:obj.n-1;
            for i = 0:obj.n-1
                obj.gf_ext_poly_eval(sigma_x, i)
                values_poly(end + 1) = obj.gf_ext_poly_eval(sigma_x, i);
            end
            
            % exponent i in index i + 1, add 1 for row index in gf_ext 
            root_indexes = find(values_poly == 1) + 1 ;
        end
        function message = decode(obj, r)
            % simplified berlekov iterative algorithm 
            sigma_x = obj.simplified_berlekamp_massey(r)
            % find error locations
            root_indexes = obj.find_roots_sigma_x(sigma_x);
            inverse_root_indexes = obj.invert_gf_ext_element(root_indexes);
            error_locations = obj.gf_ext.Power(inverse_root_indexes)
            
            r_corrected = r;
            length(error_locations);
            for i = error_locations
                if i ~= 0
                    r_corrected(i+1) = mod(r_corrected(i + 1) + 1, 2);
                end
            end
            message = r_corrected(end-obj.k+1:end);
            
   
        end
    end
end 





% % potentially encode using generator matrix, make it systemaic via imrot,
% % rref, then imrot

% classdef bch
%     properties (Access = private)
%         n; % code length
%         k; % message length
%         m; % prim poly degree
%         d; % minimum distance
%         t; % error correcting capability
%         gen_poly;  % decreasing power
%         gen_poly_deg; % degree of gen poly
%         gf_ext; % 2^m x (m-1) matrix, rows are elements of the extended gf in polynomial form
        
%     end
%     methods
%         function obj = bch(n, k,d,t,gen_poly)
%             arguments
%                 n {mustBeNumeric}
%                 k {mustBeNumeric} 
%                 d {mustBeNumeric} 
%                 t {mustBeNumeric} 
%                 gen_poly(1,:) {mustBeNumeric} % decreasing power
%             end
%             obj.n = n;
%             obj.k = k;
%             obj.m = log2(n + 1);
%             obj.d = d;
%             obj.t = t;
%             obj.gen_poly = gen_poly;
%             obj.gen_poly_deg = find(gen_poly, 1, 'last') - 1;
%             Power = [-1:2^obj.m-2]';
%             Polynomial = gftuple(Power, obj.m); 
%             Decimal = bi2de(Polynomial);
%             obj.gf_ext = table(Power, Decimal, Polynomial); % doesnt include 0
            
%         end 
%         function code_polynomial = encode(obj, msg)
%              msg_shifted = [zeros(1, obj.n -obj.k), msg]; % multiply by x^(n-k), increasing power
%              [~, rem] = gfdeconv(msg_shifted, obj.gen_poly); % get remainder
%              code_polynomial = gfadd(msg_shifted, rem, 2); 
%         end
%         function ps_exponents = gf_power_set_indexes(obj, i, max_pow)
%             i = mod(i, 2^obj.m-1);
%             ps_exponents = zeros(size(i, 2), max_pow+1);
%             % calculate power set exponents 
%             ps_exponents(:, 1) = 2; 
%             ps_exponents(:, 2) = i' + 2;
%             pows = 2:max_pow;
%             ps_exponents(:, pows+1) = mod(i' .* pows, 2^obj.m-1) + 2;
%         end 
%         function res = gf_poly_eval(obj, poly_coeffs, i) % i is for alpha^i
%             i = mod(i, 2^obj.m-1);
            
%             % calculate row indexes for each term of each power set 
%             power_set_indexes = obj.gf_power_set_indexes(i, size(poly_coeffs,2)-1); 
%             % convert to row major order for grouping (reshaping) 
%             power_set_indexes = power_set_indexes';
%             power_set_indexes = power_set_indexes(:)';
            
%             % calculate column vector [(alpha^i)^1;...; (alpha^i)^p], for
%             % all i. p is the degree of the input polynomial.
%             power_set_polynomials = obj.gf_ext.Polynomial(power_set_indexes, :);
%             % repeat for each power set
%             repeated_coeffs = repmat(poly_coeffs, 1, size(i, 2));
%             % element wise multiplication with respective coefficients 
%             power_set_polynomials = power_set_polynomials .* repeated_coeffs';
            
%             % reshape power set, grouping polynomial form of power set of 
%             % alpha^i in power_set_polynomials(:,:,i)
%             power_set_polynomials = reshape(power_set_polynomials', [obj.m, size(poly_coeffs, 2), size(i, 2)]);
            
%             % eval the polynomial, adding the terms of each set multiplied
%             % by their coefficients. summing c0*1 + c1*alpha^i + ... +
%             % cp*(alpha^i)^p for each i
%             res = mod(sum(power_set_polynomials, 2), 2);
%             res = reshape(res, [obj.m,size(i,2),1])';
%             res = res(:,:,1)

%         end 
%         function s_indexs = syndrome_indexes(obj, poly) % fix this to use matrix mult with parity check matrix in the future, for now its fine 
%             % elements [alpha, ...., alpha^2t] of gf(2^m)
%             s_indexs = zeros(1, 2*obj.t); 
%             i =  1:2*obj.t;
%             % evaluate polynomial for all elements 
%             si_poly = obj.gf_poly_eval(poly,i);
%             % convert result for each input alpha^i to decimal format.
%             si_dec = bi2de(si_poly);
%             for n = i
%                 % find row of result in gf_ext table
%                 [s_indexs(n), ~] = find(obj.gf_ext.Decimal == si_dec(n), 1);
%             end 
%             return;
%         end 
%         function message = decode(obj, r)
%             s_inds = obj.syndrome_indexes(r);
                       
            
%             newtons_coeffs = zeros(2*obj.t + 2, obj.n);
%             d = zeros(2*obj.t + 2, 1);
%             l = zeros(2*obj.t + 2, 1);
            
%             bip_tbl = table(newtons_coeffs, d, l);
%             % init values 
%             bip_tbl.newtons_coeffs(1) = 1;
%             bip_tbl.d(1) = 1;
%             bip_tbl.l(1) = 0;
%             bip_tbl.newtons_coeffs(2) = 1;
%             bip_tbl.d(2) = s_inds(1);
%             bip_tbl.l(2) = 0;
            
%             % iterate until mu = 2t
%             mu = 0;
            
            
            
% %             du = zeros(1, obj.n);
% %             I = zeros(1, obj.n);
% %             rho = zeros(1, obj.n);
% %             
% %             %
% %             % initialise variables 
% %             newtons_coeffs(1) = mod(s_inds(1), 2^obj.m-2); % row index for gf_ext table 
% %             
% %          
            
    
%         end
%     end
% end 


% % shit to fix 
% % - s2i = si^2 (makes syndrome comp much more efficient) 
% % - make it less shit, right now its a mess (potentially generate a
% % syndrome computation matrix 




% %             res = zeros(size(i,2), obj.m);
% %             for n = i
% %                 power_set_indexes = obj.gf_power_set_indexes(n, size(poly_coeffs,2)-1);
% %                 power_set_polynomials = obj.gf_ext.Polynomial(power_set_indexes, :)
% %                 res(n, :) = mod(sum(power_set_polynomials .* poly_coeffs', 1), 2)
% % 
% %             end
                        
           































% %         function res = gf_poly_eval(obj, poly, i) % i is the power of alpha 
% %             i = mod(i, 2^obj.m-1); 
% %             pow_set = zeros(size(poly, 2), obj.m); % [1, x, x^2, .... ]'
% %             
% %             % init power set 
% %             pow_set(1, :) = obj.gf_ext.Polynomial(2, :); % 1
% %             pow_set(2, :) = obj.gf_ext.Polynomial(i+2, :); % x
% %             for pow = 2:size(poly, 2)-1
% %                 n = mod(i * pow, 2^obj.m-1); % alpha^n = (alpha^i)^pow
% %                 pow_set(pow + 1, :) = obj.gf_ext.Polynomial(n+2, :); % alpha^n polynomial form 
% %             end 
% %             
% %             % dot product
% %             prod = poly' .* pow_set;
% %             
% %             % add rows for result
% %             res = mod(sum(prod, 1), 2); % add elements in gf2
% %             return;
% %         end 
%         % array of row indexes in gf_ext






