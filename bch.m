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
        function ps_inds = gf_power_set_indexes(obj, i, max_pow) % i from alpha^i
            i = mod(i, 2^obj.m-1);
            ps_inds = zeros(max_pow+1, 1);
            ps_inds(1) = 2; % 1
            ps_inds(2) = i + 2;
            pows = 2:max_pow;
            ps_inds(pows+1) = mod(i * pows, 2^obj.m-1) + 2;
        end 
        function res = gf_poly_eval(obj, poly_coeffs, i) % i is for alpha
            i = mod(i, 2^obj.m-1);
            power_set_indexes = obj.gf_power_set_indexes(i, size(poly_coeffs,2)-1);
            power_set_polynomials = obj.gf_ext.Polynomial(power_set_indexes, :);
            res = mod(sum(poly_coeffs' .* power_set_polynomials, 1), 2);
        end 
%         function res = gf_poly_eval(obj, poly, i) % i is the power of alpha 
%             i = mod(i, 2^obj.m-1); 
%             pow_set = zeros(size(poly, 2), obj.m); % [1, x, x^2, .... ]'
%             
%             % init power set 
%             pow_set(1, :) = obj.gf_ext.Polynomial(2, :); % 1
%             pow_set(2, :) = obj.gf_ext.Polynomial(i+2, :); % x
%             for pow = 2:size(poly, 2)-1
%                 n = mod(i * pow, 2^obj.m-1); % alpha^n = (alpha^i)^pow
%                 pow_set(pow + 1, :) = obj.gf_ext.Polynomial(n+2, :); % alpha^n polynomial form 
%             end 
%             
%             % dot product
%             prod = poly' .* pow_set;
%             
%             % add rows for result
%             res = mod(sum(prod, 1), 2); % add elements in gf2
%             return;
%         end 
        % array of row indexes in gf_ext
        function s_indexs = syndrome_indexes(obj, poly) % fix this to use matrix mult with parity check matrix in the future, for now its fine 
            s_indexs = zeros(1, 2*obj.t); % col i is row of si in gf_ext
            for i = 1:2*obj.t
                poly
                size(poly)
                si_poly = obj.gf_poly_eval(poly, i)
                [s_indexs(i), ~] = find(obj.gf_ext.Decimal == bi2de(si_poly), 1);
            end
            return;
        end 
    end
end 


% shit to fix 
% - s2i = si^2 (makes syndrome comp much more efficient) 
% - make it less shit, right now its a mess (potentially generate a
% syndrome computation matrix )


