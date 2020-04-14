% potentially encode using generator matrix, make it systemaic via imrot,
% rref, then imrot

classdef bch
    properties (Access = private)
        n; % code length
        k; % message length
        m;
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
            obj.gf_ext = gftuple([0:2^obj.m-2]', obj.m); % doesnt include 0
            
        end 
        function code_polynomial = encode(obj, msg)
             msg_shifted = [zeros(1, obj.n -obj.k), msg]; % multiply by x^(n-k), increasing power
             [~, rem] = gfdeconv(msg_shifted, obj.gen_poly); % get remainder
             code_polynomial = gfadd(msg_shifted, rem, 2); 
        end
        function res = gf_poly_eval(obj, poly, i) % i is the power of alpha 
            i = mod(i, 2^obj.m-1); 
            pow_set = zeros(size(poly, 2), obj.m); % [1, x, x^2, .... ]'
            
            % init power set 
            pow_set(1, :) = obj.gf_ext(1, :);
            pow_set(2, :) = obj.gf_ext(i+1, :);
            for pow = 2:size(poly, 2)-1
                n = mod(i * pow, 2^obj.m-1); % alpha^n = (alpha^i)^pow
                pow_set(pow + 1, :) = obj.gf_ext(n+1, :); % alpha^n polynomial form 
            end 
            
            % dot product
            prod = poly' .* pow_set;
            
            % add rows for result
            res = mod(sum(prod, 1), 2); % add elements in gf2
            return;
        end 
        function s = syndrome(obj, poly) % fix this to use matrix mult with parity check matrix in the future, for now its fine 
            s = zeros(2*obj.t, obj.m); % row i is si
            for i = 1:2*obj.t
                s(i, :) = obj.gf_poly_eval(poly, i);
            end
            return;
        end 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%         function res = code_poly_eval(obj, code_poly, i) % i is the power of alpha 
%             pow_set = zeros(size(code_poly,2), size(x,2)); % [1, x, x^2 ..., x^code_poly_max_degree]', converted to polynomial form 
%             % init power set
%             pow_set(1,:) = obj.gf_ext(1,:);
%             pow_set(2,:) = x;
%             disp(obj.gf_ext(1,:));
% %             for i = 1:size(code_poly,2)
% %                 x_i_poly = gf
% %             end 
%         end 
%         function s = syndrome_calc(obj, gf_poly)
%            
%             deg = find(gf_poly, 1, 'last') - 1;
%             s = zeros(2*obj.t, obj.m);
%             
%             for i = 1:2*obj.t
%                 
%             end 
%         end
        
%         
%         
%         function message = decode(obj, r)
%             message = bchdec(flip(r), obj.n, obj.k);
%             message = flip(message);
%         end
%         
%         
%         
%         
%         function rt = get_roots(obj)
%             rt = gfroots(flip(obj.gen_poly), 24, 2);
%         end
    end
end 

