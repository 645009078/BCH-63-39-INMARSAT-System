gen_poly = [1,1,1,0,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1]; % acending power
code = bch(63,39,9,4,gen_poly);
  
code.syndrome(gen_poly);



res = code.gf_poly_eval([1,1,0,0,0,0,1], 3);


syndrome = code.syndrome(gen_poly); % should be all 0

% example 

message = zeros(1,39);
message(1:20) = 1;
message


% encode 
codepoly = code.encode(message)

% check systematic 
codepoly(25:end) == message

% show parity
disp("parity bits");
codepoly(1:24)


% syndrome 
msg_synd = code.syndrome(codepoly);
disp("message syndrome: ");
msg_synd


% adding an error and showing syndrome 
disp("code poly with 1 error in position 22, syndrome");
code_poly(22) = 1;
inv_syndrome = code.syndrome(code_poly);
inv_syndrome












% checking, example 6.1 and 6.4 from textbook
gen_61 = [1,0,0,0,1,0,1,1,1];
code_61 = bch(15, 7, 5, 2, gen_61);
syndrome_example_64 = code_61.syndrome([1,0,0,0,0,0,0,0,1,0,0,0,0,0,0]);
syndrome_example_64 
% result is [alpha2, alpha4, alpha7, alpha8], shown in example 6.4







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