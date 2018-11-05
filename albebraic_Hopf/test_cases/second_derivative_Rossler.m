function DDF = second_derivative_Rossler(~,~,~)
%
% computation of DDF(y+zR), z\inB_1(0)
%
% since for Rossler the system in of order 2, the second derivative is
% constant, the inputs are not taken into account

DDF=[-2;2;2;0;-1;6;6;2;-6;2;2;-2;0;0];
