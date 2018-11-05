function plot_imag(y,varargin)
% function plot_imag(y,varargin)
%
% plot function that takes as input a complex array and plots the real and
% the imaginary part separately on two subplots

subplot(1,2,1); 
plot(real(y),varargin{:}); 
title('real part');hold on
subplot(1,2,2); plot(imag(y), varargin{:}); title('imaginary part');hold on
