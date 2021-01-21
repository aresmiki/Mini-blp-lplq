function [optW,rec,funv,Info]=min_blp_lplq(x,np,option,p,q)
%%
    % Blind deconvolution based on criterion defined by envelope spectrum
    %  code by Liu He (aremiki@163.com), July 2020
    %  Used in my PhD research at the University of SouthWest Jiaotong University.
    %
    %  This research work has been published on the Journal of Signal Processing.
    %
    % Inputs:
    %    x: 
    %       Signal to perform blind deconvolution on. 
    % 
    %    np:
    %       This is the length of the finite inpulse filter filter to 
    %       design. Investigate the performance difference using
    %       different values.
    %
    %    option: 
    %       option == 0, Analytic gradient.
    %       option != 0, Numerical gradient.
    %
    %    p and q: (OPTIONAL)
    %       The Sparsity criterion G-Lp/Lq, G-L1/L2 is the default value.
    %
    % Outputs:
    %    optW:
    %         The final 1d filter in finite impulse response format.
    %
    %
    %    rec:
    %       The input signal(s) x, filtered by the resulting our method's filter.
    %
    %    funv:
    %       the final G-Lp/Lq 
    %
    %
    %    Info:
    %       Information on the optimization process
    % 
    % 
    % Note 1:
    %    When using this code, add the minFunc file to the MATLAB environment
    %    Execute the code "addpath('./minFunc');"in this path.
    %
    % Note 2:
    %   The solution is not guaranteed to be the optimal solution 
    %   to our minimization problem, the solution is just a local 
    %   minimum and therefore a good pick.
%%
if nargin <= 3
    p=1;q=2;
end
% x-信号
%% 构造Hankel 矩阵
N=length(x);
data=[];
for i=1:floor(N-np+1)
    data(i,:)=x(i:(i+np-1));
end
NP=max(size(data));
%% 初始化权重
 % Assume initial filter as a delayed impulse
optW = zeros(size(data, 2),1); %使用高峰度的单点函数
optW(2)=1;
optW=optW./norm(optW);
%%
if option==0
Psi=fft(eye(NP,NP));    %  傅里叶正变换矩阵
[HM]=get_hilbfir_M(NP);  % Hilbert变换矩阵
% 真实包络算法
[optW,funv,~,Info]  = minFunc(@L1_L2obj_grad,optW(:), ... 
                   struct('MaxIter', 40,'Display','off','TolX',1e-4),data,Psi,HM,p,q);
else
% 数值的包络算法
[optW,funv,~,Info]  = minFunc(@L1_L2obj_grad_num,optW(:), ...
   struct('MaxIter', 400,'numDiff',1,'Display','off','TolX',1e-4),data,p,q);
end
optW=optW./norm(optW);        
rec= data*optW;
end


