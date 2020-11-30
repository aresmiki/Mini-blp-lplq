function [optW,rec,funv,Info]=min_lplq(x,np,option,p,q)
 % Blind deconvolution based on G-lplq norm
    %  code by Liu He (aremiki@163.com), January 2020
    %  Used in my PhD research at the University of SouthWest Jiaotong University.
    %   
    %
    %    Algorithm Reference:
    %    Optimized minimum generalized Lp/Lq deconvolution for recovering repetitive 
    %    impacts from a vibration mixture, Measurement, vol. 168，
    %     doi: 10.1016/j.measurement.2020.108329 
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
    % Note 1:
    %  When using this code, add the minFunc file to the MATLAB
    %  environment.
    %  Execute the code "addpath('./minFunc');"in this path.
    % 
    % Note 2:
    %   The solution is not guaranteed to be the optimal solution 
    %   to our minimization problem, the solution is just a local 
    %   minimum and therefore a good pick.
if nargin <= 3
    p=1;q=2;
end
% x-信号
% np-神经元个数
%% 构造Hankel 矩阵
N=length(x);
data=[];
for i=1:floor(N-np+1)
    data(i,:)=x(i:(i+np-1));
end
%% 初始化权重
optW = zeros(size(data, 2),1); %使用高峰度的单点函数
optW(2)=1;
optW=optW./norm(optW);
%%
if option==0
[optW,funv,~,Info]= minFunc(@L1_L2obj_grad, optW(:), ...
                   struct('MaxIter', 400,'Display','off'), data,p,q);
else
[optW,funv,~,Info]  = minFunc(@L1_L2obj, optW(:), ...
                   struct('MaxIter', 400,'numDiff',1,'Display','off'), data,p,q);
end
optW=optW./norm(optW);        
rec= data*optW;
end
%%  目标函数
function obj=L1_L2obj(w,A)
    y=A*w;
    obj=sign(log(q/p))*(norm(y,p)/norm(y,q)).^p;
end
%% 不包络归一化算法 广义Lp/Lq
function [obj,grad]=L1_L2obj_grad(w,A,p,q)
    y=A*w;
    C= sqrt(y.^2 + 1e-8);
    obj=sign(log(q/p)).*(norm(C,p)./norm(C,q)).^p;
    aJdR=sign(log(q/p)).*p*(norm(C,p)./norm(C,q)).^(p-1);
       %%    链式法则求导数 p=1,q=2
    aRdC=(norm(C,p).^(1-p).*C.^(p-1).*norm(C,q))./norm(C,q).^2-(norm(C,q).^(1-q).*C.^(q-1).*norm(C,p))./norm(C,q).^2;
    aCdy= y./C;
    grad=A'*(aJdR.*aRdC.*aCdy);  
end

