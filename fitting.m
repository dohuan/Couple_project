function out = fitting(y,A,B,C,K,D)

nt = size(y,1);
if (nargin<2)                         
    A = [NaN NaN;NaN NaN];
    B = [NaN;NaN];
    C = eye(2);
    K = zeros(2,2);
    D = zeros(2,1);
end
x0 = y(1,[1 2])';
ms = modstruc(A,B,C,D,K,x0);
th = ms2th(ms,'d',[],[],0.5);

if size(y,2)~=3
    u = ones(nt,1);
    th = pem([y u],th);
else
    th = pem(y,th);
end
% --- Regenerate the signal
u = ones(nt,1);
%e = randn(nt,2);
th = sett(th,0.5);
%out.y_fit = idsim([u e],th,x0);
%out.y_fit = idsim(u,th,x0);
out.eig = eig(th.A);
[out.y_fit,out.fitpercent,~ ] = compare(y,th);

end