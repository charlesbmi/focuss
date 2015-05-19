% Checks that x'(A'y) = y'(Ax)

%x = rand(numel(q),1);
%y = rand(sum(sum(sum(mask))),1);
%xTATy = x'*AFUN(y,2)
%yTAx = y'*AFUN(x,1)

L = 1e-2;
prev_q = rand(size(q));
x = rand(size(q));
%y = rand(size(q));
mask = rand(size(q)) < 0.5;
y = rand(sum(sum(sum(mask))) + numel(q),1);
%A = @(q) mask.*xf2kt(abs(prev_q).*q); % Sampling * Fourier * W
%AT = @(v) kt2xf(v.*mask).*abs(prev_q);

AFUN = @(x,topt) afun(x,prev_q,size(q),mask,L,topt);
Ax = AFUN(x(:),'notransp');
ATy = AFUN(y,'transp');

xTATy = x(:)'*ATy(:)
yTAx = y(:)'*Ax(:)
