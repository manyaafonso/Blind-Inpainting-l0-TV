function u=projk1d(g,lambda,niter)

tau=0.25;

ux=length(g);

pn=zeros(ux,1);

for i=1:niter
    qn = Q(Qstar(pn)-g/lambda);
    pn=(pn + tau*qn)./(1 + tau*abs(qn));
end

u = g - lambda*Qstar(pn);

function  y = Q(x)

m = length(x);

%y = zeros(m,1);

y = diff1(x);

function y = Qstar(x)

m = length(x);

y = diff2(x);

function y=diff1(x)
h=[0; 1; -1];
y=conv2c(x,h);

function y=diff2(x)
h=[1; -1; 0];
y=conv2c(x,h);
