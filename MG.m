function error2 = MG

clear all;

Ncycles = 10000;

nu1 = 3;
nu2 = 3;
omega = 0.5;

Nmin = 20;

M = 1280;
dx = 1/M;

x = (0:dx:1)';
fh = (3*x + x.^2).*exp(x);

% Compute solution by solving A u = f directly
A = (1/dx^2)*(2*diag(ones(M-1,1)) - diag(ones(M-2,1),1) - diag(ones(M-2,1),-1));
u = zeros(M+1,1);
u(2:end-1) = A\fh(2:end-1);

uh = u*0;
for icycl = 1:Ncycles
    uh = MG_rec( uh, fh, nu1, nu2, omega, Nmin);
    plot(x, u, 'r-',  x, uh, 'k-', 'LineWidth',2); grid on; drawnow
    e = u - uh;
    rh = residual(uh, fh);
    error2(icycl) = norm(rh)*sqrt(dx);
    if error2(icycl) < 1.e-3, break, end
end

end

function uh = MG_rec( uh, fh, nu1, nu2, omega, Nmin)
    uh = relax( uh, fh, nu1, omega);

    if length(uh) > Nmin
        rh = residual(uh, fh);
        r2h = restrict(rh);
        e2h = MG_rec( 0*r2h, r2h, nu1, nu2, omega, Nmin);
        eh = prolongate(e2h);
        uh = uh + eh;
    end
    
    uh = relax( uh, fh, nu2, omega);
end

function uh = relax( uh, f, nit, omega)

if nit> 0, fprintf("Smoothing: grid %d\n", length(uh)); end

M = length(uh)-1;
dx = 1/M;

uhn = uh;
for it = 1:nit
    for i = 2:M
        uhn(i) = omega*0.5*(uh(i+1) + uh(i-1) + dx^2*f(i)) + (1-omega)*uh(i);
    end
    uh = uhn;
end
end

function rh = residual( uh, fh)

M = length(uh)-1;
dx = 1/M;

rh = uh*0;
for i = 2:M
    rh(i) = fh(i) + (uh(i+1) - 2*uh(i) + uh(i-1))/dx^2;
end
end

function f2h = restrict( fh)

M = length(fh)-1;
f2h = zeros(M/2+1,1);

for i = 1:M/2+1
    f2h(i) = fh(2*i-1);
end
end

function eh = prolongate( e2h)

M = length(e2h)-1;
eh = zeros(2*M+1,1);

for i = 2:M
    eh(2*i-1) = e2h(i);
end

for i = 2:M+1
    eh(2*i-2) = 0.5*(e2h(i) + e2h(i-1));
end

end