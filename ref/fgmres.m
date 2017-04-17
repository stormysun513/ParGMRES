function [x, X] = fgmres(A, b, m, tol, maxit, P, x0)
% Flexible GMRES (FGMRES) with restart length m.
% Solves Ax = b to tolerance tol using preconditioner P and initial guess x0
% If P is a function handle, P(A) =approx= I.
% If P is a matrix, P\A =approx= I.
% Optional output X is the whole sequence of FGMRES iterates

% Written by Nick Alger 6/10/2014, released to public domain.
% Reference:
% Saad, "A flexible inner-outer preconditioned GMRES algorithm", SIAM 1993

% Parse input
N = size(b, 1);

if (isa(A, 'numeric'))
    if (size(A, 1) == N && size(A, 2) == N)
        Afct = @(u) A*u;
    else
        error(['Dimensions of A and b in Ax = b are not consistent: '...
               'A is ', num2str(size(A,1)), '-by-', num2str(size(A,2)), ', and '...
               'b is ', num2str(size(b,1)), '-by-', num2str(size(b,2)), '.']);
    end
elseif (isa(A, 'function_handle'))
    Afct = A;
else
    error('Forward operator needs to be a matrix or function handle.')
end

if (nargin < 3 || isempty(m))
    m = 50; %Default restart
end

if (nargin < 4 || isempty(tol))
    tol = 1e-6;    
end

if (nargin < 5 || isempty(maxit))
    maxit = min(N, 200);
end

if (nargin < 6 || isempty(P))
    Pfct = @(u) u;
elseif (isa(P, 'numeric'))
    Pfct = @(u) P\u;
elseif (isa(P, 'function_handle'))
    Pfct = P;
else
    error('Preconditioner needs to be a matrix or function handle.')
end

if (nargin < 7 || isempty(x0))
    x0 = zeros(size(b));
end

%Run FGMRES
X = [];
nit = 1;
while (nit <= maxit)
    disp('Starting new FGMRES cycle')
    H = zeros(m+1, m);
    r0 = b - Afct(x0);
    beta = norm(r0);
    V = zeros(N, m+1);
    V(:, 1) = (1/beta)*r0;
    Z = zeros(N, m);
    for j = 1:m
        Z(:, j) = Pfct(V(:, j));
        w = Afct(Z(:, j))
        for i = 1:j
            H(i, j) = w'*V(:, i);
            w = w - H(i, j)*V(:, i);
        end
        H(j+1, j) = norm(w);
        V(:, j+1) = (1/H(j+1, j))*w;

        e1 = zeros(j+1,1); e1(1) = 1;
        y = H(1 : j+1, 1:j)\(beta*e1);
        x = x0 + Z(:, 1:j)*y;
        if (nargout > 1)
            X = [X,x];
        end
        resnorm = norm(b-Afct(x));
        disp(['k= ', num2str(nit), ...
            ', relative residual= ', num2str(resnorm/norm(b))])
        if (resnorm < tol*norm(b))
            disp(['FGMRES converged to relative tolerance ', ...
                num2str(resnorm/norm(b)), ...
                ' at iteration ', ...
                num2str(nit)])
            return
        end

        nit = nit + 1;
        if (nit > maxit)
            break;
        end
    end
    x0 = x;
end
