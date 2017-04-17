testCase1 = load('../data/bcspwr03.mat');
A = testCase1.Problem.A;
b = zeros(size(A, 2), 1);
b(1) = 1;
m = 100;
maxit = 50000;

tol = 1e-3;
tic;
x1 = fgmres(A, b, m, tol, maxit);
toc;
disp(['Debug1: ', num2str(norm(b-A*x1, 2))]);
tic;
x2 = gmres(A, b, m, tol, maxit);
toc;
disp(['Debug2: ', num2str(norm(b-A*x2, 2))]);
