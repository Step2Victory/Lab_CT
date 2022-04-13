m = 0.127;
M = 1.206; 
J = 0.001; 
l = 0.178; 
K_f = 1.726; 
K_s = 4.487;
B_c = 5.4; 
B_p = 0.002;




pend = InvertedPendulum(m, M, J, l, K_f, K_s, B_c,  B_p);

u = @(t, x) (x(1) ^ 2);
pend = pend.addu(u);

[t, x] = ode45(@pend.linear, [0 2], [0; -0.6; 0; 0]);
[t_lin, x_lin] = ode45(@pend.nonlinear, [0 2], [0; -0.6; 0; 0]);

figure;

tiledlayout(1,2)

ax1 = nexttile;
plot(ax1, t, x(:,1), '-', t_lin, x_lin(:,1))
ylabel(ax1, 'x(t)', 'Interpreter','latex')

ax2 = nexttile;
u = [];
for i = 1: length(x) 
    u = [u, pend.u(t(i),x(i,:))];
end
plot(ax2, t, u, '-')
ylabel(ax2,'$u(t)$', 'Interpreter','latex')