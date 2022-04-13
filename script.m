m = 0.127;
M = 1.206; 
J = 0.001; 
l = 0.178; 
K_f = 1.726; 
K_s = 4.487;
B_c = 5.4; 
B_p = 0.002;




pend = InvertedPendulum(m, M, J, l, K_f, K_s, B_c,  B_p);
pend.info();
pend = pend.set_trange([0 10]);
pend = pend.set_y0([-1 -1 -1 -1]');
pend.task13([-1, -1], [-10 -10]);