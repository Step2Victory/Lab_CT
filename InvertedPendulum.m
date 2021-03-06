
classdef InvertedPendulum
    
    properties
        m
        M
        J
        l
        K_f
        K_s
        B_c
        B_p
        A
        B
        C
        u
        u_observer
        L
        observer_C
        x0
        y0
        t_range
     
    end
    properties(Constant = true)
        g = 9.8
        eps = 10e-3
    end

    methods
        function obj = InvertedPendulum(m, M, J, l, K_f, K_s, B_c, B_p)
            obj.m = m;
            obj.M = M;
            obj.J = J;
            obj.l = l;
            obj.K_f = K_f;
            obj.K_s = K_s;
            obj.B_c = B_c;
            obj.B_p = B_p;
           
            obj.x0 = [0; -0.6; 0; 0];
            obj.y0 = obj.x0;
            obj.t_range = [0 2];

            A_0 = [m + M, -m * l;
                   -m * l, (J + m * l^2)];
            A_1 = [B_c, 0;
                    0, B_p];
            A_2 = [0, 0;
                    0, -m * obj.g * l];
            obj.A = [zeros(2), eye(2);
                     -inv(A_0) * A_2, -inv(A_0) * (A_1 + [K_f * K_s, 0; 0, 0])];
            obj.B = [0; 0; inv(A_0) * [K_f; 0]];
            obj.C = obj.control_matrix(obj.A, obj.B);

            obj.u = 0;
            obj.u_observer = 0;

            obj.L = 0;
            obj.observer_C = [1, 0, 0, 0; 0, 1, 0, 0]; 
        end
        
        function info(obj)
            disp('Собственные числа:');
            disp(eigs(obj.A));
            disp('Ранг матрицы управляемости:');
            disp(rank(obj.C));
        end

        function res = Akkerman(obj, A, B, coefs_polynom) 
            o = zeros(1,length(B)); 
            o(1, length(B)) = 1; 
            res = -o * inv(obj.control_matrix(A, B)) * polyvalm(coefs_polynom, A); 
        end

        function dxdt = nonlinear(obj, t, x)
            F = obj.K_f * (obj.u(t, x) - obj.K_s * x(3));
            
            D = F - obj.m * obj.l * sin(x(2)) * (x(4) ^ 2) - obj.B_c * x(3);
            E = obj.m * obj.g * obj.l * sin(x(2)) - obj.B_p * x(4);
            
            delta = (obj.m + obj.M) * (obj.J + obj.m * obj.l^2) - obj.m^2 * obj.l^2 * cos(x(2)) ^ 2;
            delta_1 = D * (obj.J + obj.m * obj.l ^ 2) + E * obj.m * obj.l * cos(x(2));
            delta_2 = (obj.m + obj.M) * E + D * obj.m * obj.l * cos(x(2));
            
            dxdt = zeros(4,1);
            dxdt(1) = x(3);
            dxdt(2) = x(4);
            dxdt(3) = delta_1 / delta;
            dxdt(4) = delta_2 / delta;
        end
    
        function dxdt = linear(obj, t, x)
            dxdt = obj.A*x + obj.B * obj.u(t, x);
        end

        function obj = addu(obj, u)
            obj.u = u;
        end

        function obj = set_trange(obj, new_trange)
            obj.t_range = new_trange;
        end
        
        function obj = set_x0(obj, new_x0)
            obj.x0 = new_x0;
        end
        
        function obj = set_y0(obj, new_y0)
            obj.x0 = new_y0;
        end
        
        function theta = stabilize(obj, eigenvals)
            coef_polynom = [1];
            for i = 1 : length(eigenvals)
                coef_polynom = conv(coef_polynom, [1 -eigenvals(i)]);
            end
            theta = obj.Akkerman(obj.A, obj.B, coef_polynom);
        end

        function plot_res1(obj)
            [t, x] = ode45(@obj.nonlinear, obj.t_range, obj.x0);
            [t_lin, x_lin] = ode45(@obj.linear, obj.t_range, obj.x0);
            
            figure;
            tiledlayout(4,1)
            
            x = real(x);
            x_lin = real(x_lin);
            ax1 = nexttile;
            plot(ax1, t, x(:,1), '-', t_lin, x_lin(:,1))
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, x(:,2), '-', t_lin, x_lin(:,2))
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, x(:,3), '-', t_lin, x_lin(:,3))
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, x(:,4), '-', t_lin, x_lin(:,4))
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
        end
        
        function plot_res2(obj)
            z0 = [obj.x0; obj.y0];
            [t, z] = ode45(@obj.nonlinear_and_observer, obj.t_range, z0);
            [t_lin, z_lin] = ode45(@obj.linear_and_observer, obj.t_range, z0);
            
            figure;
            tiledlayout(4,1)
            
            ax1 = nexttile;
            plot(ax1,t, z(:,1), '-', t_lin, z_lin(:,1))
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, z(:,2), '-', t_lin, z_lin(:,2))
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, z(:,3), '-', t_lin, z_lin(:,3))
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, z(:,4), '-', t_lin, z_lin(:,4))
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')

            figure;
            tiledlayout(4,1)

            ax5 = nexttile;
            plot(ax5, t, z(:,5), '-', t_lin, z_lin(:,5))
            ylabel(ax5, 'Observed x(t)', 'Interpreter','latex')
            
            ax6 = nexttile;
            plot(ax6, t, z(:,6), '-', t_lin, z_lin(:,6))
            ylabel(ax6,'Observed $\varphi(t)$', 'Interpreter','latex')
            
            ax7 = nexttile;
            plot(ax7, t, z(:,7), '-', t_lin, z_lin(:,7))
            ylabel(ax7,'Observed $\dot{x(t)}$', 'Interpreter','latex')
            
            ax8 = nexttile;
            plot(ax8, t, z(:,8), '-', t_lin, z_lin(:,8))
            ylabel(ax8,'Observed $\dot{\varphi(t)}$', 'Interpreter','latex')

            figure;
            tiledlayout(2,1)
            ax9 = nexttile;
            h_lin = zeros(1, length(t_lin));
            for i = 1 : length(t_lin)
                h_lin(i) = norm(z_lin(i,1:4) - z_lin(i, 5:8));
            end
            plot(ax9, t_lin, h_lin)
            ylabel(ax9,'Невязка для линейной системы', 'Interpreter','latex')
            ax10 = nexttile;
            h = zeros(1, length(t));
            for i = 1 : length(t)
                h(i) = norm(z(i,1:4) - z(i, 5:8));
            end
            plot(ax10, t, h)
            ylabel(ax10,'Невязка для нелинейной системы', 'Interpreter','latex')

        end

        function task11a(obj, mu)
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu;
                end
            end
            theta = obj.stabilize(eigenvals);
            obj.u = @(t, x) (theta * x);

            obj.plot_res1();
        end

        function task11bc(obj, mu)
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            theta = obj.stabilize(eigenvals);
            obj.u = @(t, x) (theta * x);

            obj.plot_res1();
        end
           
        function task12ab(obj, mu, mu_observer)

            theta1 = obj.shift_eig(obj.A', -obj.observer_C(1,:)', 0, mu_observer(1));
            A_next = obj.A' - obj.observer_C(1,:)' * theta1;
            theta2 = obj.shift_eig(A_next, obj.observer_C(2,:)', 6.5418, mu_observer(2));
            theta = [theta1; theta2];
            obj.L = theta';
            
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            theta = obj.stabilize(eigenvals);
            obj.u = @(t, z) (theta * z(1:4));
            obj.u_observer = @(t, z) (theta * z(5:8));
            
            obj.plot_res2();
        end

        function task13(obj, mu, mu_observer)

            theta1 = obj.shift_eig(obj.A', -obj.observer_C(1,:)', 0, mu_observer(1));
            A_next = obj.A' - obj.observer_C(1,:)' * theta1;
            theta2 = obj.shift_eig(A_next, obj.observer_C(2,:)', 6.5418, mu_observer(2));
            theta = [theta1; theta2];
            obj.L = theta';
            eigs(obj.A - obj.L * obj.observer_C)
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            theta = obj.stabilize(eigenvals);
            obj.u = @(t, z) (theta * z(5:8));
            obj.u_observer = @(t, z) (theta * z(5:8));
            
            obj.plot_res2();

        end
       
        function task21a(obj, mu)
            h = 0.1;
            n_steps = (obj.t_range(2) - obj.t_range(1)) / h;
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu;
                end
            end
            
            x_res = [obj.x0'];
            x_lin_res = [obj.x0'];
            t_res = [0];
            t_lin_res = [0];
            u_res = [];
            u_lin_res = [];
            
            A_tmp = obj.A;
            B_tmp = obj.B;
            obj.A = expm(A_tmp * h);
            obj.B = integral(@(s) expm(A_tmp * s) * B_tmp, 0, h, 'ArrayValued', true);

            theta = obj.stabilize(exp(eigenvals * h));


            obj.A = A_tmp;
            obj.B = B_tmp;
            
            for i = 1 : n_steps

                obj.u = @(t, x) theta * x_res(end, :)';
                [t, x] = ode45(@obj.nonlinear, [(i - 1) * h, i * h], x_res(end, :)');
                x = real(x);
                x_res = [x_res; x];
                u_res = [u_res; obj.u()];
                t_res = [t_res; t];
                
                obj.u = @(t, x) theta * x_lin_res(end, :)';
                [t_lin, x_lin] = ode45(@obj.linear, [(i - 1) * h, i * h], x_lin_res(end, :)');
                x_lin = real(x_lin);
                x_lin_res = [x_lin_res; x_lin];
                u_lin_res = [u_lin_res; obj.u()];
                t_lin_res = [t_lin_res; t_lin];
                
            end
            
            t = t_res;
            t_lin = t_lin_res;
            x = x_res;
            x_lin = x_lin_res;
            
            figure;
            tiledlayout(5,1)
            
            ax1 = nexttile;
            plot(ax1,t, x(:,1), "b", t_lin, x_lin(:,1), "r")
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, x(:,2), "b", t_lin, x_lin(:,2), "r")
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, x(:,3), "b", t_lin, x_lin(:,3), "r")
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, x(:,4), "b", t_lin, x_lin(:,4), "r")
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax4_1 = nexttile;
            for i = 1 : n_steps
                plot(ax4_1, [(i - 1) * h, i * h], [u_res(i), u_res(i)], "b", [(i - 1) * h, i * h], [u_lin_res(i), u_lin_res(i)], "r")
                hold on;
            end
            ylabel(ax4_1,'$u$', 'Interpreter','latex')

        end

        function task21bc(obj, mu)
            h = 0.1;
            n_steps = (obj.t_range(2) - obj.t_range(1)) / h;
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            
            x_res = [obj.x0'];
            x_lin_res = [obj.x0'];
            t_res = [0];
            t_lin_res = [0];
            u_res = [];
            u_lin_res = [];
            
            A_tmp = obj.A;
            B_tmp = obj.B;
            obj.A = expm(A_tmp * h);
            obj.B = integral(@(s) expm(A_tmp * s) * B_tmp, 0, h, 'ArrayValued', true);

            theta = obj.stabilize(exp(eigenvals * h));


            obj.A = A_tmp;
            obj.B = B_tmp;
            
            for i = 1 : n_steps

                obj.u = @(t, x) theta * x_res(end, :)';
                [t, x] = ode45(@obj.nonlinear, [(i - 1) * h, i * h], x_res(end, :)');
                x = real(x);
                x_res = [x_res; x];
                u_res = [u_res; obj.u()];
                t_res = [t_res; t];
                
                obj.u = @(t, x) theta * x_lin_res(end, :)';
                [t_lin, x_lin] = ode45(@obj.linear, [(i - 1) * h, i * h], x_lin_res(end, :)');
                x_lin = real(x_lin);
                x_lin_res = [x_lin_res; x_lin];
                u_lin_res = [u_lin_res; obj.u()];
                t_lin_res = [t_lin_res; t_lin];
                
            end
            
            t = t_res;
            t_lin = t_lin_res;
            x = x_res;
            x_lin = x_lin_res;
            
            figure;
            tiledlayout(5,1)
            
            ax1 = nexttile;
            plot(ax1,t, x(:,1), "b", t_lin, x_lin(:,1), "r")
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, x(:,2), "b", t_lin, x_lin(:,2), "r")
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, x(:,3), "b", t_lin, x_lin(:,3), "r")
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, x(:,4), "b", t_lin, x_lin(:,4), "r")
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax4_1 = nexttile;
            for i = 1 : n_steps
                plot(ax4_1, [(i - 1) * h, i * h], [u_res(i), u_res(i)], "b", [(i - 1) * h, i * h], [u_lin_res(i), u_lin_res(i)], "r")
                hold on;
            end
            ylabel(ax4_1,'$u$', 'Interpreter','latex')

        end

        function task22ab(obj, mu, mu_observer)            
            h = 0.1;
            n_steps = (obj.t_range(2) - obj.t_range(1)) / h;
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            
            z_res = [obj.x0', obj.y0'];
            t_res = [0];
            z_lin_res = [obj.x0', obj.y0'];
            t_lin_res = [0];
            u_res = [];
            u_observer_res = [];
            u_lin_res = [];
            u_lin_observer_res = [];
            
            A_tmp = obj.A;
            B_tmp = obj.B;
            obj.A = expm(A_tmp * h);
            obj.B = integral(@(s) expm(A_tmp * s) * B_tmp, 0, h, 'ArrayValued', true);

            theta1 = obj.shift_eig(obj.A', -obj.observer_C(1,:)', 1, exp(mu_observer(1) * h));
            A_next = obj.A' - obj.observer_C(1,:)' * theta1;
            theta2 = obj.shift_eig(A_next, obj.observer_C(2,:)', exp(6.5418 * h), exp(mu_observer(2) * h));
            theta = [theta1; theta2];
            obj.L = theta';

            theta = obj.stabilize(exp(eigenvals * h));
            

            obj.A = A_tmp;
            obj.B = B_tmp;  
            
            for i = 1 : n_steps
                
                obj.u = @(t, z) (theta * z_res(end, 1:4)');
                obj.u_observer = @(t, z) (theta * z_res(end, 5:8)');
                u_res = [u_res; obj.u()];
                u_observer_res = [u_observer_res; obj.u_observer()];
                
                [t, z] = ode45(@obj.nonlinear_and_observer, [(i - 1) * h, i * h], z_res(end, :)');
                
                obj.u = @(t, z) (theta * z_lin_res(end, 1:4)');
                obj.u_observer = @(t, z) (theta * z_lin_res(end, 5:8)');
                u_lin_res = [u_lin_res; obj.u()];
                u_lin_observer_res = [u_lin_observer_res; obj.u_observer()];
                
           
                [t_lin, z_lin] = ode45(@obj.linear_and_observer, [(i - 1) * h, i * h], z_lin_res(end, :)');
                
                z_res = [z_res; z];
                t_res = [t_res; t];
                z_lin_res = [z_lin_res; z_lin];
                t_lin_res = [t_lin_res; t_lin];
                
                
            end
            
            t = t_res;
            t_lin = t_lin_res;
            z = z_res;
            z_lin = z_lin_res;
            
            figure;
            tiledlayout(5,1)
            
            ax1 = nexttile;
            plot(ax1,t, z(:,1), "b", t_lin, z_lin(:,1), "r")
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, z(:,2), "b", t_lin, z_lin(:,2), "r")
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, z(:,3), "b", t_lin, z_lin(:,3), "r")
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, z(:,4), "b", t_lin, z_lin(:,4), "r")
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax4_1 = nexttile;
            for i = 1 : n_steps
                plot(ax4_1, [(i - 1) * h, i * h], [u_res(i), u_res(i)], "b", [(i - 1) * h, i * h], [u_lin_res(i), u_lin_res(i)], "r")
                hold on;
            end
            ylabel(ax4_1,'$u$', 'Interpreter','latex')

            figure;
            tiledlayout(5,1)

            ax5 = nexttile;
            plot(ax5, t, z(:,5), "b", t_lin, z_lin(:,5), "r")
            ylabel(ax5, 'Observed x(t)', 'Interpreter','latex')
            
            ax6 = nexttile;
            plot(ax6, t, z(:,6), "b", t_lin, z_lin(:,6), "r")
            ylabel(ax6,'Observed $\varphi(t)$', 'Interpreter','latex')
            
            ax7 = nexttile;
            plot(ax7, t, z(:,7), "b", t_lin, z_lin(:,7), "r")
            ylabel(ax7,'Observed $\dot{x(t)}$', 'Interpreter','latex')
            
            ax8 = nexttile;
            plot(ax8, t, z(:,8), "b", t_lin, z_lin(:,8), "r")
            ylabel(ax8,'Observed $\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax8_1 = nexttile;
            for i = 1 : n_steps
                plot(ax8_1, [(i - 1) * h, i * h], [u_observer_res(i), u_observer_res(i)], "b", [(i - 1) * h, i * h], [u_lin_observer_res(i), u_lin_observer_res(i)], "r")
                hold on;
            end
            ylabel(ax8_1,'Observed $u$', 'Interpreter','latex')


            figure;
            tiledlayout(2,1)
            ax9 = nexttile;
            h_lin = zeros(1, length(t_lin));
            for i = 1 : length(t_lin)
                h_lin(i) = norm(z_lin(i,1:4) - z_lin(i, 5:8));
            end
            plot(ax9, t_lin, h_lin)
            ylabel(ax9,'Невязка для линейной системы', 'Interpreter','latex')
            ax10 = nexttile;
            h = zeros(1, length(t));
            for i = 1 : length(t)
                h(i) = norm(z(i,1:4) - z(i, 5:8));
            end
            plot(ax10, t, h)
            ylabel(ax10,'Невязка для нелинейной системы', 'Interpreter','latex')
        end
        
        function task23(obj, mu, mu_observer)            
            h = 0.1;
            n_steps = (obj.t_range(2) - obj.t_range(1)) / h;
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            
            z_res = [obj.x0', obj.y0'];
            t_res = [0];
            z_lin_res = [obj.x0', obj.y0'];
            t_lin_res = [0];
            u_res = [];
            u_observer_res = [];
            u_lin_res = [];
            u_lin_observer_res = [];
            
            A_tmp = obj.A;
            B_tmp = obj.B;
            obj.A = expm(A_tmp * h);
            obj.B = integral(@(s) expm(A_tmp * s) * B_tmp, 0, h, 'ArrayValued', true);

            theta1 = obj.shift_eig(obj.A', -obj.observer_C(1,:)', 1, exp(mu_observer(1) * h));
            A_next = obj.A' - obj.observer_C(1,:)' * theta1;
            theta2 = obj.shift_eig(A_next, obj.observer_C(2,:)', exp(6.5418 * h), exp(mu_observer(2) * h));
            theta = [theta1; theta2];
            obj.L = theta';

            theta = obj.stabilize(exp(eigenvals * h));
            

            obj.A = A_tmp;
            obj.B = B_tmp;  
            
            for i = 1 : n_steps
                
                obj.u = @(t, z) (theta * z_res(end, 5:8)');
                obj.u_observer = @(t, z) (theta * z_res(end, 5:8)');
                u_res = [u_res; obj.u()];
                u_observer_res = [u_observer_res; obj.u_observer()];
                
                [t, z] = ode45(@obj.nonlinear_and_observer, [(i - 1) * h, i * h], z_res(end, :)');
                
                obj.u = @(t, z) (theta * z_lin_res(end, 5:8)');
                obj.u_observer = @(t, z) (theta * z_lin_res(end, 5:8)');
                u_lin_res = [u_lin_res; obj.u()];
                u_lin_observer_res = [u_lin_observer_res; obj.u_observer()];
                
           
                [t_lin, z_lin] = ode45(@obj.linear_and_observer, [(i - 1) * h, i * h], z_lin_res(end, :)');
                
                z_res = [z_res; z];
                t_res = [t_res; t];
                z_lin_res = [z_lin_res; z_lin];
                t_lin_res = [t_lin_res; t_lin];
                
                
            end
            
            t = t_res;
            t_lin = t_lin_res;
            z = z_res;
            z_lin = z_lin_res;
            
            figure;
            tiledlayout(5,1)
            
            ax1 = nexttile;
            plot(ax1,t, z(:,1), "b", t_lin, z_lin(:,1), "r")
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, z(:,2), "b", t_lin, z_lin(:,2), "r")
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, z(:,3), "b", t_lin, z_lin(:,3), "r")
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, z(:,4), "b", t_lin, z_lin(:,4), "r")
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax4_1 = nexttile;
            for i = 1 : n_steps
                plot(ax4_1, [(i - 1) * h, i * h], [u_res(i), u_res(i)], "b", [(i - 1) * h, i * h], [u_lin_res(i), u_lin_res(i)], "r")
                hold on;
            end
            ylabel(ax4_1,'$u$', 'Interpreter','latex')

            figure;
            tiledlayout(5,1)

            ax5 = nexttile;
            plot(ax5, t, z(:,5), "b", t_lin, z_lin(:,5), "r")
            ylabel(ax5, 'Observed x(t)', 'Interpreter','latex')
            
            ax6 = nexttile;
            plot(ax6, t, z(:,6), "b", t_lin, z_lin(:,6), "r")
            ylabel(ax6,'Observed $\varphi(t)$', 'Interpreter','latex')
            
            ax7 = nexttile;
            plot(ax7, t, z(:,7), "b", t_lin, z_lin(:,7), "r")
            ylabel(ax7,'Observed $\dot{x(t)}$', 'Interpreter','latex')
            
            ax8 = nexttile;
            plot(ax8, t, z(:,8), "b", t_lin, z_lin(:,8), "r")
            ylabel(ax8,'Observed $\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax8_1 = nexttile;
            for i = 1 : n_steps
                plot(ax8_1, [(i - 1) * h, i * h], [u_observer_res(i), u_observer_res(i)], "b", [(i - 1) * h, i * h], [u_lin_observer_res(i), u_lin_observer_res(i)], "r")
                hold on;
            end
            ylabel(ax8_1,'Observed $u$', 'Interpreter','latex')


            figure;
            tiledlayout(2,1)
            ax9 = nexttile;
            h_lin = zeros(1, length(t_lin));
            for i = 1 : length(t_lin)
                h_lin(i) = norm(z_lin(i,1:4) - z_lin(i, 5:8));
            end
            plot(ax9, t_lin, h_lin)
            ylabel(ax9,'Невязка для линейной системы', 'Interpreter','latex')
            ax10 = nexttile;
            h = zeros(1, length(t));
            for i = 1 : length(t)
                h(i) = norm(z(i,1:4) - z(i, 5:8));
            end
            plot(ax10, t, h)
            ylabel(ax10,'Невязка для нелинейной системы', 'Interpreter','latex')
        end
        
        function task3(obj, Q, R)
            [X, K, L] = icare(obj.A, obj.B, diag(Q), R, 0, eye(4), 0);
            obj.u = @(t, x) (-inv(R) * obj.B' * X * x);
            [t, x] = ode45(@obj.nonlinear, obj.t_range, obj.x0);
            [t_lin, x_lin] = ode45(@obj.linear, obj.t_range, obj.x0);
            
            figure;
            tiledlayout(4,1)
            
            ax1 = nexttile;
            plot(ax1, t, x(:,1), '-', t_lin, x_lin(:,1))
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, x(:,2), '-', t_lin, x_lin(:,2))
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, x(:,3), '-', t_lin, x_lin(:,3))
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, x(:,4), '-', t_lin, x_lin(:,4))
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
        end

        function task4(obj, Ex, X, W, mu)
            
            h = 0.1;
            n_steps = (obj.t_range(2) - obj.t_range(1)) / h;
            eigenvals = eigs(obj.A);
            for i = 1 : length(eigenvals)
                if eigenvals(i) > 0
                    eigenvals(i) = mu(1);
                elseif eigenvals(i) == 0
                    eigenvals(i) = mu(2);
                end
            end
            Ad = expm(obj.A * h);
            Bd = integral(@(s) expm(obj.A * s) * obj.B, 0, h, 'ArrayValued', true);
            A_tmp = obj.A;
            B_tmp = obj.B;
            obj.A = Ad;
            obj.B = Bd;
            theta = obj.stabilize(exp(eigenvals * h));
            
            obj.A = A_tmp;
            obj.B = B_tmp;
            
            A_c = Ad + Bd * theta;
            x_plus = Ex;
            E_plus = X;
            rng('shuffle');
            x_rnd = [normrnd(Ex(1), X(1,1));normrnd(Ex(2), X(2,2));normrnd(Ex(3), X(3,3));normrnd(Ex(4), X(4,4))];
            x_res = [x_rnd'];
            t_res = [0];  
            u_res = [];
       
            for i = 1 : n_steps
                
                obj.u = @(t, x) (theta * x_plus);
                [t, x] = ode45(@obj.nonlinear, [(i - 1) * h, i * h], x_res(end, :)');
                      
                u_res = [u_res; obj.u()];
                x_res = [x_res; x];
                t_res = [t_res; t];
                
                w = [normrnd(0, W(1,1)); normrnd(0, W(2,2))];

                y = obj.observer_C * x_res(end,:)' + w;

                x_minus = A_c * x_plus;
                E_minus = obj.A * E_plus * obj.A';

                inv_E_plus = inv(E_minus) + obj.observer_C' * inv(W) * obj.observer_C;
                E_plus = inv(inv_E_plus);
                x_plus = x_minus + E_plus * obj.observer_C' * inv(W) * (y - obj.observer_C * x_minus);
                
            end
            
            x = x_res;
            t = t_res;

            figure;
            tiledlayout(5,1)
            
            ax1 = nexttile;
            plot(ax1, t, x(:,1), "b")
            ylabel(ax1, 'x(t)', 'Interpreter','latex')
            
            ax2 = nexttile;
            plot(ax2, t, x(:,2), "b")
            ylabel(ax2,'$\varphi(t)$', 'Interpreter','latex')
            
            ax3 = nexttile;
            plot(ax3, t, x(:,3), "b")
            ylabel(ax3,'$\dot{x(t)}$', 'Interpreter','latex')
            
            ax4 = nexttile;
            plot(ax4, t, x(:,4), "b")
            ylabel(ax4,'$\dot{\varphi(t)}$', 'Interpreter','latex')
            
            ax4_1 = nexttile;
            for i = 1 : n_steps
                plot(ax4_1, [(i - 1) * h, i * h], [u_res(i), u_res(i)], "b")
                hold on;
            end
            ylabel(ax4_1,'$u$', 'Interpreter','latex')
            
        end
        
        function dzdt = linear_and_observer(obj, t, z)
            dx1dt = obj.A*z(1:4) + obj.B * obj.u(t, z);
            y = obj.observer_C * z(1:4);
            dx2dt = obj.A * z(5:8) + obj.B * obj.u_observer(t,z) + obj.L * (y - obj.observer_C * z(5:8));
            dzdt = [dx1dt;dx2dt];
        end

        function dzdt = nonlinear_and_observer(obj, t, z)           
            F = obj.K_f * (obj.u(t, z) - obj.K_s * z(3));
            
            D = F - obj.m * obj.l * sin(z(2)) * (z(4) ^ 2) - obj.B_c * z(3);
            E = obj.m * obj.g * obj.l * sin(z(2)) - obj.B_p * z(4);
            
            delta = (obj.m + obj.M) * (obj.J + obj.m * obj.l^2) - obj.m^2 * obj.l^2 * cos(z(2)) ^ 2;
            delta_1 = D * (obj.J + obj.m * obj.l ^ 2) + E * obj.m * obj.l * cos(z(2));
            delta_2 = (obj.m + obj.M) * E + D * obj.m * obj.l * cos(z(2));
            
            dx1dt = zeros(4,1);
            dx1dt(1) = z(3);
            dx1dt(2) = z(4);
            dx1dt(3) = delta_1 / delta;
            dx1dt(4) = delta_2 / delta;
            y = obj.observer_C * z(1:4);
            dx2dt = obj.A * z(5:8) + obj.B * obj.u_observer(t,z) + obj.L * (y - obj.observer_C * z(5:8));
            dzdt = [dx1dt;dx2dt];
        end

        function C = control_matrix(obj, A, b)
            C = [b];
            last = b;
            for i = 2 : length(A)
                C = [C, A * last];
                last = A * last;
            end
        end
            
        function calman(obj, A, b)
            C = obj.control_matrix(A, b);
            if rank(C) == length(A)
                disp('Система полностью управляема')
            else
                bin = orth(C);
                bin = [bin, null(transpose(bin))];
                P_inv = inv(bin);
                disp(vpa(P_inv * A * bin, 4));
                disp(vpa(P_inv * b, 4));
            end
        end
        
        function theta = shift_eig(obj, A, B, mu1, mu2)
            [V, D, W] = eig(A);
            bin = [];
            for i = 1:length(V)
                if abs(D(i, i) - mu1) < obj.eps
                    bin = [bin; W(:,i)'];
                end
            end
            bin = [bin; null(bin)'];
            P_inv = bin;
            P = inv(bin);
            theta = zeros(1,length(B));
            new_A = P_inv * A * P;
            new_B = P_inv * B;
            theta(1,1) = (mu2 - new_A(1,1)) / new_B(1,1);
            theta = theta * P;
        end

    end
end
    




