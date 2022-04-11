
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
        x0
        t_range
    end
    properties(Constant = true)
        g = 9.8
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
            obj.C = control_matrix(obj.A, obj.B);
            obj.u = 0;
        end
        
        function info(obj)
            disp('Собственные числа:');
            disp(eigs(obj.A));
            disp('Ранг матрицы управляемости:');
            disp(rank(obj.C));
        end

        function res = Akkerman(obj, coefs_polynom) 
            o = zeros(1,length(obj.B)); 
            o(1, length(obj.B)) = 1; 
            res = -o * inv(obj.C) * polyvalm(coefs_polynom, obj.A); 
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
        function theta = stabilize(obj, eigenvals)
            coef_polynom = [1];
            for i = 1 : length(eigenvals)
                coef_polynom = conv(coef_polynom, [1 -eigenvals(i)]);
            end
            theta = obj.Akkerman(coef_polynom);
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
            [t, x] = ode45(@obj.nonlinear, obj.t_range, obj.x0);
            [t_lin, x_lin] = ode45(@obj.linear, obj.t_range, obj.x0);
            
            figure;
            tiledlayout(5,1)
            
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

        function task11b(obj, mu)
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
            [t, x] = ode45(@obj.nonlinear, obj.t_range, obj.x0);
            [t_lin, x_lin] = ode45(@obj.linear, obj.t_range, obj.x0);
            
            figure;
            tiledlayout(5,1)
            
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
        
        function task11c(obj, mu)
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
            [t, x] = ode45(@obj.nonlinear, obj.t_range, obj.x0);
            [t_lin, x_lin] = ode45(@obj.linear, obj.t_range, obj.x0);
            
            figure;
            tiledlayout(5,1)
            
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
    end

end
    

function C = control_matrix(A, b)
        C = [b];
        last = b;
        for i = 2 : length(A)
            C = [C, A * last];
            last = A * last;
        end
    end


