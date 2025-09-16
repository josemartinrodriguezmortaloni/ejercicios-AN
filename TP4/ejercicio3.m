function ejercicio3
  % Ejercicio 3 - M¨¦todo de Euler
  clear; clc;

  % Definir la EDO: du/dt = (t - u)/2
  f = @(t, u) 0.5 * (t - u);

  % Soluci¨®n exacta
  u_exact = @(t) 6*exp(-t/2) - 2 + t;

  % Par¨¢metros
  t0 = 0;
  tf = 10;  % Extendido para ver comportamiento
  u0 = 4;

  % Diferentes valores de Delta t
  dt_values = [0.25, 0.10, 0.01];
  colors = {'r.-', 'b.-', 'g-'};

  figure;

  for i = 1:length(dt_values)
      dt = dt_values(i);
      t = t0:dt:tf;
      N = length(t);

      % Aplicar m¨¦todo de Euler
      u = zeros(1, N);
      u(1) = u0;

      for n = 1:N-1
          u(n+1) = u(n) + dt * f(t(n), u(n));
      end

      % Calcular error en t=0.5
      idx_05 = find(abs(t - 0.5) < 1e-10);
      if ~isempty(idx_05)
          error_05 = abs(u_exact(0.5) - u(idx_05(1)));
          fprintf('dt = %.2f: u(0.5) = %.4f, Error = %.4f\n', ...
                  dt, u(idx_05(1)), error_05);
      end

      % Graficar
      if i <= 2  % Solo graficar dt=0.25 y dt=0.10
          plot(t, u, colors{i}, 'MarkerSize', 4, ...
               'DisplayName', sprintf('dt=%.2f', dt));
          hold on;
      else  % Para dt=0.01 usar l¨ªnea continua
          plot(t, u, colors{i}, 'LineWidth', 1.5, ...
               'DisplayName', sprintf('dt=%.2f', dt));
      end
  end

  % Graficar soluci¨®n exacta
  t_exact = t0:0.001:tf;
  u_ex = u_exact(t_exact);
  plot(t_exact, u_ex, 'k-', 'LineWidth', 2, ...
       'DisplayName', 'Exacta');

  % Configurar gr¨¢fico
  xlabel('t');
  ylabel('u(t)');
  title('Ejercicio 3: Comparaci¨®n de soluciones');
  legend('Location', 'best');
  grid on;
  xlim([0, 10]);
  ylim([2, 9]);

  % Gr¨¢fico del error para dt=0.01
  figure;
  dt = 0.01;
  t = t0:dt:tf;
  N = length(t);
  u = zeros(1, N);
  u(1) = u0;

  for n = 1:N-1
      u(n+1) = u(n) + dt * f(t(n), u(n));
  end

  error_func = abs(u_exact(t) - u);
  plot(t, error_func*1000, 'b-', 'LineWidth', 1.5);
  xlabel('t');
  ylabel('Error (¡Á10^{-3})');
  title('Funci¨®n error con \Delta t = 0.01');
  grid on;
  xlim([0, 10]);

  % An¨¢lisis del ETL
  fprintf('\nAn¨¢lisis del Error de Truncamiento Local:\n');
  fprintf('Segunda derivada: d^2u/dt^2 = (1-u)/4\n');
  fprintf('En t=0: d^2u/dt^2|? = (1-4)/4 = -0.75\n');
  fprintf('ETL = (dt^2/2) * (-0.75) = -0.375 * dt^2\n');
endfunction
