%% FUNCÃO CORRENTE E VORTICIDADE PARA UM ESCOAMENTO LAMINAR INCOMPRESSÍVEL 
%% NA CQTD

%--------------------------------------------------------------------------
% Código para simulação numérica do escoamento laminar incompressível 
% em uma cavidade quadrada (CQTD) utilizando a formulação de função 
% corrente e vorticidade.
%
% O código resolve as equações governantes por meio de diferenças finitas
% em uma malha bidimensional, atualizando iterativamente os campos de
% função corrente e vorticidade até a convergência.
%
% São aplicadas condições de contorno típicas do escoamento na cavidade,
% incluindo a parede superior com velocidade constante (condição de 
% escoamento de placa em movimento).
%
% Os resultados computados incluem os campos de velocidade nas direções x 
% e y, linhas de corrente, campo vetorial de velocidade e magnitude da 
% velocidade, que são visualizados por meio de gráficos.
%
% Parâmetros principais: dimensão da cavidade, propriedades do fluido 
% (densidade e viscosidade), passo de tempo, malha computacional, número 
% máximo de iterações e critério de tolerância para convergência.
%--------------------------------------------------------------------------

% Inserção dos dados de entrada
L      = 1;           % lado da cavidade;
U      = 1;           % velocidade u na parede superior da cavidade;
rho    = 1;           % densidade do fluido;
mu     = 0.01;        % viscosidade dinâmica do fluido para Re=100;
dt     = 0.001;       % tamanho do passo de tempo;
maxit  = 50000;       % máximo de iterações do método;
tol    = 1e-7;        % tolerância no método;

% Geração da malha computacional
ni    = 81;           % número de pontos da malha na direção x;
nj    = 81;           % número de pontos da malha na direção y;
h     = L/(ni-1);     % distância entre dois pontos adjacentes da malha;
x     = 0:h:L;        % primeira coordenada dos pontos da malha;
y     = 0:h:L;        % segunda coordenada dos pontos da malha;
ia    = 1:ni-2;       % coordenada i anterior;
i     = 2:ni-1;       % coordenada i corrente;
ip    = 3:ni;         % coordenada i posteior;
ja    = 1:nj-2;       % coordenada j anterior;
j     = 2:nj-1;       % coordenada j corrente;
jp    = 3:nj;         % coordenada j posterior

% Condições iniciais (a componente da i-ésima linha e j-ésima coluna da matriz corresponde ao valor da variável no ponto da malha (x_i,y_j))
Vo    = zeros(ni,nj); % vorticidade no iterado atual
Fc    = Vo;           % função corrente
Voa   = Vo;           % vorticidade no iterado anterior
u     = Vo;           % velocidade na direção x
v     = Vo;           % velocidade na direção y

% Desenvolvimento do processo iterativo
for iter = 1:maxit
    % condicões de contorno
    Vo(1:ni,nj) = -2*Fc(1:ni,nj-1)/(h^2) - U*2/h;  % parede superior;
    Vo(1:ni,1)  = -2*Fc(1:ni,2)/(h^2);             % parede inferior
    Vo(1,1:nj)  = -2*Fc(2,1:nj)/(h^2);             % parede esquerda;
    Vo(ni,1:nj) = -2*Fc(ni-1,1:nj)/(h^2);          % parede direita;

    % atualização da vorticidade
    Voa     = Vo;
    Vo(i,j) = Voa(i,j) + (-1*(Fc(i,jp)-Fc(i,ja))/(2*h) .* (Voa(ip,i)-Voa(ia,j))/(2*h) + (Fc(ip,j)-Fc(ia,j))/(2*h) .* (Voa(i,jp) -Voa(i,ja))/(2*h) + mu/rho *(Voa(ip,j)+Voa(ia,j)-4*Voa(i,j)+Voa(i,jp)+Voa(i,ja))/(h^2))*dt;

    % atualização da função corrente
    Fc(i,j) = (Vo(i,j)*h^2 + Fc(ip,j) + Fc(i,jp) + Fc(i,ja) + Fc(ia,j))/4;

    % verificação da convergência
    if iter > 10
        erro = max(max(Vo - Voa))
        if erro < tol
            break;
        end
    end
end

% Identificação da velocidade por meio da função corrente
u(2:ni-1,nj) = U;
u(i,j)       = (Fc(i,jp)-Fc(i,ja))/(2*h);
v(i,j)       = (-Fc(ip,j)+Fc(ia,j))/(2*h);


%************************************************************************%
% Plotagem dos gráficos/figuras:
%************************************************************************%

% Velocidade u 
cm = hsv(ceil(100/0.7));  
cm = flipud(cm(1:100,:));
figure(1);  
contourf(x,y,u',23,'LineColor','none');
xlabel('x');                      % eixo horizontal
ylabel('y')                       % eixo vertical
axis('equal',[0 L 0 L]);          % tamanho dos eixos
colormap(cm);                     % matriz de valores das cores do gráfico
colorbar('westoutside');          % inserir barra de cor à oeste do gráfico

% Velocidade v
cm = hsv(ceil(100/0.7));  
cm = flipud(cm(1:100,:));
figure(2);  
contourf(x,y,v',23,'LineColor','none');
xlabel('x');                      % eixo horizontal
ylabel('y')                       % eixo vertical
axis('equal',[0 L 0 L]);          % tamanho dos eixos
colormap(cm);                     % matriz de valores das cores do gráfico
colorbar('westoutside');          % inserir barra de cor à oeste do gráfico

% Linhas de corrente
N = 1000;  
xstart = max(x)*rand(N,1);  
ystart = max(y)*rand(N,1);
[X,Y]  = meshgrid(x,y);
figure(3);  
h=streamline(X,Y,u',v',xstart,ystart,[0.1, 200]);
xlabel('x');  
ylabel('y')
axis('equal',[0 L 0 L]);  
set(h,'color','k')

% Campo vetorial da velocidade
figure(4)
Len = sqrt(u.^2+v.^2+eps);
quiver(x,y,(u./Len)',(v./Len)',.4,'k-')
xlabel('x');                      % eixo horizontal
ylabel('y')                       % eixo vertical
axis('equal',[0 L 0 L]);          % tamanho dos eixos

% Magnitude da velocidade
figure(5)
Len = sqrt(u.^2+v.^2+eps);
contourf(x,y,Len',23,'LineColor','none');
xlabel('x');                      % eixo horizontal
ylabel('y')                       % eixo vertical
axis('equal',[0 L 0 L]);          % tamanho dos eixos
colormap(cm);                     % matriz de valores das cores do gráfico
colorbar('westoutside');          % inserir barra de cor à oeste do gráfico