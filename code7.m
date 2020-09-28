close all;clear;clc;
warning off;
%% r vs pz ����������

%% define parameter
n=4;  %dimension
nx = 101; ny = 101;
x = linspace(0.8,2,nx);
y = linspace(0,0.25,ny);

%% Calculation
totalTimeVal = tic;
global RelTol AbsTol calTimeVal;
if false
    lambda = zeros(nx,ny,4);
    for i=1:nx
        for j=1:ny
            RelTol = 1e-6;
            AbsTol = 1e-10;
            flag = 1;
            while flag
                try
                    y0=[0;x(i);y(j);0];
                    calTimeVal = tic;
                    lambda(i,j,:)=cal(y0,@chao_rhs_ext);
                    flag = 0;
                    disp( ['Calculating ',num2str(((i-1)*ny+j)/nx/ny*100),'%... | RelTol=',num2str(RelTol)]);
                    toc(totalTimeVal);
                catch ME
                    RelTol = RelTol * sqrt(10);
                    AbsTol = AbsTol * sqrt(10);
                    disp(['Retry... | RelTol=',num2str(RelTol)]);
                end
            end
        end
    end
    save('fig6.mat','lambda','x','y');
end
%% Plot figures
load fig6.mat;
ze = 1e-3; % less than this value Lypunov exponent counts zero
lambda = sort(lambda,3);
lambda = lambda(:,:,3:4);
N = sum(lambda>ze,3);
N = setWhite(N,x,y);
fig = figure;
set(fig, 'position', get(0,'ScreenSize')); % Fullscreen
imagesc(N','XData',x,'YData',y);
hold on;
load fig6smaller.mat;
ze = 1e-3; % less than this value Lypunov exponent counts zero
lambda = sort(lambda,3);
lambda = lambda(:,:,3:4);
N = sum(lambda>ze,3);
N = setWhite(N,x,y);
% cmap = colormap(flipud(hot)); %
cmap = colormap( [1 1 1; autumn(256)] );
imagesc(N','XData',x,'YData',y);
load fig6.mat;
axis square;
axis xy;
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
colorbar;
set(gca,'FontSize',20);
xlabel('\rho');
ylabel('p_z');
hold on;
c = ['k','b','c','g','m'];
elist = [1/32,1/36,1/48,1/128, 1/256];
for i = 1:5
    e = elist(i);
    plotEnergy(e,c(i));
end
legend;

%% plot energy curve
function plotEnergy(e,c)
r = linspace(2/(1+sqrt(1-sqrt(2*e)*4)),2/(1+sqrt(1+sqrt(2*e)*4)),200);
pz = sqrt(2*e-(1./r-1./r.^2).^2);
plot(r,pz,'LineWidth',3,'Color',c,'DisplayName',['H=1/',num2str(round(1/e))]);
end

%% set energy > 1/32 white
function N = setWhite(N,x,y)
e = 1/32;
for i = 1:size(N,1)
    for j = 1:size(N,2)
        r = x(i);
        pz = sqrt(2*e-(1./r-1./r.^2).^2);
        if real(y(j)) <1e-8
            N(i,j) = 0;
        end
        if y(j) > pz
            N(i,j) = -0.01;
        end
    end
end
r = linspace(2/(1+sqrt(1-sqrt(2*e)*4)),2/(1+sqrt(1+sqrt(2*e)*4)),200);
end

%% subfunction
function f=chao_rhs_ext(t,X)
z=X(1); r=X(2); pz=X(3); pr=X(4);
Y=[X(5),X(9),X(13),X(17);
    X(6),X(10),X(14),X(18);
    X(7),X(11),X(15),X(19);
    X(8),X(12),X(16),X(20)];
f=zeros(20,1);
f(1:4)=[ pz, pr, (3*r*z*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(5/2), -(r/(r^2 + z^2)^(3/2) - 1/r)*(1/(r^2 + z^2)^(3/2) - (3*r^2)/(r^2 + z^2)^(5/2) + 1/r^2)];
Jac =[                                                                                                                                                                          0,                                                                                                                                                                                                    0, 1, 0;
    0,                                                                                                                                                                                                    0, 0, 1;
    (3*r*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(5/2) - (9*r^2*z^2)/(r^2 + z^2)^5 - (15*r*z^2*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(7/2), (3*z*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(5/2) - (15*r^2*z*(r/(r^2 + z^2)^(3/2) - 1/r))/(r^2 + z^2)^(7/2) + (3*r*z*(1/(r^2 + z^2)^(3/2) - (3*r^2)/(r^2 + z^2)^(5/2) + 1/r^2))/(r^2 + z^2)^(5/2), 0, 0;
    ((3*z)/(r^2 + z^2)^(5/2) - (15*r^2*z)/(r^2 + z^2)^(7/2))*(r/(r^2 + z^2)^(3/2) - 1/r) + (3*r*z*(1/(r^2 + z^2)^(3/2) - (3*r^2)/(r^2 + z^2)^(5/2) + 1/r^2))/(r^2 + z^2)^(5/2),                                             (r/(r^2 + z^2)^(3/2) - 1/r)*((9*r)/(r^2 + z^2)^(5/2) - (15*r^3)/(r^2 + z^2)^(7/2) + 2/r^3) - (1/(r^2 + z^2)^(3/2) - (3*r^2)/(r^2 + z^2)^(5/2) + 1/r^2)^2, 0, 0];
f(5:20)=Jac*Y;
global calTimeVal;
time = toc(calTimeVal);
if time > 30
    throw(MException('Id:id','message'));
end
end
%% calculate lyapunov exponents
function lambda = cal(y0,odeFunc)
% Calculation
global RelTol AbsTol;
opts = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
[~,Res]=lyapunov(4,odeFunc,@(odeFunc,ts,y0)ode45(odeFunc,ts,y0,opts),0,100,10000,y0,10);
lambda = Res(end,:);
end