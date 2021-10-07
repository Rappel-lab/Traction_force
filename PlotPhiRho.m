close all;
% clc;
clear;

m      = 2^8;
n      = 2^8;
Lx     = 25;
Ly     = 25;
X      = linspace(-Lx, Lx, m+1); X = X(1:end-1); dx = 2*Lx/m;
Y      = linspace(-Ly, Ly, n+1); Y = Y(1:end-1); dy = 2*Ly/n;
[x,y]  = meshgrid(X,Y);
nu = 100;

epsilon=2; eps=epsilon;
 delta=1;
 band_width=5*dx; band_gap=5*dx;
 r0=3;
 % pst of channel
 r_floor=-3; r_ceil=3;
 % pst of band
band_floor_pst=-2;
band_ceil_pst = 20;

adh_ceil = 0;
adh_floor = 40;

eta_a0 = 1000;
gamma = 20;
proband = 2;
g = 1000;
A0 = 180;
Ma = 20;


xiM = 50;
exprhoM = 0.1;
xi0 = 0.5;

exprhoA = 0.25;
xiA = 0.1;


%% trajectory

ft = load('center_traj.txt');
fp = load('phi_profile_35.txt');
fa = load('rhoA_profile_35.txt');
fm = load('rhoM_profile_35.txt');
fux = load('ux_profile_35.txt');
fuy = load('uy_profile_35.txt');
% fp = load('phi_profile_w6.txt');
% fux = load('u_profile_w6.txt');
% fuy = load('w_profile_w6.txt');


% figure(1);  
% subplot(4,1,1);
% plot(ft(:,1),ft(:,2),'-r');title('mass center of x and z');
% subplot(4,1,2);plot(ft(:,1),ft(:,3),'--r', ft(:,1),ft(:,4),'--b');title('front and rear');ylabel('x-dir');
% subplot(4,1,3);plot(ft(:,1),ft(:,5),'k');title('cell length');ylabel('x-dir');
% subplot(4,1,4);plot(ft(:,1),ft(:,6)*dx*dy,'k');title('cell size');ylabel('size');


% figure(1);
phi=fp;
xc = fp.*x; 
xc_sin = fp/sum(fp(:)).*sin(x*pi/Lx);
xc_cos = fp/sum(fp(:)).*cos(x*pi/Lx);
mc_x = atan2(sum(xc_sin(:)), sum(xc_cos(:)))/pi*Lx;
shift_x_idx = -floor(mc_x/dx);

yc = fp.*y; 
yc_sin = fp/sum(fp(:)).*sin(y*pi/Ly);
yc_cos = fp/sum(fp(:)).*cos(y*pi/Ly);
mc_y = atan2(sum(yc_sin(:)), sum(yc_cos(:)))/pi*Ly;
shift_y_idx = -floor(mc_y/dx);

fp = circshift(fp,[shift_y_idx, shift_x_idx]);
fa = circshift(fa,[shift_y_idx, shift_x_idx]);
fm = circshift(fm,[shift_y_idx, shift_x_idx]);
fux = circshift(fux, [shift_y_idx, shift_x_idx]);
fuy = circshift(fuy, [shift_y_idx, shift_x_idx]);
 
% subplot(2,3,6);
phi_circ = contour(x,y,fp,[0.25,0.25]);

%% forces
% sum(fux(:))
% sum(fuy(:))


%% plot

% subplot(1,3,1);
% front = zeros(size(fp));
% front(x>=0) = 1;
% im_a = imagesc(X,Y,fp.*fa,'CDataMapping','scaled'); set(gca,'YDir','normal');colorbar;
% axis equal
%  axis([-25 25 -25 25]);
% %axis([-15 15 -15 15]);
% % figure(1);contour(X,Y,fa);axis equal; colorbar
% hold on; 
%  x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
% plot(x_circ,y_circ,'-r','linewidth',1); 
% set(gca,'Visible','off');
% plot([-15; -13],[15; 15],'w','linewidth',2);
% hold off
% saveas(gcf,'density.fig');


% im_m = imagesc(X,Y,fp.*fm,'CDataMapping','scaled'); set(gca,'YDir','normal');colorbar;
% axis equal
%  axis([-25 25 -25 25]);
% %axis([-15 15 -15 15]);
% % figure(1);contour(X,Y,fa);axis equal; colorbar
% hold on; 
%  x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
% plot(x_circ,y_circ,'-r','linewidth',1); 
% set(gca,'Visible','off');
% plot([-15; -13],[15; 15],'w','linewidth',2);
% hold off

figure(2);
display2Im(fp.*fa,fp.*fm);
axis equal
 axis([-25 25 -25 25]);
%axis([-15 15 -15 15]);
% figure(1);contour(X,Y,fa);axis equal; colorbar
hold on; 
 x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
plot(x_circ,y_circ,'-r','linewidth',1); 
set(gca,'Visible','off');
plot([-15; -13],[15; 15],'w','linewidth',2);
hold off
saveas(gcf,'density.fig');

% figure(3);
% imagesc(X,Y, 18*fp.*fp.*(1-fp).*(1-fp)/epsilon,'CDataMapping','scaled'); set(gca,'YDir','normal');colorbar;
% axis equal
% axis([-15 15 -15 15]);
% % figure(1);contour(X,Y,fa);axis equal; colorbar
% hold on; 
%  x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
% plot(x_circ,y_circ,'-r','linewidth',1); 
% set(gca,'Visible','off');
% plot([-9; -7],[9; 9],'w','linewidth',2);
% hold off

figure(4);
% fm = 2*fa - fp;
abs_vel = sqrt(fux.^2+fuy.^2).*fp;
fa(fp<0.5)=0;fm(fp<0.5)=0;

% imagesc(X, Y, fa-fm,'CDataMapping','scaled');
%caxis([0, max(abs_vel(:))]);
% caxis([-1,1]);
im1 = fa.*fp;
im2 = fm.*fp;
im1 = im1-min(im1(:));
% im1 = im1/max(im1(:));
im2 = im2-min(im2(:));
% im2 = im2/max(im2(:));
I = zeros([size(im1),3]);
I(:,:,1) = im1;
I(:,:,2) = im2;
%     image(I);
imagesc(X, Y, I); set(gca,'YDir','normal');
caxis([0,1]);

hold on
% ux = fux.*fp; uy = fuy.*fp;
ux = fux; uy = fuy;
sep = 15;
ux(fp<0.65)=0; uy(fp<0.65)=0;
quiver(X(1:sep:end), Y(1:sep:end), ux(1:sep:end,1:sep:end), uy(1:sep:end,1:sep:end),'w','AutoScale','on','MaxHeadSize',2,'linewidth',0.5); axis equal; 
set(gca,'YDir','normal');
% annotation('arrow',...
%         'Color', cmap(angled,:),...
%         'headStyle','cback1','HeadLength',50,'HeadWidth',headWidth);
axis([-20 20 -20 20]);

% phi_circ = contour(x,y,fp,[0.5,0.5]);
% x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
plot(x_circ,y_circ,'-w','linewidth',1); 


set(gca,'Visible','off');
plot([-15; -10],[15; 15],'w','linewidth',2);
colormap jet
hold off
saveas(gcf,'force_density_flow.fig');

figure(5);
% fm = 2*fa - fp;
% abs_vel = sqrt(fux.^2+fuy.^2).*fp;
% fa(fp<0.5)=0;fm(fp<0.5)=0;
imagesc(X, Y, abs_vel,'CDataMapping','scaled');
caxis([0, 0.1]);
% caxis([-1,1]);
hold on
% ux = fux.*fp; uy = fuy.*fp;
% sep = 15;
% normalize the velocity vector
ux_n = ux./abs_vel;
uy_n = uy./abs_vel;
ux_n(isnan(ux_n)) = 0;
uy_n(isnan(uy_n)) = 0;

ux_n(fp<0.5)=0; uy_n(fp<0.5)=0;
quiver(X(1:sep:end), Y(1:sep:end), ux_n(1:sep:end,1:sep:end), uy_n(1:sep:end,1:sep:end),'w','AutoScale','on','MaxHeadSize',2,'linewidth',0.5); axis equal; 
set(gca,'YDir','normal');
% annotation('arrow',...
%         'Color', cmap(angled,:),...
%         'headStyle','cback1','HeadLength',50,'HeadWidth',headWidth);
axis([-20 20 -20 20]);

% phi_circ = contour(x,y,fp,[0.5,0.5]);
% x_circ=phi_circ(1,2:end); y_circ=phi_circ(2,2:end);
plot(x_circ,y_circ,'-w','linewidth',1); 


set(gca,'Visible','off');
plot([-15; -10],[15; 15],'w','linewidth',2);
colormap jet
hold off
saveas(gcf,'force_stress_flow.fig');

%% calculate velocity
% figure(5);
% plot(ft(:,1),ft(:,2));
% t_itv=(ft(:,1)> 200 & ft(:,1)< 300);
% t=ft(t_itv,1);
% mx=ft(t_itv,2);
% % figure(4); plot(t,mx);
% p=polyfit(t,mx,1);
% vel=p(1);

