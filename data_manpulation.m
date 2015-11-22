x= reshape(x,[25,75]) ;
y= reshape(y,[25,75]) ;
density= reshape(density,[25,75]);
densityu= reshape(densityu,[25,75]);
densityv= reshape(densityv,[25,75]);
densityw= reshape(densityw,[25,75]);
energy= reshape(energy,[25,75]);

for i = 1:25
    for j =	1:75	
        u(i,j) = densityu(i,j)/density(i,j) ; 
        v(i,j) = densityv(i,j)/density(i,j) ; 
        w(i,j) = densityw(i,j)/density(i,j) ; 
        pressure(i,j) = 0.4*(energy(i,j) - 0.5*density(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)));
        temperature(i,j) = pressure(i,j) /(287.14*density(i,j)) ; 
        velocity(i,j) = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j)+w(i,j)*w(i,j)) ;
        mach(i,j) = velocity(i,j) / sqrt(1.4*287.14*temperature(i,j)) ;
		
		% non_dimensinalisation
		non_dim_x(i,j) = x(i,j)/(y(25,75)-y(1,1));
		non_dim_y(i,j) = y(i,j)/(y(25,75)-y(1,1));
		non_dim_den(i,j) = density(i,j)/density(1,1);
		non_dim_vel(i,j) = velocity(i,j)/velocity(1,1);
		non_dim_press(i,j) = pressure(i,j)/pressure(1,1);
		non_dim_temp(i,j) = temperature(i,j)/temperature(1,1) ;

    end
end
% % plotting
% k=20 ;
% % for k=1:75 ;
% 	if rem(k,10) ==0 
% 		plot(non_dim_x(k,:),non_dim_den  (k,:),'-o'); hold on ;
% 		plot(non_dim_x(k,:),non_dim_vel  (k,:),'-o'); hold on ;
% 		plot(non_dim_x(k,:),non_dim_press(k,:),'-o'); hold on ;
% 		plot(non_dim_x(k,:),non_dim_temp (k,:),'-o'); hold on ;
% 		plot(non_dim_x(k,:),mach         (k,:),'-o'); hold on ;
% 	end
% % end

% h = surf(non_dim_x,non_dim_y) ;

h = surf(non_dim_x,non_dim_y,mach) ;
% h = surf(non_dim_x,non_dim_y,non_dim_den) ;
% h = surf(non_dim_x,non_dim_y,non_dim_press) ;
% h = surf(non_dim_x,non_dim_y,non_dim_temp) ;

colormap hsv
colorbar
hcb=colorbar
shading interp;

t = title(hcb,'\bf {Mach}') ;
xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf Mach');
title(' \bf Mach Number (M), flow over a bump')

% t = title(hcb,'$\frac{\rho}{\rho_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf \rho / \rho_{\infty}');
% title(' \bf Nondimensional density(\rho/ \rho_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{V}{V_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf V/V{\infty}');
% title(' \bf Nondimensional velocity(V/V_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{P}{P_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf P/P{\infty}');
% title(' \bf Nondimensional pressure(P/P_{\infty}), flow over a bump')

% t = title(hcb,'$\frac{T}{T_{\infty}}$') ;
% set(t,'Interpreter','Latex');
% xlabel('\bf{X(m)}'); ylabel('\bf y(m)'); zlabel('\bf T/T{\infty}');
% title(' \bf Nondimensional temperature (T/T_{\infty}), flow over a bump')