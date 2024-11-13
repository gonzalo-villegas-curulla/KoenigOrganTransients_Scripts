% Parameters
A = 1;                  % Maximum amplitude
T = 1e-1;                  % Desired period
omega = 2 * pi / T;     % Angular frequency
alpha = 10;              % Damping factor for amplitude modulation
t = linspace(0, 50*T, 1000);  

% Signal definition
s = 2 * A * cos(omega * t) .* cosh(alpha * t);


lambda1 = alpha + i*omega;
lambda2 = -alpha + i*omega;
s1 = A*exp(lambda1*t);
s2 = A*exp(lambda2*t);


figure(1); clf;
plot(t, s);


figure(2);clf; 
% plot(t,real(s1) );
hold on;
plot(t, real( A*exp(-alpha*t).*exp(-i*omega*t)));
% plot(t,s);

grid on;

%%



figure(3); clf;

s21 = ( COMP(2,1:hL) )  ;
s22 = ( envel_second );

Ptar = max(envel_second);%mean(envel_second(end- fix(fs*0.5): end));
idx20 = find( envel_second/Ptar<0.2,1,'last');
idx80 = find( envel_second/Ptar>0.8, 1, 'first');


subplot(211);
plot(log(time), log(s21) );
hold on;
plot(log(time), log(s22)  ) ;

plot(log([1,1]*time(idx20)),[-5,4 ],'--k');
plot(log([1,1]*time(idx80)),[-5,4 ],'--k');




% xlim([0 0.4]);
% xlim([0 2]);

%%
tt = [0:1/fs:20];
ff = 3;
sig = sin(ff*tt);

figure;
plot(log(tt), log(1e-7+ (sig) ) );

% plot( (tt),  (1e-7+ (sin(1*tt)))  );

