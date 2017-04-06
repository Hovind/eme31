Pp = 4; % [1] Nombre de paires de pôles
Rs = 0.18; % [Omega] Résistance stator vue en Park
Ld = 1.15e-3; % [H] Inductaice stator de Park directe
Lq = 3.31e-3; % [H] Inductaice stator de Park quadratique
1.3; % [T] Induction ŕemanente des aimants
Phim = 200e-3; % [Wb] Flux des aimants vu par l'eisamble des spires d'une phase

RC1 = 50.0;
Wmbase = 350;
Emax = 800.0;
J = 800e-6;

Wv = 70; % [rad/s] Pulsation des tensions
E = 12; % [V] Tension continue d'alimentation
k = 0.3; % [?] 

D = 100e-3;
T = 1e-5;
N = round(D / T);

% Allocation de mémoire
n = 4;
phi = zeros(2, N);
vout = zeros(3, N);
vpark = zeros(2, N);
iout = zeros(3, N);
ipark = zeros(2, N);
y = zeros(3, N);
t = 0:T:D-T;
[a wm phid phiq] = deal(0, 0, Phim, 0);

P = @(a) ...
	sqrt(2/3)*[cos(a) cos(a - 2/3*pi) cos(a - 4/3*pi);
	-sin(a) -sin(a - 2/3*pi) -sin(a - 4/3*pi)];

Pinv = @(a) ...
    sqrt(2/3)*[cos(a) -sin(a);
    cos(a-2/3*pi) -sin(a-2/3*pi);
    cos(a-4/3*pi) -sin(a-4/3*pi)];

for i=2:N,
    theta = Wv * t(i) + pi;
	vo = E * sign(sin([theta theta - 2/3*pi theta - 4/3*pi]'));

	vp = P(a) * vo;
	[vd vq] = deal(vp(1), vp(2));
	
	id = 1 / Ld * (phid - Phim);
	iq = 1 / Lq * phiq;
    ip = [id iq]';
    io = Pinv(a) * ip;

	we = Pp * wm;
	c = Pp * (phid*iq - phiq*id);
	cr = k * we;

	Da = we;
	Dwm = 1/J * (c - cr);
	Dphid = vd + we * phiq - Rs * id;
	Dphiq = vq - we * phid - Rs * iq;
	
    a = a + Da * T;
    wm = wm + Dwm * T;
    phid = phid + Dphid * T;
    phiq = phiq + Dphiq * T;
    
    phi(:,i) = [phid phiq]';
    vout(:,i) = [vo(1) + 5 vo(2) vo(3) - 5]';
    vpark(:,i) = vp;
    iout(:,i) = io;
    ipark(:,i) = ip;
    y(:,i) = [c we a]';
end

%plot shit;
figure(1);
plot(t, vout);
legend('Vo_1', 'Vo_2', 'Vo_3');

figure(2);
plot(t, vpark);
legend('V_d', 'V_q');

figure(3);
plot(t, phi);
legend('phi_d', 'phi_q');

figure(4);
plot(t, ipark);
legend('I_d', 'I_q');

figure(5);
plot(t, iout);
legend('I_1', 'I_2', 'I_3');

figure(6);
plot(t, y);
legend('C', 'We', 'A');