Pp = 4; % [1] Nombre de paires de pôles
Rs = 0.18; % [Omega] Résistance stator vue en Park
Ld = 1.15e-3; % [H] Inductaice stator de Park directe
Lq = 3.13e-3; % [H] Inductaice stator de Park quadratique
1.3; % [T] Induction ŕemanente des aimants
200e-3; % [Wb] Flux des aimants vu par l'eisamble des spires d'une phase

RC1 = 50.0;
Wmbase = 350;
Emax = 800.0;
J = 800e-6;

Wv = 70; % [rad/s] Pulsation des tensions
E = 12; % [V] Tension continue d'alimentation
k = 0.3; % [?] 

D = 100e-6;
T = 1e-6;
N = D / T;

% Allocation de mémoire
n = 4;
x = zeros(n, N);
t = zeros(1, N);
[phid phiq wm a] = deal(0);
phim = 3; % WHAT IS THIS?

P = @(a) ...
	[cos(a) cos(a - 2/3*pi) cos(a + 2/3*pi);
	sin(a) sin(a - 2/3*pi) sin(a + 2/3*pi)];

for i=2:N,
	vo = E * sign([a a - 2/3*pi a - 4/3*pi]');

	vp = P(a) * vo;
	vd = vp(1);
	vq = vp(2);

	
	id = 1 / Ld * (phid - phim);
	iq = 1 / Lq * phiq;

	we = Pp * wm;
	Cr = k * we;
	C = Pp * (phid*iq - phiq*id);

	Da = we;
	Dwm = 1/J * C - Cr;
	Dphid = vd + we * phiq - Rs * id;
	Dphiq = vq - we * phid - Rs * iq;
	

	dx = [Da Dwm Dphid Dphiq]';
	x(:,i) = x(i-1) + T * dx;
	t(i) = t(i-1) + T;
end

%plot shit;
