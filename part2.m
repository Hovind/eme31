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

Wv = 220; % [rad/s] Pulsation des tensions
E = 50; % [V] Tension continue d'alimentation
k = 0.2; % [?]
Fd = 10e3;
Tp = 100e-6;

D = 20e-3;
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

Ma = 0;

P = @(a) ...
	sqrt(2/3)*[cos(a) cos(a - 2/3*pi) cos(a - 4/3*pi);
	-sin(a) -sin(a - 2/3*pi) -sin(a - 4/3*pi)];

Pinv = @(a) ...
    sqrt(2/3)*[cos(a) -sin(a);
    cos(a-2/3*pi) -sin(a-2/3*pi);
    cos(a-4/3*pi) -sin(a-4/3*pi)];

for i=2:N,
    % Dérivation
    
    % Inputs
    Mio = io;
    Map = Ma;
    Ma = a;
    
    % Controle
    Mip = P(Ma) * Mio;
    [Mid Miq] = deal(Mip(1), Mip(2));
    
    Mwe = (Ma - Map) / T;
    
    Mphid = Phim + Ld * Mid;
    Mphiq = Lq * Miq;
    
    Ephid = Rphid - Mphid;
    Ephiq = Rphiq - Mphiq;
    
    Rvd = SEphid * Gi + Ephid * Gp - Mphid * Mwe;
    Rvq = SEphiq * Gi + Ephiq * Gp + Mphiq * Mwe;
    Rvp = [Rvd Rvq]';
    
    Rvo = Pinv(Ma) * Rvp;
    
    theta = Wv * t(i) + pi;
	vo = E * sign(sin([theta theta - 2/3*pi theta - 4/3*pi]'));

    % Puissance
	vp = P(a) * vo;
	[vd vq] = deal(vp(1), vp(2));
	
	id = 1 / Ld * (phid - Phim);
	iq = 1 / Lq * phiq;
    ip = [id iq]';
    io = Pinv(a) * ip;

	we = Pp * wm;
	c = Pp * (phid*iq - phiq*id);
	cr = k * we;

    % Controle
    DSEphid = Ephid;
    DSEphiq = Ephiq;
    
    % Puissance
	Da = we;
	Dwm = 1/J * (c - cr);
	Dphid = vd + we * phiq - Rs * id;
	Dphiq = vq - we * phid - Rs * iq;
	
    % Intégration
    % Controle
    SEphid = SEphid + DSEphid * T;
    SEphiq = SEphiq + DSEphiq * T;
    
    %Puissance
    a = a + Da * T;
    wm = wm + Dwm * T;
    phid = phid + Dphid * T;
    phiq = phiq + Dphiq * T;
    
    % Sauveguarde résultat pour les courbes
    phi(:,i) = [phid phiq]';
    vout(:,i) = [vo(1) + 5 vo(2) vo(3) - 5]';
    vpark(:,i) = vp;
    iout(:,i) = io;
    ipark(:,i) = ip;
    y(:,i) = [c we a]';
end

%plot shit
figure(1);
plot(t, vout);
title('Vo');
legend('Vo_1', 'Vo_2', 'Vo_3');
saveas(gcf,'Barchart', '.pdf')

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