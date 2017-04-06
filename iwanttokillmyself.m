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

D = 20e-3;
Tc = 100e-6;
Tr = 1e-6;
N = round(D / Tc);
M = round(Tc / Tr);

% Allocation de mémoire
n = 4;
phi = zeros(2, N);
rvout = zeros(3, N);
vout = zeros(3, N);
vpark = zeros(2, N);
iout = zeros(3, N);
ipark = zeros(2, N);
y = zeros(3, N);
tout = 0:Tr:D-Tr;
[a wm phid phiq] = deal(0, 0, Phim, 0);

P = @(a) ...
	sqrt(2/3)*[cos(a) cos(a - 2/3*pi) cos(a - 4/3*pi);
	-sin(a) -sin(a - 2/3*pi) -sin(a - 4/3*pi)];

Pinv = @(a) ...
    sqrt(2/3)*[cos(a) -sin(a);
    cos(a-2/3*pi) -sin(a-2/3*pi);
    cos(a-4/3*pi) -sin(a-4/3*pi)];

PWM = @(a, t, tr, tf) ...
    -E+2*E*((tr < t) .* (t < tf)); 

for i=1:N*M,
    t = mod(tout(i), Tc); 
    if t == 0
        % Controle
        theta = Wv * tout(i) + pi;

        Rvo = E * sin([theta theta - 2/3*pi theta - 4/3*pi]');
        tr = Tc * (1 - Rvo / E) / 4;
        tf = Tc * (3 + Rvo / E) / 4;
    end
    
    
    vo = PWM(Rvo, t, tr, tf);
    % Tension de Park
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
	
    a = a + Da * Tr;
    wm = wm + Dwm * Tr;
    phid = phid + Dphid * Tr;
    phiq = phiq + Dphiq * Tr;
    
    phi(:,i) = [phid phiq]';
    rvout(:,i) = Rvo;
    vout(:,i) = [vo(1) + 5 vo(2) vo(3) - 5]';
    vpark(:,i) = vp;
    iout(:,i) = io;
    ipark(:,i) = ip;
    y(:,i) = [c we a]';
end

%plot shit;
figure(1);
plot(tout, vout);
legend('Vo_1', 'Vo_2', 'Vo_3');

figure(2);
plot(tout, vpark);
legend('V_d', 'V_q');

figure(3);
plot(tout, phi);
legend('phi_d', 'phi_q');

figure(4);
plot(tout, ipark);
legend('I_d', 'I_q');

figure(5);
plot(tout, iout);
legend('I_1', 'I_2', 'I_3');

figure(6);
plot(tout, y);
legend('C', 'We', 'A');


figure(1);
plot(tout, rvout);
legend('RVo_1', 'RVo_2', 'RVo_3');
