%{
function zmpp = calculate_zmpp(r, c, w, t, d, v, p, yv)
    zmpp = -((w / (c * p) * (t + 8 * d * v / (3 * pi())) / yv) * 1i) * r * c;
end

function yv = calculate_yv(j,d,kv)
    yv = (besselj(2,kv*d/2)/besselj(0,kv*d/2));
end

function kv = calculate_kv(w, p, u)
    kv = sqrt(-(w * p / u) * 1i);
end

function result = calculate_v(p)
    % Define the coefficients
    a = [1, -1.4092, 0, 0.33818, 0, 0.06793, -0.02287, 0.063015, -0.01614];
    
    % Calculate the sum
    result = 0;
    for n = 0:8
        result = result + a(n + 1) * (sqrt(p))^n;
    end
end

function zac = calculate_zac(r,c,w,D)
    zac = -(r*c*1/tan(w*D/c))*1i;
end

function [za,zb,zc] = calculate_dyTransform(z1,z2,z3)
    za = (z2*z3)/(z1+z2+z3);
    zb = (z2*z1)/(z1+z2+z3);
    zc = (z3*z1)/(z1+z2+z3);
end

function zTotal = calculate_zTotal(za,zb,zc,z1,z2)
    zTotal = (z1+zb)*(z2+zc)/(z1+z2+zb+zc) + za;
end

%% Initial Values

f = 200:2000;
wf = 2*f*pi();
absCoef = f*0;

% Air properties
dd = 0.25; % air cavity depth
r0 = 1.204; % air density (kg m^-3)
c0 = 343.1; % sound speed (m/s)
u0 = 1.825*10^(-5); % dynamic viscosity (kg (m s)^-1)

% Initial values for MPP1, MPP2, MPP3
t = [0.002 0.002 0.001] ; % panel thickness (m)
d = [0.001 0.002 0.0001] ; % perforation diameter (m)
p = [0.01 0.04 0.01]; % perforation porosity (%)
n = 1; % perforation number

%% Input (Optional)

%{
for i = 1:3
    prompt = {'Enter value for t:', 'Enter value for d:', 'Enter value for p:'};
    dlgtitle = ['Input MPP-', num2str(i)];
    dims = [1 50];
    definput = {num2str(t(i)), num2str(d(i)), num2str(p(i))};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    t(i) = str2double(answer{1});
    d(i) = str2double(answer{2});
    p(i) = str2double(answer{3});
    fprintf('For MPP-%d:\n', i);
    fprintf('t: %f\n', t(i));
    fprintf('d: %f\n', d(i));
    fprintf('p: %f\n\n', p(i));
end
%}

%% Absorption Coefficient Prediction Across Frequency Range
% 2 Primary Mpps, 1 Secondary Mpp, 2 Air cavities
% f = 0.2-2 kHz

z = zeros(1,5);
a = 0;

for i = wf
    % mpp impedance
    for j = 1:3
        v = calculate_v(p(j));
        kv = calculate_kv(i,r0,u0);
        yv = calculate_yv(d(j),kv);
        z(j) = calculate_zmpp(r0,c0,i,t(j),d(j),p(j),v,yv);
    end
    
    % air cavity impedance
    for j = 1:2
        z(j+3) = calculate_zac(r0,c0,i,dd);
    end
    
    % delta-wye transform
    [z(3),z(4),z(5)] = calculate_dyTransform(z(3),z(4),z(5));

    % total impedance
    zTotal = calculate_zTotal(z(1),z(2),z(4),z(5),z(3));

    % reflection coefficient
    rfl = (zTotal-r0*c0)/(zTotal+r0*c0);

    % absorption coefficient
    a = a + 1;
    absCoef(a) = 1-abs(rfl)^2;
end
disp(a);
% Plot
plot (f,absCoef);
grid on;
%}