%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Code to find out the position of Satellites using Orbital Parameters for time duration of 24 hours
% Satellite Constellation
% Ref: https://kr.mathworks.com/matlabcentral/fileexchange/65845-satellite-constellation
% Input: 6 orbital elements of satellites
% Output: 24 hours satellites constellation
close all; clc;
%% Parameters
input_data = readmatrix('Starlink constellation.csv');

RE = 6378;                          % Earth's radius                        [km]
muE = 3.986004418e+14;              % Earth gravitational parameter         [m^3/sec^2]

% 6 orbital elements
eccen=input_data(:,1).';            % ECCENTRICITY
RAANangle=input_data(:,2).';        % RAAN ANGLE
sa = input_data(:,3).';             % SEMI-MAJAR AXIS PARAMETER
arg=input_data(:,4).';              % ARGUMENT OF PERIGEE
Meananom=input_data(:,5).';         % MEAN ANOMALY
orbinc=input_data(:,6).';           % INCLINATION

numofsat=1:1:1584;                  % # of SATELLITES
%% Satellites' Constellation
parfor m=1:length(numofsat)
    Mo=Meananom(1,m);                    % Mean anomaly of satellite
    RAAN  = RAANangle(1,m);              % RAAN                          [rad]
    w     = arg(1,m);                    % Argument of perigee           [rad]
    inc   = orbinc(1,m);                 % inclination                   [rad]
    ecc   = eccen(1,m);                  % eccentricity
    a     = sa(1,m);                     % semi-major axis               [m]
    n     = sqrt(muE/a^3);               % Mean Motion                   [rad/s]
    t=24*60*60;                          % time span (24 hrs)
    J=1:600:t;                           % 10 minutes
    xs=zeros(1,72);                      
    ys=zeros(1,72);                      
    zs=zeros(1,72);                      
        for k=1:length(J)
        Mi=Mo+n*J(k);
        Tol = 0.001;
        E1=Mi;
        E2=E1-((E1-ecc*(sin(E1))-Mi)/(1-ecc*(cos(E1)))); % NEWTON-RAPHSON APPROACH
        error=100;
                while error > Tol
                E1=E2;
                E2=E1-((E1-ecc*(sin(E1))-Mi)/(1-ecc*(cos(E1))));
                error=abs(E2-E1);
                end
         r=a*(1-ecc*cos(E2)); % RADIUS
         b=2*atan(sqrt((1-ecc)/(1+ecc))*tan(E2/2)); % TRUE ANOMALY
         u = w+b;
        X=r*cos(u); % IN-PLANE
        Y=r*sin(u); % IN-PLANE
        xs(:,k) = X*cos(RAAN)-Y*cos(inc)*sin(RAAN);    % ECEF x-coordinate SAT         [m]
        ys(:,k) = X*sin(RAAN)+Y*cos(inc)*cos(RAAN);    % ECEF y-coordinate SAT         [m]
        zs(:,k) = Y*sin(inc);                          % ECEF z-coordinate SAT         [m]
        end
Xs(m,:)=[xs];
Ys(m,:)=[ys];
Zs(m,:)=[zs];
x_points(m,1) = xs(1,1);
y_points(m,1) = ys(1,1);
z_points(m,1) = zs(1,1);
end
plot3(Xs,Ys,Zs,'.b', MarkerSize=1) % SATELLITE CONSTELLATION
hold on
scatter3(x_points,y_points,z_points,'ob', 'filled') % SATELLITES
%% Plot Earth6
% earth image file
image_file = 'https://images.alphacoders.com/590/590923.jpg';
set(gca, 'NextPlot','add', 'Visible','off'); % Turn off the normal axes

[X,Y,Z] = sphere;
X = RE*X; Y = RE*Y; Z = RE*Z;

GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0
globe = surf(X, Y, -Z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
cdata = imread(image_file);
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha);