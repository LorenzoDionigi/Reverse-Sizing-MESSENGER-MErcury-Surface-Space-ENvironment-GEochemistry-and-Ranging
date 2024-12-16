function [r_eci,v_eci] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
p=a.*(1-e.^2);
r=p./(1+e.*cos(th));

r_pf=r*[cos(th); sin(th); 0];
v_pf=sqrt(mu./p)*[-sin(th); e+cos(th); 0];

R_OM=[cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_om=[cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

T=[R_om*R_i*R_OM]';

r_eci=T*r_pf;
v_eci=T*v_pf;

