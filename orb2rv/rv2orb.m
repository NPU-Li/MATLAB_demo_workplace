%----------------------- Begin Code Sequence -----------------------------%
% Purpose:                                                                %
% Convert a given set of state vectors in ECI reference frame to orbital  %
% elements.                                                               %
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
% Outputs:
%---------                                                                %
%a                      [1 x N]                         Semi-Major Axis
%                                                       (km)
%
%eMag                   [1 x N]                         Eccentricity
%                                                       (unitless)
%
%i                      [1 x N]                         inclination
%                                                       (radians)
%
%O                      [1 x N]                         Right Ascention of
%                                                       the ascending node
%                                                       (radians)
%
%o                      [1 x N]                         Argument of perigee
%                                                       (radians)
%
%M                      [1 x N]                         Mean Anomaly
%                                                       (radians)
%
%truLon                 [1 x N]                         True Longitude
%                                                       (radians)
%
%argLat                 [1 x N]                         Argument of Latitude
%                                                       (radians)
%
%lonPer                 [1 x N]                         Longitude of Periapse
%                                                       (radians)
%
%p                      [1 x N]                         Semilatus Rectum
%                                                       (km)
%%----------------------- 开始代码序列 -----------------------------%
% 目的：                                                               %
% 将给定的一组ECI参考系中的状态向量转换为轨道根数。                    %
%-------------------------------------------------------------------------%
%                                                                         %
% 输入：                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         在ECI坐标系中的
%                                                       位置矢量
%
%v_ECI                  [3 x N]                         在ECI坐标系中的
%                                                       速度矢量
%
%mu                     double                          引力常数
%                                                       如未指定，默认为地球引力常数
%
% 输出：                                                                 %
%---------                                                                %
%a                      [1 x N]                         半长轴 (km)
%
%eMag                   [1 x N]                         偏心率 (无量纲)
%
%i                      [1 x N]                         轨道倾角 (弧度)
%
%O                      [1 x N]                         升交点赤经 (弧度)
%
%o                      [1 x N]                         近地点幅角 (弧度)
%
%M                      [1 x N]                         平近点角 (弧度)
%
%truLon                 [1 x N]                         真经度 (弧度)
%
%argLat                 [1 x N]                         纬度幅角 (弧度)
%
%lonPer                 [1 x N]                         近地点经度 (弧度)
%
%p                      [1 x N]                         半通径 (km)
%
% References:
%-------------
%Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
%
% Function Dependencies:
%------------------
%None
%------------------------------------------------------------------       %
% Programed by Darin Koblick  03-04-2012                                  %
% Updated to address circular equatorial orbits       12/12/2013          %
%------------------------------------------------------------------       %
function [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(r,v,mu)
if ~exist('mu','var');  t = getConst(); mu = t.Earth.Mu; end
%Specific angular momentum
h = cross(r,v);
n = cross(repmat([0;0;1],[1,size(r,2)]),h); nMag = sqrt(sum(n.^2,1));
vMag = sqrt(sum(v.^2,1)); 
rMag = sqrt(sum(r.^2,1)); 
hMag = sqrt(sum(h.^2,1));
e = (1./mu).*(bsxfun(@times,(vMag.^2 - mu./rMag),r) - bsxfun(@times,dot(r,v),v)); 
eMag = sqrt(sum(e.^2,1));
zeta = (vMag.^2)./2 - mu./rMag;
%Special Procedure when we have a parabolic orbit
idx = eMag ~= 1;
a = NaN(size(eMag));
p = NaN(size(eMag));
if any(idx)
    a(idx) = -mu./(2.*zeta(idx)); 
    p = a(idx).*(1-eMag(idx).^2); 
else
    a(idx) = Inf; 
    p(idx) = (hMag(idx).^2)./mu; 
end
%Compute the angles
i = acos(h(3,:)./hMag); 
O = acos(n(1,:)./nMag);
o = acos(dot(n,e)./(nMag.*eMag));
nu = acos(dot(e,r)./(eMag.*rMag));
lonPer = acos(e(1,:)./eMag);
argLat = acos(dot(n,r)./(nMag.*rMag));
truLon = acos(r(1,:)./rMag);
%Account for those cases where satellite is in circular orbit
         O(n(1,:) == 0) = 0;
       o(dot(n,e) == 0) = 0;
    lonPer(e(1,:) == 0) = 0;
      nu(dot(e,r) == 0) = 0;
  argLat(dot(n,r) == 0) = 0;
%Apply Quadrant Checks to All Determined Angles
idx = n(2,:) < 0; if any(idx);  O(idx) = 2*pi - O(idx);  end
idx = e(3,:) < 0; if any(idx); o(idx) = 2*pi - o(idx); end
idx = dot(r,v) < 0; if any(idx); nu(idx) = 2*pi - nu(idx); end
idx = e(2,:) < 0; if any(idx); lonPer(idx) = 2*pi-lonPer(idx);  end
idx = r(3,:) < 0; if any(idx); argLat(idx) = 2*pi - argLat(idx); end
idx = r(2,:) < 0; if any(idx); truLon(idx) = 2*pi - truLon(idx); end
end
