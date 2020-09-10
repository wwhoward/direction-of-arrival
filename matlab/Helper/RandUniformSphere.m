function [azi,ele] = RandUniformSphere(parSize)

% Description:

% Creates a vector of random angles for a sphere that are uniformly distributed along the sphere (not concentarated at the
% poles)

% Sources:
% https://www.mathworks.com/matlabcentral/answers/289416-x-y-z-sphere-except-uniform-or-random-distributed-points
% http://corysimon.github.io/articles/uniformdistn-on-sphere/

% Parameters:

% parSize: [1x2] vector containing [Nsig,Npass];

u = randn(3,prod(parSize));
u = u./sqrt(sum(u.^2,1));
[azi,ele] = cart2sph(u(1,:),u(2,:),u(3,:));
azi = reshape(azi,[parSize(1),parSize(2)]);
ele = reshape(ele,[parSize(1),parSize(2)]);
ele = ele + pi/2;

end

