function problem = randPoseProblem(N)
% generate a random camera pose estimation problem
% N: number of keypoints

% define camera field of view (FOV)
FOV = deg2rad(150);
minDepth = 2;
maxDepth = 10;

% generate N 3D keypoints inside the camera field of view (FOV)
kpts3d = zeros(3,N);
kpts2d = zeros(2,N);
for i = 1:N
    alpha   = -FOV/2 + FOV * rand;
    beta    = -pi + 2*pi*rand;
    v       = [sin(alpha)*cos(beta); sin(alpha)*sin(beta); cos(alpha)];
    v       = v / v(3);
    depth   = minDepth + (maxDepth - minDepth) * rand;

    kpts3d(:,i)  = v * depth;
    kpts2d(:,i)  = v(1:2);
end

% generate groundtruth rotation and translation
R = quat2rotm(randrot);
t = randn(3,1);

% obtain 3D keypoints in world coordinate
kpts3d = R'*(kpts3d - t);

% generate 2D balls enclosing the groundtruth 2D keypoints
bound = 0.2;
centers = zeros(2,N);
radii   = zeros(1,N);
for i = 1:N
    [center, radius] = enclosingBall(kpts2d(:,i),bound);
    centers(:,i) = center;
    radii(i)     = radius;
end

problem.N       = N;
problem.kpts3d  = kpts3d;
problem.kpts2d  = kpts2d;
problem.R       = R;
problem.t       = t;
problem.centers = centers;
problem.radii   = radii;
end


function [center, radius] = enclosingBall(x,bound)
dim     = length(x);
radius  = bound*rand;
v       = randn(dim,1);
v       = v/norm(v);
center  = x + (radius * rand) * v;
% the distance between center and x is (radius * rand), therefore, a ball
% centered around "center" with radius "radius" must enclose x
end