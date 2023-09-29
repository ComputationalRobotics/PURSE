classdef PURSE
    properties
        quadCon
        linCon
        centers
        radii
        kpts3d
        intrinsics
    end

    methods
        %% Constructor
        function obj = PURSE(mu,r,Y,P)
            % construct a purse from 2D balls
            % mu: 2 x N matrix, each column is a 2D point
            % r: 1 x N matrix, each column is the radius of a ball
            % Y: 3 x N matrix, each column is a 3D keypoint
            % P: 3 x 3 camera intrinsics, default P is identity
            if nargin < 4
                P = eye(3);
            end
            obj.intrinsics       = P;
            obj.centers          = mu;
            obj.radii            = r;
            obj.kpts3d           = Y;
            N                    = size(Y,2);
            A                    = zeros(12,12,N);
            b                    = zeros(12,N);
            for i = 1:N
                radius      = r(i);
                center      = mu(:,i);
                Lambda      = (1 / radius)^2 * eye(2);
                Yi          = Y(:,i);
                U           = [kron(Yi',P), P];
                u1          = U(1,:)';
                u2          = U(2,:)';
                u3          = U(3,:)';
                b(:,i)      = u3;
                Umu         = [u1-center(1)*u3, u2-center(2)*u3];
                A(:,:,i)    = Umu*Lambda*Umu' - u3*u3';
            end
            obj.quadCon = A;
            obj.linCon  = b;
        end

        %% normalize 2D keypoints (remove camera intrinsics)
        function y_norm = calibrate_kpts2d(obj,y)
            % y of shape (2,N)
            n_kpts = size(y,2);
            y      = [y;ones(1,n_kpts)];
            y_norm = obj.intrinsics \ y;
            y_norm = y_norm(1:2,:);
            % return y_norm also in shape (2,N)
        end
        
        %% Check membership
        function isvalid = checkMembership(obj,R,t)
            N = size(obj.centers,2);
            s = [R(:);t];
            isvalid = true;
            for i = 1:N
                bi = obj.linCon(:,i);
                Ai = squeeze(obj.quadCon(:,:,i));
                if bi'*s <= 0
                    isvalid = false;
                    break
                end
                if s'*Ai*s > 0
                    isvalid = false;
                    break
                end
            end
        end

        %% RANSAG
        function [y,Y] = sample3pairs(obj)
            N   = size(obj.centers,2);
            ids = randsample(N, 3);
            Y   = obj.kpts3d(:,ids);
            mu  = obj.centers(:,ids);
            r   = obj.radii(:,ids);
            v   = normc(randn(2,3));
            y   = mu + (r .* rand(1,3)) .* v;
        end

        function [Rset,tset,Ravg,tavg] = ransag(obj,T)
            % T: maximum number of iterations
            Rset = {};
            tset = {};
            Ravg = eye(3);
            tavg = zeros(3,1);
            Rsum = zeros(3,3);
            tsum = zeros(3,1);

            for i = 1:T
                [y,Y] = obj.sample3pairs();
                % remove camera intrinsics from 2D keypoints
                y     = obj.calibrate_kpts2d(y);
                [R,t] = solveP3P(y,Y);
                if ~isempty(R)
                    for j = 1:size(t,2)
                        inside = obj.checkMembership(R(:,:,j),t(:,j));
                        if inside
                            Rsum = Rsum + R(:,:,j);
                            tsum = tsum + t(:,j);
                            Rset = [Rset;{R(:,:,j)}];
                            tset = [tset;{t(:,j)}];
                        end
                    end
                end
            end

            n = length(Rset);
            
            if any(isnan(Rsum)) | any(isinf(Rsum))
                idx1 = find(isnan(Rsum));
                idx2 = find(isinf(Rsum));

                Rsum(idx1) = 0;
                Rsum(idx2) = 0;
            end
            
            if n > 0
                tavg = tsum / n;
                Ravg = project2SO3(Rsum);
            end
        end

    end


end

