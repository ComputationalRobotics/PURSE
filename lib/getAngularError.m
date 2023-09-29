
function rotError = getAngularError(R_gt,R_est)
% R_est can be a 3x3 matrix or a 3x3xN matrix
R_gt = squeeze(R_gt);
assert(isequal(size(R_gt),[3 3]));
assert(size(R_est, 1)==3 && size(R_est, 2)==3);

rotError = abs(acos((tensorprod(eye(3), tensorprod(R_gt', R_est, 2, 1), [1, 2], [1, 2])-1) / 2 ));
rotError = rad2deg( rotError );

end