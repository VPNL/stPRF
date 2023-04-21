function inNum = nearZerotoZero(inNum,tol)

if nargin == 1
    tol = 1e-2;
end
inNum( abs(inNum) < tol) = 0;

end