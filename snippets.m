
% sphereHarm:
% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
% trapSphereR: integration of region R defined by tt and pp
% \int_R Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})

% convert -delay 50 -loop 0 slep*.png slep.gif
% convert -delay 30 -quality 100 aus_* movie.mov
% check=bsxfun(@rdivide,D*w,w) % check eigenvector works as expected

setenv('PATH', [getenv('PATH') ';D:\Perl\bin']);

setenv('PATH', [getenv('PATH') ':/usr/local/bin:/usr/local/bin']);