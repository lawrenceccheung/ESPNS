function calcerror(root, kmin, kmax)
%
%

kvec = kmin:kmax;

% These are the exact values for k = 1 to +3
exactu = [1, 0.25, 1.0/18.0];

% make the exact root list
exactroot=zeros(1,abs(kmin)+1);
exactroot=[exactroot, exactu(1:abs(kmax))];
%disp(exactroot);

N=abs(kmin)+abs(kmax)+1;
iroot=1;

fprintf('%2s %35s   %12s\n', 'k',  'Computed u','Abs(err)')
fprintf('%2s %35s   %12s\n', '--', '----------','--------')
for i=1:N,
    err = abs(root(1,i)-exactroot(i));
    if abs(root(1,i)) < 1.0e-6, 
        fprintf('%+i %18.8e %+16.8ei    %12.6e\n', kvec(i), real(root(1,i)), imag(root(1,i)), err);
    else
        fprintf('%+i %18.13f %+16.13fi    %12.6e\n', kvec(i), real(root(1,i)), imag(root(1,i)), err);
    end
end


