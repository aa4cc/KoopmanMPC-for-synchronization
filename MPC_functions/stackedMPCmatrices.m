function [Ab_bar, Ab, Bb_bar, Bb] = stackedMPCmatrices(A,B,N)
    n = size(A,1);
    m = size(B,2);
    
    Ab_bar = zeros((N+1)*n + N*m, n);
    Ab_bar(1:n,:) = eye(n,n);
    for i = 2:N+1
        Ab_bar((i-1)*n+1:i*n,:) = Ab_bar((i-2)*n+1:(i-1)*n,:)*A;
    end
    

    
    
    Bb_bar = zeros((N+1)*n + N*m, N*m);
    for i = 2:N+1
        Bb_bar((i-1)*n+1:i*n,:) = A * Bb_bar((i-2)*n+1:(i-1)*n,:);
        Bb_bar((i-1)*n+1:n*i,(i-2)*m+1:m*(i-1)) = B;
    end
    Bb_bar((N+1)*n +1:end, :) = eye(N*m);
    
%     %first block rows of Ab and Bb corresponding to x0 removed
    Ab = Ab_bar(n+1:end-N*m,:);
    Bb = Bb_bar(n+1:end-N*m,:);    
    
%     % Debug
%     Ab_bar = Ab_bar(n+1:end,:);
%     Bb_bar = Bb_bar(n+1:end,:); 
end