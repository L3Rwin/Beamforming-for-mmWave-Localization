function a=getResponse(N,phi)
    ii=0:N-1;
    a=exp(-1j*pi*ii'*phi)/sqrt(N);