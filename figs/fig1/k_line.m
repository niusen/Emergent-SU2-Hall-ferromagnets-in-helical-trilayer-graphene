function kline=k_line(k1,k2,L)
kx=linspace(k1(1),k2(1),L);
ky=linspace(k1(2),k2(2),L);

kline=[kx;ky]';
end