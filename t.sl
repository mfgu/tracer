variable R, a, sa, ca, lam, d, cb, b, sb, b2, ds, dd, L, S, rp, di, dw;

R = 44.3E2;
a = [1.5:4.0:0.01]*PI/180;
sa = sin(a);
ca = cos(a);
lam = 15.0E-8;
d = 0.1/2400.0;
cb = ca - lam/d;
b = acos(cb);
sb = sin(b);
b2 = -4.5e-3*R;
ds = 55E-4;
dd = 25E-4;

L = 1.0/(1.0/(2*R*sb) - b2/(R*cb));
S = (sa+sb-(ca-cb)*cb/sb)/(R*sa*sa) - (sa*sa+(ca-cb)*(ca-cb))/(L*sa*sa);
S = 1.0/S;
a *= 180.0/PI;
b *= 180.0/PI;
rp = (ca-cb)/sqrt(sa*sa*ds*ds/(S*S) + sb*sb*dd*dd/(L*L));

di = ds*(sa/sb)*(L/S);

dw = d*sb/L;
