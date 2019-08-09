variable rp, n, d, r, s, a, e, i, x1, x2, x3, x4, x5;
variable rmax, smax, amax, emax, dmax, k, lo, hi, h;

rp = [500.0, 750.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0];
n = length(rp);
d = Array_Type[n];
r = @d;
s = @d;
a = @d;
e = @d;
rmax = Double_Type[n];
smax = @rmax;
amax = @rmax;
emax = @rmax;
dmax = @rmax;

for (i = 0; i < n; i++) {
  (x1, x2, x3, x4, x5) = readcol(sprintf("torg%d.dmp", i+1), 1, 2, 4, 5, 6);
  d[i] = x1;
  r[i] = x3;
  s[i] = x2;
  a[i] = x5;
  e[i] = x4;
  k = where(e[i] == max(e[i]));
  k = k[0];
  dmax[i] = d[i][k];
  emax[i] = e[i][k];
  rmax[i] = r[i][k];
  smax[i] = s[i][k];
  amax[i] = a[i][k];
}

k = open_plot("eff-length.ps/ps");
yrange(3e-7, 3e-5);
xrange(10.0,500);
ylog;
xlabel("Length (cm)");
ylabel("Efficiency");
for (i = 0; i < n; i++) {
  color(1);
  if (i == 0) plot(d[i], e[i]);
  else oplot(d[i], e[i]);
  color(1);
  xylabel(d[i][0], log10(e[i][0]), sprintf("%d",int(rp[i])), 90, 1);
}
close_plot(k);

k = open_plot("max-eff.ps/ps");
ylin;
xlabel("Resolving Power");
ylabel("Max. Efficiency");
xrange(400,2100);
yrange(1e-6,3e-5);
plot(rp, emax);
close_plot(k);

k = open_plot("radius.ps/ps");
ylabel("Radius of Curvature (m)");
xlabel("Resolving Power");
xrange(400,2100);
yrange(15.0, 60.0);
plot(rp, rmax*1e-2);
close_plot(k);

k = open_plot("length.ps/ps");
ylabel("Total Length (cm)");
xlabel("Resolving Power");
xrange(400,2100);
yrange(50.0, 350);
plot(rp, dmax);
close_plot(k);

k = open_plot("image.ps/ps");
(x1,x2) = readcol("image_torg.15",1,2);
xlabel("X ("+latex2pg("\\mu")+"m)");
ylabel("Y (cm)");
xrange(-75,75);
yrange(-1.0, 1.0);
(lo,hi) = linear_grid(-150,150,256);
h = histogram(x1, lo, hi);
pointstyle(-6);
connect_points(0);
plot(x1, x2);
close_plot(k);

k = open_plot("spec.ps/ps");
ylabel("Counts");
xlabel("X ("+latex2pg("\\mu")+"m)");
xrange(-75,75);
yrange(0,max(h)*1.2);
pointstyle(-1);
connect_points(1);
hplot(lo, hi, h);
close_plot(k);

k = open_plot("eff-lambda.ps/ps");
connect_points(1);
pointstyle(-1);
(x1,x2) = readcol("grasp_torg.r",1,2);
ylabel("Efficiency");
xlabel("Wavelength ("+latex2pg("\\A")+")");
xrange(5,50);
yrange(1e-6,3e-6);
plot(x1, x2);
close_plot(k);
