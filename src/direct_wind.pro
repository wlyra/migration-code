
Mass=0.7
Lx=1.

nx=500
lnr=grange(-1,2,nx) & r=10^lnr

Rin=10. ; AU

x_in=R-Rin

a = -4.3822599312953403E-01
b = -1.0658387115013401E-01
c = 5.6994636363474982E-01
d = 1.0732277075017336E-02
e = -1.3180959703632333E-01
f = -1.3228570869396541E+00

y = (a * exp(b*x_in) + c * exp(d*x_in) + e * exp(f*x_in))

i=where(x_in lt 0) & y[i]=0.

store=deriv(R,y)

sigmaw=store/R*exp(-((R-Rin)/57.)^10)

tmp=R*sigmaw & dr=r[1:nx-1]-r[0:nx-2]
tmp2=total(tmp[1:nx-1]*dr)
sigmaw_normalized=sigmaw/tmp2

mdot_msunyr = 4.8d-9 * Mass^(-0.148) * LX^(1.14) ; for Mass in solar units and Lx in units of 10^30 erg/s

!p.charsize=2

swind_direct = sigmaw_normalized * mdot_msunyr

plot,r,swind_direct

end
