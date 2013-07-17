
Mass=0.7
Lx=1.

a=0.15138
b=-1.2182
c=3.4046
d=-3.5717
e=-0.32762
f=3.6064
g=-2.4918



nx=500
lnRau=grange(-1,2,nx) & Rau=10^lnRau
Mdot=replicate(0.,nx) + 1e-20
R=0.85*(Mass/1)^(-1)*Rau

rlog10 = alog10(R)

for i=1,nx-1 do begin
    if (R[i] le 100 and R[i] gt 0.7) then begin
        Mdot[i]=10.^(a*rlog10[i]^6+$
                     b*rlog10[i]^5+$
                     c*rlog10[i]^4+$
                     d*rlog10[i]^3+$
                     e*rlog10[i]^2+$
                     f*rlog10[i]  +$
                     g                 )
     endif
end

;plot,r,mdot,/xlog,/ylog,yr=[1e-5,1e1],xr=[1e-1,100],xs=3,ys=3

v=alog(10.)

wind=(Mdot/(2*!pi)*((6*a*alog(R)^5)/(R^2*v^7) + $
5*b*alog(R)^4./(R^2*v^6) + 4*c*alog(R)^3./(R^2*v^5) + $
3*d*alog(R)^2./(R^2*v^4)+2*e*alog(R)/(R^2*v^3) + $
f/(R^2*v^2)))*exp(-(R/100)^10)


;plot,r,wind ;,/xlog,/ylog,yr=[1e-8,1e-2],xr=[1e-1,100],xs=3,ys=3

a=2*!pi*R*wind 
ddr=r[1:nx-1] - r[0:nx-2]
Mdot=total(a[1:nx-1]*ddr)

mdot_msunyr = 6.25d-9 * Mass^(-0.068) * LX^(1.14) ; for Mass in solar units and Lx in units of 10^30 erg/s

!p.charsize=2

swind = wind / mdot * mdot_msunyr

plot,rau,swind

end
