function get_wind,wind,mwind_input,AU,msun,yr
;
    ;common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
  common grid,rr,rrr,dr,dr1,dr2
;
  nx=n_elements(rr)
;
    case wind of 
       'central': begin
          mwind=mwind_input*(msun/yr)
          rmax=30*AU;rr[nx-1]         ; 30*AU
          rg=5*AU               ; this can be better calculated GM/c^2, or similar
          ;den=2*!dpi*(rmax-rg)*rr
          den=2*!dpi*(rr/rg)^2.5 * AU^2 
          il=where(rr le rg,nl) & ig=where(rr gt rg,ng)
          swind=dblarr(nx)
          if (nl ne 0) then swind[il]=0.
          if (ng ne 0) then swind[ig]=mwind/den[ig]
       end
;       
       'owen': begin
          a=0.15138 & b=-1.2182  & c=3.4046
          d=-3.5717 & e=-0.32762 & f=3.6064
          g=-2.4918
;
          Mass=0.7
          Lx=1.
;
          r_aux=0.85*(Mass/1)^(-1)*rr/AU ; rr in AU, mass in solar masses
          rlog10 = alog10(r_aux)
          Mdot_wind=dblarr(nx)
;
          for i=1,nx-1 do begin
             if (r_aux[i] le 100 and r_aux[i] gt 0.7) then begin
                Mdot_wind[i]=10.^(a*rlog10[i]^6+$
                                  b*rlog10[i]^5+$
                                  c*rlog10[i]^4+$
                                  d*rlog10[i]^3+$
                                  e*rlog10[i]^2+$
                                  f*rlog10[i]  +$
                                  g                 )
             endif
          endfor
;
          v=alog(10.)
;
          wind=(Mdot_wind/(2*!pi)*((6*a*alog(r_aux)^5)/(r_aux^2*v^7) + $
                              5*b*alog(r_aux)^4./(r_aux^2*v^6) + 4*c*alog(r_aux)^3./(r_aux^2*v^5) + $
                              3*d*alog(r_aux)^2./(r_aux^2*v^4)+2*e*alog(r_aux)/(r_aux^2*v^3) + $
                              f/(r_aux^2*v^2)))*exp(-(r_aux/100)^10)
;
          a=2*!pi*r_aux*wind
          ddr=r_aux[1:nx-1] - r_aux[0:nx-2]
          Mdot_norm=total(a[1:nx-1]*ddr)

          mdot_msunyr = 6.25d-9 * Mass^(-0.068) * LX^(1.14) ; for Mass in solar units 
                                                            ; and Lx in units of 10^30 erg/s
          swind = wind / mdot_norm * mdot_msunyr
;
;  add direct wind
;
;          if (ldirect_wind) then begin
;             Rin=10. ; AU                                                                           
;
;             x_in=R-Rin
;
;a = -4.3822599312953403E-01
;b = -1.0658387115013401E-01
;c = 5.6994636363474982E-01
;d = 1.0732277075017336E-02
;e = -1.3180959703632333E-01
;f = -1.3228570869396541E+00
;
;y = (a * exp(b*x_in) + c * exp(d*x_in) + e * exp(f*x_in))
;
;i=where(x_in lt 0) & y[i]=0.
;
;store=deriv(R,y)
;
;sigmaw=store/R*exp(-((R-Rin)/57.)^10)
;
;tmp=R*sigmaw & dr=r[1:nx-1]-r[0:nx-2]
;tmp2=total(tmp[1:nx-1]*dr)
;sigmaw_normalized=sigmaw/tmp2
;
;mdot_msunyr = 4.8d-9 * Mass^(-0.148) * LX^(1.14) ; for Mass in solar
;units and Lx in units of 10^30 erg/s                                                                  ;  
;
;!p.charsize=2
;
;swind_direct = sigmaw_normalized * mdot_msunyr
;          endif
;          plot,rr/AU,swind
;          stop
       end
    endcase
;
    return,swind
;
end
