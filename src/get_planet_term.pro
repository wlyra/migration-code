function get_planet_term,sigma,tmid,cp,gamma,q,sqrtgm,ap,omega
;
    ;common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
  common grid,rr,rrr,dr,dr1,dr2
;
  nx=n_elements(rr)
  rr1=1./rr
;
  cs2=tmid*cp*(gamma-1)
  H=sqrt(cs2*gamma)/omega
;
  coef=sign(ap-rr) * q^2 * sqrtgm * rr1
  n_smooth=10.
  num = rr < ap 
;
; The original term was max(abs(rr-ap),H)
;
  tmp = abs(rr-ap)
  den = (tmp^n_smooth + H^n_smooth)^(1./n_smooth)
;
  arg = sqrt(rr)*sigma*(num/den)^4.
  atmp=[[0,0,0],arg,[0,0,0]]
  dexp=der(atmp)
;
  temp = coef * dexp
  planet_term = temp - total(temp*rr)/total(rr)

  return,planet_term
;
end
