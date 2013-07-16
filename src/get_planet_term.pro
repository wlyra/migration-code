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
  n=10
  arg = sqrt(rr)*sigma*( (rr < ap) / ( abs(rr-ap)^n + H^n )^(1./n) )^4
  atmp=[[0,0,0],arg,[0,0,0]]
  dexp=der(atmp)
;
  temp = coef * dexp
  planet_term = temp - total(temp*rr)/total(rr)

  return,planet_term
;
end
