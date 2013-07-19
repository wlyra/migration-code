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
  temp_pt = coef * dexp

  AU = 1.49d13 ; cm/AU
  m0 = min(temp_pt,mid)
  m1 = max(temp_pt[0:mid],p1)
  m2 = max(temp_pt[mid:nx-1],p2)
  p2 = p2 + mid
  sum = total(temp_pt); not mass, but scales assuming narrow nonzero region of temp_pt
  hw = (p2-p1)/2. ; half-width of parabola, in terms of index
  x_zero = (p1+p2)/2. ; x-coordinate of the parabola's vertex
  para = dblarr(nx)
  x_dom = dindgen(nx) ; x domain of the function
  ptemp = (x_dom - (x_zero + hw))*(x_dom - (x_zero - hw))
  para = para < 3.*sum/(4.*hw^3) * ptemp

  planet_term = temp_pt + para
  ; print,total(temp_pt),total(para) ; different by ~3% right now :/

  return,planet_term
;
end