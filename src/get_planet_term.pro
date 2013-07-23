function get_planet_term,sigma,timid,cp,gamma,q,sqrtgm,ap,omega
;
    ;common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
  common grid,rr,rrr,dr,dr1,dr2
;
  nx=n_elements(rr)
  rr1=1./rr
;
  cs2=timid*cp*(gamma-1)
  H=sqrt(cs2*gamma)/omega
;
  coef=sign(ap-rr) * q^2 * sqrtgm * rr1
  num = rr < ap 
;
; The original term was max(abs(rr-ap),H)
;
  tmp = abs(rr-ap)
  ; den = (tmp^n_smooth + H^n_smooth)^(1./n_smooth)
  den = tmp > H
;
  arg = sqrt(rr)*sigma*(num/den)^4.
  atmp=[[0,0,0],arg,[0,0,0]]
  dexp=der(atmp)
;
  temp_pt = coef * dexp

  AU = 1.49d13 ; cm/AU
  m0 = min(temp_pt,imid)
  m1 = max(temp_pt[0:imid],ip1)
  m2 = max(temp_pt[imid:nx-1],ip2)
  ip2 = ip2 + imid
  mass = total(temp_pt*rr); proportional to mass
  hw = (ip2-ip1)/2. ; half-width of parabola, in terms of index
  x_zero = (ip1+ip2)/2. ; x-coordinate of the parabola's vertex
  para = dblarr(nx)
  x_dom = dindgen(nx) ; x domain of the function
  ptemp = (x_dom - (x_zero + hw))*(x_dom - (x_zero - hw))
  scale = mass / total(abs(ptemp[ip1:ip2]*rr[ip1:ip2]))
  para[ip1:ip2] = scale * ptemp[ip1:ip2] ; everything else remains zero

  temp2_pt = temp_pt + para
  planet_term = smooth(temp2_pt, 3, /EDGE_TRUNCATE)
  ; smoothing causes some nonzero mass contribution, but it is small

  return,planet_term
;
end