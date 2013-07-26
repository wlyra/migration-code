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
  temp_pt = coef * dexp ; temporary planet term

  ; This term should have zero mass contribution and should be nonzero only
  ; near the planet's radius, so a parabola with the negative of its mass
  ; contribution and roots at the analytic term's peaks is added to it to
  ; achieve these properties.
  m0 = min(temp_pt,imid) ; This line's only purpose is to assign to imid
  m1 = max(temp_pt[0:imid],ip1) ; This line assigns to ip1
  m2 = max(temp_pt[imid:nx-1],ip2) ; This line assigns to ip2
  ip2 = ip2 + imid ; Corrects index, the temp_pt slice starts at zero not imid
  mass = total(temp_pt*rr); proportional to mass, actually
  hw = (ip2-ip1)/2. ; half-width of parabola, in terms of index
  x_zero = (ip1+ip2)/2. ; x-coordinate of the parabola's vertex
  para = dblarr(nx) ; initializes to all zeros with double precision
  x_dom = dindgen(nx) ; x domain of the parabola function, double precision
  ptemp = (x_dom - (x_zero + hw))*(x_dom - (x_zero - hw)) ; Creates a parabola
  ; with the correct roots, and arbitrary area between
  scale = mass / total(abs(ptemp[ip1:ip2]*rr[ip1:ip2])) ; scales parabola area
  para[ip1:ip2] = scale * ptemp[ip1:ip2] ; sets values between root, everything
  ; else remains zero

  temp2_pt = temp_pt + para
  planet_term = smooth(temp2_pt, 3, /EDGE_TRUNCATE) ; smooth temp2_pt
  ; smoothing causes some nonzero mass contribution, but it is small

  return,planet_term
;
end