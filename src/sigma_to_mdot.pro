function sigma_to_mdot,sigma
;
  common coeffs,c1,c2,c3,cprime
  common aux,aux_mx,aux_nx
;
  lgsigma1=(c1-cprime)/2.1 
  lgsigma2=(c3-c2)/0.9     
;
  lgsigma=alog10(sigma)
;
  lgmdot_thin = c1 +     lgsigma
  lgmdot_int  = c2 + 2.0*lgsigma
  lgmdot_thick= c3 + 1.1*lgsigma
;
  ithin =where( lgsigma le lgsigma1                         ,nthin) 
  iint  =where((lgsigma gt lgsigma1)and(lgsigma lt lgsigma2),nint)
  ithick=where( lgsigma ge lgsigma2                         ,nthick) 
;
  if (nthin ne 0)  then aux_nx[ithin] =lgmdot_thin[ithin]
  if (nint ne 0)   then aux_nx[iint]  =lgmdot_int[iint]
  if (nthick ne 0) then aux_nx[ithick]=lgmdot_thick[ithick]
;
  mdot=10^aux_nx
  return,mdot
;
end
