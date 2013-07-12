function mdot_to_sigma,mdot
;
  common coeffs,c1,c2,c3,cprime
  common aux,aux_mx,aux_nx
;
;fits of Papaloizou & Terquem 1999
;
  lgmdot1 = (3.1*c1 - cprime)/2.1 
  lgmdot2 = (2*c3-1.1*c2)/0.9     
;
  lgmdot=alog10(mdot)
;
  lgsigma_thin =  lgmdot - c1
  lgsigma_int  = (lgmdot - c2)/2.0
  lgsigma_thick= (lgmdot - c3)/1.1
;
  ithin =where( lgmdot le lgmdot1                       ,nthin) 
  iint  =where((lgmdot gt lgmdot1)and(lgmdot lt lgmdot2),nint)
  ithick=where( lgmdot ge lgmdot2                       ,nthick)
;
  if (nthin ne 0)  then aux_nx[ithin] =lgsigma_thin[ithin]
  if (nint ne 0)   then aux_nx[iint]  =lgsigma_int[iint]
  if (nthick ne 0) then aux_nx[ithick]=lgsigma_thick[ithick]
;
  sigma=10^aux_nx
  return, sigma
;
end
