function get_coeff,alpha
;
  common grid,rr,rrr,dr,dr1,dr2
;
  logc1 = 0.9360636 + 0.1195816*alog10(alpha) + $
       (0.0233002 - 0.0061733*alog10(alpha))*alog10(rr)
  c1=10^logc1
  logc3 = 0.7782080 + 0.0545617*alog10(alpha) + $
         (0.0366565 - 0.0019087*alog10(alpha))*alog10(rr)
  c3=10^logc3
  cprime = 16.0897161 + 2.0665*alog10(alpha)
  c2 = (1.1*c1+cprime)/2.1
;
  coeff={c1:c1,c2:c2,c3:c3,cprime:cprime}
  return,coeff
;
end
