function get_tau,tt,sigma,omega
;
; from Bell et al. 97
;
  common teffs_for_tau,T1,T2,T3,T4,T5,T6,T7,T8,T9
  common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
;
  kk=0. ; just initialize it
;
  if (TT le T1) then begin
      k=2d-4 & a=0 & b= 2.1  & kk=k*tt^b
  endif else $
  if ((TT gt T1) and (TT le T2)) then begin
      k=3.   & a=0 & b=-0.01 & kk=k*tt^b
  endif else $
  if ((TT gt T2) and (TT le T3)) then begin
      k=0.01 & a=0 & b= 1.1  & kk=k*tt^b
  endif else $
  if ((TT gt T3) and (TT le T4)) then begin
      k=5d4  & a=0 & b=-1.5  & kk=k*tt^b
  endif else $
  if ((TT gt T4) and (TT le T5)) then begin
      k=0.1  & a=0 &  b= 0.7 & kk=k*tt^b
  endif else $
  if ((TT gt T5) and (TT le T6)) then begin
      k=2d15 & a=0 & b=-5.2  & kk=k*tt^b
  endif else $
  if ((TT gt T6) and (TT le T7)) then begin
      k=0.02 & a=0 & b= 0.8  & kk=k*tt^b
  endif else $
  if ((TT gt T7) and (TT le T8)) then begin
      logk=81.3010 & a=1. & b=-24.
      H=sqrt(TT*cp*(gamma-1))/omega
      rho=sigma/(2*H)
      logkk=logk+a*alog10(rho)+b*alog10(TT)
      kk=10^(logkk)
      k=1d33
  endif else $
  if ((TT gt T8) and (TT le T9)) then begin
      k=1d-8 & a=2./3 & b=3.
      H=sqrt(TT*cp*(gamma-1))/omega
      rho=sigma/(2*H)
      kk=k*rho^a*tt^b
  endif else begin
      print,'get_tau: temperature higher than ',T9
      stop
  endelse
;
  ctau={kappa_cte:k,b_exp:b,kappa:kk,a_exp:a}
;
  return,ctau
;
end
