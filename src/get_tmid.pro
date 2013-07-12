function get_tmid,sigma,omega,mdot,alpha,told,label
;
  common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
  common root,method
  common iteration,maxit
;
; we know the temperature always decreases, so use the last 
; temperature as starting point
; 
  right=told+epsi
;
; the temperature should not change much between timesteps
;
  if (label eq 'start') then begin
      left=tb
  endif else begin
      left=.95*told > tb
      ;left=told > tb
  endelse
;
; make sure we bracketed the root
;  
  x1=left & x2=right 
  if (x1 eq x2) then begin
      print,'bad initial range for newton-raphson'
      stop
  endif
;      
  func=get_phi(x1,sigma,omega,mdot,alpha) & fl=func.f
  func=get_phi(x2,sigma,omega,mdot,alpha) & fh=func.f
;
  for j=0,maxit-1 do begin
      if (fl*fh lt 0) then goto,root_bracketed_okay
      factor=1.6
      if (abs(fl) lt abs(fh)) then begin
          x1 = x1-factor*(x2-x1) > tb
          func=get_phi(x1,sigma,omega,mdot,alpha) & fl=func.f
          ;print,'reduce x1 to bracket root ij,x1,x2,fl,fh',j,x1,x2,fl,fh
      endif else begin
          x2 = x2+factor*(x2-x1)
          func=get_phi(x2,sigma,omega,mdot,alpha) & fh=func.f
          ;print,'increase x2 to bracket root ij,x1,x2,fl,fh',j,x1,x2,fl,fh
      endelse
  endfor
root_bracketed_okay:

;
  case method of
      'bisection': temperature=bisect(left,right,sigma,omega,mdot,alpha)
      'newton':    temperature=newton_raphson(left,right,sigma,omega,mdot,alpha,fl,fh)
  endcase
;
  return,temperature
;
end
