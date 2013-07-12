function get_phi,temp,sigma,omega,mdot,alpha
;
  common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
  common root,method
;
  ;try the nagakawa fit
;
  ctau=get_tau(temp,sigma,omega) 
  kappa_cte=ctau.kappa_cte 
  b_exp=ctau.b_exp & a_exp=ctau.a_exp
  kappa=ctau.kappa & tau=.5*kappa*sigma & tau1=1./tau
  taueff=0.375*tau + 0.43301270 + .25*tau1
  edot=3./(4*!pi)*mdot*omega^2
  tb1=temp^(b_exp-1)
  tb2=temp^(b_exp-2)
;
  temp1=1./temp
;
  case norm of 
      'terquem': begin
;
; Papaloizou & Terquem 1999
;
          aa=2*stbz*(temp^4/taueff-tb^4)
          bb=-.5*edot
          cc=-9*alpha*omega*temp*cp*(gamma-1)*taueff^0.25/(1.6d-3*temp)
          phi=aa+bb+cc
          dphi=0. & d2phi=0.
; not calculate for newton-raphson
          if (method eq 'newton') then begin
             print,"Cannot use newton-raphson and terquem at the same time. Slow "
             print,"convergence. Use norm='nakagawa' instead"
             stop
          endif
      end
      'nakagawa': begin
;
; Nakamoto & Nakagawa 1994
; 
          phi  = 2*stbz*(temp^4-tb^4)-taueff*edot
; analytical derivative
          dtaudt=(b_exp-.5*a_exp)*tau*temp1
          dphi  =  8*stbz*temp^3 - .25*edot*dtaudt*(1.5-tau1^2)
;
          dtau2dt=(b_exp-.5*a_exp)*(b_exp-.5*a_exp-1)*tau*temp1^2
          d2phi = 24*stbz*temp^2 - .25*edot*dtau2dt*(1.5-tau1^2) $
            + .5*edot*tau1^3*dtaudt^2
      end
  endcase
;
  func={f:phi,df:dphi,d2f:d2phi}
  return,func
;
end
