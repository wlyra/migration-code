function get_torque,sigma,tmid,q,GM,ap,omegap,eos
;
    common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
    common grid,rr,rrr,dr,dr1,dr2
    common dim,nx,mx,l1,l2
    common aux,aux_mx,aux_nx
    common gridtype,gridtype
;
; use linear interpolation
;
    if (gridtype eq 'linear') then begin
        ddr1=mean(dr1)
        i1=floor((ap-rr[0])*ddr1)
    endif else begin
        dlnr=alog(rr[1])-alog(rr[0])
        i1=floor((alog(ap)-alog(rr[0]))/dlnr)
        ddr1=dr1[i1]
    endelse
;   
    deltar=ap-rr[i1]
    i2=i1+1
;
    tmidp = tmid[i1]+( tmid[i2]- tmid[i1])*deltar*ddr1
    sigmap=sigma[i1]+(sigma[i2]-sigma[i1])*deltar*ddr1
;
    aux_mx[l1:l2]=tmid  & aux_mx=update_bounds(aux_mx,'outflow')
    gtmid=der(aux_mx)
;
    aux_mx[l1:l2]=sigma & aux_mx=update_bounds(aux_mx,'outflow')
    gsigma=der(aux_mx)
;
    gtmidp = gtmid[i1]+( gtmid[i2]- gtmid[i1])*deltar*ddr1
    gsigmap=gsigma[i1]+(gsigma[i2]-gsigma[i1])*deltar*ddr1
;
    beta = -ap/tmidp*gtmidp
    alpha= -ap/sigmap*gsigmap
;
    torque_iso=-(2.5-0.5*beta-0.1*alpha) - 1.4*beta + 1.1*(1.5-alpha)
    ksi=beta-(gamma-1)*alpha
    torque_adi=gamma1*(-2.5-1.7*beta + 0.1*alpha + 1.1*(1.5-alpha) $
                       + 7.9*ksi*gamma1)
;
    case eos of 
;
        'iso':torque_norm=torque_iso
;
        'adi':torque_norm=torque_adi
;
        'blend': begin
;
; Blend the iso and adiabatic torques with the radiation number
; For that we need to calculate the radiation number - tcool/tdyn
; The cooling time is the timescale of E/Edot
;
;     tcool=E/Edot 
;
; Edot is the cooling, div(F). Integrating it in volume and applying 
; Gauss' theorem, it becomes a surface density 
;    
;     tcool=Int E dv/ Int F da 
; 
; Assuming both E and F to be constant 
;
;     tcool=(E*V)/(F*A)
; 
; We have E=cv*rho*T, V=4/3*pi*H^3, A=4*pi*H^2 and 
; F=stbz*Teff^4=sigma*Tmid^4/tau, so 
;
;    tcool=cv*Sigma*tau/(6*stbz*tmid^3)
;    tdyn=2*!dpi/omegap
;    rn = radiation_number = tcool/tdyn
;
            ctaup=get_tau(tmidp,sigmap,omegap) & kappa=ctaup.kappa
            taup=.5*kappa*sigmap
            taueff=0.375*taup + 0.43301270 + .25/taup
            rn=cv*sigmap*taueff*omegap/(12*!dpi*stbz*tmidp^3)
            torque_norm=(rn^2*torque_adi+torque_iso)/(rn+1)^2
        end
    endcase
;
; These torques are normalized. The normalization factor is 
;
;   torque01=q^2/h2*sigmap*omegap^2*ap^4
;
; The following substitutions simplify the matter, and eliminate ap
;
;;omega2p=GM/ap^3
;;scale_height_sq=cs2p/omega2p
;;h2=scale_height_sq/ap^2
;h2=cs2*ap/GM
;
    cs2p=tmidp*cp*(gamma-1) 
    torque0=q^2/cs2p*sigmap*GM^2
;
    torque=torque_norm*torque0
;
    return,torque
;
end
