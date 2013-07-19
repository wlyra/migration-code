pro main
;
    common boolean,true,false
    common grid,rr,rrr,dr,dr1,dr2
    common gridtype,gridtype
    common dim,nx,mx,l1,l2
    common order,nghost,order
    common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
    common teffs_for_tau,T1,T2,T3,T4,T5,T6,T7,T8,T9
    common coeffs,c1,c2,c3,cprime
    common aux,aux_mx,aux_nx
    common root,method
    common iteration,maxit
;
    false=0 & true=1
    maxit=100
;
; Read the input file
;
    exist_input=file_test('start.in')
    if (exist_input eq 0) then begin
;
; There is no input file. Print an 
; informative error message.
;
    noinputfile
;
    endif else begin
      openr, 1, 'start.in'
      line=''
      while (not eof(1)) do begin
        readf, 1, line
        if (execute(line) ne 1) then $
          message, 'There was a problem with start.in', /INF
      endwhile
      close, 1
    endelse
;
; Define switches
;
    if (evolve_density eq 'no') then begin
        ldensity=false
    endif else if (evolve_density eq 'yes') then begin 
        ldensity=true
    endif else begin
        print,"use evolve_density='yes' or 'no' in start.in"
    endelse
;
    if (evolve_temperature eq 'no') then begin
        ltemperature=false
    endif else if (evolve_temperature eq 'yes') then begin 
        ltemperature=true
    endif else begin
        print,"use evolve_temperature='yes' or 'no' in start.in"
    endelse
;
    if (evolve_planet eq 'no') then begin
        lplanet=false
    endif else if (evolve_planet eq 'yes') then begin 
        lplanet=true
    endif else begin
        print,"use evolve_planet='yes' or 'no' in start.in"
    endelse
;
    if (update_timestep eq 'no') then begin
        lupdate_timestep=false
    endif else if (update_timestep eq 'yes') then begin 
        lupdate_timestep=true
    endif else begin
        print,"use update_timestep='yes' or 'no' in start.in"
    endelse
;
    if (restart eq 'no') then begin
        lrestart=false
    endif else if (restart eq 'yes') then begin 
        lrestart=true
    endif else begin
        print,"use restart='yes' or 'no' in start.in"
    endelse
;
    if (print_to_screen eq 'no') then begin
        lprint_to_screen=false
    endif else if (print_to_screen eq 'yes') then begin
        lprint_to_screen=true
    endif else begin
        print,"use print_to_screen='yes' or 'no' in start.in"
    endelse
;
    if (default_plots eq 'no') then begin
        ldefault_plots=false
    endif else if (default_plots eq 'yes') then begin
        ldefault_plots=true
    endif else begin
        print,"use default_plots='yes' or 'no' in start.in"
    endelse

;
    window,retain=2,ysize=900,xsize=700
    !p.multi=[0,1,4]
    !p.charsize=3.
;
; coefficients for 3rd order Runge-Kutta
;
    alpha_ts   = [0.   , -5./9.  ,-153./128.]
    beta_ts    = [1./3., 15./16. ,   8./15. ]
;
; constants
;
    msun=1d33                   ;g
    mearth=6d27                 ;g
    AU=1.49d13                  ;cm
    r0=AU
    yr=31556926.                ;s
    Myr=1d6*yr
    GNewton=6.67d-8
    Rgas=8.314d7
    stefan_boltzmann=5.6704d-5  ; erg cm-2 s-1 K-4
    dsigma_mass = [] ; originally nothing
    swind_mass = [] ; originally nothing
    term1_mass = [] ; originally nothing
;
; grid
;
    if (not((order eq 2) or (order eq 6))) then begin
      print,'the order has to be 2 or 6, you are using ',order
      stop
    endif
    nghost=order/2
    mx=nx+2*nghost & l1=nghost & l2=mx-1-nghost
    tmp=dblarr(mx) 
;
    if (gridtype eq 'linear') then begin
        rr=grange(r_int,r_ext,nx)*AU & dr=rr[1]-rr[0] 
        tmp[l1:l2]=rr
        for i=1,nghost do begin
            tmp[l1-i]=tmp[l1]-i*dr
            tmp[l2+i]=tmp[l2]+i*dr
        endfor
        dr_tmp=dr
        dr=replicate(dr_tmp,nx)
        rrr=tmp
    endif else begin
        lnrr=grange(alog(r_int*AU),alog(r_ext*AU),nx) & dlnr=lnrr[1]-lnrr[0]
        tmp[l1:l2]=lnrr
        for i=1,nghost do begin
            tmp[l1-i]=tmp[l1]-i*dlnr
            tmp[l2+i]=tmp[l2]+i*dlnr
        endfor
        rr=exp(lnrr) & rrr=exp(tmp) & dr=rrr[l1:l2]-rrr[l1-1:l2-1]
    endelse
    rr1=1./rr & dr1=1./dr  & dr2=dr1^2
    r_ref=AU
;
; write to data directory, test first if the directory exists
;
    exist_datadir=file_test('data',/directory)
    if (exist_datadir eq 0) then begin
        spawn,'mkdir data'
    endif else begin
; it exists and it is a new run, so delete the old VARS to avoid
; confusion
        if (lrestart eq false) then $
          spawn,'rm -rf ./data/VAR*sav'
    endelse
;
; write grid to a file it is not yet written
; 
    exist_gridfile=file_test("./data/grid.sav")
    if (exist_gridfile eq 0) then save,rr,r_ref,dr,filename='./data/grid.sav'
;
; auxiliaries
;
    aux_nx=dblarr(nx) & aux_mx=dblarr(mx)
;
; fix the mass accretion rate
;
    mdot_init = mdot_input*(msun/yr) 
    mdot=replicate(mdot_init,nx)
;
; Get the coefficients of Papaloizou & terquem
; 
    coeff=get_coeff(alpha)
    c1=coeff.c1 & c2=coeff.c2 & c3=coeff.c3 & cprime=coeff.cprime
;
; use these coeffcients to fix the density
;
    sigma=dblarr(nx)
    sigma=mdot_to_sigma(mdot)
    sigma_init=sigma
;
; set also the viscosity, and a density floor
;
    nu=mdot/(3*!dpi*sigma)
    density_floor=1d-7
;
; sigma_nu=mdot/(3*pi) is a handy quantity
;
    sigma_nu=dblarr(mx)
;
; Calculate temperature at surface, 
; following (Paps & Terquem 1999)
;
    GM=GNewton*Msun
    sqrtgm=sqrt(GM)
    omega=sqrtgm/rr^1.5
    gamma=1.4
    gamma1=1./gamma
    mmol=2.4 
    Rgasmu=Rgas/mmol
    cp=gamma*Rgasmu/(gamma-1)
    cv=cp/gamma
;
; shortcuts
;
    stbz=stefan_boltzmann
    ;epsi=4*temperature_precision^2
    epsi=2*temperature_precision
    one_over_three_pi=1./(3*!dpi)
    tsave1_myr=1./(tsave*myr)
    myr1=1./myr
    msun1=1./msun
    msun_yr=msun/yr
    msun_yr1=1./msun_yr
    mearth1=1./mearth
    r_ref1=1./r_ref
;
; opacities - have to define the temperature where change occur
;
    T1=132. & T2=170. & T3=375. & T4=390.
    T5=580. & T6=680. & T7=960. & T8=1570. & T9=3730.
    tau=dblarr(nx)
;
; go directly to tmid without calculating teff
;
    tmid=dblarr(nx) & tb=background_temperature
    tmid_old=replicate(.9*T9,nx)
    for i=0,nx-1 do begin
        tmid[i]=get_tmid(sigma[i],omega[i],mdot[i],alpha,tmid_old[i],'start')
    endfor
;    
    tmid_orig=tmid
    ;!p.multi=0
    ;plot,rr/au,tmid,xs=3,ys=3,/ylog,/xlog
    ;stop
;
; wind
;
    swind = get_wind(wind,mwind_input,AU,msun,yr)

    mwind=mwind_input*(msun/yr)
    rmax=rr[nx-1]               ; 30*AU
    rg=5*AU         ; this can be better calculated GM/c^2, or similar
    den=2*!dpi*(rmax-rg)*rr
    il=where(rr le rg,nl) & ig=where(rr gt rg,ng)
    swind=dblarr(nx)
    if (nl ne 0) then swind[il]=0.
    if (ng ne 0) then swind[ig]=mwind/den[ig]
;
; planet
;
    ap=planet_position*AU & dap=ap*0.
    mp=planet_mass*mearth
    q=mp/msun
    print,q; originally asdfasdf
;
; temperature where the transition to type 2 migration occurs
; computed as H = r_hills define. The quantities do not change in time
;
    ;temperature_type2 = GM/(rr*cp*(gamma-1)) * (q/3.)^(2./3)
;
; diagnostics
;
    ndim=itmax/itdiagnos+1
;
    time=dblarr(ndim)
    mass=dblarr(ndim)
    np=n_elements(mp) 
    position=dblarr(ndim,np)
    tmax=tmax_myr*myr
    tcheck=0.
    t=0.
;
;  timestep
;
    courant=0.4
    dt_density=1d33 & dt_planet=1d33
    if (ldensity) then dt_density=courant*min(dr^2/(3*nu))
    if (lplanet)  then begin
    dap=dblarr(np)
      for ip=0,np-1 do begin
        omegap=sqrtgm*ap[ip]^(-1.5)
        torque=get_torque(sigma,tmid,q[ip],GM,ap[ip],omegap,eos)
        dap[ip] = 2.*torque/(mp[ip]*ap[ip]*omegap)
      endfor
      dt_planet=courant*min(dr/abs(dap))
    endif  
    tmp = dt_planet < dt_density
    dt_max = 1000.*yr ;minimum timestep
    dt=tmp < dt_max
    dt_beta_ts=dt*beta_ts
;
    print,'dt=',dt/yr,'yr'
    print,'entering evolution'
;
    ic=0
;
; Write the VAR0 file with initial condition
;
    if (lrestart eq false) then begin
        isave=0
        save,t,sigma,tmid,ap,filename='./data/VAR0.sav'
        print,'time=',t*myr1,' Myr. Wrote VAR0 to disk.'
    endif else begin
        isave=isave_input
        restore,filename='./data/VAR'+strtrim(isave,2)+'.sav'
        tcheck=fix(t*tsave1_myr)
    endelse
;
; Write time-series file. Check first if the file exists. 
; If so, rename the old.
;
    exist_timeseries=file_test("./data/time_series.dat")
    if (exist_timeseries eq 1) then begin
        spawn,'mv ./data/time_series.dat ./data/time_series_old.dat'
    endif 
;
    openw,lun,'./data/time_series.dat',/get_lun
    printf,lun,'#it----time--dt--rhomax----rhomin----mass----TTmax----TTmin----mdot----ap'
    free_lun,lun
    loadct,5
;
; Start computations
;
    print,'#it----time--dt--rhomax----rhomin----mass----TTmax----TTmin----mdot----ap'
    for it=0L,itmax-1 do begin
; 
; recompute the timestep if density changes
;
        if (lupdate_timestep eq true) then begin
          if (ldensity eq true) then begin
              mdot=sigma_to_mdot(sigma)
              nu=mdot/(3.*!dpi*sigma)
              dt_density=courant*min(dr^2/(3*nu))
          endif else begin
              dt_density=1d33
          endelse
;
          if (lplanet eq true) then begin
              dt_planet=courant*min(dr/abs(dap))
          endif else begin
              dt_planet=1d33
          endelse
          tmp=dt_density < dt_planet
          dt=tmp < dt_max
          ;if (tmp gt dt_max) then begin
          ;    print,'using dt_max',tmp/yr,dt_max/yr
          ;endif
          dt_beta_ts=dt*beta_ts
      endif

;
; subtimestepping
;
      for itsub=0,2 do begin
;
; evolution
;
        if (itsub eq 0) then begin
          dsigma=0.
          dap=replicate(0.,np)
          ds=0
        endif else begin
          dsigma=alpha_ts[itsub]*dsigma
          dap=alpha_ts[itsub]*dap
          ds=alpha_ts[itsub]*ds
        endelse
        ds=ds+1.
;
; Get the new mdot, temperature, and viscosity (for the
; other timesteps)
;
        if ((ldensity eq true) or (ltemperature eq true)) then $
          mdot=sigma_to_mdot(sigma)
;        
        if (ltemperature eq true) then begin
          tmid_old=tmid ;upper limit of the bisection
          for i=0,nx-1 do begin
              tmid[i]=get_tmid(sigma[i],omega[i],mdot[i],alpha,tmid_old[i],'run')
          endfor
      endif
;
; Evolution density. First take Sigma*nu 
;
        if (ldensity eq true) then begin 
          sigma_nu[l1:l2]=mdot*one_over_three_pi
          sigma_nu=update_bounds(sigma_nu,'outflow')
;
; derivatives grad(Sigma*nu) and Laplace(Sigma*nu)
;
          del2sigmanu = der2(sigma_nu)
          gsigmanu    =  der(sigma_nu)       
;
          planet_term=0
          for ip=0,np-1 do begin
             tmp_planet = get_planet_term(sigma,tmid,cp,gamma,q[ip],sqrtgm,ap[ip],omega)
             planet_term=planet_term + tmp_planet
          endfor
;
; evolution equation with sigma*nu as dependent variable
;
          dsigma = 3*del2sigmanu + 4.5*rr1*gsigmanu - swind + planet_term
          ; plot,rr*r_ref1,dsigma
          ; dsigma_mass = [dsigma_mass, [2*!dpi*total(dsigma*rr*dr)] ]
          ; swind_mass = [swind_mass, [2*!dpi*total(-swind*rr*dr)] ]
          ; term1_mass = [term1_mass, [2*!dpi*total((dsigma+swind-planet_term)*rr*dr)] ]
        endif
;
; planet evolution (Lagrangian)
;
        if (lplanet eq true) then begin
          for ip=0,np-1 do begin
            omegap=sqrtgm*ap[ip]^(-1.5)
            torque=get_torque(sigma,tmid,q[ip],GM,ap[ip],omegap,eos)
            dap[ip] = 2.*torque/(mp[ip]*ap[ip]*omegap)
          endfor
        endif
;
; advance runge-kutta
;
        if (ldensity eq true) then sigma = sigma + dt_beta_ts[itsub]*dsigma        
        if (lplanet  eq true) then ap    = ap    + dt_beta_ts[itsub]*dap
;
; advance time
;
        t = t + dt_beta_ts[itsub]*ds
;
; apply density floor
;
        if (ldensity eq true) then begin
          i=where(sigma le 0,nii) 
          if (nii ne 0) then sigma[i]=density_floor
        endif
        if (lplanet eq true) then begin
          i=where(ap le rr[1],nii) 
          if (nii ne 0) then ap[i]=rr[1]
          i=where(ap ge rr[nx-1],nii) 
          if (nii ne 0) then ap[i]=rr[nx-1]
        endif


      endfor                      ; end subtimesteps
;
; diagnostics
;
      if ((it mod itdiagnos) eq 0) then begin
        time[ic]=t
        mass[ic]=2*!dpi*total(sigma*rr*dr)
        position[ic,*]=ap
        rhomax=max(sigma)
        rhomin=min(sigma)
        TTmax=max(tmid)
        TTmin=min(tmid)
        mdotm=mean(mdot)
;
; print to screen
;
        
        if (lprint_to_screen) then begin
          line="print,it,time[ic]*myr1,dt*myr1,rhomax,rhomin,mass[ic]*msun1,ttmax,ttmin,mdotm*msun_yr1"
          for ip=0,np-1 do begin
              line=line+",ap["+strtrim(ip,2)+"]*r_ref1"
          endfor
          if (execute(line) ne 1) then begin
              print,line
              print,'There was a problem with printing to screen'
              stop
          endif
        endif
;        print,it,time[ic]*myr1,rhomax,rhomin,mass[ic]*msun1,ttmax,ttmin,mdotm*msun_yr1,ap[0]*r_ref1,ap[1]*r_ref1,ap[2]*r_ref1
;
; print to file
;
        openw,lun,'./data/time_series.dat',/get_lun,/append
        line="printf,lun,it,time[ic]*myr1,rhomax,rhomin,mass[ic]*msun1,ttmax,ttmin,mdotm*msun_yr1"
        lformat="format='(i8,e17.8,e17.8,e17.8,e17.8,e17.8,e17.8,e17.8"
        for ip=0,np-1 do begin
            line=line+",ap["+strtrim(ip,2)+"]*r_ref1"
            lformat=lformat+",e17.8"
        endfor
        lformat=lformat+")'"
        lline=line+','+lformat
        if (execute(lline) ne 1) then begin
            print,lline
            print, 'There was a problem with printing to time_series.dat'
            stop
        endif
        free_lun, lun
;
        if (time[ic] gt tmax) then begin
          print,'Maximum time exceeded'
          goto, out_of_loop
        endif
;
        ;if (ldensity eq true) then begin
        ;  if (mass[ic] lt mearth) then begin
        ;    print,'Total mass less than an Earth mass'
        ;    goto, out_of_loop
        ;  endif
        ;endif
;        
        ;if (ltemperature eq true) then begin 
        ;  if (min(tmid-temperature_type2) le 0) then begin
        ;    print,'Type I migration should start - stop computations.'
        ;    goto, out_of_loop
        ;  endif
        ;endif
;
; increase counter for on-the-fly plotting
;
        ic=ic+1
;
    endif
;
; write the variables to the disk at a given frequency
;
    tcheck_old=tcheck
    tcheck=fix(t*tsave1_myr) ; t/tsave  - when it changes, time to write to disk
    if (tcheck ne tcheck_old) then begin
        isave=isave+1
        save,t,sigma,tmid,ap,filename='./data/VAR'+strtrim(isave,2)+'.sav'
        if (lprint_to_screen) then begin
          print,'time=',t*myr1,' Myr. Wrote VAR'+strtrim(isave,2)+' to disk.'
        endif
    endif

      if ((it mod itplot) eq 0) then begin
;
          erase
          if (ldefault_plots) then begin
              plot,rr*r_ref1,sigma,ys=1,$
                title='t='+strtrim(t*myr1,2)+' Myr',$
                yr=[.00004,500.],xtitle='!8r!x';,/ylog,xs=3,/xlog; originally .00004,5000.
              oplot,rr*r_ref1,sigma_init,li=1

              plot,rr*r_ref1,tmid,xs=3,$
                title='temperature',xtitle='!8r!x',/ylog,yr=[9,2000],ys=1,/xlog
              oplot,rr*r_ref1,tmid_orig,li=1

              plot,time[0:ic-1]*Myr1,100*mass[0:ic-1]/mass[0],ys=3,$
                title='Disk mass - Remaining percentage',xtitle='time (Myr)',xr=[0,tmax_myr],yr=[0,100]

              plot,time[0:ic-1]*Myr1,position[0:ic-1,0]*r_ref1,$
                title='planet position',xtitle='time (myr)',xr=[0,tmax_myr],$
                yr=[0.1,30],ys=1,/ylog,/nodata
              for ip=0,np-1 do begin
                  if (mp[ip]/mearth eq 1.) then begin
                      cor=100
                  endif else $
                  if (mp[ip]/mearth eq 10) then begin
                      cor=150
                  endif else begin
                      cor=50
                  endelse
                  oplot,time[0:ic-1]*Myr1,position[0:ic-1,ip]*r_ref1,color=cor
              endfor
          endif else begin
              ; Alternate plots

              plot,rr*r_ref1,sigma,ys=1,$
                title='t='+strtrim(t*myr1,2)+' Myr',$
                yr=[.00004,500.],xtitle='!8r!x';,/ylog,xs=3,/xlog; originally .00004,5000.
              oplot,rr*r_ref1,sigma_init,li=1

              ; plot,time[0:ic-1]*Myr1,dsigma_mass,xr=[0,tmax_myr]
              
              ; plot,time[0:ic-1]*Myr1,swind_mass,xr=[0,tmax_myr]

              ; plot,time[0:ic-1]*Myr1,term1_mass,xr=[0,tmax_myr],title='term1_mass vs. time'
              
              plot,rr*r_ref1,dsigma,title='dsigma vs. radius'

              plot,rr*r_ref1,planet_term,title='planet_term vs. radius'

              plot,time[0:ic-1]*Myr1,100*mass[0:ic-1]/mass[0],ys=3,$
                title='Disk mass - Remaining percentage',xtitle='time (Myr)',xr=[0,tmax_myr],yr=[0,100]
          endelse

;        wait,0.1

        endif
;
    endfor
out_of_loop:

    print,'done'
    !p.multi=0
end
