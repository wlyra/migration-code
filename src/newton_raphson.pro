function newton_raphson,left,right,sigma,omega,mdot,alpha,fl,fh
;
    common thermo,stbz,tb,cp,cv,gamma,gamma1,epsi,norm
    common iteration,maxit
;
    it=0
;
    x1=left & x2=right 
;
    if (fl eq 0) then begin
        out=x1
        goto,found_root
    endif
;
    if (fh eq 0) then begin
        out=x2
        goto,found_root
    endif
;    
    if (fl lt 0) then begin
;Orient the search so that f(x1)<0
        xl=x1 & xh=x2
    endif else begin
        xh=x1 & xl=x2
    endelse

    rts=0.5*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    func=get_phi(rts,sigma,omega,mdot,alpha)
    f=func.f & df=func.df & d2f=func.d2f
;    
    for j=0,maxit-1 do begin
        bla1=(rts-xh)*df-f
        bla2=(rts-xl)*df-f
;
        if ((bla1*bla2 gt 0) or (abs(2.*f) gt abs(dxold*df))) then begin
;                  
            dxold=dx
            dx=0.5*(xh-xl)
            rts=xl+dx
            if (xl eq rts) then begin
                out=xl
                goto,found_root
            endif
        endif else begin  
            dxold=dx
;            
            a=1.2 < 1-f*d2f/(2.*df^2)
            halley_factor=0.8 > a
;
            dx=f/(df*halley_factor)
            temp=rts
            rts=rts-dx
            if (temp eq rts)  then begin
                out=rts
                goto,found_root
            endif
        endelse
        if (abs(dx) lt epsi) then begin
            out=rts
            goto,found_root
        endif
        func=get_phi(rts,sigma,omega,mdot,alpha)
        f=func.f & df=func.df & d2f=func.d2f
        if (f lt 0) then begin
            xl=rts
        endif else begin
            xh=rts
        endelse
    endfor

    print,'exceed maximum of iterations'
    print,'left,right=',xl,xh
    stop
;
found_root:    
;
    ;print,j
    return,out
;
end
