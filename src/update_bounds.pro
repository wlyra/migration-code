function update_bounds,f,boundary
;
  common grid,rr,rrr,dr,dr1,dr2
  common dim,nx,mx,l1,l2
  common order,nghost,order
  common gridtype,gridtype
;
  case boundary of 
;
      'outflow': begin
          for i=1,nghost do begin
              f[l1-i]=f[l1] 
              f[l2+i]=f[l2]
          endfor
       end
;
      'zero-torque': begin
          for i=1,nghost do begin
              f[l1-i]=0.
              f[l2+i]=0.
          endfor
      end
;
      'reflective': begin
          for i=1,nghost do begin
              f[l1-i]=f[l1+i]
              f[l2+i]=f[l2-i]
          endfor
      end
;
      'powerlaw': begin
          if (order ne 2) then begin
              print,'powerlaw boundaries only work for 2nd order derivatives'
              stop
          endif  
          if (gridtype ne 'linear') then begin
              print,'powerlaw boundary not coded for nonlinear grid'
              stop
          endif
;           
          x=rrr
          k1=alog(f[2]/f[1])/alog(x[2]/x[1])
          f[0]=f[1]*(x[0]/x[1])^k1
;
          k2=alog(f[mx-3]/f[mx-2])/alog(x[mx-3]/x[mx-2])
          f[mx-1]=f[mx-2]*(x[mx-1]/x[mx-2])^k2
      end
;
  endcase
;
  return,f
;
end
