pro panelplots,cu=cu,ipanel=p,lont=lont,latt=latt,w0=w0,zoom=zoom,swcorner=swc,size=size,xran=xran,yran=yran $
              ,mxdis=mxdis,block=block,dev=dev,smooth=smooth,raw=raw,profi=profi,xcuqi=cuqi,dblock=dblock $
              ,dprofi=dprofi,corr=corr,aniso=aniso,pdev=pdev,sdev=sdev,super=super,nvar=nvar $
              ,wedge=wedge,nodes=nodes,wedgo=wedgo,nodos=nodos,dnodes=dnodes,fnodes=fnodes $
              ;  processing
              ,nsm=sm $
              ;  overplotting  
              ,xl=x $
              ;  control
              ,stop=stop

if not keyword_set(w0) then w0 =1

if keyword_set(xran) then begin
   xr=xran
endif
if keyword_set(yran) then begin
   yr=yran
endif
if keyword_set(zoom) then begin
   xr=[zoom(0) , zoom(2) ]
   yr=[zoom(1) , zoom(3) ]
endif
if keyword_set(zoom) then begin
   xr=[zoom(0) , zoom(2) ]
   yr=[zoom(1) , zoom(3) ]
endif
if keyword_set(swc) then begin
   xr=[swc(0),swc(0)+1.2*size]
   yr=[swc(1),swc(1)+size]
endif

if keyword_set(x) then begin
   ox=where( x.xspk ge xr(0) and x.xspk le xr(1) $ 
         and x.yspk ge yr(0) and x.yspk le yr(1) $
         and x.panel eq p )
endif

cell_fill=0
Title$='Title '
unit$ ='unit'

if keyword_set(raw) then begin
   f=cu.raw
   lev=(findgen(16)-7.999999)*300.
   title$='Raw Topography' & unit$='m'
endif
if keyword_set(smooth) then begin
   f=cu.smooth
   lev=(findgen(16)-7.99999)*200.
   title$='Smooth Topography' & unit$='m'
endif
if keyword_set(dev) then begin
   f=cu.dev
   lev=(findgen(16)-7.99999)*200.
   title$='Raw Topography minus Smooth' & unit$='m'
endif
if keyword_set(pdev) then begin
   f=cu.pdev
   lev=(findgen(16)-7.99999)*200.
   title$='+ve (Raw Topography minus Smooth)' & unit$='m'
endif
if keyword_set(profi) then begin
   f=cu.profi
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Profi' " & unit$='m'
endif
if keyword_set(dprofi) then begin
   f=cu.profi-cu.dev
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Profi' - Dev " & unit$='m'
endif
if keyword_set(mxdis) then begin
   f=cu.mxdis
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Mxdis' " & unit$='m'
endif
if keyword_set(block) then begin
   f=cu.block
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Block' " & unit$='m'
endif
if keyword_set(dblock) then begin
   f=cu.block-cu.pdev
   lev=(findgen(16)-7.99999)*200.
endif
if keyword_set(fnodes) then begin
   f=cu.nodes
   g=cu.dev
   oo=where( g lt 0.)
   f(oo)=g(oo)
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Nodes' fixed over Dev<0 " & unit$='m'
endif
if keyword_set(nodes) then begin
   f=cu.nodes
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Nodes' " & unit$='m'
endif
if keyword_set(nodos) then begin
   f=cu.nodos
   lev=(findgen(16)-7.99999)*200.
   title$=" 'NodOs' (restricted) " & unit$='m'
endif
if keyword_set(wedge) then begin
   f=cu.wedge
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Wedge' " & unit$='m'
endif
if keyword_set(wedgo) then begin
   f=cu.wedgo
   lev=(findgen(16)-7.99999)*200.
   title$=" 'WedgO' (restricted) " & unit$='m'
endif
if keyword_set(dnodes) then begin
   f=cu.nodes-cu.dev
   lev=(findgen(16)-7.99999)*200.
   title$=" 'Nodes'-Dev " & unit$='m'
endif

if keyword_set(aniso) then begin
   f=cu.aniso
   lev=(findgen(16)-1)/16.
endif

if keyword_set(sdev) then begin
   f=cu.sdev
   lev=(findgen(16)-7.99999)*200.
   cell_fill=1
endif
if keyword_set(corr) then begin
   f=cu.corr
   lev=0.4+(0.6/16)*findgen(16)
   cell_fill=1
endif
if keyword_set(cuqi) then begin
   f=cu.cuqi
   lev=findgen(16)-1
endif

sm$=''
if keyword_set(sm) then begin
   f=smoothxyp(f,sm,sm)
   sm$='  - smoothed -'+padstring( sm )
endif


;Plotting ...
rex
!p.charsize=1.75

window,re=2,xs=1200,ys=900,w0
wset,w0
amwgct
contour,f(*,*,p-1),lev=lev,c_colo=indgen(16),/fill,/xst,/yst,pos=[.075,.1,.8,.9], xtit='Cell #',ytit='Cell #',xr=xr,yr=yr,cell=cell_fill

if keyword_set(mxdis) or keyword_set(super) then begin
   contour,f(*,*,p-1),lev=lev,c_colo=indgen(16),/noer,/xst,/yst,pos=[.075,.1,.8,.9],xr=xr,yr=yr,c_thick=1
   contour,f(*,*,p-1),lev=[.1,1] ,/noer,/xst,/yst,pos=[.075,.1,.8,.9],xr=xr,yr=yr,c_thick=1
endif

if keyword_set(x) then begin 
   oplot,x.xspk(ox)-1, x.yspk(ox)-1 ,ps=1,syms=.5
endif

xcolorbar,pos=[.83,.2,.85,.8],clev=lev,unit=unit$,labsz=1.5,uns=1.5
xyouts,/norm,align=.5,.4,.92,Title$+sm$,size=3.
xyouts,/norm,align=.5,.5,.02, cu.rf,size=1.3







return
end
