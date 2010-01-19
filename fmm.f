      subroutine bs2(n0,n1,xi,yi,ui,vi,n2,n3,xj,yj,gj,sj,mp,lm)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      real tim(2),tfmm(9)

      complex*16  a(mpmax,nbne),ao(mpmax,nbne), b(mpmax,nbne,lmax)
      complex*16 ss(mpmax,mpmax),sr(mpmax,mpmax),rr(mpmax,mpmax)
      complex*16 aa(mpmax),bb(mpmax),c(mpmax),d(mpmax)
      complex*16 invt,invt2,t,ziic,zjjc

      dimension xi(npmax),yi(npmax),ui(npmax),vi(npmax)
      dimension xj(npmax),yj(npmax),gj(npmax),sj(npmax)

      dimension xdi(npmax),ydi(npmax),udi(npmax),vdi(npmax)
      dimension xdj(npmax),ydj(npmax),gdj(npmax),sdj(npmax)

      dimension ndi(npmax,2),ndbi(npmax),ndj(npmax,2),ndbj(npmax)
      dimension nbi(npmax),nfi(nbmax),nei(nbmax),nfio(nbmax),neio(nbmax)
      dimension nbj(npmax),nfj(nbmax),nej(nbmax),nfjo(nbmax),nejo(nbmax)
      dimension nii(nbmax),njj(nbmax),nc(3),nfn(9)

      call cpu_time(tic)
      do i = 1,9
       tfmm(i) = 0
      end do
      pi = 2*acos(0.0)
      xmax = -1e10
      xmin = 1e10
      ymax = -1e10
      ymin = 1e10

      do i = n0,n1
       xmax = max(xi(i),xmax)
       xmin = min(xi(i),xmin)
       ymax = max(yi(i),ymax)
       ymin = min(yi(i),ymin)
      end do

      do i = n2,n3
       xmax = max(xj(i),xmax)
       xmin = min(xj(i),xmin)
       ymax = max(yj(i),ymax)
       ymin = min(yj(i),ymin)
      end do

      rdx = (xmax-xmin)*(1+1e-5)
      rdy = (ymax-ymin)*(1+1e-5)
      rd = max(rdx,rdy)
      nmax = max(n1-n0+1,n3-n2+1)

      lev = lm
      lml = 2

      call boxdata(n0,n1,xi,yi,rd,lev,nbi,ndbi,ndi,nei,nfi,
     *    lbi,rb,xmin,ymin)
      call boxdata(n2,n3,xj,yj,rd,lev,nbj,ndbj,ndj,nej,nfj,
     *    lbj,rb,xmin,ymin)

      n = 0
      do i = n0,n1
       n = n+1
       xdi(n) = xi(ndbi(n))
       ydi(n) = yi(ndbi(n))
      end do

      n = 0
      do i = n2,n3
       n = n+1
       xdj(n) = xj(ndbj(n))
       ydj(n) = yj(ndbj(n))
       gdj(n) = gj(ndbj(n))
       sdj(n) = sj(ndbj(n))
      end do

      toc = tic
      call cpu_time(tic)
      tfmm(1) = tfmm(1)+tic-toc

c Step 1. S-expansion

      do i = 1,lbj
       call boxc(nfj(i),2,nc)
       xjc = xmin+(nc(1)+.5)*rb
       yjc = ymin+(nc(2)+.5)*rb
       do j = 1,mp
        bb(j) = 0
       end do

       do j = ndj(i,1),ndj(i,2)
        xjjc = xdj(j)-xjc
        yjjc = ydj(j)-yjc
        zjjc = dcmplx(xjjc,yjjc)
        do k = 1,mp
         bb(k) = bb(k)+gdj(j)*zjjc**(k-1)
        end do
       end do

       do j = 1,mp
        b(j,i,lev) = bb(j)
       end do
      end do

      toc = tic
      call cpu_time(tic)
      tfmm(2) = tfmm(2)+tic-toc

c Step 2. S|S-translation

      if(lm.gt.lml)then
       do lev = lm-1,lml,-1
        do i = 1,4**(lev+1)
         nejo(i) = nej(i)
        end do
        do i = 1,lbj
         nfjo(i) = nfj(i)
        end do
        lbjo = lbj
        call boxdata(n2,n3,xj,yj,rd,lev,nbj,ndbj,ndj,nej,nfj,
     *    lbj,rb,xmin,ymin)
        do i = 1,lbj
         call boxc(nfj(i),2,nc)
         xjc = (nc(1)+.5)*rb
         yjc = (nc(2)+.5)*rb
         do j = 1,mp
          bb(j) = dcmplx(0.0,0.0)
         end do
         do j = 1,4
          nfjc = nfj(i)*4+j-1
          do k = 1,lbjo
           if(nfjo(k).eq.nfjc)then
            call boxc(nfjc,2,nc)
            xjco = (nc(1)+.5)*rb/2
            yjco = (nc(2)+.5)*rb/2
            t = dcmplx(xjc-xjco,yjc-yjco)

            ss(1,1) = 1
            nminus = 1
            do l = 1,mp
             c(l) = t**(l-1)
            end do
            do n = 2,mp
             nminus = -nminus
             ss(n,1) = 0
             ss(1,n) = nminus*c(n)
             mminus = nminus
             d(1) = 1
             do m = 2,n-1
              d(m) = real(n-m+1)/(m-1)*d(m-1)
              mminus = -mminus
              ss(n,m) = 0
              ss(m,n) = mminus*c(n-m+1)*d(m)
             end do
             ss(n,n) = ss(n-1,n-1)
            end do 

            do l = 1,mp
             do m = 1,mp
              bb(l) = bb(l)+ss(m,l)*b(m,nejo(nfjc+1),lev+1)
             end do 
            end do 

           end if 
          end do
         end do 
         do j = 1,mp
          b(j,i,lev) = bb(j)
         end do
        end do 
       end do 
       lev = lml
      end if
      call boxdata(n0,n1,xi,yi,rd,lev,nbi,ndbi,ndi,nei,nfi,
     *    lbi,rb,xmin,ymin)

      toc = tic
      call cpu_time(tic)
      tfmm(3) = tfmm(3)+tic-toc

c Step 3. S|R-translation

      do i = 1,lbi
       call boxc(nfi(i),2,nc)
       xic = (nc(1)+.5)*rb
       yic = (nc(2)+.5)*rb
       do j = 1,mp
        aa(j) = 0
       end do
       call e3b(nfi(i),nfj,lbj,nii,li,lev)
       do j = 1,li
        jj = nii(j)
        call boxc(jj,2,nc)
        xjc = (nc(1)+.5)*rb
        yjc = (nc(2)+.5)*rb
        t = dcmplx(xic-xjc,yic-yjc)

        invt = 1/t
        sr(1,1) = invt
        invt2 = invt*invt
        nminus = 1
        do k = 1,mp
         c(k) = invt**k
        end do
        do n = 2,mp
         nminus = -nminus
         sr(n,1) = c(n)
         sr(1,n) = nminus*c(n)
         d(1) = c(n)
         mminus = nminus
         do m = 2,n-1
          mminus = -mminus
          d(m) = -(1+real(n-1)/(m-1))*d(m-1)*invt
          sr(n,m) = d(m)
          sr(m,n) = mminus*d(m)
         end do
         sr(n,n) = -2*invt2*(2-1.0/(n-1))*sr(n-1,n-1)
        end do

        do k = 1,mp
         do l = 1,mp
          aa(k) = aa(k)+sr(l,k)*b(l,nej(jj+1),lev)
         end do 
        end do 

       end do 
       do j = 1,mp
        a(j,i) = aa(j)
       end do
      end do

      toc = tic
      call cpu_time(tic)
      tfmm(4) = tfmm(4)+tic-toc

c Step 4. R|R-translation

      if(lm.gt.lml)then
       do lev = lml+1,lm

c R|R-translation

        do i = 1,4**(lev-1)
         neio(i) = nei(i)
        end do
        do i = 1,lbi
         nfio(i) = nfi(i)
         do j = 1,mp
          ao(j,i) = a(j,i)
         end do
        end do
        lbio = lbi
        call boxdata(n0,n1,xi,yi,rd,lev,nbi,ndbi,ndi,nei,nfi,
     *    lbi,rb,xmin,ymin)
        do i = 1,lbi
         call boxc(nfi(i),2,nc)
         xic = (nc(1)+.5)*rb
         yic = (nc(2)+.5)*rb
         nfip = nfi(i)/4
         call boxc(nfip,2,nc)
         xico = (nc(1)+.5)*rb*2
         yico = (nc(2)+.5)*rb*2
         t = dcmplx(xic-xico,yic-yico)

         rr(1,1) = 1
         do j = 1,mp
          c(j) = t**(j-1)
         end do
         do m = 2,mp
          rr(m,1) = c(m)
          rr(1,m) = 0
          d(1) = 1
          do n = 2,m-1
           d(n) = real(m-n+1)/(n-1)*d(n-1)
           rr(n,m) = 0
           rr(m,n) = c(m-n+1)*d(n)
          end do 
          rr(m,m) = rr(m-1,m-1)
         end do

         do j = 1,mp
          t = 0
          do k = 1,mp
           t = t+rr(k,j)*ao(k,neio(nfip+1))
          end do
          a(j,i) = t
         end do

        end do
        call boxdata(n2,n3,xj,yj,rd,lev,nbj,ndbj,ndj,nej,nfj,
     *    lbj,rb,xmin,ymin)

        toc = tic
        call cpu_time(tic)
        tfmm(5) = tfmm(5)+tic-toc

c S|R-translation

        do i = 1,lbi
         call boxc(nfi(i),2,nc)
         xic = (nc(1)+.5)*rb
         yic = (nc(2)+.5)*rb
         do j = 1,mp
          aa(j) = a(j,i)
         end do
         call e4b(nfi(i),nfj,lbj,nii,li,lev)
         do j = 1,li
          jj = nii(j)
          call boxc(jj,2,nc)
          xjc = (nc(1)+.5)*rb
          yjc = (nc(2)+.5)*rb
          t = dcmplx(xic-xjc,yic-yjc)

          invt = 1/t
          sr(1,1) = invt
          invt2 = invt*invt
          nminus = 1
          do k = 1,mp
           c(k) = invt**k
          end do
          do n = 2,mp
           nminus = -nminus
           sr(n,1) = c(n)
           sr(1,n) = nminus*c(n)
           d(1) = c(n)
           mminus = nminus
           do m = 2,n-1
            mminus = -mminus
            d(m) = -(1+real(n-1)/(m-1))*d(m-1)*invt
            sr(n,m) = d(m)
            sr(m,n) = mminus*d(m)
           end do
           sr(n,n) = -2*invt2*(2-1.0/(n-1))*sr(n-1,n-1)
          end do

          do k = 1,mp
           do l = 1,mp
            aa(k) = aa(k)+sr(l,k)*b(l,nej(jj+1),lev)
           end do
          end do

         end do
         do j = 1,mp
          a(j,i) = aa(j)
         end do
        end do

        toc = tic
        call cpu_time(tic)
        tfmm(6) = tfmm(6)+tic-toc

       end do
       lev = lm
      end if

c Step 5. R-expansion and Final summation

      do i = 1,lbi
       call boxc(nfi(i),2,nc)
       xic = xmin+(nc(1)+.5)*rb
       yic = ymin+(nc(2)+.5)*rb
       call e2b(n2,n3,nbj,nfi(i),nfj,lbj,nfn,ln,lev)
       lj = 0
       do j = 1,ln
        jn = nej(nfn(j)+1)
        do k = ndj(jn,1),ndj(jn,2)
         lj = lj+1
         njj(lj) = k
        end do
       end do

       toc = tic
       call cpu_time(tic)
       tfmm(7) = tfmm(7)+tic-toc

       do j = ndi(i,1),ndi(i,2)
        xiic = xdi(j)-xic
        yiic = ydi(j)-yic
        ziic = dcmplx(xiic,yiic)

        t = 0
        do k = 1,mp
         t = t+a(k,i)*ziic**(k-1)
         uidu = -0.5/pi*imag(t)
         vidu = -0.5/pi*dble(t)
        end do

        toc = tic
        call cpu_time(tic)
        tfmm(8) = tfmm(8)+tic-toc

        do k = 1,lj
         kk = njj(k)
         xij = xdi(j)-xdj(kk)
         yij = ydi(j)-ydj(kk)
         rij = xij**2+yij**2
         if(rij.lt.eps)then
          rij = eps
         end if
         sij = sdj(kk)**2

         cutoff = 1/rij*(1-exp(-rij/2/sij))
c         cutoff = (rij**2+3*rij*sij+4*sij**2)/(rij+sij)**3
c         cutoff = (rij+2*sij)/(rij+sij)**2
c         cutoff = 1/(rij+sij)
         uidu = uidu+0.5/pi*gdj(kk)*yij*cutoff
         vidu = vidu-0.5/pi*gdj(kk)*xij*cutoff
        end do

        toc = tic
        call cpu_time(tic)
        tfmm(9) = tfmm(9)+tic-toc

        udi(j) = uidu
        vdi(j) = vidu
       end do
      end do
      n = 0
      do i = n0,n1
       n = n+1
       ui(ndbi(n)) = udi(n)
       vi(ndbi(n)) = vdi(n)
      end do
      print*,'init   : ',tfmm(1)
      print*,'index  : ',tfmm(7)
      print*,'p2p    : ',tfmm(9)
      print*,'p2m    : ',tfmm(2)
      print*,'m2m    : ',tfmm(3)
      print*,'m2l    : ',tfmm(4)+tfmm(6)
      print*,'l2l    : ',tfmm(5)
      print*,'l2p    : ',tfmm(8)
      print*,'------------------'
      print*,'total  : ',tfmm(1)+tfmm(2)+tfmm(3)+tfmm(4)+tfmm(5)+tfmm(6)
     *    +tfmm(7)+tfmm(8)+tfmm(9)
      return
      end

c Bookmarking box data structure

      subroutine boxdata(n0,n1,x,y,rd,lev,nb,ndb,nd,ne,nf,
     *    lb,rb,xmin,ymin)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension x(npmax),y(npmax)
      dimension nb(npmax),nbd(npmax),ndb(npmax),nd(npmax,2)
      dimension ne(nbmax),nf(nbmax),nxs(3,npmax)

      kl = 2**lev
      rb = rd/kl
      n = 0
      do i = n0,n1
       n = n+1
       nxs(1,n) = (x(i)-xmin)/rb
       nxs(2,n) = (y(i)-ymin)/rb
      end do

      call boxn(n,nxs,2,lev,nb)
      do i = 1,n
       nbd(i) = nb(i)
       ndb(i) = i+n0-1
      end do
      call sort(n,nbd,ndb)

      lb = 0
      nbo = -1
      do i = 1,n
       if(nbd(i).ne.nbo)then
        lb = lb+1
        ne(nbd(i)+1) = lb
        nf(lb) = nbd(i)
        nd(lb,1) = i
        if(lb.ge.2)then
         nd(lb-1,2) = i-1
        end if
        nbo = nbd(i)
       end if
      end do
      nd(lb,2) = n
      return
      end

c Retrive index of E3 boxes

      subroutine e3b(nf,nff,lb,nii,li,lev)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension nff(nbmax),nii(nbmax),nc(3),nfn(9)
      dimension nxs(3,npmax)

      kl = 2**lev
      call boxc(nf,2,nc)
      call neighbor(nc,kl,nxs,n)
      call boxn(n,nxs,2,lev,nfn)
      li = 0
      do i = 1,lb
       jj = 0
       do j = 1,n
        if(nff(i).eq.nfn(j))then
         jj = jj+1
        end if
       end do
       if(jj.eq.0)then
        li = li+1
        nii(li) = nff(i)
       end if
      end do
      return
      end

c Retrive index of E4 boxes

      subroutine e4b(nf,nff,lb,nii,li,lev)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension nff(nbmax),nii(nbmax),nc(3),nfn(9),nfpn(9)
      dimension nxs(3,npmax)

      kl = 2**lev
      nfp = nf/4
      call boxc(nfp,2,nc)
      call neighbor(nc,2*kl,nxs,np)
      call boxn(np,nxs,2,lev+1,nfpn)
      call boxc(nf,2,nc)
      call neighbor(nc,kl,nxs,n)
      call boxn(n,nxs,2,lev,nfn)
      li = 0
      do i = 1,np
       do j = 1,4
        nfnpn = nfpn(i)*4+j-1
        ii = 0
        do k = 1,n
         if(nfn(k).eq.nfnpn)then
          ii = ii+1
         end if
        end do
        if(ii.eq.0)then
         do k = 1,lb
          if(nff(k).eq.nfnpn)then
           li = li+1
           nii(li) = nfnpn
          end if
         end do
        end if
       end do
      end do
      return
      end

c Retrive index of E2 boxes

      subroutine e2b(n0,n1,nb,nf,nff,lb,nfn,n,lev)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension nb(npmax),nff(nbmax),nc(3),nfn(9),nfnn(9)
      dimension nxs(3,npmax)

      kl = 2**lev
      call boxc(nf,2,nc)
      call neighbor(nc,kl,nxs,nn)
      call boxn(nn,nxs,2,lev,nfnn)
      n = 0
      do i = 1,nn
       do j = 1,lb
        if(nff(j).eq.nfnn(i))then
         n = n+1
         nfn(n) = nfnn(i)
        end if
       end do
      end do
      return
      end

c Retrive box index from coordinate numbers

      subroutine boxn(n,nxs,ndi,lev,nb)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension nb(npmax),nxs(3,npmax)

      do i = 1,n
       nb(i) = 0
      end do
      do i = 1,lev
       do j = 1,n
        nbx = mod(nxs(1,j),2)
        nxs(1,j) = nxs(1,j)/2
        nb(j) = nb(j)+nbx*2**(ndi*(i-1)+1)
       end do

       if(ndi.ge.2)then
        do j = 1,n
         nby = mod(nxs(2,j),2)
         nxs(2,j) = nxs(2,j)/2
         nb(j) = nb(j)+nby*2**(ndi*(i-1))
        end do
       end if

       if(ndi.eq.3)then
        do j = 1,n 
         nbz = mod(nxs(3,j),2)
         nxs(3,j) = nxs(3,j)/2
         nb(j) = nb(j)+nbz*2**(ndi*(i-1)+2)
        end do
       end if
      end do
      return
      end

c Retrive coordinate numbers from box index

      subroutine boxc(nf,ndi,nc)
      implicit real*8(a-h,o-z)

      dimension nc(3)

      do i = 1,3
       nc(i) = 0
      end do
      nb = nf
      k = 0
      i = 0
      do while(nb.ne.0)
       j = ndi-k
       nc(j) = nc(j)+mod(nb,2)*2**i
       nb = nb/2
       k = mod(k+1,ndi)
       if(k.eq.0)then
        i = i+1
       end if
      end do

      if(ndi.eq.3)then
       du = nc(1)
       nc(1) = nc(2)
       nc(2) = nc(3)
       nc(3) = du
      end if
      return
      end

c Retrive coordinate numbers of neighbor boxes

      subroutine neighbor(nc,kl,nxs,n)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension nc(3),nxs(3,npmax)

      if(nc(1).eq.0)then
       nx = 2
       nxi = 0
      elseif(nc(1).eq.kl-1)then
       nx = 2
       nxi = kl-2
      else
       nx = 3
       nxi = nc(1)-1
      end if

      if(nc(2).eq.0)then
       ny = 2
       nyi = 0
      elseif(nc(2).eq.kl-1)then
       ny = 2
       nyi = kl-2
      else
       ny = 3
       nyi = nc(2)-1
      end if

      n = 0
      do i = 1,nx
       do j = 1,ny
        n = n+1
        nxs(1,n) = nxi+i-1
        nxs(2,n) = nyi+j-1
       end do
      end do
      return
      end

c Quicksort

      subroutine sort(n,na,nb)
      implicit real*8(a-h,o-z)

      include 'memory.f'

      dimension na(npmax),nb(npmax),istack(50)

      jstack = 0
      m = 7
      l = 1
      ir = n
    1 if(ir-l.lt.m)then
       do j = l+1,ir
        naa = na(j)
        nbb = nb(j)
        do i = j-1,1,-1
         if(na(i).le.naa) goto 2
         na(i+1) = na(i)
         nb(i+1) = nb(i)
        end do
        i = 0
    2   na(i+1) = naa
        nb(i+1) = nbb
       end do
       if(jstack.eq.0) return
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack-2
      else
       k = (l+ir)/2
       nd = na(k)
       na(k) = na(l+1)
       na(l+1) = nd
       nd = nb(k)
       nb(k) = nb(l+1)
       nb(l+1) = nd

       if(na(l+1).gt.na(ir))then
        nd = na(l+1)
        na(l+1) = na(ir)
        na(ir) = nd
        nd = nb(l+1)
        nb(l+1) = nb(ir)
        nb(ir) = nd
       end if

       if(na(l).gt.na(ir))then
        nd = na(l)
        na(l) = na(ir)
        na(ir) = nd
        nd = nb(l)
        nb(l) = nb(ir)
        nb(ir) = nd
       end if

       if(na(l+1).gt.na(l))then
        nd = na(l+1)
        na(l+1) = na(l)
        na(l) = nd
        nd = nb(l+1)
        nb(l+1) = nb(l)
        nb(l) = nd
       end if

       i = l+1 
       j = ir
       naa = na(l)
       nbb = nb(l)
    3  continue
       i = i+1
       if(na(i).lt.naa) goto 3
    4  continue
       j = j-1
       if(na(j).gt.naa) goto 4
       if(j.lt.i) goto 5
       nd = na(i)
       na(i) = na(j)
       na(j) = nd
       nd = nb(i)
       nb(i) = nb(j)
       nb(j) = nd
       goto 3
    5  na(l) = na(j)
       na(j) = naa
       nb(l) = nb(j)
       nb(j) = nbb
       jstack = jstack+2
       if(jstack.gt.50) pause 'NSTACK too small in sort'
       if(ir-i+1.ge.j-1)then
        istack(jstack) = ir
        istack(jstack-1) = i
        ir = j-1
       else
        istack(jstack) = j-1
        istack(jstack-1) = l
        l = i
       end if
      end if 
      goto 1
      return
      end
