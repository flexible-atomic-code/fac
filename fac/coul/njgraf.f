c*****************************************************************
c===============================================  n j g r a f  =====
c*****************************************************************
c*****************************************************************
c
c   section 4          njgraf 
c
c     ******************************
      subroutine njgraf(recup,fail)
      implicit real*8(a-h,o-z)
c     ******************************
c
c
c  ***this is the main program.it handles  all the analysis of the
c  ***recoupling coefficient without referring explicitly to the values
c  ***of angular momenta which are in j1(j),except for zero in case free
c  ***=.false. .like njsym it prepares arrays of arguments for phase
c  ***factors,(2*j+1) factors and 6j coefficients to be computed in
c  ***gensum,which can also be called separately when only the numerical
c  ***values of angular momenta change.these variable angular momenta 
c  ***should be declared free(j)=.true.,so that the formula prepared for
c  ***gensum should be correct when  j1 is not zero.
c  ***fail will be true when the recoupling coefficient is zero because
c  ***of unsatisfied delta or other similar causes.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
      parameter(mfact=100)
c
      logical fail,find,tabs,cut,free,sumvar
c
      integer arrow,arr,tab1
c
      character*6 name,namsub 
c
      common/nam/namsub
c
      common/const/i6c,i7c,i8c,i9c,idel,iwc
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
      common/cutdig/cut
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),
     + free(mangm)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
c     common block /facts / used to store ln(i!) in the new racah
c     routine written by stan scott.
c
      common /facts / gam(mfact)
c
      common/dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc
c
      common/sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),
     + jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),
     + k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),
     + jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)
c
c
      data name/'njgraf'/
c
c
c ***building up of the unstructured graph
c
      fail=.false.
      j6c=0
      j7c=0
      j8c=0
      j9c=0
      jwc=0
      jdel=0
      i6c=1
      i7c=1
      i8c=1
      i9c=1
      idel=1
      iwc=1
      call setdim
      nfin=0
      cut=.false.
      call settab(fail)
      m=m+1
      if(fail)go to 7
      m=m-1
      jf=0
      jf1=0
c
c   ***locating and handling of zeros
c
      call zero(jf1,jf,fail)
      if(fail)go to 7
      mp=m
      if(nbtr.eq.0)go to 6
      jump=1
c
c   ***building of a flat diagram out of the unstructured graph.
c   ***there may be several flat diagrams out of the original
c   ***graph,in case of possible cuts.then the flat diagrams
c   ***will have free ends.
c
    1 call diagrm(jump)
      nfin=max0(0,nfree-2)
c
      if(nfin.ne.0) then
        jump=3
c
c  ***handling of free ends if a cut was found
c
        call cutnl(fail)
        if(fail)go to 7
      else
       jump=2
       if(nfree .eq. 1) then
          call cut1l(fail)
          if(fail)go to 7
       else
         if (nfree .gt. 1) then
           call cut2l(fail)
           if(fail)go to 7
         endif
       endif
      endif
c
      nbtr=nbtr+nfin
      if(nbtr.ne.0)cut=.true. 
c
c  ***analysis of the flat diagram.
c  ***closed circuits of increasing order nc are searched,analysed,and
c  ***taken out of the flat diagram,thus reducing the number of nodes,
c  ***nbnode.
  
c
      nc=0
   10 nc=nc+1
      call search(find)
      if(.not.find)go to 10
      ncp=nc-2
      jpol=0
      if(m.eq.mp.and.nc.gt.3)call setdim
      if(ipartl.gt.2)call polygn(jpol)
      go to (11,12,13,14),nc
   11 call lolpop(fail)
      if(fail)go to 7
      go to 15
   12 call bubble(jpol,fail)
      if(fail)go to 7
      go to 15
   13 call triang(fail)
      if(fail)go to 7
      go to 15
   14 call square
   15 nbnode=nbnode-2
      if(nbnode.eq.0)go to 9
      ifirst=ih(1)
      ilast=ih(nbnode)
c
c **printj is an all purpose printing subroutine called from many places
c
      call printj(namsub,8)
      if(nbnode.eq.nfin)go to 9
      nc=ncp
c
c  ***proceed to other circuits of order nc-1
c
      go to 10
    9 if(nbtr.eq.0)go to 6
      if(jump.eq.3)call ordtri 
c
c ***at this stage,the flat diagram has been reduced to nodes
c ***involving free ends.proceed to build other flat diagrams
c ***if necessary.
c
      go to 1
c
c  ***all parts of the original graph have been reduced.
c
    7 recup=0.
      m=m-1
      return
    6 call printj(name,0)
c
c
c  ***preparation of the results,and separation in several sums
c  *** if cuts have been detected,also in the flat diagram itself
c
      call sprate(m)
      m=m-1
c
c
c  ***gensum computes the numerical value of the recoupling 
c  ***coefficient.
c
      call gensum(recup)
c
      return
      end 
c
c
c     ****************************
      subroutine bubble(jpol,fail)
      implicit real*8(a-h,o-z)
c     ****************************
c
c  ***reduces a circuit of order 2,giving delta function and phase
c  ***factors.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,sumvar
c
      integer arr,tab1
c
      character*6 name,namsub 
c
      common/nam/namsub
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'bubble'/
c
      namsub=name
      k2=2
      k23=3
      i1=1
      i2=1
      it1=npoint(1) 
      it2=npoint(2) 
c
      if(it2.eq.ilast) then
        if(it1.ne.ifirst) then
          it2=it1
          it1=ilast 
        endif
          i1=-1
          k23=2
          i2=2
      endif
c
      call phase(it1,jdiag,m4trd)
      k=iabs((3*arr(it2,1)+2*arr(it2,2)+arr(it2,3))/2)+1
      if(k.ne.4)call phase2(jdiag(it2,k))
      if(nbnode.eq.2)return
      il1=il(it2)+i1
      it=ih(il1)
      arr(it,k23)=arr(it1,k23)
      l=jdiag(it1,k23)
      l1=jdiag(it,k23)
      jdiag(it,k23)=l
c
c
      if(jpol.ne.1) then
        call delta(l,l1,fail) 
        if(fail)return
      else
        mp=mp-1
        kw(2,jwc)=l 
        j6(j6c-1)=l 
        j6(j6c)=l
        if(k.eq.2)j8(j8c)=l
      endif
c
      tab1(l,i2)=it 
c
      if(it1.ne.ilast)then
        if(it2.eq.ilast) then 
          tab1(l,1)=ih(2)
          il1=2
          k2=1
        endif
c
      do 5 i=il1,nbnode
        it=ih(i)
        il(it)=i-k2 
        ih(i-k2)=it 
    5 continue
c
      endif
c
    6 j9(j9c+1)=l
      j9c=j9c+2
      j9(j9c)=l
c
      return
      end 
c
c
c     **********************
      subroutine change(l,k)
      implicit real*8(a-h,o-z)
c     **********************
c
c     exchanges the free ends in either first or last triad of jdiag. 
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm)
c
      integer arr,tab1
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      call phase(l,jdiag,m4trd)
      jp=jdiag(l,k) 
      jdiag(l,k)=jdiag(l,1)
      jdiag(l,1)=jp 
      jar=arr(l,k)
      arr(l,k)=arr(l,1)
      arr(l,1)=jar
c
      return
      end 
c
c
c     *****************************************
      subroutine chvar(jp,nbc,kbc,jt,jinv,nsum)
      implicit real*8(a-h,o-z)
c     *****************************************
c
c
c   ***change the order of summation variable to be able to perform
c   ***separately the summations in gensum.
c
      logical jt(nsum)
      dimension jp(nbc),jinv(nsum)
c
c
      kb=kbc+1
c
      if(kb.le.nbc) then
        do 1 i=kb,nbc
          jk=jp(i)
          if(jt(jk))then
            kbc=kbc+1
            jp(i)=jp(kbc)
            jp(kbc)=jinv(jk)
          endif
    1   continue
      endif
c
      return
      end 
c
c
c     **********************
      subroutine cut1l(fail)
      implicit real*8(a-h,o-z)
c     **********************
c
c  ***cut on one line,that was left as a free end in jdiag.puts
c  ***corresponding delta in j23.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,sumvar,free
c
      integer arr,tab1
c
      character*6 name
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
      data name/'cut1l '/
c
c
      it=itfree(1)
      j0=jdiag(it,1)
      call delta (j0,m,fail)
      if(fail)go to 2
      call delta(jdiag(it,3),jdiag(it,2),fail)
      if(fail)go to 2
      jdiag(it+1,3)=jdiag(it,3)
c
      if(arr(it,2) .eq. arr(it,3))then
        arr(it+1,3)=1
        arr(it-1,2)=-1
      else
      if (arr(it,2) .lt. arr(it,3))then 
        arr(it+1,3)=-1
        arr(it-1,2)=1
      endif
      endif
c
      j9c=j9c+1
      j9(j9c)=jdiag(it,3)
      j=2 
      call zero(j,j0,fail)
      if(fail)go to 2
      il1=il(it+1)
c
      do 1 i=il1,nbnode
        it=ih(i)
        ilp=i-1
        il(it)=ilp
        ih(ilp)=it
    1 continue
c
      nbnode=nbnode-1
c
    2 call printj(name,mtriad)
      return
      end 
c
c
c     **********************
      subroutine cut2l(fail)
      implicit real*8(a-h,o-z)
c     **********************
c
c  ***cut on two lines that were left as free ends in jdiag.puts
c  ***corresponding delta in j23.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,tabs,sumvar
c
      integer arr,tab1,arrow
c
      character*6 name
c
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'cut2l '/
c
c
      it1=itfree(1) 
      it2=itfree(2) 
      jt1=jdiag(it1,1)
      jt2=jdiag(it2,1)
      call delta(jt1,jt2,fail)
      if(fail) go to 1
      if(arr(it1,1).eq.arr(it2,1))call phase2(jt1)
      arr(it2,1)=-arr(it1,1)
      jdiag(it2,1)=jt1
      tab1(jt1,2)=it2
      j9(j9c+1)=jt1 
      j9c=j9c+2
      j9(j9c)=jt1
      call otherj(0,jt1,l1,lc1,k1)
      call otherj(0,jt2,l2,lc2,k2)
      j23(l2,lc2)=jt1
      line(jt1,k1)=l2
      lcol(jt1,k1)=lc2
      arrow(l2,lc2)=-arrow(l1,lc1)
c
    1 call printj(name,mtriad)
c
      return
      end 
c
c
c     **********************
      subroutine cutnl(fail)
      implicit real*8(a-h,o-z)
c     **********************
c
c  ***this subroutine  examines the case where there are more than
c  ***two free ends,but they are contiguous,so that the graph can
c  ***be cut without destroying the flat structure.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      integer arrow,arr,tab1
c
      logical tabs,sumvar,fail
c
      character*6 name
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/keep/jkp(2,3),jarr(2,3),it2,it3,it5
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'cutnl '/
c
c
      ntf=itfree(nfree)-itfree(1)
      if(ntf.gt.nfree)go to 8 
      it2=itfree(1) 
      it3=itfree(nfree)
      it1=it2-1
      it4=it3+1
c
      if(ntf.ne.nfree)then
        jt=jdiag(it2,3)
        call delta(jt,jdiag(it3,2),fail)
c
        if(fail)go to 8
c
        if(arr(it2,3).eq.arr(it3,2))then
          call phase2(jt)
          arr(it2,3)=-arr(it2,3)
          arr(it1,2)=-arr(it1,2)
        endif
c
        jdiag(it3,2)=jt
        jdiag(it4,3)=jt
        j9(j9c+1)=jt
        j9c=j9c+2
        j9(j9c)=jt
        nbtr=nbtr+nfree
        it5=0
       go to 6
      endif
      nfr=0
c
      do 3 it5=it2,it3
        nfr=nfr+1
        if(itfree(nfr).gt.it5)go to 4
    3 continue
c
    4 jkp(1,1)=jdiag(it5,1)
      jarr(1,1)=-arr(it5,1)
      jkp(1,2)=jdiag(it2,3)
      jarr(1,2)=-arr(it2,3)
      jkp(1,3)=jdiag(it3,2)
      jarr(1,3)=-arr(it3,2)
c
      do 5 j=1,3
        jkp(2,j)=jdiag(it5,j) 
        jarr(2,j)=arr(it5,j)
    5 continue
c
      jdiag(it5,2)=jdiag(it3,2)
      arr(it5,2)=arr(it3,2)
      jdiag(it5,3)=jdiag(it2,3)
      arr(it5,3)=arr(it2,3)
      ilp=il(it2)
      il(it5)=ilp
      ih(ilp)=it5
      nbtr=nbtr+nfree+2
      call phase(it5,jdiag,m4trd)
      k=iabs((3*arr(it5,1)+2*arr(it5,2)+arr(it5,3))/2+1)
      if(k.ne.4) call phase2(jdiag(it5,k))
    6 il1=il(it4)
c
      do 7 i=il1,nbnode
        it=ih(i)
        ilp=i-nfree 
        il(it)=ilp
        ih(ilp)=it
    7 continue
c
      nbnode=nbnode-nfree
      nfin=0
c
    8 call printj(name,8)
c
      return
      end 
c
c
c     ****************************
      subroutine delta(ja,jb,fail)
      implicit real*8(a-h,o-z)
c     ****************************
c
c  ***test for delta(ja,jb).if they are summation variables,the second
c  ***is changed into the first everywhere.if they are fixed,their
c  ***value is checked,and fail put to .true. if they differ.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,cut,sumvar,free
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/cutdig/cut
      common/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6
      common/dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
 1000 format(/2x,'from delta',2x,'ja=',i2,l2,5x,'jb=',i2,l2)
c
c
c
      if(ibug3.eq.1)print 1000,ja,sumvar(ja),jb,sumvar(jb)
      if(sumvar(ja).and.sumvar(jb))go to 2
      if(free(ja).or.free(jb)) then
        jdel=jdel+1 
        ldel(jdel,1)=ja
        ldel(jdel,2)=jb
        sumvar(ja)=.false.
        sumvar(jb)=.false.
        return
      endif
c
      if(j1(ja).ne.j1(jb))fail=.true.
      cut=.true.
      return
c
    2 if(j6c.ne.j6cc) then
        j61=j6cc+1
c
        do 3 i=j61,j6c
          if(j6(i).eq.jb)j6(i)=ja
    3   continue
c
      endif
c
      if(j7c.ne.j7cc) then
        j71=j7cc+1
c
        do 5 i=j71,j7c
          if(j7(i).eq.jb)j7(i)=ja
    5   continue
      endif
c
      if(j8c.ne.j8cc) then
        j81=j8cc+1
c
        do 7 i=j81,j8c
          if(j8(i).eq.jb)j8(i)=ja
    7   continue
      endif
c
      if(j9c.ne.j9cc) then
        j91=j9cc+1
c
        do 9 i=j91,j9c
          if(j9(i).eq.jb)j9(i)=ja
    9   continue
      endif
c
      if(jwc.ne.jwcc) then
       jw1=jwcc+1
c
        do 14 i=jw1,jwc
         do 13 j=1,6
           if(kw(j,i).eq.jb)kw(j,i)=ja
   13    continue
   14   continue
      endif
c
      if(jdel.ne.jdelc) then
        jdel1=jdelc+1
c
        do 17 i=jdel1,jdel
          do 16 j=1,2
            if(ldel(i,j).eq.jb)ldel(i,j)=ja
   16     continue
   17   continue
c
        sumvar(jb)=.false.
      endif
c
      return
      end 
c
c
c     *********************** 
      subroutine diagrm(jump) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c     this subroutine builds up a flat diagram from the triads j23 and
c  ***places them in jdiag.arrows are in arr (integer).the diagram is 
c  ***built so as to maximize the number of triads involved,within  a 
c  ***one-step-forward-check process.if the diagram does not
c  ***include all the nbtr triads,it will have 'free ends'.jdiag has
c  ***dimension double that of j23,because the path may proceed either
c  ***way.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical tabs,free,sumvar
c
      integer arr,tab1,arrow
c
      character*6 name
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/build/ial(m4trd),if1,if2,node
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
      data name/'diagrm'/
c
c
c
c  ***initialization
c
      if(jump .gt. 2) go to 17
      if(jump .lt. 2) nb=0
    1 nb=nb+1
      if(tabs(nb))go to 1
      node=nbtr
      ilast=nbtr
c
      do 2 j=1,3
        jdiag(node,j)=j23(nb,j)
        arr(node,j)=arrow(nb,j)
    2 continue
c
      tabs(nb)=.true.
c
      do 15 i=1,mp
        ial(i)=0
   15 continue
c
      if1=jdiag(node,1)
      if2=jdiag(node,3)
      ial(if1)=1
      ial(if2)=1
   17 ntime=0
      i1=1
      k1=1
      k2=2
      k3=3
    3 jb=jdiag(node,k2)
      call otherj(0,jb,l,lc,kp)
      call neibor(lc,l1,l2)
c  november 22 1989 check consistency of triads
      if(tabs(l))stop' building diagram impossible '
      call way(l,l1,l2,ich,nd)
      node=node+i1
      tabs(l)=.true.
      jdiag(node,k3)=j23(l,lc)
      arr(node,k3)=arrow(l,lc)
      ict=ich*i1
c
      if (ich .le. 0) then
        lp=l1
        l1=l2
        l2=lp
      endif
c
      if (ict .le. 0) call phase(l,j23,m2trd)
      jdiag(node,k1)=j23(l,l1)
      arr(node,k1)=arrow(l,l1)
      jdiag(node,k2)=j23(l,l2)
      arr(node,k2)=arrow(l,l2)
      j=j23(l,l1)
      ial(j)=ial(j)+1
      j=j23(l,l2)
      ial(j)=ial(j)+1
      if(nd.lt.1)go to 3
      ntime=ntime+1 
      ilast=max0(node,ilast)
      ifirst=min0(node,nbtr)
      nbp=ial(if1)+ial(if2)
      if (nbp .gt. 3 .or. ntime .gt. 1) then
        nbnode=ilast-ifirst+1 
        nbtr=nbtr-nbnode
c
c  ***definition of free ends and other quantities.
c
        call intab
        call printj(name,mtriad)
        go to 50
      endif
c
      if (nbp .gt. 2) then
        if (ial(if1) .le. ial(if2)) then
          jt=jdiag(nbtr,1)
          jar=arr(nbtr,1)
          jdiag(nbtr,1)=jdiag(nbtr,3)
          arr(nbtr,1)=arr(nbtr,3)
          jdiag(nbtr,3)=jt
          arr(nbtr,3)=jar
          call phase(nbtr,jdiag,m4trd)
      endif
      endif
c
      node=nbtr
      i1=-1
      k2=3
      k3=2
      go to 3
c
   50 return
      end 
c
c
c     **********************
      subroutine dracah(rac)
      implicit real*8(a-h,o-z)
c     **********************
c
c
c  ***subroutine to calculate racah coefficients
c  ***the arguments i,j,k,l,m,n should be twice their actual value
c  ***works for integer and half-integer values of angular momenta.
c  ***the routine makes use of the gam array, thus subroutine factt
c  ***must be called before this routine is used. 
c  ***written by n s scott.
c
      parameter(mfact=100)
c
      common /facts / gam(mfact)
      common/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6
      common/jrac/i,j,k,l,m,n
c
      data zero,one,two/0.D0,1.D0,2.D0/
c
c
      j1=i+j+m
      j2=k+l+m
      j3=i+k+n
      j4=j+l+n
      if((2*max0(i,j,m)-j1).gt.0.or.mod(j1,2).ne.0) go to 2 
      if((2*max0(k,l,m)-j2).gt.0.or.mod(j2,2).ne.0) go to 2 
      if((2*max0(i,k,n)-j3).gt.0.or.mod(j3,2).ne.0) go to 2 
      if((2*max0(j,l,n)-j4).gt.0.or.mod(j4,2).ne.0) go to 2 
      go to 1
   2  rac=zero
      return
c
   1  continue
      j1=j1/2
      j2=j2/2
      j3=j3/2
      j4=j4/2
      j5=(i+j+k+l)/2
      j6=(i+l+m+n)/2
      j7=(j+k+m+n)/2
      numin=max0(j1,j2,j3,j4)+1
      numax=min0(j5,j6,j7)+1
      rac=one
c
      if(numin.eq.numax)go to 4
      numin=numin+1 
c
      do 3 ki=numax,numin,-1
        xnom=ki*(j5-ki+2)*(j6-ki+2)*(j7-ki+2)
        xdnom=(ki-1-j1)*(ki-1-j2)*(ki-1-j3)*(ki-1-j4)
        rac=one-rac*xnom/xdnom
   3  continue
c
      numin=numin-1 
   4  rac=rac*((-one  )**(j5+numin+1))*exp((gam(numin+1)-gam(numin-j1)
     * -gam(numin  -j2)-gam(numin  -j3)-gam(numin  -j4)-gam(j5+2-numin)
     * -gam(j6+2-numin)-gam(j7+2-numin))+((gam(j1+1-i)+gam(j1+1-j)
     * +gam(j1+1-m)-gam(j1+2)+gam(j2+1-k)+gam(j2+1-l)+gam(j2+1-m)
     * -gam(j2+2)+gam(j3+1-i)+gam(j3+1-k)+gam(j3+1-n)-gam(j3+2)
     * +gam(j4+1-j)+gam(j4+1-l)+gam(j4+1-n)-gam(j4+2))/two  ))
c
  6   return
      end 
c
c
c     ****************
      subroutine factt
      implicit real*8(a-h,o-z)
c     ****************
c
c
c  ***calculates the logs of factorials required by the racah
c  ***coefficient routine dracah
c  *** written by n.s. scott
c
      parameter(mfact=100)
      common /facts / gam(mfact)
      data thirty,one,two/30.D0,1.D0,2.D0/
c
c
      gam(1)=one
      gam(2)=one
      x=two
c
      do 10 i=3,30
        gam(i)=gam(i-1)*x
        x=x+one
   10 continue
c
      do 20 i=1,30
        gam(i)=log(gam(i))
  20  continue
c
      x=thirty
c
      do 30 i=31,mfact
        gam(i)=gam(i-1)+log(x)
        x=x+one
  30  continue
c
      return
      end 
c
c     ************************
      subroutine gensum(recup)
      implicit real*8(a-h,o-z)
c     ************************
c
c
c
c  ***carries out the summation over coefficients defined by the arrays
c  ***the arrays j6,j7,j8,ldel and jw to give recup
c  ***the entry is either made from njgraf or directly assuming that the
c  ***arrays j6,...,jw have already been determined by a previous
c  ***entry to njgraf and that the summation is required for another set
c  ***of j values defined by the array j1
c
c  ***recup is the recoupling coefficient
c
c  ***subroutine called: dracah
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
      parameter(mfact=100)
c
      logical ldiag,noel,free,sumvar
c
      common/graph/j12(4,mtriad,mtriad)
      common/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6
c
      common/inform/ ird , ipd ,
     +               imene , iiene , icene ,
     +               imradd, iiradd, icradd,
     +               imcolx, iicolx, iccolx,
     +               imradr, iiradr, icradr,
     +               imcoli, iicoli, iccoli,
     +               imaug , iiaug , icaug
c
      common /facts / gam(mfact)
      common/jrac/ist(6)
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
      common/sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),
     + jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),
     + k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),
     + jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),jw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      dimension mat(mtriad,mtriad),jmnp(5),jmxp(5),noel(mtriad),
     + maxlp(mtriad),jsum2(mtriad),jsum3(mtriad),
     + jsum(2,m6j),jwtest(m6j),wstor(m6j),ipair(2,2),ldiag(mtriad)
c
      dimension xj1(mangm)
c
      dimension jlow(2),jhig(2)
c
      logical already
c
      data zero,one,two /0.D0,1.D0,2.D0/
      data epsil/1.D-10/
      data mxcsvr/4/
c
c  ***format statements used in gensum
c
  302 format('   sum nr.',i3)
  303 format(' no summation. recoupling coefficient=',g15.8)
  304 format(' recoupling coefficient=',g15.8)
  305 format(6f5.1,10x,g15.8)
  315 format(i6,2x,6f5.1,10x,g15.8)
  306 format(' number of independent sums:',i3)
  307 format(' sum nr.',i2,' sum value=',g15.8,' recup=',g15.8)
  308 format(' fail in gensum at 310')
  309 format('   not involving summation variable')
  400 format(//' print out from subroutine gensum'//' values of angular
     + momenta in *real* format'/,(14f5.1))
  401 format(/' racah w functions(6j)'/' arguments in *real* format'
     +, 18x,'value')
c
c
c
c   ***evaluates all terms in j6,j7,j8,j9,ldel,and jw which do not invol
c   ***a summation.the result is stored in recup and iastor
c
      if(ibug3.eq.1) then
c
        do 139 i=1,m
          xj1(i)=(j1(i)-one)/two
  139   continue
c
        write(ipd ,400) (xj1(i),i=1,m)
        write(ipd ,306) nlsum
        write(ipd ,401)
      endif
c
      mm=m+1
      j1(mm)=1
c
c  ***test delta functions
c
      j1(mm)=1
      if(jdel.le.0) go to 180
c
      do 141 i=1,jdel
        i1=ldel(i,1)
        i2=ldel(i,2)
        if(i1.gt.mm.or.i2.gt.mm)then
          if(i1.gt.mm)j1(i1)=j1(i2)
          if(i2.gt.mm)j1(i2)=j1(i1)
        else
        if(j1(i1).ne.j1(i2)) then
          recup=zero
          return
        endif
        endif
  141 continue
c
  180 recup=one
      if(jwc.ne.0)then
c
c  ***multiply recup by all racah coefficients which do not involve a
c  ***summation
c
      if(ibug3.eq.1) write(ipd ,309)
c
      do 7 i=1,jwc
        if(inv6j(i).gt.0)go to 7
       do 3 j=1,6
         i1=jw(j,i)
         ist(j) = j1(i1) - 1
    3  continue
c
       call dracah(x1)
       if(ibug3.eq.1) write(ipd ,305) (xj1(jw(k,i)),k=1,6),x1
       recup = recup*x1
c
    7 continue
c
      endif
c
      sqr=1.0
c
      if(j6c.ne.0) then
        do 12 i=1,j6c
          i1=j6(i)
          sqr=sqr*j1(i1)
   12   continue
      endif
c
      spr=1.0
c
      if(j9c.ne.0) then
        do 144 i=1,j9c
          i1=j9(i)
          spr=spr*j1(i1)
  144   continue
      endif
c
      recup=recup*sqrt(sqr/spr)
      if(abs(recup).lt.epsil)go to 145
      iastor = 0
c
      if(j7c.ne.0) then
        do 17 i=1,j7c
          i1=j7(i)
          iastor = iastor + j1(i1) -1
   17   continue
      endif
c
      if(j8c.ne.0) then
        do 22 i=1,j8c
          i1=j8(i)
          iastor = iastor +2*(j1(i1)-1)
   22   continue
      endif
c
      if(nlsum.le.0) then
        iastor=iastor/2
c
c  ***no summation involved.end of computation
c
        stor1=one
        stor=one
        if(mod(iastor,2).eq.1)recup=-recup
        if(ibug3.eq.1) write(ipd ,303) recup
        return
c
      endif
c
c
c  ***evaluation of the part involving summations.
c
c
      nfs=0
      jwr=0
      j6f=0
      j7f=0
      j8f=0
      j9f=0
      nps=0
   25 nps=nps+1
      if(ibug3.eq.1)write(ipd ,302) nps
c
c
c   *** loop on the disconnected summations
c
c
      ias=0
      nsum=nbj(nps)-nfs
      jwrd=nb6j(nps)-jwr
      j6cp=k6cp(nps)
      j7cp=k7cp(nps)
      j8cp=k8cp(nps)
      j9cp=k9cp(nps)
c
c
c     ***the range of values of each summation variable is
c     ***defined by establishing a matrix of the links between
c     ***variables.mat(i,j) contains:
c        i=j    number of possible values of i due to triangular
c               relations with non-variables,i.e. constants.
c       i.gt.j  number of links between i and j through constants
c       i.lt.j  value of the constant,if the above is 1.if not,
c               these values are srored in j12(l,i,j) where there
c               is room for mxcsvr such values (l.le.4)
c
c
      do 52 i=1,nsum
       do 152 j=1,nsum
         mat(i,j)=0
  152  continue
   52 continue
c
      do 66 i1=1,nsum
        i1t=i1+nfs
        i2=jsum6(i1t)
       do 65 i3=1,i2
         i=jsum5(i1t,i3)
         j=jsum4(i1t,i3)
         go to (54,55,56,57,58,59),j
c
c  ***the rows of the ipair arrays give limits of summation imposed
c
c
   54    ipair(1,1) = jword(2,i)
         ipair(1,2) = jword(5,i)
         ipair(2,1) = jword(3,i)
         ipair(2,2) = jword(6,i)
         go to 60
c
   55    ipair(1,1) = jword(1,i)
         ipair(1,2) = jword(5,i)
         ipair(2,1) = jword(4,i)
         ipair(2,2) = jword(6,i)
         go to 60
c
   56    ipair(1,1) = jword(1,i)
         ipair(1,2) = jword(6,i)
         ipair(2,1) = jword(4,i)
         ipair(2,2) = jword(5,i)
         go to 60
c
   57    ipair(1,1) = jword(2,i)
         ipair(1,2) = jword(6,i)
         ipair(2,1) = jword(3,i)
         ipair(2,2) = jword(5,i)
         go to 60
c
   58    ipair(1,1)= jword(1,i)
         ipair(1,2) = jword(2,i)
         ipair(2,1) = jword(3,i)
         ipair(2,2) = jword(4,i)
         go to 60
c
   59    ipair(1,1) = jword(1,i)
         ipair(1,2) = jword(3,i)
         ipair(2,1) = jword(2,i)
         ipair(2,2) = jword(4,i)
c
   60    do 63 i4=1,2
           km=0
          do 62 i5=1,2
            if(ipair(i4,i5).gt.mp)km=km+1
   62     continue
c
          jj1=ipair(i4,1)
          jj2=ipair(i4,2)
          if(km .eq. 1) go to 67
          if(km .gt. 1) go to 63
c
c  ***one variable linked to two constants.fix the diagonal mat(i,i)
c
          jt1=j1(jj1)-1
          jt2=j1(jj2)-1
          jmin=iabs(jt1-jt2)
          jmax=jt1+jt2
c
          if(mat(i1,i1) .gt. 1) then
c
c  ***if there are several couples of constants ,take the more
c  ***stringent combination
c
            jmin=max0(jmin,jsum(1,i1))
            jmax=min0(jmax,jsum(2,i1))
            if(jmax.ge.jmin)then
             jsum(1,i1)=jmin
             jsum(2,i1)=jmax
             mat(i1,i1)=(jmax-jmin)/2+1
             go to 63
            else
             recup=zero
             go to 110
            endif
          else
          if(mat(i1,i1) .lt. 1) then
c
c  ***first time
c
            mat(i1,i1)=(jmax-jmin)/2+1
            jsum(1,i1)=jmin
            jsum(2,i1)=jmax
c
          endif
          endif
c
          go to 63
c
c
c
c
c  ***one variable linked to one constant and one variable  non diagonal
c  ***element
c
   67     jt1=min0(jj1,jj2)
          jt2=max0(jj1,jj2)-mp
          if(jt2.gt.i1)go to 63
          jt4=j1(jt1)-1
          k=mat(i1,jt2)
          if(k.eq.0)go to 107
c
      do 71 ll=1,k
        if(jt4.eq.j12(ll,jt2,i1))go to 63
   71 continue
c
  107 k=k+1
      if(k.gt.mxcsvr)go to 63
      mat(i1,jt2)=k
      j12(k,jt2,i1)=jt4
c
   63 continue
   65 continue
   66 continue
c
c  ***reduce the diagonal elements by taking into account the non
c  ***diagonal elements,and keep the latter only if needed
c
  150 ichan=0
c
      do 74 i=1,nsum
        noel(i)=.true.
        i1=i-1
        jlow(1) = 1
        jhig(1) = i1
        jlow(2) = i+1
        jhig(2) = nsum
        do 720 ireduc = 1,2
        if(i1.eq.0)go to 170
       do 72  j=jlow(ireduc),jhig(ireduc)
       if(ireduc.eq.1)then
         ik1=i
         ik2=j
       else
         ik1 = j
         ik2 = i
       endif
         if(mat(ik1,ik2).eq.0 .or. mat(j,j) .eq. 0) go to 72
      jmin1=0
      jmax1=1000
      k=mat(ik1,ik2)
c
      do 203 l1=1,k
c
        l3=mat(j,j)
        jj1=jsum(1,j)
        jnd=j12(l1,ik2,ik1)
        jmin=1000
        jmax=0
        jmnp(l1)=0
        jmxp(l1)=1000
c
      do 204 l2=1,l3
c
        jmn=iabs(jnd-jj1)
        jmx=jnd+jj1
        jmin=min0(jmn,jmin)
        jmax=max0(jmx,jmax)
        jmnp(l1)=max0(jmn,jmnp(l1))
        jmxp(l1)=min0(jmx,jmxp(l1))
        jj1=jj1+2
c
  204 continue
c
      jmin1=max0(jmin1,jmin)
      jmax1=min0(jmax1,jmax)
c
  203 continue
c
      if(mat(i,i).eq.0) then
        jsum(1,i)=jmin1
        jsum(2,i)=jmax1
        mat(i,i)=(jmax1-jmin1)/2+1
        ichan=ichan+1
        go to 206
      endif
c
      if(jsum(1,i).lt.jmin1) then
        jsum(1,i)=jmin1
        ichan=ichan+1
      endif
c
      if(jsum(2,i).gt.jmax1) then
        jsum(2,i)=jmax1
        ichan=ichan+1
      endif
c
  206 k1=0
c
      do 207 l1=1,k
        if(jmnp(l1).le.jsum(1,i).and.jmxp(l1).ge.jsum(2,i))go to 207
        k1=k1+1
        j12(k1,ik2,ik1)=j12(l1,ik2,ik1)
  207 continue
c
      if(k1.ne.k) then
        mat(ik1,ik2)=k1
        ichan=ichan+1
      endif
c
      mat(ik2,ik1)=j12(1,ik2,ik1)
      if(ireduc.eq.1)noel(i)=.false.
   72  continue
c
c
  170 if(i.eq.nsum)go to 74
720   continue
c
c
   74 continue
c
      if(ichan.ne.0)go to 150
c
c
c
c
c
c  ***carry out the summations.
c
  220 do 230 i=1,nsum
        jsum3(i)=1
        ldiag(i)=.false.
        if(mat(i,i).eq.1)ldiag(i)=.true.
  230 continue
c
      do 231 i=1,jwrd
        jwtest(i)=1
  231 continue
c
      stor=zero
      stor1=one
      nolp=0
      ip=1
      ipold=ip
  240 nolp=nolp+1
c
c
c  ***find the range of jsum2(nolp)
c  ***nolp is the index  of the summation variable
c
      jmin=jsum(1,nolp)
      jmax=jsum(2,nolp)
      if(noel(nolp))go to 241
      no1=nolp-1
c
      do 242 nj=1,no1
        if(mat(nolp,nj) .eq. 1) then
          jj1=mat(nj,nolp)
          jj2=jsum2(nj)
          jmin=max0(jmin,iabs(jj2-jj1))
          jmax=min0(jmax,jj1+jj2)
        else
        if(mat(nolp,nj) .gt. 1) then
          k=mat(nolp,nj)
          jj2=jsum2(nj)
c
         do 245 i=1,k
          jj1=j12(i,nj,nolp)
          jmin=max0(jmin,iabs(jj2-jj1))
          jmax=min0(jmax,jj1+jj2)
  245    continue
c
        endif
        endif
c
  242 continue
c
  241 jsum2(nolp)=jmin
      maxlp(nolp)=jmax
      if(ldiag(nolp))jsum3(nolp)=0
      if(nolp.lt.nsum)go to 240
c
      already=.false.
      ipold=ip
      do 260 jj=jmin,jmax,2
        jsum2(nsum)=jj
c
c
c  ***determine which racah coefficients need re-evaluating and
c  ***set jwtest appropriately
c
      do 114 j=ip,nsum
        if(jsum3(j).le.0) go to 114
        i2=jsum6(j)
c
      do 113 i1=1,i2
          i3=jsum5(j,i1)
          jwtest(i3)=1
c
  113 continue
c
  114 continue
c
      do 98 j=1,jwrd
        if(jwtest(j).eq.0)go to 98
        jwj=j+jwr
c
      do 90 i=1,6
        if(jword(i,jwj).le.mp) then
          i1=jword(i,jwj)
          ist(i) = j1(i1) - 1
        else
          i1=jword(i,jwj)-mp-nfs
          ist(i) = jsum2(i1)
        endif
   90 continue
c
      call dracah(x1)
      wstor(j)=x1
      if(ibug3.eq.1) then
        do 99 i=1,6
          xj1(i)=ist(i)/two
   99   continue
c
        write (ipd,315) jwj,(xj1(i), i=1,6),x1
      endif
   98 continue
c
c
c  ***form product of racah coefficients,(2j+1) factors and (-1)
c  ***factors in stor1
c
      do 126 i=1,jwrd
        stor1 = stor1*wstor(i)
  126 continue
c
c  ***iastor contains the power of (-1)which is common to all terms
c
      ix2 = 0
      xij6cp=1.
      if(j6cp.ne.j6f) then
        jb=j6f+1
c
        do 128 i=jb,j6cp
          i1=j6p(i)-nfs
          xij6cp=xij6cp*(jsum2(i1)+1.)
  128   continue
      endif
c
      if(j9cp.ne.j9f) then
        jb=j9f+1
c
        do 147 i=jb,j9cp
          i1=j9p(i)-nfs
          xij6cp=xij6cp/(jsum2(i1)+1.)
  147   continue
      endif
c
      stor1 = stor1*sqrt(xij6cp)
c
      if(j7cp.ne.j7f) then
        jb=j7f+1
c
        do 131 i=jb,j7cp
          i1=j7p(i)-nfs
          ix2 = ix2 + jsum2(i1)
  131   continue
      endif
c
      if(j8cp.ne.j8f) then
        jb=j8f+1
c
        do 134 i=jb,j8cp
          i1=j8p(i)-nfs
          ix2 = ix2 + 2*(jsum2(i1))
  134   continue
      endif
c
      if(mod(ix2,2).eq.1) then
        ias=-1
        ix2=ix2+1
      endif
c
      ix2 = ix2/2
c
c
c  ***add term into stor and reset stor1 to 1 ready for next term
c
      if (mod(ix2,2) .eq. 1) stor1 = -stor1
      stor = stor + stor1
      stor1=one
      nsum1 =nsum-1
      already=.true.
      if(nsum1.eq.0)go to 260
c
      do 261 ik=1,nsum1
        jsum3(ik)=0
  261 continue
c
      do 262 ik=1,jwrd
        jwtest(ik)=0
  262 continue
c
  260 continue
c
  250 nolp=nolp-1
c
      if(nolp.ne.0) then
        if(ldiag(nolp))go to 250
        jsum3(nolp)=1
        jsum2(nolp)=jsum2(nolp)+2
        if(jsum2(nolp).gt.maxlp(nolp))go to 250
        ip=nolp
        if(.not.already)ip=min(ipold,ip)
c
c
c    ***proceed to next variable
c
        go to 240
c
      endif
c
      recup=recup*stor
      if(ibug3.eq.1) write(ipd ,307) nps,stor,recup
      if(abs(recup).lt.epsil)go to 145
      jwr=jwrd+jwr
      nfs=nsum+nfs
      j6f=j6cp
      j7f=j7cp
      j8f=j8cp
      j9f=j9cp
      iastor=iastor+ias
c
c
c  ***proceed to next sum
c
      if(nps.lt.nlsum)go to 25
      iastor=iastor/2
      if(mod(iastor,2).ne.0)recup=-recup
      if(ibug3.eq.1) write(ipd ,304) recup
  110 return
c
c  ***no summations. check that there are no inconsistencies. then
c  ***multiply by (-1) factor and exit
c
  145 recup=zero
      return
      end
c
c     ****************
      subroutine intab
      implicit real*8(a-h,o-z)
c     ****************
c
c  ***this subroutine called at the end of diagrm,fixes the arrays ih 
c  ***and il-so to speak hardware and logical addresses of triads in
c  ***jdiag.also determines the number of free ends nfree and their
c  ***location itfree.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm)
c
      logical free
c
      integer arr,tab1
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/build/ial(m4trd),if1,if2,node
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
c
c
      do 1 i=1,m
    1   ial(i)=1
c
      do 3 i=ifirst,ilast
        j=jdiag(i,1)
        k=ial(j)
        tab1(j,k)=i 
        ial(j)=k+1
    3 continue
c
      ifr=ifirst-1
c
      do 4 i=ifirst,ilast
        it=i-ifr
        il(i)=it
        ih(it)=i
    4 continue
c
      j=jdiag(ifirst,3)
      k=ial(j)
      if(k .gt. 1) tab1(j,2)=tab1(j,1)
      tab1(j,1)=ifirst
      ial(j)=3
      j=jdiag(ilast,2)
      tab1(j,2)=ilast
      ial(j)=3
      nfree=0
c
      do 7 i=ifirst,ilast
        j=jdiag(i,1)
        if(ial(j).ne.3)then
          nfree=nfree+1
          itt=ilast+nfree
          tab1(j,2)=itt
          il(itt)=nfree*1000
          itfree(nfree)=i
        endif
    7 continue
c
      return
      end 
c
c
c     *********************** 
      subroutine lolpop(fail) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c  ***reduces a loop with one line and one node in the flat graph.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
c
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),
     + free(mangm)
c
      logical fail,sumvar
c
      integer arr,tab1
c
      character*6 name,namsub 
c
      dimension kp(3),ks(3)
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
      logical free
c
      common/nam/namsub
c
      data name/'lolpop'/
      data kp/2,3,1/
      data ks/0,1,-1/
c
c
c
      namsub=name
      i1=npoint(1)
      k3=2
      if(i1.eq.ilast)k3=3
      l=jdiag(i1,k3)
      call delta(l,mp,fail)
      if(fail)return
      k=kp(k3)
      if(arr(i1,k).lt.0)call phase2(jdiag(i1,k))
      k1=ks(k3)
      il1=il(i1)+k1 
      i2=ih(il1)
      l1=jdiag(i2,1)
      call delta(l1,jdiag(i2,k3),fail)
      if(fail)return
      if(arr(i2,k3).eq.k1)call phase2(l1)
      il2=il(i2)+k1 
      i3=ih(il2)
      k2=k3+k1
      jdiag(i3,k2)=l1
      arr(i3,k2)=arr(i2,1)
      j9c=j9c+1
      j9(j9c)=l1
      j6c=j6c+1
      j6(j6c)=jdiag(i1,1)
      if(k3.eq.3)return
c
      do 1 i=3,nbnode
        it=ih(i)
        ilp=i-2
        il(it)=ilp
        ih(ilp)=it
    1 continue
c
      return
      end 
c
c
c     ***************************
      subroutine neibor(lc,l1,l2)
      implicit real*8(a-h,o-z)
c     ***************************
c
c  ***gives the positions of the other two arguments in the triad.
c
c
      if (lc .lt. 2) then
        l1=2
        l2=3
      else
      if (lc .eq. 2) then
        l1=3
        l2=1
      else
        l1=1
        l2=2
      endif
      endif
      return
      end 
c
c
c     ****************
      subroutine ordtri
      implicit real*8(a-h,o-z)
c     ****************
c
c  ***this subroutine orders the triads which were left with free ends
c  ***as consequence of cutting,so that the new graph will start there.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar,tabs
c
      integer arrow,arr,tab1
c
      character*6 name
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/build/ial(m4trd),if1,if2,node
      common/keep/jkp(2,3),jarr(2,3),it2,it3,it5
c
      data name/'ordtri '/
c
c
      do 10 i=1,mp
        ial(i)=0
   10 continue
c
      if(nfin.ne.0)then
        nbt1=nbtr-1 
        nbt=nbt1+nfin
        nbtt=nbt+1
        nb=0
        go to 31
      endif
c
      nf=nbtr-itfree(1)
c
      if(it5.eq.0) then
        nbt1=nbtr-1 
        n0=0
        nft=nfree
        isw=2
        go to 100
      endif
c
      nft=it5-it2
      nm=nft+nbtr+1 
      nbt1=nbtr
c
      do 21 j=1,3
        jdiag(nbtr,j)=jkp(1,j)
        arr(nbtr,j)=jarr(1,j) 
c  June 25  1989  (avi)
c        jdiag(nm,j)=jkp(2,j)
c        arr(nm,j)=jarr(2,j)
   21 continue
c
      n0=0
      isw= 1
      go to 100
c
   22 n0=nft
c June 25 1989 (avi)
      do 211 j=1,3
        jdiag(nm,j)=jkp(2,j)
        arr(nm,j)=jarr(2,j)
 211  continue
c June 25 1989 (avi)
c      nbt1=nbt1+n0
      nbt1=nbt1+1
      nft=it3-it5
      isw= 3
      go to 100
c
c June 25 1989 (avi)
c   24 nft=nft+1
   24 nbt1=k-nft
c
   23 node=nbt1+nft 
      call change(node,2)
      go to 40
c
c
   31 do 35 i=1,nbnode
        i1=ih(i)
        if(il(i1).gt.ilast)go to 35
        i2=nbt1+i
        if(i1.gt.nbtt)go to 33
        if(i1.eq.i2)go to 32
        if(il(i2).le.nbnode)go to 35
c
   33 do 34 j=1,3
        jdiag(i2,j)=jdiag(i1,j)
        arr(i2,j)=arr(i1,j)
   34 continue
c
      il(i1)=ilast+i
   32 nb=nb+1
      il(i2)=0
c
   35 continue
c
      if(nb.ne.nfin)go to 31
      node=nbt
   40 if1=jdiag(nbtr,1)
      if2=jdiag(nbtr,3)
c
      do 51 i=nbtr,node
       do 50 k=1,3
         j=jdiag(i,k)
         ial(j)=ial(j)+1
   50  continue
   51 continue
c
      ilast=node
      call printj(name,8)
c
      return
c
  100 if(nf.le.0)then
        nfr=n0
        i1=1
      else
        nfr=nft+1
        i1=-1
      endif
c
      do 4 i=1,nft
        ik=nfr+i1*i 
        it=itfree(ik)
        k=nbt1+ik
c
      do 3 j=1,3
        jdiag(k,j)=jdiag(it,j)
        arr(k,j)=arr(it,j)
    3 continue
c
    4 continue
c
      go to (22,23,24),isw
      end 
c
c
c     ********************************* 
      subroutine otherj(lin,j,lo,lco,k) 
      implicit real*8(a-h,o-z)
c     ********************************* 
c
c
c  ***gives the other triad where a given j occurs and its position.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm)
c
      logical tabs
      integer arrow 
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
c
      lo=line(j,1)
      if(lo.eq.lin.or.tabs(lo)) then
        k=1
        lo=line(j,2)
        lco=lcol(j,2)
      else
        k=2
        lco=lcol(j,1)
      endif
c
      return
      end 
c
c
c     ***************************
      subroutine phase(l,jm,ndim)
      implicit real*8(a-h,o-z)
c     ***************************
c
c  ***phase factor arising from non-cyclic permutation of arguments in
c  ***triad l.jm may be either j23 or jdiag.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar
      dimension jm(ndim,3)
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
c
      j7(j7c+1)=jm(l,1)
      j7(j7c+2)=jm(l,2)
      j7c=j7c+3
      j7(j7c)=jm(l,3)
c
      return
      end 
c
c
c     ********************
      subroutine phase2(j)
      implicit real*8(a-h,o-z)
c     ********************
c
c  ***adds a phase factor (-1)**2j
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
c
      j8c=j8c+1
      j8(j8c)=j
c
      return
      end 
c
c
c     *********************** 
      subroutine polygn(jpol) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c  ***this routine reduces a circuit of arbitrary order nc.it exchanges
c  ***nodes on the flat diagram until the distance on the axis between
c  ***nodes equeals one.each exchange introduces a summation variable 
c  ***and a 6j symbol.the circuit has a maximum of npart=2 disconnected
c  ***parts on the axis.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      integer arr,tab1
      logical sumvar
      character*6 name
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'polygn'/
c
c
      nc1=nc+1
      nc2=nc
      nbc=ipartl-2
c
   10 do 8 i=1,nbc
        it2=npoint(nc1-i)
        it1=npoint(nc2-i)
        jb=jdiag(it1,1)
        jc=jdiag(it2,1)
        jdiag(it1,1)=jc
        jdiag(it2,1)=jb
        jar=arr(it1,1)
        arr(it1,1)=arr(it2,1) 
        arr(it2,1)=jar
        je=jdiag(it1,2)
        mp=mp+1
        sumvar(mp)=.true.
        jdiag(it1,2)=mp
        jdiag(it2,3)=mp
c
        if(tab1(jb,1) .eq. it1) then
          tab1(jb,1)=it2
        else
          tab1(jb,2)=it2
        endif
c
        if(tab1(jc,1) .eq. it2) then
          tab1(jc,1)=it1
        else
          tab1(jc,2)=it1
        endif
c
        if(arr(it1,2).le.0)then
          call phase2(je)
          arr(it1,2)=1
          arr(it2,3)=-1
        endif
c
        jwc=jwc+1
        kw(1,jwc)=jb
        kw(2,jwc)=mp
        kw(3,jwc)=je
        kw(4,jwc)=jc
        kw(5,jwc)=jdiag(it2,2)
        kw(6,jwc)=jdiag(it1,3)
        j6(j6c+1)=mp
        j6c=j6c+2
        j6(j6c)=mp
    8 continue
c
      nc=nc-nbc
c
      if(nc .gt. 4) then
        nbc=iparts-2
        nc1=iparts+1
        nc2=iparts
        go to 10
      endif
c
      if(npart .ne. 1) then
        npoint(3)=npoint(nc1) 
        npoint(4)=npoint(nc1+1)
      endif
c
      if(nc.eq.2)jpol=1
      call printj(name,msum)
c
      return
      end 
c
c
c     ***************************
      subroutine printj(names,jp)
      implicit real*8(a-h,o-z)
c     ***************************
c
c  ***this subroutine prints intermediate results in standard form from
c  ***wherever it is called.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
      parameter(mtab=30)
c
      integer arr,tab1,arrow
c
      logical tabs,sumvar,free
c
      character im,ip,is(3)
      character*4 i6,i7,i8,i9,ij1
      character*6 names,nsettb
      character*8 iblank,ifree,ifr
c
      dimension ix(6),jtab(mtab,3)
c
      equivalence(i6c,ix(1))
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/const/i6c,i7c,i8c,i9c,idel,iwc
      common/zer/nzero,jzero(m6j)
      common/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      data iblank,ifree,ip,im/'        ','free end','+','-'/
      data nsettb/'settab'/
      data i6,i7,i8,i9,ij1/'i6= ','i7= ','i8= ','i9= ','j1= '/
c
 1000 format(/10x,'nbnode=',i3,10x,'nbtr=',i3,10x,'nfin=',i3,
     + /10x,'ifirst=',i3,10x,'ilast=',i3,9x,'nfree=',i3)
 1001 format(//7x,'il',3x,'ih',14x,'jdiag'//)
 1002 format(28x,3(a1,2x))
 1003 format(7x,i2,3x,i2,2x,a8,2x,3i3/) 
 1004 format(/5x,'tab1'/)
 1005 format(4(i3,1h),2x,i3,i5,5x))
 1006 format(/2x,'sumvar=',15(i3,l1))
 1010 format(//10x,'j23',10x,'nbtr1=',i3//)
 1012 format(18x,3(a1,2x))
 1013 format(i9,i5,2x,3i3/)
 1014 format(/3x,'j  l1 k1  l2 k2')
 1015 format(4(i4,1h),i3,i3,i4,i3))
 1020 format(/3x,a4,3x,3(20i3/))
 1021 format(/3x,'delta=',7(i5,i3))
 1022 format(/3x,'kw(arg. of 6j)',6i3)
 1030 format(//2x,'nc=',i2,4x,'npart=',i2,4x,'ipartl=',i2,4x,
     +  'iparts=',i2,4x,'icross=',i2,4x,/2x,'npoint=',20i3) 
 1040 format(//2x,'nzero=',i2,5x,12(i4,1h),i3))
 1050 format(///3x,'print out after calling subroutine ',a7)
c
c
c
      if(ibug3.ne.1)return
      print 1050,names
c
c  ***initialise variables
c
c
      jump=jp
      if(jump.eq.0)then
c
        do 9 i=1,6
          ix(i)=1
    9   continue
c
        print 1020,ij1,(j1(i),i=1,m)
      endif
c
      if(jump.lt.8)go to 20
      print 1000,nbnode,nbtr,nfin,ifirst,ilast,nfree
      jump=jump-8
      print 1001
      k=0 
c
      do 1 i=1,nbnode
        it=ih(i)
        ifr=iblank
        jt=jdiag(it,1)
c
        if(tab1(jt,2).ne.it.or.jt.eq.jdiag(ifirst,3))then
          k=k+1
          jtab(k,1)=jt
          jtab(k,2)=tab1(jt,1)
          jtab(k,3)=tab1(jt,2)
        endif
c
      if(tab1(jt,2).gt.ilast)ifr=ifree
c
      do 2 j=1,3
        is(j)=ip
        if(arr(it,j).lt.1)is(j)=im
    2 continue
c
      print 1002,(is(j),j=1,3)
      print 1003,il(it),it,ifr,(jdiag(it,j),j=1,3)
c
    1 continue
c
      print 1004
      ntime=0
      jt=jdiag(ifirst,3)
      if(jt.ne.jdiag(ilast,2))then
       if(tab1(jt,2).lt.1000)go to 5
      endif
    4 k=k+1
      jtab(k,1)=jt
      jtab(k,2)=tab1(jt,1)
      jtab(k,3)=tab1(jt,2)
    5 ntime=ntime+1 
c
      if(ntime.ne.2) then
        jt=jdiag(ilast,2)
        if(tab1(jt,2).eq.1000)go to 4
      endif
c
      print 1005,((jtab(i,j),j=1,3),i=1,k)
      print 1006,(i,sumvar(i),i=1,mp)
   20 if(jump.lt.4)go to 30
      jump=jump-4
      nbtr1=2*n-2
      print 1010,nbtr1
      k=0 
c
      do 11 i=1,nbtr1
        if(tabs(i))go to 11
        k=k+1
c
      do 12 j=1,3
        is(j)=ip
        if(arrow(i,j).lt.1)is(j)=im
   12 continue
c
      print 1012,(is(j),j=1,3)
      print 1013,k,i,(j23(i,j),j=1,3)
c
   11 continue
c
      print 1014
      mm=m
      if(names.ne.nsettb) mm=m-1
      print 1015,(i,(line(i,j),lcol(i,j),j=1,2),i=1,mm)
c
   30 if(jump.ge.2)then
        jump=jump-2 
        print 1030,nc,npart,ipartl,iparts,icross,(npoint(i),i=1,nc)
      endif
c
      if(jump.ge.1) print 1040,nzero,(i,jzero(i),i=1,nzero) 
      if(j6c.ge.i6c)print 1020,i6,(j6(i),i=i6c,j6c)
      if(j7c.ge.i7c)print 1020,i7,(j7(i),i=i7c,j7c)
      if(j8c.ge.i8c)print 1020,i8,(j8(i),i=i8c,j8c)
      if(j9c.ge.i9c)print 1020,i9,(j9(i),i=i9c,j9c)
      if(jdel.ge.idel)print 1021,((ldel(i,j),j=1,2),i=idel,jdel)
      if(jwc.ge.iwc)print 1022,((kw(j,i),j=1,6),i=iwc,jwc)
      i6c=j6c+1
      i7c=j7c+1
      i8c=j8c+1
      i9c=j9c+1
      idel=jdel+1
      iwc=jwc+1
      return
      end 
c
c
c     *********************** 
      subroutine search(find) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm)
c
      logical find
      integer arr,tab1
      character*6 name
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      data name/'search'/
c
 1000 format(' error in search.i,i1,i2,i3,npart,ipart,nc=',7i5)
c
c  ***this subroutine locates circuits or loops of order nc.npoint(nc)
c  ***are the indices of the points(triads) pertaining to the first
c  ***such loop found.
c  ***npart is the number of separate parts(groups of contiguous points)
c  ***on the axis of the flat graph.iparts is the number of points in 
c  ***the smallest part.ipartl is the number of points in the largest 
c  ***part.
c  ***this subroutine finds all the possible loops of order 3 and 4.for
c  ***nc.ge.5,it looks for only those who are partitionned in npart.le.2
c  ***which can eventually reduce tp a loop of order 4 without breaking
c  ***the basic structure of the flat graph. icross=-1,if lines cross 
c-------------------------------------------------------------------- 
c
c
c  ***initialization
c
      find=.false.
      ncm1=nc-1
      ncm=nc-2
      icross=0
c
c  ***first are treated two cases that do not involve do loops
c  ***1.one isolated point,either the first or the last
c
      npart=1
      ipartl=nc-1
      iparts=1
c
c  ***a.first
c
      i1=ifirst
      k3=3
      k2=2
  200 ja=jdiag(i1,1)
      jc=jdiag(i1,k3)
c
      if(ja.eq.jc) then
        if(nc.gt.1) go to 800 
        npoint(1)=i1
        go to 900
      endif
c
      i2=tab1(ja,k2)
      i3=tab1(jc,k2)
c
      if(iabs(il(i3)-il(i2))-ncm .lt. 0) go to 800
c
      if(iabs(il(i3)-il(i2))-ncm .gt. 0) then
c
c  ***b.last
c
        if(i1.ne.ifirst) go to 250
        i1 = ilast
        k3=2
        k2=1
        go to 200
      endif
c
      ic=1
      npoint(ic)= i1
      i20=min0(i2,i3)
      i21=il(i20)
      i31=i21+ncm1
c
      do 203 ii=i21,i31
        ic=ic+1
        npoint(ic)=ih(ii)
  203 continue
c
      if(nc.le.2) then
        if(jdiag(ifirst,1).ne.jdiag(ilast,1))call phase(i1,jdiag,m4trd)
        go to 900
      endif
c
      if(i1.ne.ilast) then
        it=i2
        jt=jdiag(ilast,2)
        k4=2
        i4=ilast
      else
        it=i3
        jt=jdiag(ifirst,3)
        k4=3
        i4=ifirst
      endif
c
      if(it.eq.i20)call phase(i1,jdiag,m4trd)
      if(jt.eq.ja.or.jt.eq.jc)call change(i4,k4)
      go to 900
c
c  ***2.two isolated points,first and last.
c
  250 if(nc.eq.1)return
      if(nc.le.3) go to 100
      ipartl=nc-2
      iparts=1
      i1=ifirst
      i2=ilast
      ja=jdiag(i1,1)
      jb=jdiag(i1,3)
c
      if(tab1(ja,2).ne.i2) then
        ja=jdiag(i1,3)
        jb=jdiag(i1,1)
        if(tab1(ja,2).ne.i2) go to 100
      endif
c
      if(ja.eq.jdiag(i2,1)) then
        jc=jdiag(i2,2)
      else
        jc=jdiag(ilast,1)
      endif
c
      i3=tab1(jb,2) 
      i4=tab1(jc,1) 
      idist=il(i4)-il(i3)
c
      if(iabs(idist)-(ncm-1) .lt. 0) go to 800
      if(iabs(idist)-(ncm-1) .eq. 0) then
        npoint(1)= ilast
        npoint(2) = ifirst
        icross=isign(1,idist) 
        ic=2
        i20=min0(i3,i4)
        i21=il(i20) 
        i31=i21+ncm 
c
        do 261 ii=i21,i31
          ic=ic+1
          npoint(ic) = ih(ii) 
  261   continue
c
        if(ja.eq.jdiag(ifirst,1))call change(ifirst,3)
        if(ja.eq.jdiag(ilast,1))call change(ilast,2)
        go to 900
      endif
c
c  ***first general case:all points in one group
c
  100 npart=1
      iparts=0
      ipartl=nc
      k3=1
c
      do 101 in=1,nbnode
        i=ih(in)
  108   ja=jdiag(i,k3)
        if(i.ne.tab1(ja,2))then
          i2=tab1(ja,2)
c
          if(il(i2)-in-ncm1 .lt. 0) go to 800
          if(il(i2)-in-ncm1 .eq. 0)then 
            i21=il(i2)
            ic=0
c
            do 103 ii=in,i21
             ic=ic+1
             npoint(ic)=ih(ii)
  103       continue
c
            if(ja.eq.jdiag(ifirst,3))call change(ifirst,3)
            if(ja.eq.jdiag(ilast,2))call change(ilast,2)
            go to 900
        endif
      endif
c
        if(in.eq.1) then
          if(k3.ne.3) then
            k3=3
            go to 108
          else
            k3=1
          endif
        endif
c
  101 continue
c
c  ***search did not find loop nc.le.3
c
      if(nc.le.3) return
c
c  ***general case of loop partitionned in 2 groups.do loop 
c  ***on iparts
c
      npart=2
      nc2=nc/2
      k3=1
      k2=1
c
      do 400 ips=2,nc2
        jps=ips-1
        nbn=nbnode-jps
c
      do 301 i1=1,nbn
        i=ih(i1)
        i2=ih(i1+jps)
  302   ja=jdiag(i,k3)
        jd=jdiag(i2,k2)
c
        if(i.eq.tab1(ja,1)) then
          ii2=tab1(jd,2)
          ii1=tab1(ja,2)
        else
          ii1=tab1(ja,1)
          ii2=tab1(jd,1)
        endif
c
        idist=il(ii1)-il(ii2) 
c
        if(iabs(idist)-(ncm-jps) .lt. 0) go to 800
        if(iabs(idist)-(ncm-jps) .gt. 0) go to 320
  306   icross=isign(1,idist) 
        ic=0
        i21=il(i2)
c
      do 310 ii=i1,i21
        ic=ic+1
        npoint(ic)=ih(ii)
  310 continue
c
      i20=min0(ii1,ii2)
      i30=max0(ii1,ii2)
      i21=il(i20)
      i31=il(i30)
c
      do 311 ii=i21,i31
        ic=ic+1
        npoint(ic)=ih(ii)
  311 continue
c
      iparts=ips
      ipartl=nc-ips 
      if(jdiag(ifirst,3).eq.ja.or.jdiag(ifirst,3).eq.jd)call
     +change(ifirst,3)
      if(jdiag(ilast,2).eq.ja.or.jdiag(ilast,2).eq.jd)call
     +change(ilast,2)
      go to 900
c
  320 if(i1.eq.1) then
        if(k3.eq.3) then
          k3=1
          go to 301 
        else
          k3=3
          go to 302 
      endif
      endif
c
      if(i2.eq.ilast)then
       if(k2.ne.2) then
        k2=2
        go to 302
       endif
      endif
c
  301 continue
  400 continue
c
c  ***search did not find circuit of order nc
c
      return
c
c  ***loop found
c
  900 find=.true.
      call printj(name,msum)
c
      return
c
c  ***error printout
c
  800 print 1000,i,i1,i2,i3,npart,iparts,nc
      stop
c
      end 
c
c
c     *****************
      subroutine setdim
      implicit real*8(a-h,o-z)
c     *****************
c
c  ***set dimensions of arrays.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc
c
c
c
      jwcc=jwc
      jdelc=jdel
      j6cc=j6c
      j7cc=j7c
      j8cc=j8c
      j9cc=j9c
c
      return
      end 
c
c
c     *********************** 
      subroutine settab(fail) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c  ***builds up the unstructured graph
c  ***sets the array j23,containing the two lists of original triads
c  ***j2 and j3,and the corresponding arrows on the angular momenta
c  ***lines.also establishes the numerical and phase factors connecting
c  ***recoupling coefficient and graphs,according to yutsis,levinson and
c  ***vanagas.for this purpose determines the total j
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,tabs,free,sumvar
c
      integer arrow 
c
      character*6 name
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
      common/build/ial(m4trd),if1,if2,node
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'settab'/
c
c
c
c
      ipr=n-1
      nbtr=ipr+ipr
c
      do 4 i=1,ipr
       do 5 j=1,2
         j23(i,j)=j2(i,j)
         arrow(i,j)=1
    5  continue
       tabs(i)=.false.
       j23(i,3)=j2(i,3)
       arrow(i,3)=-1
    4 continue
c
      ipr1=ipr+1
c
      do 7 i=ipr1,nbtr
        ii=i-ipr
      do 6 j=1,2
        j23(i,j)=j3(ii,j)
        arrow(i,j)=-1
    6 continue
      tabs(i)=.false.
      j23(i,3)=j3(ii,3)
      arrow(i,3)=1
    7 continue
c
      do 11 j=1,nbtr
        j8(j)=j23(j,1)
   11 continue
c
      j8c=nbtr+ipr
      nb1=nbtr+1
c
      do 12 j=nb1,j8c
        i=j-ipr
        j8(j)=j23(i,3)
   12 continue
c
      j6c=nbtr
c
      do 13 j=1,j6c 
        j6(j)=j23(j,3)
   13 continue
c
      do 10 i=1,m
        sumvar(i)=.false.
        ial(i)=1
   10 continue
c
      do 9 i=1,nbtr 
       do 8 j=1,3
         ji=j23(i,j)
         k=ial(ji)
         line(ji,k)=i
         lcol(ji,k)=j
         ial(ji)=k+1
    8  continue
    9 continue
c
      it=0
c
      do 18 i=1,nbtr
        jt=j23(i,3) 
c
        if(ial(jt).eq.3) then 
          call otherj(i,jt,l,lc,k)
          if(lc.eq.3)go to 19 
          go to 18
        endif
c
        if(it.eq.1) then
          call delta(jt1,jt,fail)
          if(fail)go to 20
          k=line(jt,1)
          kc=lcol(jt,1)
          line(jt1,2)=k
          lcol(jt1,2)=kc
          line(jt,2)=line(jt1,1)
          lcol(jt,2)=lcol(jt1,1)
          j23(k,kc)=jt1
          ial(jt)=1 
          go to 19
        endif
c
        jt1=jt
        it=1
c
   18 continue
c
   19 j9(j9c+1)=jt
      j9c=j9c+2
      j9(j9c)=jt
c
   20 call printj(name,4)
c
      return
      end 
c
c
c     ********************
      subroutine sprate(m)
      implicit real*8(a-h,o-z)
c     ********************
c
c  ***this subroutine prepares the information to be transfered to
c  ***gensum for numerical evaluation.the common blocks /graph/ and
c  ***/tree/ are used as working memory,and their previous content
c  ***is destroyed. 
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sum6j,t6j,jt,js,sumvar,cut
c
      character*4 name
c
      common/graph/jtem4(mtriad,m6j),jtem5(mtriad,m6j),jtem6(mtriad), 
     +  nsum6j(m6j),j6sum(m6j),idum1(44)
c
      common/tree/sum6j(m6j),t6j(m6j),jt(mtriad),js(mtriad),
     +  inver(mangm),jnsum(mtriad),jinv(mtriad),n6jn(m6j),in6j(m6j),
     +  jsumt(m6j,6),idum(101)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),jw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/cutdig/cut
      common/dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc
c
      common/sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),
     + jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),
     + k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),
     + jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)
c
c
      data kflw,kfl6,kfl7,kfl8,kfl9/m6j,m3mngm,m3mngm,m3mngm,mangmp/
      data kfls/mtriad/
c
 1000 format(2x,'dimension error for  ',a4,i5,' is out of allowed ', 
     +'range',i3)
c
c
c
c  ***test that array dimensions have not been exceeded.
c
      if(jwc.gt.kflw)then
        nmx=kflw
        npx=jwc
        name='kflw' 
      else
      if(j6c.gt.kfl6)then
        nmx=kfl6
        npx=j6c
        name='kfl6' 
      else
      if(j7c.gt.kfl7)then
        nmx=kfl7
        npx=j7c
        name='kfl7' 
      else
      if(j8c.gt.kfl8)then
        nmx=kfl8
        npx=j8c
        name='kfl8' 
      else
      if(j9c.le.kfl9)go to 54 
        nmx=kfl9
        npx=j9c
        name='kfl9' 
      endif
      endif
      endif
      endif
c
   60 print 1000,name,npx,nmx 
      stop
c
c  ***determination of effective summation variables and their
c  ***relationships with 6j coefficients.
c
c
   54 do 2 i=1,jwc
        inv6j(i)=0
        sum6j(i)=.false.
    2 continue
c
      nsum=0
      nlsum=0
      if(mp.eq.m)return
      m1=m+1
c
      do 1 i=m1,mp
        if(sumvar(i))then
          nsum=nsum+1
          jsum6(nsum)=0
          inver(i)=nsum
        endif
    1 continue
c
      if(nsum.eq.0)return
c
      if(nsum.gt.kfls)then
        nmx=kfls
        npx=nsum
        name='nsum' 
        go to 60
      endif
c
      kt=0
c
      do 4 i=1,jwc
       do 5 j=1,6
         ik=jw(j,i) 
         if(.not.sumvar(ik))go to 5
c
         if(.not.sum6j(i))then
           sum6j(i)=.true.
           kt=kt+1
           j6sum(kt)=0
           nsum6j(kt)=i
           inv6j(i)=kt
         endif
c
         isk=inver(ik)
         i2=jsum6(isk)+1
         jsum6(isk)=i2
         jsum4(isk,i2)=j
         jsum5(isk,i2)=kt
         i3=j6sum(kt)+1
         j6sum(kt)=i3
         jsumt(kt,i3)=isk
    5 continue
    4 continue
c
      call var(j6,j6p,j6c,j6cp,j6cc,sumvar,mp,m,inver)
      call var(j7,j7p,j7c,j7cp,j7cc,sumvar,mp,m,inver)
      call var(j8,j8p,j8c,j8cp,j8cc,sumvar,mp,m,inver)
      call var(j9,j9p,j9c,j9cp,j9cc,sumvar,mp,m,inver)
c
      if(.not. cut)then
        nlsum=1
        nbj(1)=nsum 
        nb6j(1)=kt
        k6cp(1)=j6cp
        k7cp(1)=j7cp
        k8cp(1)=j8cp
        k9cp(1)=j9cp
c
        do 21 i=1,kt
          i1=nsum6j(i)
        do 22 j=1,6 
          jword(j,i)=jw(j,i1) 
   22   continue
   21   continue
  
        do 80 i=1,nsum
          isu=jsum6(i)
         do 81 j=1,isu
           i1=jsum5(i,j)
           j1=jsum4(i,j)
           jword(j1,i1)=mp+i
   81    continue
   80   continue
c
       return
      endif
c
c  ***separation of variables and sums in case a cut was detected.
c
      k6c=0
      k7c=0
      k8c=0
      k9c=0
      nj=0
      n6j=0
c
      do 9 i=1,kt
        t6j(i)=.false.
    9 continue
c
      do 7 i=1,nsum 
        jt(i)=.false.
        js(i)=.false.
    7 continue
c
      j=1 
c
   10 nj=nj+1
      jnsum(nj)=j
      jinv(j)=nj
      jt(j)=.true.
   18 js(j)=.true.
      js6=jsum6(j)
c
      do 11 i=1,js6 
        i6j=jsum5(j,i)
c
        if(.not.t6j(i6j))then 
          t6j(i6j)=.true.
          n6j=n6j+1 
          n6jn(n6j)=nsum6j(i6j)
          in6j(i6j)=n6j
        endif
c
        j6j=j6sum(i6j)
c
      do 12 k=1,j6j 
        jk=jsumt(i6j,k)
        if(.not.jt(jk))then
          nj=nj+1
          jnsum(nj)=jk
          jinv(jk)=nj
          jt(jk)=.true.
        endif
   12 continue
c
   11 continue
c
      do 13 jj=1,nsum
        j=jj
        if(.not.js(jj) .and. jt(jj))go to 18
   13 continue
c
      nlsum=nlsum+1 
      nbj(nlsum)=nj 
      nb6j(nlsum)=n6j
c
      if(j6cp.ne.0)call chvar(j6p,j6cp,k6c,jt,jinv,nsum)
      k6cp(nlsum)=k6c
      if(j7cp.ne.0)call chvar(j7p,j7cp,k7c,jt,jinv,nsum)
      k7cp(nlsum)=k7c
      if(j8cp.ne.0)call chvar(j8p,j8cp,k8c,jt,jinv,nsum)
      k8cp(nlsum)=k8c
      if(j9cp.ne.0)call chvar(j9p,j9cp,k9c,jt,jinv,nsum)
      k9cp(nlsum)=k9c
c
      if(nj.ne.nsum)then
        do 16 jj=1,nsum
          j=jj
          if(.not.jt(jj))go to 10
   16   continue
      endif
c
      do 26 i=1,kt
        i1=n6jn(i)
       do 27 j=1,6
        jword(j,i)=jw(j,i1)
   27  continue
   26 continue
c
      do 28 i=1,nsum
        ik=jnsum(i) 
        i2=jsum6(ik)
        jtem6(i)=i2 
       do 29 j=1,i2 
         jtem4(i,j)=jsum4(ik,j)
         k=jsum5(ik,j)
         jtem5(i,j)=in6j(k)
   29  continue
   28 continue
c
      do 40 i=1,nsum
        i2=jtem6(i) 
        jsum6(i)=i2 
       do 41 j=1,i2 
         i1=jtem5(i,j)
         j1=jtem4(i,j)
         jsum4(i,j)=j1
         jsum5(i,j)=i1
         jword(j1,i1)=i+mp
   41  continue
   40 continue
c
      return
      end 
c
c
c     *****************
      subroutine square
      implicit real*8(a-h,o-z)
c     *****************
c
c  ***reduces a circuit of order 4 in the two cases which are left
c  ***over by polygn,namely two disconnected groups of two points
c  ***and one group of two points plus the two ends of the axis.in
c  ***the latter, the end of the axis is transferred to the beginning.
c  ***in this process,one summation variable and two 6j symbols are
c  ***introduced.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar
      integer arr,tab1
      character*6 name,namsub 
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/nam/namsub
c
      data name/'square'/
c
c
c
      namsub=name
      mp=mp+1
      sumvar(mp)=.true.
      k=1 
      it1=npoint(1) 
      it2=npoint(2) 
c
      if(icross.eq.1)then
        it3=npoint(3)
        it4=npoint(4)
        k23=3
        k32=2
      else
        it3=npoint(4)
        it4=npoint(3)
        k23=2
        k32=3
      endif
c
      l4=jdiag(it2,1)
c
      if(arr(it2,1).le.0)then 
        call phase2(l4)
        arr(it2,1)=1
        arr(it3,1)=-1
      endif
c
      l2=jdiag(it1,1)
      if(arr(it1,1).gt.0)call phase2(l2)
      jwc=jwc+1
      kw(1,jwc)=l4
      kw(2,jwc)=l2
      kw(3,jwc)=jdiag(it2,2)
      jj1=jdiag(it1,3)
      kw(4,jwc)=jj1 
      kw(5,jwc)=mp
      kw(6,jwc)=jdiag(it1,2)
      if(arr(it1,2).lt.0)call phase2(jdiag(it1,2))
      jwc=jwc+1
      kw(1,jwc)=l4
      kw(2,jwc)=l2
      jj3=jdiag(it3,k23)
      jj2=jdiag(it4,k32)
      kw(3,jwc)=jj3 
      kw(4,jwc)=jj2 
      kw(5,jwc)=mp
      kw(6,jwc)=jdiag(it3,k32)
      if(arr(it3,k32).lt.0)call phase2(jdiag(it3,k32))
      j6(j6c+1)=mp
      j6c=j6c+2
      j6(j6c)=mp
c
      if(npart.eq.1) then
        itmin=it2
        itmax=it3
        itl=max0(it3,it4)
        ith=ilast
      else
        itmin=min0(it2,it3)
        itmax=max0(it2,it3)
        itmn=min0(it1,it4)
        itmx=max0(it1,it4)
        itl=max0(itmin,itmn)
        ith=min0(itmax,itmx)
      endif
c
      tab1(mp,1)=itmin
      tab1(mp,2)=itmax
      jdiag(it2,1)=mp
      jdiag(it3,1)=mp
      jdiag(it2,3)=jj1
      arr(it2,3)=arr(it1,3)
      jdiag(it3,k32)=jj2
      arr(it3,k32)=arr(it4,k32)
c
      if(icross .eq. 1) then
        j7(j7c+1)=l2
        j7(j7c+2)=l4
        call phase2(l4)
        j7c=j7c+3
        j7(j7c)=mp
      else
        call phase2(jj2)
      endif
c
      itll=il(itl)
      if(npart.eq.1.and.icross.eq.1)itll=itll+1
      ithl=il(ith)
      if(npart.eq.2.and.icross.ne.1)ithl=ithl-1
c
    5 do 6 i=itll,ithl
        it=ih(i)
        ilp=i-k
        il(it)=ilp
        ih(ilp)=it
    6 continue
c
      if(ithl.ne.nbnode) then 
        itll=ithl+2 
        if(itll.le.nbnode) then
          ithl=nbnode
          k=2
          go to 5
        endif
      endif
c
      if(npart.ne.2)then
        tab1(jj1,1)=ih(1)
        tab1(jj1,2)=ih(nbnode-2)
      endif
c
      return
      end 
c
c
c     **************************************
      subroutine trdel(jj1,jj2,jj3,nbn,fail)
      implicit real*8(a-h,o-z)
c     **************************************
c
c  ***test for triangular delta.if not satisfied fail=.true.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,sumvar,cut,free
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/cutdig/cut
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
c
c
c
      if(sumvar(jj1).or.sumvar(jj2).or.sumvar(jj3))return
      if(nbn.gt.4)cut=.true.
      if(.not.free(jj1).and..not.free(jj2).and..not.free(jj3))then
        i1=j1(jj1)
        i2=j1(jj2)
        i3=j1(jj3)
        if(i1.lt.(iabs(i2-i3)+1).or.i1.gt.(i2+i3-1))fail=.true.
      endif
c
      return
      end 
c
c
c     *********************** 
      subroutine triang(fail) 
      implicit real*8(a-h,o-z)
c     *********************** 
c
c  ***reduces a triangle having one apex at either end of the axis of 
c  ***the flat diagram.
c  ***this introduces one 6j symbol and some phase factors .
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,sumvar
      integer arr,tab1
      character*6 name,namsub 
c
      common/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),
     + ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,ipartl,npart,
     + icross,nfree,itfree(m6j),nfin,nc,idummy(18)
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      common/nam/namsub
c
      data name/'triang'/
c
c
c
      namsub=name
      it1=npoint(1) 
      it2=npoint(2) 
      it3=npoint(3) 
      jwc=jwc+1
      kw(1,jwc)=jdiag(it3,2)
      kw(2,jwc)=jdiag(it2,3)
      kw(3,jwc)=jdiag(it3,1)
      if(arr(it3,1).gt.0)call phase2(kw(3,jwc))
      kw(4,jwc)=jdiag(it2,1)
      if(arr(it2,1).lt.0)call phase2(kw(4,jwc))
      k23=3
      if(it1.eq.ifirst)k23=2
      kw(5,jwc)=jdiag(it1,k23)
      kw(6,jwc)=jdiag(it3,3)
      call trdel(kw(1,jwc),kw(2,jwc),kw(5,jwc),nbnode,fail) 
      if(fail)go to 15
      if(arr(it3,3).gt.0)call phase2(kw(6,jwc))
      jt1=kw(5,jwc) 
      jdiag(it3,1)=jt1
      jdiag(it3,3)=kw(2,jwc)
      arr(it3,1)=arr(it1,k23) 
      arr(it3,3)=arr(it2,3)
c
      if(it1.ne.ifirst)then
        tab1(jt1,1)=it3
        tab1(jt1,2)=ih(nbnode-1)
        k12=1
      else
        tab1(jt1,1)=ih(2)
        tab1(jt1,2)=it3
        k12=2
      endif
c
      il3=il(it3)
c
      if(it1.ne.ilast)then
        il2=il(it2)-1
c
        do 2 i=2,il2
          it=ih(i)
          ilp=i-1
          il(it)=ilp
          ih(ilp)=it
    2   continue
      endif
c
      do 1 i=il3,nbnode
        it=ih(i)
        ilp=i-k12
        il(it)=ilp
        ih(ilp)=it
    1 continue
c
   15 return
      end 
c
c
c     ***************************************************** 
      subroutine var(jn,jns,jnc,jnsc,jbc,sumvar,mp,m,inver) 
      implicit real*8(a-h,o-z)
c     ***************************************************** 
c
c  ***test for variable character and put in jns if yes,and jn now
c  ***contains 0.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical sumvar(mp)
      dimension jn(jnc),jns(mangmp),inver(mp)
c
c
c
      jnsc=0
      if(jbc.ne.jnc)then
        jbbc=jbc+1
c
        do 1 i=jbbc,jnc
          i1=jn(i)
          if(sumvar(i1))then
            jnsc=jnsc+1
            j=inver(i1)
            jns(jnsc)=j
            jn(i)=m 
          endif
    1 continue
      endif
c
      return
      end 
c
c
c     ******************************
      subroutine way(l,ka,kb,ich,nb)
      implicit real*8(a-h,o-z)
c     ******************************
c
c  ***tests one step forward  if the way is free.first and second
c  ***arguments are interchanged or not according to ich=-1,or +1
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm)
c
      logical tabs
      integer arrow 
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/build/ial(m4trd),if1,if2,node
c
c
      k1=j23(l,ka)
      k2=j23(l,kb)
      nb=ial(k1)+ial(k2)-1
      if(nb)3,2,8
    2 nb1=ial(k1)-ial(k2)
      if(nb1)9,8,8
    3 call otherj(l,k1,l1,lc1,la)
      call otherj(l,k2,l2,lc2,lb)
      call neibor(lc1,i1,i2)
      call neibor(lc2,i3,i4)
      ji1=j23(l1,i1)
      ji2=j23(l1,i2)
      ji3=j23(l2,i3)
      ji4=j23(l2,i4)
      ia=ial(ji1)+ial(ji2)
      ib=ial(ji3)+ial(ji4)
      nbp=ib+ia+1
      nbm=ib-ia
      go to (8,4,5,4,6),nbp
    4 if(nbm)9,8,8
    5 if(nbm)9,6,8
    6 if(ji3.eq.if1.or.ji3.eq.if2.or.ji4.eq.if1.or.ji4.eq.if2)go to 9 
    8 ich=1
      go to 10
    9 ich=-1
   10 return
      end 
c
c
c     **************************
      subroutine zero(j,jz,fail)
      implicit real*8(a-h,o-z)
c     **************************
c
c  ***suppresses one line and two nodes of the unstructured graph
c  ***introduces  zeros in the triads j23.as a consequence the other
c  ***two arguments of the triad are put equal.if there was already
c  ***a zero in the triad which is changed,it is a special case.
c
c
      parameter(mangm=60,mtriad=12,m2trd=2*mtriad,m4trd=4*mtriad)
      parameter(m6j=20,msum=10,m3mngm=3*mangm,mangmp=2*(mangm/3))
c
      logical fail,tabs,free,sumvar,cut,nocut
c
      integer arrow 
c
      character*6 name
c
      common/zer/nzero,jzero(m6j)
      common/cutdig/cut
      common/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)
      common/keep/jkp(2,3),jarr(2,3),it2,it3,it5
      common/build/ial(m4trd),if1,if2,node
c
      common/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),
     +  lcol(mangm,2),tabs(m2trd),nbtr
c
      common/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),
     + j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp
c
      data name/'zero  '/
c
c
c
c
c
      nocut=.false. 
      nzero=0
c
      if(j.ge.1)then
       call otherj(0,jz,lin,lc,k1)
        i=nzero
        go to 8
      endif
c
      do 11 i=1,m
        if(j1(i).ne.1.or.free(i).or.ial(i).le.1)go to 11
        nzero=nzero+1
        jzero(nzero)=i
   11 continue
c
      nocut=.true.
      m=m+1
      j1(m)=1
      sumvar(m)=.false.
      free(m)=.false.
      if(nzero.eq.0)go to 7
      call printj(name,1)
      i=0 
    1 i=i+1
      jz=jzero(i)
      j=0 
   13 j=j+1
      lin=line(jz,j)
      if(tabs(lin))go to 2
      lc=lcol(jz,j) 
    8 call neibor(lc,l1,l2)
      jj1=j23(lin,l1)
      jj2=j23(lin,l2)
c
      if(jj1.eq.jj2)then
        j6c=j6c+1
        j6(j6c)=jj1 
        go to 10
      endif
c
      call delta(jj1,jj2,fail)
      if(fail)go to 7
c
      if(j1(jj1).ne.1.and.j1(jj2).ne.1)go to 15
      if(j1(jj1) .lt. j1(jj2))go to 15
      if(j1(jj1) .gt. j1(jj2))go to 19
c
      if(nzero.ne.0)then
        do 17 jjx=i,nzero
          jjz=jzero(jjx)
          if(jj1 .eq. jjz)go to 15
   18     if(jj2 .eq. jjz)go to 19
   17   continue
      endif
c
      go to 15
c
   19 jjz=jj2
      jj2=jj1
      jj1=jjz
c
   15 call otherj(lin,jj1,lo1,lco1,k1)
      call otherj(lin,jj2,lo2,lco2,k2)
      j9c=j9c+1
      j9(j9c)=jj1
      j23(lo2,lco2)=jj1
      line(jj1,k1)=lo2
      lcol(jj1,k1)=lco2
c
   10 if(arrow(lin,l1) .lt. arrow(lin,l2)) then
        call phase2(jj1)
      else
      if(arrow(lin,l1) .eq. arrow(lin,l2)) then
        arrow(lo1,lco1)=1
        arrow(lo2,lco2)=-1
      endif
      endif
c
      tabs(lin)=.true.
      nbtr=nbtr-1
      if(nbtr.eq.0)go to 7
c      if(lo1.ne.lo2)go to 2
c      l=6-lco1-lco2 
c      jt=j23(lo1,l) 
c      if(j1(jt).eq.1.and..not.free(jt))go to 2
c      call delta(jt,m,fail)
c      if(fail)go to 7
c      call neibor(l,l1,l2)
c      jtf=j23(lo1,l1)
c      if(arrow(lo1,l1) .lt. arrow(lo1,l2))call phase2(jtf)
c      j6c=j6c+1
c      j6(j6c)=jtf
c      nbtr=nbtr-1
c      tabs(lo1)=.true.
c      call otherj(lo1,jt,lin,lc,k)
c      go to 8
c     november 22 1989  
      if(lo1.eq.lo2)then
        l=6-lco1-lco2 
        jt=j23(lo1,l) 
        if(j1(jt).eq.1.and..not.free(jt))go to 2
        call delta(jt,m,fail)
        if(fail)go to 7
        nzero=nzero+1
        jzero(nzero)=jt
      end if
    2 if(j.eq.1)go to 13
c
      if (nbtr .ne. 0) then
        if(i.lt.nzero) go to 1
      endif
c
    7 call printj(name,4)
      if(nocut)cut=.false.
c
      return
      end 

