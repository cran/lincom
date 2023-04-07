      subroutine gridsch(mark,n,p,n1,sn,grdpt,coef,sp,
     &     ccf,ord,wk)
      integer n,p,n1,grdpt,ccf(p),ord(n1)
      double precision mark(n,p),sn,coef(0:p),sp,wk(n)

      integer r,i
      double precision mesh

      r=n1-ceiling(dble(n1)*sn)+1
      mesh=1.0d0/dble(grdpt)
      sp=-1.0d0

      do i=1,n1
         ord(i)=i
      enddo

      call gcore(mark,n,p,n1,coef,sp,mesh,r,grdpt,1,
     &     ccf,ord,wk)

      end

      recursive subroutine gcore(mark,n,p,n1,coef,sp,mesh,r,l1,q,
     &     ccf,ord,wk)
      integer n,p,n1,r,l1,q,ccf(p),ord(n1)
      double precision mark(n,p),coef(0:p),sp,mesh,wk(n)

      integer i,j,k
      double precision csp

      if(q .gt. p) then
      
         do i=1,n
            wk(i)=0.0d0
            do j=1,p
               wk(i)=wk(i)+dble(ccf(j))*mesh*mark(i,j)
            enddo
         enddo

         do i=2,n1
            j=i
            do while(j .gt. 1 .and.
     &           wk(ord(min(j-1,r))) .gt. wk(ord(j)))
               k=ord(j)
               ord(j)=ord(min(j-1,r))
               ord(min(j-1,r))=k
               j=min(j-1,r)
            enddo
         enddo
         k=ord(r)
         csp=0.0d0
         do i=n1+1,n
            if(wk(i) .lt. wk(k)) csp=csp+1.0d0
         enddo
         csp=csp/dble(n-n1)
         if(csp .gt. sp) then
            sp=csp
            coef(0)=wk(k)
            do i=1,p
               coef(i)=dble(ccf(i))*mesh
            enddo
         endif
         return
      endif

      k=1
      if(q .eq. p .and. l1 .gt. 0) k=2*l1
      do i=-l1,l1,k
         ccf(q)=i
         call gcore(mark,n,p,n1,coef,sp,mesh,r,l1-abs(i),q+1,
     &        ccf,ord,wk)
      enddo

      end
