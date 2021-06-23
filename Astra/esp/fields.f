           subroutine field_c

         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'

         common /com_mf/ bpol(nrp,ntp),btor(nrp,ntp),
     &                   br(nrp,ntp),bz(nrp,ntp),
     &                   rc(nrp,ntp),zc(nrp,ntp)

!magnetic field at cells

         do i=1,iplas-1
          do j=2,nt-1

           bp_ij=(psi(i+1,j)-psi(i,j))/st(i,j)
           bp_ij1=(psi(i+1,j+1)-psi(i,j+1))/st(i,j+1)

         if(i.ne.1) then
           bp_2= (
     &           bp_ij**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j) )+
     &           bp_ij1**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j) )
     &           )/vol(i,j)
         else
           bp_2= (
     &           bp_ij**2*(                     vol2(i,j)/sin2(i,j) )+
     &           bp_ij1**2*( vol3(i,j)/sin3(i,j)                     )
     &           )/vol(i,j)
         endif
           bpol_ij=dsqrt(bp_2)
           btor_ij=f(i)*s(i,j)/vol(i,j)

           bpol(i,j)=bpol_ij
           btor(i,j)=btor_ij
           btot_ij=dsqrt(bp_2+btor_ij**2)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          r12=(r1+r2)*0.5d0
          r23=(r3+r2)*0.5d0
          r34=(r3+r4)*0.5d0
          r14=(r1+r4)*0.5d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)
          z0=(z1+z2+z3+z4)*0.25d0
        
          z12=(z1+z2)*0.5d0
          z23=(z3+z2)*0.5d0
          z34=(z3+z4)*0.5d0
          z14=(z1+z4)*0.5d0

          dr=r34-r12 
          dz=z34-z12 
          dr=dr/dsqrt(dr**2+dz**2)
          dz=dz/dsqrt(dr**2+dz**2)

          br(i,j)=bpol_ij*dr
          bz(i,j)=bpol_ij*dz

          rc(i,j)=r0
          zc(i,j)=z0

          enddo
         enddo

         do i=1,iplas-1
           bpol(i,1)=bpol(i,nt-1)
           btor(i,1)=btor(i,nt-1)
           br(i,1)=br(i,nt-1)
           bz(i,1)=bz(i,nt-1)
           rc(i,1)=rc(i,nt-1)
           zc(i,1)=zc(i,nt-1)
         enddo

        open(1,file='fields.wr',form='formatted')
         write(1,*) iplas-1,nt-1
         write(1,*) ((rc(i,j),i=1,iplas-1),j=1,nt-1)
         write(1,*) ((zc(i,j),i=1,iplas-1),j=1,nt-1)
         write(1,*) ((br(i,j),i=1,iplas-1),j=1,nt-1)
         write(1,*) ((bz(i,j),i=1,iplas-1),j=1,nt-1)
         write(1,*) ((btor(i,j),i=1,iplas-1),j=1,nt-1)
         write(1,*) ((bpol(i,j),i=1,iplas-1),j=1,nt-1)
        close(1)

           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


