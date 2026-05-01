
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2014, All rights reserved
!
! Code is under development, and the source code documentation and 
!  user's guide are presently in a preliminary form. The authors 
!  therefore request that users consult with the code authors, 
!  particularly before publishing ELITE results.   We also request 
!  that users report any errors in the documentation or code to the 
!  authors, and share any modifications made to the code or accompanying 
!  idl routines with the authors.  Please request permission from
!  P.B. Snyder or H.R. Wilson before distributing the code, as we wish 
!  to keep track of all users and changes to the code, and maintain a 
!  master version of the code to prevent unnecessary branching.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program elite

! main routine to solve for non-local plasma stability and
!  growth rate
!-----------------------------------------------------------------------
      use elite_data, only: init,vacuum,nxinterp,outfile,n1term,verbose
      implicit none
      integer ios,i

      call init                 !read in parameters from data file runname.in
      call surfcalc(1)   ! call surfcalc for reference surface
      call mercy(1)
      if (n1term == 2) then
         if (verbose .ge. 5) write(*,*) 'using new nustuff including 1/n terms, n1term=',n1term
         if (verbose .ge. 3) write(outfile,*) 'using new nustuff including 1/n terms'
         call nustuff(1)       ! setup matrix for ixinterp=1
      else
         if (verbose .ge. 2) write(*,*) 'using old nustuff, n1term=',n1term
         if (verbose .ge. 2) write(outfile,*) 'using old nustuff, n1term=',n1term
         call oldnustuff(1)
      endif
      if(vacuum) call rdvacfile
      do i=2,nxinterp
         call surfcalc(i)
         call mercy(i)
         if (n1term == 2) then
            call nustuff(i)
         else
            call oldnustuff(i)
         endif
      enddo

      if (verbose .ge. 5) write (outfile,*) 'call solveit'
      call solveit
!      close(65)

end program elite
