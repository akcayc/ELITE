
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2005, All rights reserved
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

program eliteeq

!----------------------------------------------------------------------
!  Equilibrium module for the ELITE stability code.
!    This version of ELITE reads in an equilibrium
!    from an EFIT file and maps it to the R,Z,B_p format
!    required by the main body of the ELITE, or it directly
!    reads equilibrium information from dskbal or dskgato
!    files produced by TOQ or CORSICA, or elite equilibrium 
!    files produced by HELENA or SCENE.  
!       P.B. Snyder and H.R. Wilson
!         building upon earlier code by R.L. Miller and H.R. Wilson

       use eliteeq_data, only: init,outfile,verbose

       call init
       call gendat
       call wrtdsk
       call testeq
       if (verbose .gt. 3) close(outfile)
!      stop
     
end program eliteeq
