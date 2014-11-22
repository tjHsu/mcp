       program ljmd2Dvv

       implicit real*8(a-h,o-z)

c **** Specify parameters:
c **** Number of trials
       parameter(Ncycle=5000, Nshots = 100, Nprint=Ncycle/Nshots)
c **** Temperature to set the initial velocities
       parameter(Tinit=2.5 )
c **** Box size
       parameter(bx=15.0d0, by=15.0d0)
c **** Number of particles
       parameter(Npart=100)
c **** Time step
       parameter(dt = 0.001, dt2=dt/2, dtsq2 = dt2*dt)
c **** Random number generator's seed
       parameter(iseed=1234567)

c **** Set up arrays 
c **** Coordinates
       dimension x(2,Npart) 
c **** Velocities 
       dimension v(2,Npart) 
c **** Forces (accelerations) 
       dimension f(2,Npart) 

       rtemp = sqrt(Tinit)

c **** set counters
c **** Running Averages    
c **** Energy
       avE  = 0.0d0
       avE2 = 0.0d0
c **** Potential Energy
       avP  = 0.0d0
       avP2 = 0.0d0
c **** Shifted Potential Energy
       avPc  = 0.0d0
       avPc2 = 0.0d0
c **** Kinetic Energy
       avK  = 0.0d0
       avK2 = 0.0d0

c *** set seed for the random numbers generator (notice, that the name for srand
c *** may be different on different architectures )
       call srand(iseed)

c **** Set up the initial lattice configuration, 

c***       Nlat =  int(sqrt(real(Npart)))
       Nlat=10    
       dlatx = bx/Nlat 
       dlaty = by/Nlat 

       ipart = 0
       do ix = 1,Nlat
       do iy = 1,Nlat
        ipart=ipart+1
        x(1,ipart) = ix * dlatx
        x(2,ipart) = iy * dlaty
c *** Initial velocities according to Maxwell 
c *** from the Gaussian Distribution
        v(1,ipart) =(rand()+rand()+rand()+
     >               rand()+rand()+rand()+
     >               rand()+rand()+rand()+
     >               rand()+rand()+rand())-6.
        v(2,ipart) =(rand()+rand()+rand()+
     >               rand()+rand()+rand()+
     >               rand()+rand()+rand()+
     >               rand()+rand()+rand())-6.
        v(1,ipart)=v(1,ipart)*rtemp
        v(2,ipart)=v(2,ipart)*rtemp
       enddo
       enddo

c *** Get rid of the center of mass motion
      px = 0.d0
      py = 0.d0
      do ipart = 1, Npart
       px = px + v(1,ipart)
       py = py + v(2,ipart)
      enddo
      px = px /Npart
      py = py /Npart

      do ipart = 1, Npart
        v(1,ipart) = v(1,ipart) - px
        v(2,ipart) = v(2,ipart) - py
      enddo

c *** Start MD cycle 

      time = 0.0d0

c *** open trajectory file and save initial coordinates
c *** Make it's status new, so it will not overwrite 
c *** the previous

      open(unit=7,file="ljmd2Dvv.xyz",status="unknown")

c *** Save the snapshot (each ends with the empty record)
      write(7,11) time
11    format(e12.5)
      do ipart = 1, Npart
      write(7,12) x(1,ipart), x(2,ipart) 
12    format(2(1x,e12.5))
      enddo
      write(7,*)
 

c *** Calculate the current forces to start
       call energy(x,f,P,w,Pc,bx,by)

       do icycle = 1, Ncycle

c *** Advance coordinates and velocities
       do ipart = 1, Npart
         x(1,ipart) = x(1,ipart) + dt* v(1,ipart) + dtsq2* f(1,ipart)
         x(2,ipart) = x(2,ipart) + dt* v(2,ipart) + dtsq2* f(2,ipart)
         v(1,ipart) = v(1,ipart) + dt2* f(1,ipart)
         v(2,ipart) = v(2,ipart) + dt2* f(2,ipart)
       end do 
c **** Get new energy and forces
       call energy(x,f,P,w,Pc,bx,by)

c **** Complete advancing velocities (required?) and get Kinetic energy 
       Ekin = 0.d0
       do ipart =1, Npart
         v(1,ipart) = v(1,ipart)+dt2 * f(1,ipart)
         v(2,ipart) = v(2,ipart)+dt2 * f(2,ipart)
         Ekin = Ekin +v(1,ipart)**2 + v(2,ipart)**2
       end do
	   Ekin = Ekin*0.5
       E = Pc + Ekin
       Press = w + 2.d0/2.d0 * Ekin
       Temp = Ekin / Npart
c **** increment averages
c **** Total energy
         avE = avE + E 
         avE2= avE2+ E * E
c **** Potential Energy
         avP = avP + P 
         avP2= avP2+ P * P
c **** Shifted Potential Energy
         avPc = avPc + Pc 
         avPc2= avPc2+ Pc * Pc
c **** Kinetic Energy
         avK = avK + Ekin 
         avK2= avK2+ Ekin * Ekin 
c **** Pressure 
         avPress = avPress + Press 
         avPress2= avPress2+ Press * Press
c **** Temperature 
         avTemp = avTemp + Temp 
         avTemp2= avTemp2+ Temp * Temp
         
         
         time = time + dt

c **** Check if it's time to save print out 
         if(mod(icycle,Nprint).eq.0) then
c *** Save the snapshot (each ends with the empty record)
           write(7,13) time
13         format(e12.5)
           do ipart = 1, Npart
             write(7,14) x(1,ipart), x(2,ipart) 
14           format(2(1x,e12.5))
           enddo
           write(7,*)
 
c **** Print out instant values
             write(6,1) time, Ekin, Pc, E
1            format(" time = ",e12.5, " Instant Kin.energy=",e12.5
     >    ," Pot.E = ",e12.5," Total = ",e12.5)

c **** Running averages
             Erunav = avE / icycle 
             E2runav = avE2 / icycle 
             E2runav = sqrt((E2runav-Erunav*Erunav))
             Prunav = avP / icycle 
             P2runav = avP2 / icycle 
             P2runav = sqrt((P2runav-Prunav*Prunav))
             Pcrunav = avPc / icycle 
             Pc2runav = avPc2 / icycle 
             Pc2runav = sqrt((Pc2runav-Pcrunav*Pcrunav))
             eKrunav = avK / icycle 
             eK2runav = avK2 / icycle 
             eK2runav = sqrt((eK2runav-eKrunav*eKrunav))
             Pressrunav = avPress / icycle 
             Press2runav = avPress2 / icycle 
             Press2runav = sqrt((Press2runav-Pressrunav*Pressrunav))
             Temprunav = avTemp / icycle 
             Temp2runav = avTemp2 / icycle 
             Temp2runav = sqrt((Temp2runav-Temprunav*Temprunav))
             write(6,2) Erunav, E2runav, 
     >                  Prunav,P2runav,
     >         Pcrunav,Pc2runav,
     >         eKrunav,eK2runav,
     >         Temprunav,Temp2runav,
     >         Pressrunav,Press2runav
2      format(
     >       " Av. Tot. En.=",e12.5," M.sq.Dev=",e12.5/
     >       " Av. Pot. En.=",e12.5," M.sq.Dev=",e12.5/  
     >       " Av. Shft.En.=",e12.5," M.sq.Dev=",e12.5/
     >       " Av. Kin. En.=",e12.5," M.sq.Dev=",e12.5/  
     >       " Temperature =",e12.5," M.sq.Dev=",e12.5/
     >       " Pressure    =",e12.5," M.sq.Dev=",e12.5/)
         endif
c **** End loop over the MD cycles
       enddo

c *** That's all folks
c ***  close the trajectory file
        close(7)

        end

       subroutine energy(x,f,P,w,Pc,bx,by)
       implicit real*8(a-h,o-z)
c **** Cut-off radius
       parameter(rcut=3, r2cut=rcut*rcut  )
c **** Number of particles
       parameter(Npart=100)  
c ***
c **** Set up arrays
c **** Coordinates
       dimension x(2,Npart)
c **** Forces (accelerations)
       dimension f(2,Npart)

       bxinv = 1.0d0/bx
       byinv = 1.0d0/by
       ncut = 0
       w = 0.d0
       P   = 0.0d0
       do i = 1,Npart
        f(1,i) = 0.d0
        f(2,i) = 0.d0
       enddo

       do i = 1,Npart-1
       fxi = f(1,i)
       fyi = f(2,i)
       x1i=x(1,i)
       x2i=x(2,i)
       do j = i+1,Npart
        dx = x1i-x(1,j) 
        dy = x2i-x(2,j) 
c **** Here are the Periodic Boundary Conditions
        dx=dx - bx * anint(dx * bxinv)
        dy=dy - by * anint(dy * byinv)
        r2ij=dx*dx+dy*dy
c ****  Cut off. Go to the next j if too far away
        if(r2ij.gt.r2cut) goto 99
        sr2=1.0d0/r2ij
        sr6 =sr2 * sr2 * sr2
        sr12=sr6 * sr6
        vij = sr12 - sr6 
        P = P + vij
        wij = vij + sr12
        w = w + wij
        fij = wij * sr2 
        fxij = fij * dx
        fyij = fij * dy
        fxi = fxi + fxij         
        fyi = fyi + fyij
        f(1,j) = f(1,j) - fxij
        f(2,j) = f(2,j) - fyij
        ncut = ncut + 1       
99      continue
       enddo
       f(1,i) = fxi
       f(2,i) = fyi
       enddo
c **** Calculate shifted potential
       sr2 = 1.0d0 / r2cut
       sr6 = sr2 * sr2 * sr2
       sr12=sr6 * sr6
       vij = sr12 - sr6
       Pc = P - ncut * vij

c **** Multiply force, virial and energy by factors

       do i = 1,Npart
        f(1,i) = f(1,i) * 24.0d0
        f(2,i) = f(2,i) * 24.0d0
       enddo
       P = P* 4.0d0
       Pc= Pc*4.0d0
       w = w *24.0d0 / 2.0d0


       return
       end
