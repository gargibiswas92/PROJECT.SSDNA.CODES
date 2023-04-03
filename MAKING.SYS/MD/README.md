# Instructions for MD:  
1. For each of the length of DNA we run simulations keeping the ends fixed, also by applying force at its both ends or at some intermediate point.  
2. In both the cases codes are more or less same except for the LD.f, where we employed the condition for keping it either fixed or not fixed.  
3. The MDWrapper.prefs is a bit different where we apply force.  
4. Also, when applying multiple forces at multiple points on the DNA we changed the external_force.f program as well.  
5. Each time we changed anything in any of the fortran programs inside the MD folder we comiled it using the command:  
   make -f MDmake_check (This version of MDmake we used).
