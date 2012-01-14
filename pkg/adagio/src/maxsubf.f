C Maximal contiguous subset problem
      subroutine maxsubf(x, n, s, i1, i2)
      integer n, i1, i2, j1, j2
      double precision x(n)
      double precision s, ss

      ss = 0
      j1 = 1
      j2 = 1

      do 3000 i = 1,n 
      if (ss .gt. -x(i)) then
        ss = ss + x(i)
        j2 = i
        if (ss .gt. s) then
          s  = ss
          i1 = j1
          i2 = j2
        endif
      else
        ss = 0
        j1 = i+1
        j2 = i+1
      endif
 3000 continue

      return
      end
