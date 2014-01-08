

      
      subroutine bilin(nx,ny, nnew, pNew, x0, y0, dx, dy, dat)
      implicit double precision (a-h,o-z)
      double precision pNew(nnew,3), dat(nx,ny), xs(nx),ys(ny)
c      write(99,*) nx,ny
      do i = 1,nx
        xs(i) = x0+(i-1)*dx 
      enddo
      do i = 1,ny
        ys(i) = y0+(ny-i+1)*dy
      enddo
     
      do i = 1,nnew
c        write(*,*), i, nnew, nx, ny
        x = pNew(i,1)
        y = pNew(i,2)
        ix = floor((x-x0)/dx)+1
        iy = ny-floor((y-y0)/dy)

        write(*,'(2I7,8f12.5)') ix,iy,x,y,x0,y0,ys(1),ys(ny),
     1                       xs(1),xs(nx)
        x1 = xs(ix)
c        write(*,*) ix,ix+1,x1
        x2 = xs(ix+1)
c        write(*,*) ix+1,iy,x2
        y1 = ys(iy)
c        write(*,*) iy,iy+1,y2
        y2 = ys(iy+1)
c        write(*,*) ix,ix+1,iy,dx,dy
        f1 = (x2-x)/dx *dat(ix,iy) + (x-x1)/dx * dat(ix+1,iy)
c        write(*,*) f1
        f2 = (x2-x)/dx *dat(ix,iy+1) + (x-x1)/dx * dat(ix+1,iy+1)
c        write(*,*) f2
        fnew = (y-y2)/dy*f1 + (y1-y)/dx*f2
c        write(*,*) fnew
        pNew(i,3) = fnew
        write(*,'(3I5,8f10.4)')i,ix,iy,x,y,y1,y2,dat(ix,iy),f1,f2,fnew
      enddo
      
      end
   