
 character*20 :: a='true'
 
 select case(a)
  case('true');
  write(0,*) a.eq.'true'
 end select
 end
 
 