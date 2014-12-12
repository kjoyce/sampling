roll_until_boxcar = function(last_roll,rolls) {
  this_roll = runif(1)<=1/6
  if( this_roll & last_roll ) {
    return(rolls+1)
  }
  else {
    roll_until_boxcar(this_roll,rolls+1) 
  }
}

(mean(replicate(10000,roll_until_boxcar(0,0))))
