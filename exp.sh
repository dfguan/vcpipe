# this gives you an example how to run parallel lsf jobs, and wait for all of them end 

sleeps() {
	sleep $1
}

export -f sleeps

for i in 1 2 3
do
	let t="30 * $i"
	bsub -J$i -K sleeps $t & # if & is not appended here, the command will run serially, -K wait for job finished 
done
wait

echo done
