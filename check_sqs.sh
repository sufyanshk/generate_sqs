#set time counter=0
t=0

#write in seconds how much time do you want your MCSQS to run
sqs_time=172800

while [ $t -lt $sqs_time ]
do
	#wait for 1 second in every loop
	sleep 1
	if grep -q Perfect_match 'bestcorr.out'; then
		#creat stopsqs file if the Perfect_match is found in bestcorr.out file
		touch stopsqs 
		#go out of the loop
		break
	else
		#else increase the timer
		t=$(( $t+1 ))
	fi
done
if [ ! -f "stopsqs" ]; then
	touch stopsqs
fi
