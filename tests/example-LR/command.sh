frequency=5.43320000
templatefile=i0-template
file=i3-nl-shg
simulatedtime=83 # in fs
periodicsteps=9  # odd number please

# From now on, I should not need to edit
sed "s/REPLACEFREQ/$frequency/g" $templatefile > $file
period=`python3 -c "import math;HBAR_eVfs=6.58211899E-1;print(HBAR_eVfs*2*math.pi/$frequency)"`
echo 'Frequency [eV] = ' $frequency
echo 'Period [fs]    = ' $period
text='           fs    # [NL] Simulation Time'

for i in `seq 1 1 $periodicsteps`; do 
	t=`python3 -c "print($simulatedtime+$period/$periodicsteps*($i-0))"`; 
	sed "s/NLtime=.*/NLtime=$t$text/g" $file > $file\_$i
done

postperiod_i=`echo 1+$periodicsteps|bc`
postperiod_f=`echo 3+$periodicsteps|bc`
for i in `seq $postperiod_i 1 $postperiod_f`; do 
	t=`python3 -c "print($simulatedtime+$period*($i-$periodicsteps+0.25))"`; 
	sed "s/NLtime=.*/NLtime=$t$text/g" $file > $file\_$i
done
