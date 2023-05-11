#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include <Remove Points>

// # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # # 	               # # # # # # # # # # # # 
// # # # # # # # 	S P I K E S    # # # # # # # # # # # # 
// # # # # # # # 	               # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

Threadsafe FUNCTION/DF ReturnSpikeTimeWave(WAVE InWave,[VARIABLE V_Threshold, VARIABLE MinISI, VARIABLE RateThresh4Onset, VARIABLE Bool_RatePositive])
// this returns spike times in the time unit of the wave
// spike time is time of crossing Threshold, which defaults to 0 mV

// the diversity of spike shapes requires a detection criterion that is more sophisticated than a 
// simple crossing of a threshold from below
// main problem was noise in very slowly falling spikes

// the optional parameters work as follows:
// V_Threshold  - crossing this from below makes a sample a candidate for spike time
// 					in the absence of other optional parameters thats the only criterion and the default for the
// 					threshold value is 0 mV
// MinISI  	    - a candidate spike time is rejected if it occurs sooner than MinISI (seconds) after the last accepted spike
// Bool_RatePositive  	    - a candidate spike time is rejected if it occurs
// 					in a region (2 pnts before + 2pnts after) in which the rate of rise is on average negative 

// InWave			// contains V_mem vs time

IF (ParamIsDefault(V_Threshold ) )
	V_Threshold = 0
ENDIF
IF (ParamIsDefault(MinISI ) )
	MinISI = 0.0002 //(5 kHz instantaneous rate)
ENDIF
IF (ParamIsDefault(Bool_RatePositive ) )
	Bool_RatePositive = 0
ENDIF
IF (ParamIsDefault(RateThresh4Onset))
	RateThresh4Onset = 30 // V/s
ENDIF
// contruct name of spike time wave 
// if there is any "_" in the inWave's name, replace the ending after the last "_"
// with "ST":
		VARIABLE		N_underScore=ItemsInList(NameOfWave(InWave),"_")
		STRING		OutName

		IF (N_underScore >0)
			OutName =StringFromList(N_underScore-1, NameOfWave(InWave),"_")
			// now it contains the name AFTER the last underscore; next, find the pos of that last underscore
			// and take the name until (excl.) this underscore
			OutName=NameOfWave(InWave)[0,FindListItem(OutName, NameOfWave(InWave)+"_","_")-1]
		ELSE	
			OutName=NameOfWave(InWave)+"_"	// take full name if no underscores are found 
														// end hence the composition of the original
														// names is unclear
		ENDIF

		OutName=OutName+"ST"
		
	// for multithreading folder return empty folder in case error occured
	DFREF dfWRONG= NewFreeDataFolder()
	// for multithreading folder return correct folder
	DFREF dfFree= NewFreeDataFolder()
	SetDataFolder dfFree		

// quick check: any threshold crossings?
	Wavestats/Q/M=1 InWave
	IF (V_max<V_Threshold)
		MAKE/O/N=(0) $OutName
		Print "no threshold crossings in " +OutName
		Return dfWRONG
	ENDIF


// quick health check: dV/dt should not be larger than say 2000 V/s
	Differentiate InWave/D=Dummy;
	IF (wavemax(Dummy)>4000 ||  wavemin(Dummy)<-4000 ) 
// there are probably artefacts
		WAVESTATS/Q Dummy
		KillWaves/Z Dummy
		print 0,"Presumed artefact in "+GetWavesDataFolder(InWave,2)
		print  GetWavesDataFolder(InWave,2), V_min,"@",V_minloc,V_Max,"@",V_maxloc
		Return dfWRONG
	ENDIF
	KillWaves/Z Dummy
	FindLevels /DEST=$OutName /EDGE=1 /M=(MinISI) /Q  InWave, V_Threshold 
	// use P, so the result comes in index rather than x-scaling
	WAVE		ST=$OutName
	VARIABLE		dt=DimDelta(InWave,0)
	IF (Bool_RatePositive)
	// test whether average trend from 3 pnts before to 3 pnts after spike candidate time 
	// is upward. This is equivalent to just testing whether the point 2 later is above the point 2 earlier
		VARIABLE tt,nst=DimSize(ST,0)
		FOR (tt=0;tt<nst;tt++)
			IF (!(InWave[x2pnt(Inwave,ST[tt])+4]> InWave[x2pnt(Inwave,ST[tt])-4]) )
				DeletePoints tt,1,ST
				nst--
			ENDIF
		ENDFOR
	ENDIF
	// now call spikestatistics
	// this is optional, but makes sense to do right now
	// if there are any spikes at all
	IF (DimSize(ST,0)>0)
		SpikeStats(InWave, ST,RateThresh=RateThresh4Onset)
	ENDIF
	//make a note about the threshold that was used:
	Note ST,"Threshold:"+num2str(V_Threshold)+" "+Waveunits(InWave,1)
	return dfFree
END

FUNCTION ReturnSpikeTimes(WAVE InWave,[VARIABLE V_Threshold, VARIABLE MinISI, VARIABLE RateThresh4Onset, VARIABLE Bool_RatePositive])
// this returns spike times in the time unit of the wave
// spike time is time of crossing Threshold, which defaults to 0 mV

// the diversity of spike shapes requires a detection criterion that is more sophisticated than a 
// simple crossing of a threshold from below
// main problem was noise in very slowly falling spikes

// the optional parameters work as follows:
// V_Threshold  - crossing this from below makes a sample a candidate for spike time
// 					in the absence of other optional parameters thats the only criterion and the default for the
// 					threshold value is 0 mV
// MinISI  	    - a candidate spike time is rejected if it occurs sooner than MinISI (seconds) after the last accepted spike
// Bool_RatePositive  	    - a candidate spike time is rejected if it occurs
// 					in a region (2 pnts before + 2pnts after) in which the rate of rise is on average negative 

// InWave			// contains V_mem vs time

IF (ParamIsDefault(V_Threshold ) )
	V_Threshold = 0
ENDIF
IF (ParamIsDefault(MinISI ) )
	MinISI = 0.0002 //(5 kHz instantaneous rate)
ENDIF
IF (ParamIsDefault(Bool_RatePositive ) )
	Bool_RatePositive = 0
ENDIF
IF (ParamIsDefault(RateThresh4Onset))
	RateThresh4Onset = 30 // V/s
ENDIF
// contruct name of spike time wave 
// if there is any "_" in the inWave's name, replace the ending after the last "_"
// with "ST":
		VARIABLE		N_underScore=ItemsInList(NameOfWave(InWave),"_")
		STRING		OutName

		IF (N_underScore >0)
			OutName =StringFromList(N_underScore-1, NameOfWave(InWave),"_")
			// now it contains the name AFTER the last underscore; next, find the pos of that last underscore
			// and take the name until (excl.) this underscore
			OutName=NameOfWave(InWave)[0,FindListItem(OutName, NameOfWave(InWave)+"_","_")-1]
		ELSE	
			OutName=NameOfWave(InWave)+"_"	// take full name if no underscores are found 
														// end hence the composition of the original
														// names is unclear
		ENDIF
		OutName=OutName+"ST"
		
		
// quick check: any threshold crossings?
	Wavestats/Q/M=1 InWave
	IF (V_max<V_Threshold)
		MAKE/O/N=(0) $OutName
		Return 0
	ENDIF


// quick health check: dV/dt should not be larger than say 2000 V/s
	Differentiate InWave/D=Dummy;
	IF (wavemax(Dummy)>4000 ||  wavemin(Dummy)<-4000 ) 
// there are probably artefacts
		WAVESTATS/Q Dummy
		KillWaves/Z Dummy
		print 0,"Presumed artefact in "+GetWavesDataFolder(InWave,2)
		print  GetWavesDataFolder(InWave,2), V_min,"@",V_minloc,V_Max,"@",V_maxloc
		Return 0
	ENDIF
	KillWaves/Z Dummy
	FindLevels /DEST=$OutName /EDGE=1 /M=(MinISI) /Q  InWave, V_Threshold 
	// use P, so the result comes in index rather than x-scaling
	WAVE		ST=$OutName
	VARIABLE		dt=DimDelta(InWave,0)
	IF (Bool_RatePositive)
	// test whether average trend from 3 pnts before to 3 pnts after spike candidate time 
	// is upward. This is equivalent to just testing whether the point 2 later is above the point 2 earlier
		VARIABLE tt,nst=DimSize(ST,0)
		FOR (tt=0;tt<nst;tt++)
			IF (!(InWave[x2pnt(Inwave,ST[tt])+4]> InWave[x2pnt(Inwave,ST[tt])-4]) )
				DeletePoints tt,1,ST
				nst--
			ENDIF
		ENDFOR
	ENDIF
	// now call spikestatistics
	// this is optional, but makes sense to do right now
	// if there are any spikes at all
	IF (DimSize(ST,0)>0)
		SpikeStats(InWave, ST,RateThresh=RateThresh4Onset)
	ENDIF
	//make a note about the threshold that was used:
	Note ST,"Threshold:"+num2str(V_Threshold)+" "+Waveunits(InWave,1)
	return DimSize(ST,0)
END

Threadsafe FUNCTION SpikeStats(InWave_V, InWave_ST,[RateThresh])
// call it like this SpikeStats(InWave_V, InWave_ST,RateThresh=30)
WAVE	InWave_V, InWave_ST		// voltage and spike-time wave
VARIABLE	RateThresh				// threshold value for dV/dt to define onset
// returns threshold voltage and onset rapidness for every spike

// take spike time and voltage trace (assuming units are idetical (seconds) and volts
// 1. health test to check whether the spike times are plausible given the voltage trace
// - the spike times are used, intead of a spike definition based on the voltage tracce
//   to allow for a subset of spikes to be used (instead of all)
// 2. append to the wave note of the spike time the coefficient of variation
// 3. go back from each spike to find onset as defined by the rate threshold variable
//    - if not given, it defaults to 30 V/s
// 4. create seven waves with an entry for every spike in ST
	//		one for threshold voltage "_THR" 
	//		one for the onset rapidness, the phase slope at the threshold value
	//		one for the inter-spike interval (first spike entry is NaN)
	// 	one for the maximal voltage
	// 	one for the maximal rate of rise (dV/dt)
	// 	one for the voltage at which maximal dV/dt occurs
	//		one for full width at half maximum (half way up between threshold and peak potential)
	// 	one for the minimal voltage reached after a spike (after hyperpolarization)
	// 	one for the time of time at which the RateThresh was crossed
	
// 5. append to Wavenote of ST average values for most of the above properties plus 
//		the local variability of the spike train plus  the coefficient of variation of the ISI

		VARIABLE		dt=DimDelta(InWave_V,0)
		
// health test: Spike times should define a narrow set of voltages in the voltage trace
// (assuming that a VOLTAGE threshold rather than a rate threshold was used as definition
		Duplicate /O InWave_ST, Dummy
		WAVE 	dummy
		dummy[]=InWave_V[x2pnt(InWave_V,InWave_ST[p])]
		WAVESTATS/Q dummy
		VARIABLE	expectedVjitter=1500*dt/2*2 // dV/dt max * dt/2  * safety factor of 1.5
		IF ( (V_max-V_min) > expectedVjitter)
			Print 0,"The range of threshold voltages related to the spike times\ris "+num2str((V_max-V_min)*1000)+" mV.\rThat seems too big, aborting..."
			KillWaves dummy
			Return -1
		
		ENDIF

	// decompose name of inwaves to create new names: remove ending after last"_"
		VARIABLE		N_underScore=ItemsInList(NameOfWave(InWave_V),"_")
		STRING		SlopeName, ThreshName, VmaxName, mxRateName, minRateName, VofMaxRateName, ISIName,FWHMName, AHPName, AHPdelName, ThreshCrossTime
		IF (N_underScore >0)
			slopeName =StringFromList(N_underScore-1, NameOfWave(InWave_V),"_")
			// now it contains the name AFTER the last underscore; next, find the pos of that last underscore
			// and take the name until (excl.) this underscore
			SlopeName=NameOfWave(InWave_V)[0,FindListItem(SlopeName, NameOfWave(InWave_V)+"_","_")-1]
		ELSE	
			SlopeName=NameOfWave(InWave_V)	// take full name if no underscores are found 
														// end hence the composition of the original
														// names is unclear
		ENDIF
		ThreshName=SlopeName+"VTHR"
		VmaxName=SlopeName+"Vmx"
		mxRateName=SlopeName+"mxRt"
		VofMaxRateName=SlopeName+"VmxRt"
		minRateName=SlopeName+"minRt"
		FWHMName=SlopeName+"FWHM"
		ISIName=SlopeName+"ISI"
		AHPName=SlopeName+"AHP"
		AHPdelName=SlopeName+"AHPdel"
		ThreshCrossTime=SlopeName+"SonsetT"
		SlopeName=SlopeName+"RPD"

		Duplicate /O InWave_ST, $SlopeName, $ThreshName, $VmaxName, $mxRateName, $VofMaxRateName, $minRateName, $ISIName, $FWHMName, $AHPName,$AHPdelName, $ThreshCrossTime
		WAVE	slope=$SlopeName
		WAVE	thresh=$ThreshName
		WAVE	Vmx=$VmaxName
		WAVE	mxRt=$mxRateName
		WAVE	VmxRt=$VofMaxRateName
		WAVE	minRt=$minRateName
		WAVE	FWHM=$FWHMName
		WAVE	ISI=$ISIName  
		WAVE	AHP=$AHPName  
		WAVE	AHPdel=$AHPdelName  
		WAVE	STX=$ThreshCrossTime  
		slope	=NaN
		thresh=NaN
		Vmx	=NaN
		mxRt	=NaN
		minRt	=NaN
		VmxRt	=NaN
		FWHM	=NaN
		AHP 	=NaN
		AHPdel 	=NaN
		STX 	=NaN
		SetScale d 0,0,"1/s", slope
		SetScale d 0,0,"V", thresh, VmxRt, Vmx, AHP
		SetScale d 0,0,"V/s", mxRt, minRt
		SetScale d 0,0,"s", ISI
		SetScale d 0,0,"s", FWHM, AHPdel
		SetScale d 0,0,"s", STX
		// make sure all waves that were created have index scaling
		SetScale/P x,0,1, slope,thresh, VmxRt, Vmx,ISI, mxRt, minRt, FWHM, AHP, AHPdel, STX

// get local variation coefficient and append to wave note
		Note/K 	InWave_ST
		Note/NOCR InWave_ST, "local variation of spike times (LV):"
		Note/NOCR InWave_ST, num2str(ReturnLocalCV(InWave_ST))
// 10/2016 aNNe:  decided to append coefficient of variation as well.
		Note InWave_ST, "coef. of variation of spike times (CV):"
		Note/NOCR InWave_ST, num2str(ReturnCV(InWave_ST))

// find threshold as defined by rate threshold 
		IF (ParamIsDefault(RateThresh))
	           	RateThresh = 30 //(V/s)
		ENDIF

		VARIABLE		nSpikes=DimSize(InWave_ST,0), sp, sPnt// Pnt number in Voltage wave where spike was detected
		VARIABLE	 	nPntsBefore, nPntsAfter, upsamplefac=max(1,round(100000*dt))// stretch of data before and after that is copied for analysis
																			// factor of how much the temporal sampling is increased
		VARIABLE		nPntsInDummy, SlopeMaxWasApplied=0, SlopeMax=3000000
		FOR (sp=0; sp< nSpikes; sp +=1)
			
			// copy a stretch of voltages around the detection point
			// use a max of 2 ms before but no more than the distance to the previous spike
			// to define how much time after the spike is taken, determine 
			// i) where threshold (defined by rate) is reached again (from above) 
			// ii) when next spike comes 
			
			// if threshold is reached again, take 5 ms after this (if next spike appears later) 
			// else take time halfway to next spike 
			
			// in order to find rate threshold, first take shorter part of spike and study it in detail --> 2 ms before, 2 ms after 
			
			IF (sp >0)	// not the first spike)
				nPntsBefore=min(0.002, InWave_ST[sp]-InWave_ST[sp-1])/dt
			ELSE
				nPntsBefore=0.002/dt
			ENDIF
			IF (sp+1 < nSpikes)	// not the last spike)
				nPntsAfter=min(0.002, InWave_ST[sp+1]-InWave_ST[sp])/dt
			ELSE
				nPntsAfter=0.002/dt
			ENDIF
			
			// calculate slope
			// check for drop below threshold slope (going backwards from voltage threshold crossing

			
			sPnt=x2pnt(InWave_V,InWave_ST[sp])
			DUPLICATE/O/R=[sPnt-nPntsBefore,sPnt-1+nPntsAfter] InWave_V, dummy
			// there is trouble if the first spike occurs less than 2 ms after voltage trace start
			
			// at this point upsample dummy and do all further analysis on the upsampled versions
			// currently a sinc reconstruction is used, this is 
			nPntsInDummy=DimSize(dummy,0)
			Duplicate/O dummy, dummy_US
			Resample/up=(upsamplefac) dummy_US	// this does a sinc reconstruction
			//			Interpolate2/T=2/N=((nPntsInDummy-1)*upsamplefac+1)/E=2/Y=dummy_US dummy

			WAVE	dummy_US			// upsampled version covering the exact same time
			// create rate
			// it turns out that differentiating the upsampled data gives a more strongly oscillating/fluctuating rate
			// it is better to first differentiate the original and then apply upsampling
			
			Differentiate /METH=0 dummy  /D=dummy_dV 
			Resample/up=(upsamplefac) dummy_dV
			WAVE	dummy_dV
			SetScale d 0,0,"V/s", dummy_dV
			
			
			// now find dV/dt max and voltage at which it occurs
			WAVESTATS/Q/M=1 dummy_dV
			mxRt[sp]	=	V_max
			VmxRt[sp]=	dummy_US[V_maxRowloc]
			
			// start searching for threshold BEFORE the peak location
			  
			FindLevel /EDGE=1 /P/Q/R=[V_maxRowloc,0] dummy_dV,RateThresh
			IF (V_Flag==0)	// level crossing found
			   // V_LevelX now holds the interpolated point (index) of threshold crossing	
			  		   
			   thresh[sp]=dummy_US[floor(V_LevelX)]+(dummy_US[ceil(V_LevelX)]-dummy_US[floor(V_LevelX)])*(V_LevelX-floor(V_LevelX))
				// linear interpolation finds estimate of voltage at which rate threshold was crossed
				
				// to get slope, take differentiated voltage and determine the slope d (dV/dt) /dV
				// at the time of the rate threshold crossing
				Differentiate /METH=0 dummy_dV	/X=dummy_US /D=dummy_phsSlp
				slope[sp]=dummy_phsSlp[V_LevelX]
				
				// replace unbelievable slopes with NaN
				IF (slope[sp]<0)
					slope[sp]=NaN
				ENDIF
				IF (slope[sp]>SlopeMax) // carefully! hard coded
					slope[sp]=NaN
					SlopeMaxWasApplied+=1
				ENDIF
				
				// get the time at which this occurs
				STX[sp]= pnt2x(dummy_dV,floor(V_LevelX)) + dt/upsamplefac*(V_LevelX-floor(V_LevelX))

				
			ELSE
				thresh[sp]=NaN
				slope[sp]=NaN
				STX[sp]= NaN
			ENDIF
			
			
			// now find peak voltage
			WAVESTATS/Q/M=1 dummy_US
			Vmx[sp]	=	V_max
			
			
			Variable HM, FWHMval=0	//half maximum and full width at half maximum
			// now use the existing thresh and Vmx values to define "half max" and then find the duration between up and down crossing
			// the parameters used are ok for cortical pyramidal cells at and above room temp
			// they might not work for other cases and in particular not for very small inter spike intervals
			VARIABLE 	minAmp		// max voltage to be considered for after hyper polarization minimum finding
			VARIABLE 	tAfter		// max time after spike to be taken into account
			VARIABLE 	tEnd= pnt2x(Inwave_V,DimSize(InWave_V,0)-1)		// end o input trave (time of last point
			VARIABLE	tRecross	// time at which threshold is crossed again
			IF (thresh[sp]!=thresh[sp]) // this is the case if threshold was could not be determined and is set to NaN
				FWHMval = NaN
				minAmp = -0.02		// -20 mV as a default
			ELSE 
				//define half maximum to search full width at half max later
				HM=thresh[sp]+0.5*(Vmx[sp]-thresh[sp])

				// at this point, first determine, where threshold (defined by rate) is reached again (from above) 
				IF (sp+1 < nSpikes)	// not the last spike)
					tAfter=InWave_ST[sp+1]-InWave_ST[sp]
				ELSE // last spike
					tAfter=tEND-InWave_ST[sp]
				ENDIF	
				// now search downward threshold crossing in this time range
				FindLevel /EDGE=2/Q/R=(InWave_ST[sp],InWave_ST[sp]+tAfter) InWave_V,thresh[sp]
				IF (V_Flag) // no crossing found
					IF (sp+1 < nSpikes)	// not the last spike)
						tAfter=0.5*(InWave_ST[sp+1]-InWave_ST[sp]) // half time to next spike
					ELSE // last spike
						tAfter=pnt2x(Inwave_V,DimSize(InWave_V,0)-1)-InWave_ST[sp] // all the time to end of wave
					ENDIF
					tRecross=NaN
				ELSE //  crossing found 
					tRecross=V_LevelX
					// time range of interest is time to reach threshold again plus 5 ms
					tAfter=tRecross-InWave_ST[sp]
					IF (sp==nSpikes-1) // last spike
						tAfter+=min(0.005,tEND-InWave_ST[sp]) // if possible 5 ms, otherwise to the end of wave 
					ELSE // not last spike
						tAfter+=min(0.005,InWave_ST[sp+1]-InWave_ST[sp]-0.001) // if possible 5 ms, otherwise until 1 ms before next spike
					ENDIF
				ENDIF
				// create another upsampled spike waveform and analyse minRate, FWHM
				
				IF (sp >0)	// not the first spike)
					nPntsBefore=min(0.002, InWave_ST[sp]-InWave_ST[sp-1])/dt
				ELSE
					nPntsBefore=0.002/dt
				ENDIF
				nPntsAfter=tAfter/dt
				
				
				// calculate slope
				
				sPnt=x2pnt(InWave_V,InWave_ST[sp])
				DUPLICATE/O/R=[sPnt-nPntsBefore,sPnt-1+nPntsAfter] InWave_V, dummy_US
				nPntsInDummy=DimSize(dummy_US,0)
				Differentiate /METH=0 dummy_US  /D=dummy_dV 
				Resample/up=(upsamplefac) dummy_US	// this does a sinc reconstruction
				//			Interpolate2/T=2/N=((nPntsInDummy-1)*upsamplefac+1)/E=2/Y=dummy_US dummy
	
				WAVE	dummy_US			// upsampled version covering the exact same time
				// create rate
				// it turns out that differentiating the upsampled data gives a more strongly oscillating/fluctuating rate
				// it is better to first differentiate the original and then apply upsampling
				
				Resample/up=(upsamplefac) dummy_dV
				WAVE	dummy_dV
				SetScale d 0,0,"V/s", dummy_dV
				
				
				// now find dV/dt min 
				WAVESTATS/Q/M=1 dummy_dV
				minRt[sp]	=	V_min
				
				FindLevel /EDGE=1 /Q/R=(InWave_ST[sp]-0.0015,InWave_ST[sp]+0.0015) dummy_US,HM
				IF (V_Flag==0)
					 FWHMval=0 -	V_LevelX
					 FindLevel /EDGE=2 /Q/R=(V_LevelX+0.00001,V_LevelX+tAfter) dummy_US,HM
					 IF (V_Flag==0)
					 	FWHMval += V_LevelX
					 ELSE 
						FWHMval=NaN
					 ENDIF
				ELSE 
					FWHMval=NaN
				ENDIF
				FWHM[sp]=FWHMval
				minAmp= thresh[sp]
			ENDIF
			
			// finally look for a minimum after the peak and deeper than the threshold
			// go through the original data
			
			// take same time interval as just before, but from raw data
			// only need to do this, if threshold was crossed again, i.e. if tRecross!=Nan
			IF (tRecross==tRecross)
				DUPLICATE/O/R=[sPnt-nPntsBefore,sPnt-1+nPntsAfter] InWave_V, dummy
				Differentiate /METH=0 dummy  /D=dummy_dV 

				//average with window of 0.02 ms
				VARIABLE BoxWidth=round(2e-4/dt/2)*2+1
			
				Smooth/B BoxWidth, dummy_dV;
			
				//Find Zero crossing in dV after threshold was reached again from above
				FindLevel/Edge=1/Q/R=(tRecross,tEnd) dummy_dV,0
				IF (V_flag==0)
					AHPdel[sp]= V_LevelX-tRecross
					AHP[sp]= dummy[x2pnt(dummy,V_LevelX)]
				ENDIF
			ENDIF // if it bever recrossed, just don't write anything into AHP
			

				
			 
//			Variable endOfAHPSearch
//			IF (sp< nSpikes-1)
//				endOfAHPSearch=x2pnt(InWave_V,InWave_ST[sp+1])
//			ELSE	
//			   endOfAHPSearch=DimSize(Inwave_V,0)-1
//			ENDIF
//			FindPeak /B=(max(10,0.0001/dt))/M=(minAmp) /N /P/Q/R=[sPnt,endOfAHPSearch] InWave_V	
//			IF (V_flag==0)	// minimum found
//				AHP[sp]= V_PeakVal
//			ENDIF
			KillWaves dummy,dummy_US,dummy_dV,dummy_phsSlp
		ENDFOR // loop over all spikes
		
		IF (SlopeMaxWasApplied)
			Print "in "+SlopeName+" a hard-coded cutoff of "+num2str(slopemax)+" 1/s was applied"+num2istr(SlopeMaxWasApplied)+" times."
		ENDIF
		// now create Inter-spike-interval, first entry will be NaN as there is no ISI before the first spike
		// and the idea is to be able to plot all properties against each other and hence the first spike 
		// also needs an entry here
		
		ISI[0]=NaN
		IF (nSpikes>1)
			ISI[1,nSpikes-1]-=InWave_ST[p-1]
		ENDIF
		VARIABLE 	Reporter
		
		IF (mean(thresh)!=mean(thresh) ) // contains NaNs
			Duplicate/O thresh, dummy
			thdsfRemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(thresh)
		ENDIF
		Note		 InWave_ST, "average threshold (rate crosses "+num2str(RateThresh)+" V/s):"
		Note/NOCR InWave_ST, num2str(Reporter) + " V"
		
		IF (mean(slope)!=mean(slope) ) // contains NaNs
			Duplicate/O slope, dummy
			thdsfRemoveNaNs(dummy)
			Reporter=Median(dummy)
		ELSE
			Reporter=Median(slope)
		ENDIF
		
		Note		 InWave_ST, "median phase slope at threshold:"
		Note/NOCR InWave_ST, num2str(Reporter) + " 1/s"
		
		Note		 InWave_ST, "average peak potential:"
		Note/NOCR InWave_ST, num2str(mean(Vmx)) + " V"
		Note		 InWave_ST, "average minimal rate of rise:"
		Note/NOCR InWave_ST, num2str(mean(minRt)) + " V/s"
		Note		 InWave_ST, "average peak rate of rise (PRR):"
		Note/NOCR InWave_ST, num2str(mean(mxRt)) + " V/s"
		Note		 InWave_ST, "average voltage at PRR:"
		Note/NOCR InWave_ST, num2str(mean(VmxRt)) + " V"
		Note		 InWave_ST, "average ISI:"
		Note/NOCR InWave_ST, num2str(mean(ISI,1,nSpikes-1)) + " s"
		Note		 InWave_ST, "average firing rate:"
		Note/NOCR InWave_ST, num2str(DimSize(InWave_ST,0)/(Dimsize(Inwave_V,0)*dt)) + " Hz"
		
		IF (mean(FWHM)!=mean(FWHM) ) // contains NaNs
			Duplicate/O FWHM, dummy
			thdsfRemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(FWHM)
		ENDIF
				
		Note		InWave_ST, "average FWHM:"
		Note/NOCR	InWave_ST, num2str(Reporter)  + " s"
		
		IF (mean(AHP)!=mean(AHP) ) // contains NaNs
			Duplicate/O AHP, dummy
			thdsfRemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(AHP)
		ENDIF
		Note	 	InWave_ST, "average AHP:"
		Note/NOCR 	InWave_ST, num2str(Reporter) + " V"

		IF (mean(AHPdel)!=mean(AHPdel) ) // contains NaNs
			Duplicate/O AHPdel, dummy
			thdsfRemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(AHPdel)
		ENDIF
			Note		 InWave_ST, "average AHP:"
		Note/NOCR	InWave_ST, num2str(Reporter) + " V"
		
		
		// copy the wave note of the spike time wave to the notes of all derived waves
		note/K/NOCR slope, note(InWave_ST)
		note/K/NOCR thresh , note(InWave_ST)
		note/K/NOCR Vmx, note(InWave_ST)
		note/K/NOCR mxRt, note(InWave_ST)
		note/K/NOCR minRt, note(InWave_ST)
		note/K/NOCR VmxRt, note(InWave_ST)
		note/K/NOCR FWHM, note(InWave_ST)
		note/K/NOCR ISI, note(InWave_ST)  
		note/K/NOCR AHP, note(InWave_ST)  
		note/K/NOCR AHPdel, note(InWave_ST)  
		note/K/NOCR STX, note(InWave_ST)
		
KillWaves/Z Dummy
END 


Threadsafe FUNCTION ReturnLocalCV(WAVE InWave)
// this returns a single number, the average local variation 
// according to equation 2 in 
// http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000433
// Shinomoto et al. 2009
// 
// InWave // contains spike times!
		
		VARIABLE	LV, kk, nS=DimSize(InWave,0) // number of spikes
		LV=0
		FOR (kk=0; kk<nS-2; kk+=1)
			LV+= ( ((InWave[kk+1]-InWave[kk])-(InWave[kk+2]-InWave[kk+1]))/(InWave[kk+2]-InWave[kk]) )^2
			
		ENDFOR 
			LV*=3/(nS-2)		// ns is number of spikes, the original paper used number of ISI, hence here we divide by Nspikes-2
									// while Shinomoto divides by nIntervals-1
			Return LV
END		//ReturnLocalCV

Threadsafe FUNCTION ReturnCV(WAVE InWave)
// this returns a single number, the coefficient of variation of the inter-spike intervals

// InWave // contains spike times!
IF (DimSize(InWave,0) < 2)
	Return 0
ENDIF 
		VARIABLE	nS=DimSize(InWave,0) // number of spikes
		Duplicate/O InWave aNNe_W_ISI
		WAVE	ISI=aNNe_W_ISI
		ISI[1,nS-1]-=InWave[p-1]	// from spike times to inter-spike intervals
		DeletePoints 0,1,ISI			// kill first spike time
		WAVESTATS/Q ISI
		KillWaves/Z aNNe_W_ISI
		
		Return V_sdev/v_avg
		
END		//ReturnCV

FUNCTION CollectSpikePropertiesByFolder()
// assumes that waves "~_ST" exist (spike times) and that they 
// contain a wave note with info about spike properties
// this function collects those info for all_ST waves in a folder
// they are averaged
// in the higher folder, one wave per property is created that holds the 
// averages of all folders

	DFREF		topDf=GetDataFolderDFR()
	
	STRING		FolderListe=ReplaceString(",", StringByKey("FOLDERS", DataFolderDir(1,topDf)+",",  ":", ";") ,";") 
	FolderListe= RemoveFromList("Packages;IGNORE;", FolderListe, ";")
	FolderListe=SortList(FolderListe,";",16)
	VARIABLE	runde,numFolders=ItemsInList(FolderListe)
	STRING ST_Liste
	STRING BaseName= "ByFolder_"
// create the waves that hold the average properties
	STRING ThreshName=BaseName+"VTHR"
	STRING VmaxName=BaseName+"Vmx"
	STRING mxRateName=BaseName+"mxRt"
	STRING minRateName=BaseName+"minRt"
	STRING VofMaxRateName=BaseName+"VmxRt"
	STRING FWHMName=BaseName+"FWHM"
	STRING ISIName=BaseName+"ISI"
	STRING AHPName=BaseName+"AHP"
	STRING AHPdelName=BaseName+"AHPdel"
	STRING ThreshCrossTime=BaseName+"SonsetT"
	STRING SlopeName=BaseName+"RPD"
	STRING LVName= Basename +"LV"
		MAKE/N=(numFolders) /O $SlopeName, $ThreshName, $VmaxName, $minRateName, $mxRateName, $VofMaxRateName, $ISIName, $FWHMName, $AHPName, $AHPdelName, $ThreshCrossTime, $LVName
		WAVE	slope=$SlopeName
		WAVE	thresh=$ThreshName
		WAVE	Vmx=$VmaxName
		WAVE	minRt=$minRateName
		WAVE	mxRt=$mxRateName
		WAVE	VmxRt=$VofMaxRateName
		WAVE	FWHM=$FWHMName
		WAVE	ISI=$ISIName  
		WAVE	AHP=$AHPName  
		WAVE	AHPdel=$AHPdelName  
		WAVE	LV=$LVName  
		slope	=NaN
		thresh=NaN
		Vmx	=NaN
		mxRt	=NaN
		VmxRt	=NaN
		FWHM	=NaN
		AHP 	=NaN
		LV		=NaN
		ISI	=NaN
		SetScale d 0,0,"1/s", slope
		SetScale d 0,0,"V", thresh, VmxRt, Vmx, AHP
		SetScale d 0,0,"V/s", mxRt, minRt
		SetScale d 0,0,"s", ISI
		SetScale d 0,0,"s", FWHM, AHPdel
		SetScale d 0,0,"", LV
		
		// make sure all waves that were created have index scaling
		SetScale/P x,0,1, slope,thresh, VmxRt, Vmx,ISI, mxRt, FWHM, AHP, LV

		
	VARIABLE	cST, nST
	VARIABLE	LVavg, THRavg, SLOPEavg, VMXavg, MINRTavg, MXRTavg, VMXRTavg, ISIavg, FWHMavg, AHPavg, AHPdelavg
	STRING	STListe, Notiz	
	STRING	LVstring, THRstring, SLOPEstring, VMXstring, MINRTstring, MXRTstring, VMXRTstring, ISIstring, FWHMstring, AHPstring, AHPdelstring
	STRING	currFolderName
	

	FOR (runde=0; runde< numFolders; runde+=1)
		currFolderName=StringFromList(runde,FolderListe,";")
		SetDataFolder $(currFolderName)
		LVavg=0
		THRavg=0
		SLOPEavg=0
		VMXavg=0		 
		MINRTavg=0
		MXRTavg=0
		VMXRTavg=0
		ISIavg=0
		FWHMavg=0
		AHPavg=0
		AHPdelavg=0

		
		STListe = WaveList("*_ST",";","")
		nST=ItemsInList(STListe,";")
				
		// now loop over all st waves, extract all the info
		FOR (cST=0; cST<nST; cST+=1)
			WAVE	currST=$(StringFromList(cST,STListe,";"))
			Notiz = note(currST)
			
			LVstring		= StringByKey("(LV)",Notiz,":"," ")
			THRstring	= StringByKey("V/s)",Notiz,":"," ")
			SLOPEstring	= StringByKey("threshold",Notiz,":"," ")
			VMXstring	= StringByKey("potential",Notiz,":"," ")
			MINRTstring	= StringByKey("rise",Notiz,":"," ")
			MXRTstring	= StringByKey("(PRR)",Notiz,":"," ")
			VMXRTstring	= StringByKey("PRR",Notiz,":"," ")
			ISIstring	= StringByKey("ISI",Notiz,":"," ")
			FWHMstring	= StringByKey("FWHM",Notiz,":"," ")
			AHPstring	= StringByKey("AHP",Notiz,":"," ")
			AHPdelstring= StringByKey("delay",Notiz,":"," ")
			
			// cumulate numbers
			LVavg+=str2num(LVstring)
			THRavg+=str2num(THRstring)
			SLOPEavg+=str2num(SLOPEstring)
			VMXavg+=str2num(VMXstring)
			MINRTavg+=str2num(MINRTstring)
			MXRTavg+=str2num(MXRTstring)
			VMXRTavg+=str2num(VMXRTstring)
			ISIavg+=str2num(ISIstring)
			FWHMavg+=str2num(FWHMstring)
			AHPavg+=str2num(AHPstring)
			AHPdelavg+=str2num(AHPdelstring)
	
		ENDFOR // loop across waves in a given folder
		slope	[runde]=SLOPEavg/nST
		thresh[runde]=THRavg/nST
		Vmx	[runde]=VMXavg/nST
		mxRt	[runde]=MXRTavg/nST
		minRt	[runde]=MINRTavg/nST
		VmxRt	[runde]=VMXRTavg/nST
		FWHM	[runde]=FWHMavg/nST
		AHP 	[runde]=AHPavg/nST
		AHPdel	[runde]=AHPdelavg/nST
		ISI 	[runde]=ISIavg/nST
		LV		[runde]=LVavg/nST
		
		SetDimLabel 0,runde,$(Cln2Prg(currFolderName)),	slope,thresh,Vmx,mxRt,minRt,VmxRt,FWHM,AHP, AHPdel,LV,ISI	


		SetDataFolder topDf	
	ENDFOR // loop across folders
END

Function/S Cln2Prg(String dimLabel)
	dimLabel = ReplaceString(":", dimLabel, "§")
	return dimLabel
End

Function/S Prg2Cln(String dimLabel)
	dimLabel = ReplaceString("§", dimLabel, ":")
	return dimLabel
End


FUNCTION AnaIDProtocol(STRING Prefix)
// start in folder containing voltage (*_V) and current (_I) traces
// only work with waves starting with prefix

VARIABLE useSameTime4All=1 // switches between stimulus detection for all at once or for individual stimuli
	STRING	V_list=WaveList(Prefix+"*_V",";","")
	STRING	I_list=WaveList(Prefix+"*_I",";","")
	// healthcheck
	Variable Num_V=ItemsInList(V_list,";")
	Variable Num_I=ItemsInList(I_list,";")
	IF (Num_I != Num_V)
		DoAlert 0,"Number of voltage and current traces differs ("+num2istr(Num_V)+" vs "+num2istr(Num_I)+"."
		Num_V=min(num_V,num_I)	// continue with the minimal number
	ENDIF
	IF (Num_V<1)
		DoAlert 0,"Did not find any waves matching the prefix "+Prefix+" in "+GetDataFolder(1)+"."
		Num_V=min(num_V,num_I)	// continue with the minimal number
	ENDIF
	
	// attempt to obtain a few key parameters automatically
	// assume that two large transients in the stimulus wave betray the stimulus onset and offset
	// 
	VARIABLE BaselineCurrent
	VARIABLE rr, thresh, min_dT
	MAKE/FREE/N=(max(1,(1-useSameTime4All)*num_V)) W_OnsetT,W_OffsetT

	WAVE stim = $(StringFromList(0,I_list,";")) // getting the first stimulus to understand how many points 
												// these waves have
	VARIABLE	npntsInStim=DimSize(stim,0)

IF (useSameTime4All)

	// make stim time estimation more robust by averaging across all stimuli diffs
	// then use the same times for all
	
	FOR (rr=0; rr<num_V; rr++)
		WAVE stim = $(StringFromList(rr,I_list,";"))
		IF (rr==0)								// first round, establish waves
			Duplicate/O /Free stim averageStim_D
			averageStim_D=0
					
			Duplicate/O/Free stim curr_I		// use differentiation to check where the biggest jump occurs
			
			Differentiate curr_I
			
			averageStim_D+= abs(curr_I)

		ELSEIF (DimSize(stim,0)!=npntsInStim) 		// health check: number of points should be consistent

			DoAlert 0,"Warning: stimulus "+StringFromList(rr,I_list,";")+" has "+num2istr(DimSize(stim,0))+"points. First stimulus had "+num2istr(npntsInStim)+" points. Stopping the evaluation!"
			rr=num_V 
		ELSE // normal run
		
			Duplicate/O/Free stim curr_I		// use differentiation to check where the biggest jump occurs
			
			Differentiate curr_I
			
			averageStim_D+= abs(curr_I)
		ENDIF
	ENDFOR
	averageStim_D /= num_V
	
	WAVESTATS/Q/M=1 averageStim_D
	thresh = V_max/2
	min_dt=DimDelta(averageStim_D,0)*max(5,DimSize(averageStim_D,0)/100)
	Findlevels/Q/D=OnOffTimes/M=(min_dT) averageStim_D, thresh
	// this uses the assumption that between onset and offset are at least 	5 points 
	// or even 1 % of the total duration
	// now the newly created wave OnOffTimes has ideally two entries 
	// which point to onset and offset time of the stimulus
	IF(DimSize(onOffTimes,0)!=2)
		DoAlert 0,"Check out stimulus time estimation"
	ENDIF
	W_OnsetT[0]=onOffTimes[0]
	W_OffsetT[0]=onOffTimes[1]


ELSE // get a separate time for onset and offset for each trial
	
	FOR (rr=0; rr<num_V; rr++)
		WAVE stim = $(StringFromList(rr,I_list,";"))
		Duplicate/O/Free stim curr_I
		Differentiate curr_I
		curr_I = abs(curr_I)
		WAVESTATS/Q/M=1 curr_I
		thresh = V_max/2
		min_dt=DimDelta(curr_I,0)*max(5,DimSize(curr_I,0)/100)
		Findlevels/Q/D=OnOffTimes/M=(min_dT) curr_I, thresh
		// this uses the assumption that between onset and offset are at least 	5 points 
		// or even 1 % of the total duration
		// now the newly created wave OnOffTimes has ideally two entries 
		// which point to onset and offset time of the stimulus
		W_OnsetT[rr]=onOffTimes[0]
		W_OffsetT[rr]=onOffTimes[1]
	ENDFOR
ENDIF
	
	VARIABLE onsetT=median(W_OnsetT)
	VARIABLE StimDuration=median(W_OffsetT)-median(W_OnsetT)
	
	printf "estimated onset %g, duration %g \r",onsetT,StimDuration
	STRING FSL_Name, IR_Name, R_Name, Stim_Name
	
	FSL_Name	=Prefix+"FirstSpikeLatency"
	IR_Name	=Prefix+"InitialRate"
	R_Name		=Prefix+"Rate" 
	Stim_Name	=Prefix+"StimAmp"
	Make/O/N=(num_V) $FSL_Name, $IR_Name, $R_Name, $Stim_Name 
	
	WAVE 	FirstSpikeLate=$FSL_Name
	WAVE	InitRate=$IR_Name
	WAVE	Rate=$R_Name
	WAVE	StimAmp=$Stim_Name
	
	FirstSpikeLate=NaN
	InitRate = 0
	VARIABLE	BaseLineStim	
	FOR (rr=0; rr<num_V; rr++)
		WAVE stim = $(StringFromList(rr,I_list,";"))
		onsetT=	W_OnsetT[min(rr,DimSize(W_OnsetT,0)-1)]
		StimDuration = W_OffsetT[min(rr,DimSize(W_OnsetT,0)-1)] - onsetT
		BaseLineStim+= mean(stim,DimOffset(stim,0),OnsetT-DimDelta(stim,0))  // collect baseline across all stimuli, then average
		StimAmp[rr] = mean(stim,OnsetT+DimDelta(stim,0),OnsetT+StimDuration-DimDelta(stim,0))
	ENDFOR
	BaseLineStim/=num_V
	StimAmp-=BaseLineStim	// stimamp is relative to baseline
	STRING W_Name
	VARIABLE	nSpikes, firstIndex, lastIndex, offsetT
	FOR (rr=0; rr<num_V; rr++)
		W_Name = StringFromList(rr,V_list,";")
		WAVE resp = $(W_Name)
		nSpikes = ReturnSpikeTimes(resp,V_Threshold=0, MinISI=0.002, Bool_RatePositive=1)
		IF (nSpikes > 0)
			WAVE ST=$(W_Name[0,strlen(W_Name)-2]+"ST")
			Duplicate/FREE/O ST,dummy	// a wave to manipulate to remove spikes outside stim
			InsertPoints 0,1,dummy
			dummy[0]=DimOffset(resp,0)		// introduced a fake spike right at time onset
			dummy[DimSize(dummy,0)]={rightx(resp)}		// introduced a fake spike at the end of time
			// now FindLevel will work routinely

			
			// only consider spikes during stimulus
			onsetT=	W_OnsetT[min(rr,DimSize(W_OnsetT,0)-1)]
			FindLevel/EDGE=1/P/Q dummy, onsetT
			IF (V_flag) // no crossing found
				Abort "This should not be: could not find spike time"
			ENDIF
			firstIndex=ceil(V_LevelX)


			offsetT = W_OffsetT[min(rr,DimSize(W_OnsetT,0)-1)]
			FindLevel/EDGE=1/P/Q dummy, offsetT
			IF (V_flag) // no crossing found
				Abort "This should not be: could not find spike time"
			ENDIF
			LastIndex=floor(V_LevelX)
			IF (LastIndex >0) // any spikes before stimulus offset
				FirstSpikeLate[rr] = dummy[firstIndex]-onsetT
				IF (LastIndex > FirstIndex) // more than 1 spike in stimulus 
					InitRate[rr] = 1/( dummy[FirstIndex+1]-dummy[FirstIndex])
				ENDIF
				Rate[rr]=(LastIndex-FirstIndex+1) / (offsetT-onSetT)
			ELSE
				Rate[rr]=0			
			ENDIF
		ENDIF
	ENDFOR
	IF ((stringmatch("s",WaveUnits(stim,0))) || (stringmatch("S",WaveUnits(stim,-1)))|| (stringmatch("SEC",WaveUnits(stim,-1)))|| (stringmatch("sec",WaveUnits(stim,-1))) )
		SetScale/P d,0,0,"Hz", InitRate, Rate
		SetScale/P d,0,0,"s",FirstSpikeLate
	ENDIF
	SetScale/P d,0,0,WaveUnits(stim,-1), StimAmp
	KillWaves W_OnsetT,W_OffsetT
	Note/K Rate,"nSpikes during stim divided by stim duration\r"
	Note/K InitRate,"inverse of first ISI\r"

	Note Rate, I_list
	Note InitRate, I_list
	Note FirstSpikeLate, I_list
	Note StimAmp, I_list

END // AnaIDProtocol