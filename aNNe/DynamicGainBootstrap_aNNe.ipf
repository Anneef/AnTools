#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include "DynamicGainAnalysis_aNNe"

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # # 	                         # # # # # # # # # # # # 
// # # # # # # #  STA-based Gain Bootstrap   # # # # # # # # # # # # 
// # # # # # # # 	                         # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


FUNCTION Bootstrap_STA_DerivedGain(nRounds,DoNoiseFloor[, earliestST, latestST, SpikeTimeSuffix])
VARIABLE		nRounds	//number of bootstrap rounds
// This calls all the necessary procedures for confidence interval bootstrap
// the waves AC_avg_scaled and FreqPoints have to exist before this can be run!
VARIABLE		DoNoiseFloor		// if DoNoiseFloor> 0, the spike times get shifted
										//and the outcome represents the noise floor rather than the confidence interval

VARIABLE		earliestST,latestST	// those are boundaries of possible spike times; needed in order to properly shift spike times, i.e. only needed for NoiseFloor 
	// if not provided they will be approximated by the time of the first and last spike

STRING		SpikeTimeSuffix // standard is "ST", but for special cases, it is usefull to define otherwise, for instance STwC ("with Criteria")


IF ( (ParamIsDefault(earliestST)) || (ParamIsDefault(LatestST))) 
	earliestST=1000
	latestST=-1000
	// values with earliest > latest will indicate to later stages that the values should not be used
ENDIF

IF (ParamIsDefault(SpikeTimeSuffix))
	SpikeTimeSuffix = "ST"
ENDIF


//		IF (exists("AC_avg_scaled" )!=1)
//			DoAlert 0,"AC_avg_scaled not found"
//			Return -1
//		ELSE
//			WAVE AC=AC_avg_scaled
//		ENDIF 
		AvgACFromList(ListDataByNote("AC", "spikerate",-100, 1000),"Local_AC_avg_scaled");
		WAVE Local_AC_avg_scaled
		// a local, updated version of the average autocorrelation is computed
			

		IF (exists("FreqPoints" )!=1)
			DoAlert 0,"FreqPoints not found"
			Return -1
		ELSE
			WAVE Freqs = FreqPoints
		ENDIF 
		
		PrepareST_List(SpikeTimeSuffix)
			// creates the IndexWave "BST_SpikeIndicies"
			// and the CountWave "BST_SpikeCount"
			// also creates 2 WAVES of WaveREFERENCES, which hold all STWaves and all IWaves
			// which are used, in the sequence in which they are listed in the 

		WAVE IndexWave = BST_SpikeIndicies
		WAVE	CountWave= BST_SpikeCount
		WAVE/WAVE STReferences=BST_STReferences
		WAVE/WAVE IReferences=BST_IReferences
		
		Shuffle_BST(IndexWave,CountWave, nRounds)
		
		WAVE IdxMatrix	= BST_IdxMx
		// next create the matrix that holds all STAs from Bootstrap
		// #rows taken from AC_avg_scaled
		
		Duplicate/O Local_AC_avg_scaled BST_MTX
		Redimension/D/N=(-1, nRounds) BST_MTX
		WAVE BSTResult = BST_MTX		// holds bootstrap STAs
		BSTResult=0
		
		// now this will take time
		FillBST_MultiThread(BSTResult, CountWave, IdxMatrix, STReferences, IReferences, DoNoiseFloor, earliestST, latestST)

		
		// collect complex gains from all those STAs, but first some preparations
		
		VARIABLE		nFreqs=DimSize(Freqs,0)
		MAKE/O/D/N=(nFreqs,nRounds) BST_GainMtx
		WAVE	Magnitudes= BST_GainMtx		// will hold the Magnitudes
		MAKE/O/D/N=(nFreqs,nRounds) BST_PhsMtx
		WAVE	Phases= BST_PhsMtx		// will hold the phases
		
Print "Bootstrap completed. Starting the gain calculation"
		SplitBeforeFFT(Local_AC_avg_scaled,0)
		WAVE	Local_AC_avg_scaled_splt		// created by the last command
		FFT/OUT=1 /DEST=Local_AC_avg_scaled_splt_FFT Local_AC_avg_scaled_splt // needed to calculate gains
		WAVE/C Local_AC_avg_scaled_splt_FFT
		// get overall spike rate from note of the Bootstrap STA
		STRING Notiz=Note(BSTResult)
		STRING		cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
					cellFolderList = ReplaceString(",", cellFolderList, ";")
					

		STRING Note_String=cellFolderList+"\rAVG_SD:"+StringByKey("AVG_SD", Notiz,":" ,"\r")
	
		// now read out the number of spikes that went into the STA
		
		// now read ou the duration of the trial for that STA
			VARIABLE	totduration=str2num(StringByKey("totalduration", Notiz,":" ,"\r"))
			Note_String+="\rtotalduration:"+StringByKey("totalduration", Notiz,":" ,"\r")
			VARIABLE	totnSpikes=str2num(StringByKey("total#spikes", Notiz,":" ,"\r"))
			Note_String+="\rtotal#spikes:"+StringByKey("total#spikes", Notiz,":" ,"\r")

		VARIABLE globalrate=totnSpikes/totduration
		//		globalrate=4.56144&& possiblz give a fixed global ate to fix issue
		IF (globalrate!=globalrate)
			Abort "Everything went ok but I cannot determine global rate"
		ENDIF
//
//		
		VARIABLE	rr
		// create Gain from Local_AC_avg and all the different bootstrap STAs
		// for Multithreading, create a wave to hold all wave references
		MAKE/O/WAVE/N=(nRounds) W_holder_cmplx
		MUltiThread/NT = (ThreadProcessorCount -2) W_holder_cmplx[]=CreateFilteredGain(BSTResult,p, Local_AC_avg_scaled_splt_FFT,DoNoiseFloor,globalRate,FreqPoints=Freqs)
		
		
		FOR (rr=0; rr< nRounds; rr+=1)
			WAVE/C	W_G_cmplxFlt=W_holder_cmplx[rr]
			Magnitudes[][rr]	=real(r2polar(W_G_cmplxFlt[p]))
			Phases[][rr]		=imag(r2polar(W_G_cmplxFlt[p]))
		ENDFOR
		KillWaves/Z W_D,W_G, W_D_splt, W_D_splt_FFT,  W_G_cmplxFlt
		STRING	NameOfResult
		// sort bootstrap magnitudes to identify the required percentiles
		// this is done for each frequency point separately and independent from any other point
		
		// but first
		IF (DoNoiseFloor == 0) // if it is not noise floor, then also create confidence intervals for phases
			MAKE/O/D/N=(nFreqs,3) Phs_CI 
			WAVE Phs_CI
			// also need a temporary wave that holds phases and corresponding magnitudes in two columns as input for statscircularMeans
			MAKE/O/N=(nRounds,2) W_tempphsmag
			FOR (rr=0; rr< nFreqs; rr+=1)		// go through frequency points one by one
				W_tempphsmag[][0]=Phases[rr][p]
				W_tempphsmag[][1]=Magnitudes[rr][p]
				StatsCircularMeans/Q/ALPH=0.05/CI W_tempphsmag
				WAVE	W_CircularMeans
				// now W_CircularMeans holds the average and confidence interval boundaries in the repectively labelled rows:
				Phs_CI[rr][0]=W_CircularMeans[%CI_t1]
				Phs_CI[rr][1]=W_CircularMeans[%tBar]
				Phs_CI[rr][2]=W_CircularMeans[%CI_t2]
				// this gives strange values for mean when the noise is extremely large, the CIs seem to be true CI of phases, not CI around average vector(phase, mag)
			ENDFOR
				// for unwrapping try 
				// -2*Pi
				// Unwrap 2*Pi, wave
		ENDIF		

		
		Variable upperLimit, lowerLimit
		IF (DoNoiseFloor >0)
			NameOfResult="Noise_Floor"
			upperLimit=0.90
			lowerLimit=0.95
		ELSE
			NameOfResult="Mag_CI"
			upperlimit=0.975
			lowerLimit=0.025
		ENDIF
		Duplicate/D/O Magnitudes, $NameOfResult 
 		WAVE Conf_Int=$NameOfResult
		FOR (rr=0; rr< nFreqs; rr+=1)		// go through frequency points one by one
			Duplicate/O/R=[rr][0,nRounds-1] Conf_Int Dum 	// make temporary copy of all the magnitude entries 
			Redimension/N=(nRounds,0) Dum                	// create column vector out of the 1 x nRounds matrix
			Sort Dum, Dum										// those entries are now sorted
			Conf_Int[rr][0,nRounds-1]= Dum[q]				// write over the unsorted values
		ENDFOR
		// now each column contains sorted entries of magnitudes
		
 		WAVE Conf_Int=$NameOfResult
		FOR (rr=0; rr< nFreqs; rr+=1)
			Duplicate/O/R=[rr][0,nRounds-1] Conf_Int Dum
			Redimension/N=(nRounds,0) Dum
			Sort Dum, Dum
			Conf_Int[rr][0,nRounds-1]= Dum[q]
		ENDFOR

	

		KillWaves W_holder_cmplx, Magnitudes, phases
		KillWaves/Z Dum, STReferences, IReferences//
		
		
		// get the avg of the bootstrap distribution
		MatrixOP/O Boot_avg=averagecols(Conf_Int^t)^t
		
		// sort out the 5 and 95% and the median
		Conf_Int[][0]=Conf_Int[p][nRounds*lowerLimit]
		Conf_Int[][2]=Conf_Int[p][(nRounds-1)*upperLimit]
		Conf_Int[][1]=(Conf_Int[p][floor(nRounds/2)]+Conf_Int[p][ceil(nRounds/2-1)])/2
		Redimension/N=(-1,3) Conf_Int
		Note/K Conf_Int, Note_String
END


FUNCTION PrepareST_List(ST_Suff)
// this creates two waves to hold 
// BST_SpikeIndicies 	- of all spike time waves here are the indicies of usable spikes listed
// BST_SpikeCount  		- of all spike time waves this holds the wave names (in DimLabel) and it holds
//							  the number of usable spike times (the once too close to begin / end cannot be used

// this has to be called ONCE before ANY Bootstrap
STRING			ST_Suff					// suffix identifying spike times, without leading "_"
	
	
	VARIABLE	range=1		// duration of STA (in seconds)
	
	STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
					cellFolderList = ReplaceString(",", cellFolderList, ";")
						cellFolderList= RemoveFromList("Packages;IGNORE;", cellFolderList, ";")

	VARIABLE		ii, nSubs=ItemsInList(cellFolderList,";")
	STRING			currentDFPath = GetDataFolder(1)		
					//nsubs=1 // hard wired when going through individual subfolders
	
	IF (nSubs==0)	// there are no subfolders, work inside the current folder
		cellFolderList = currentDFPath+";"
	ENDIF
	VARIABLE		kk, nSTs
	DFREF			rootFolderRf=GetDataFolderDFR()
	STRING			STList, currIName, currSTName
	VARIABLE		minST, maxST, count
	
	
	// make the list of indices of spike times 
	// also make a wave that holds the number of usable spikes per spike time wave
	
	MAKE/O/N=0 BST_SpikeIndicies
	WAVE BST_SpikeIndicies = BST_SpikeIndicies
	MAKE/O/N=0 BST_SpikeCount
	WAVE	BST_SpikeCount= BST_SpikeCount
	// now make waves to hold references of the STWaves and IWaves 
	// these are only needed for the Multithreading
	MAKE/O/WAVE/N=0 BST_STReferences
	MAKE/O/WAVE/N=0 BST_IReferences
	WAVE/WAVE BST_STReferences
	WAVE/WAVE BST_IReferences
	// go through all subfolders, get spike time waves and put the usable spikes in a list
	
	FOR (ii=0;ii<nSubs;ii+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(ii, cellFolderList,";")
		SetDataFolder $(FolderName)
		StList=WaveList("*_"+ST_Suff,";","")
		nSTs=ItemsInList(STList,";")
		FOR (kk=0; kk< nSTs; kk+=1)
			currSTName = StringFromList(kk, STList,";")
			WAVE currST=$currSTName
			currIName=RemoveEnding(currSTName,ST_Suff) + "I" 
			WAVE	currI=$currIName
			
			// only needed to find out the duration and with this the range of spike times 
			// that can be used for the STA
			minST=DimOffset(currI, 0)+range/2
			maxST=DimOffset(currI, 0)-range/2+DimDelta(currI, 0)*Dimsize(currI, 0)
			
			FindLevel /Q /P currST, minST
			IF (!V_flag)	// level found
				minST=ceil(V_LevelX)
			ELSEIF (currST[0]> minST)
				minST=0
			ELSE
				DoAlert 0,"no spike time > "+num2str(minST)+" found in "+GetDataFolder(0)+currSTName
				minST=0
			ENDIF
			
			FindLevel /Q/P/R=[DimSize(currST,0)-1,0] currST , maxST
			IF (!V_flag)	// level found
				maxST=floor(V_LevelX)
			ELSEIF (currST[DimSize(currST,0)-1]< maxST)
				maxST=DimSize(currST,0)-1
			ELSE
				DoAlert 0,"no spike time < "+num2str(maxST)+" found in "+GetDataFolder(0)+currSTName
				maxST=DimSize(currST,0)-1
			ENDIF
			BST_SpikeCount[DimSize(BST_SpikeCount,0)]={maxST-minST+1}
			SetDimLabel 0,DimSize(BST_SpikeCount,0)-1,$(Cln2Prg(GetDataFolder(0)+":"+currSTName)) BST_SpikeCount
			count = DimSize(BST_SpikeIndicies,0)
			Redimension/N=(count+maxST-minST+1) BST_SpikeIndicies
			BST_SpikeIndicies[count,count+maxST-minST]=minST+p-count
			BST_STReferences[DimSize(BST_STReferences,0)]={currST}
			WAVE	currI=$currIName
			BST_IReferences[DimSize(BST_IReferences,0)]={currI}

		ENDFOR
	ENDFOR
	SetDataFolder rootFolderRf
	Printf "Found %g spikes in %g sweeps in %g folders\r",DimSize(BST_SpikeIndicies,0), DimSize(BST_SpikeCount,0), nsubs


END
			
	
FUNCTION Shuffle_BST(IndexWave,CountWave, nBSTrounds)
WAVE	IndexWave,CountWave
VARIABLE	nBSTrounds	// number of bootstrap rounds

// CountWave contains  the names (incl.path) of all ST waves
// in the dim label and the respective number of spikes that 
// can be used from each of those spike time waves

// IndexWave is simply a concatenation of the indicies of the spikes that are usable

// now we need to create a matrix with indicies such that there is a random subset 
// of the indicies of every wave in the different columns of that matrix
// the subsets have to be balanced, i.e. every spike index is occuring BSTNumber times
// in total

	MAKE/O/U/W/N=(DimSize(IndexWave,0),nBSTrounds) BST_IdxMx  // 16 bit unsigned integer, 
																		//only indicies up to 2^16-1 are allowed
	
	WAVE	BST_IdxMx
	BST_IdxMx=NaN
	
	VARIABLE		nWaves=DimSize(CountWave,0), ww, length, previous
	previous=0		// number of entries in IndexWave previous to current ST wave
	
	FOR (ww=0; ww< nWaves; ww+=1)		// goes over all the waves containing spike times (ST waves)
		length=CountWave[ww]				// the number of spikes in the ST wave currently accessed
		MAKE/O/N=(length*nBSTrounds) W_dummy, W_rand	// create two waves that hold (nBSTrounds x current spike time count) numbers
																	// i.e. each usefull spike in the current ST Wave will appear nBSTrounds times
																	
																	
		W_dummy[]=IndexWave[previous+mod(p,length)]	// fill with the indicies,
																	// that apply to the current ST wave
																	// over and over (BSTnumber times).
																	// 	the offset "previous" assures
																	// that we start with the correct spike index
																	// the mod function assures that we keep going over the same indicies
																	// again and again (nBSTrounds times)
		
		W_rand = enoise(1,2)									// random number sequence
		SORT W_rand , W_dummy								// now the indicies are randomly shuffled
		// next, we fill the wave of BST indicies
		BST_IdxMx[previous,previous+length-1][0,nBSTrounds-1] = W_dummy[p-previous+q*length]
		
		// the shuffling here is using the same number of spikes for a given wave in each BST round. 
		// only which spikes are used differs between rounds, but the number of spikes stays the same or a given ST wave. 
		// It is the number that is laid out in the count wave
		
		previous+=length	
	ENDFOR

	KillWaves/Z W_dummy, W_rand
END		// Shuffle_BST



Threadsafe	FUNCTION/DF	Bootstrap_MultiThread(CountWave,IdxMatrix, currentBSTround, STWaveReferences, IWaveReferences, DoNoiseFloor,earliestST, latestST)
WAVE			CountWave, IdxMatrix
WAVE/WAVE	STWaveReferences, IWaveReferences
VARIABLE		currentBSTround			// the current bootstrap round
VARIABLE		DoNoiseFloor,earliestST, latestST // if DoNoiseFloor the spike times are shifted (by more than one orrelation time)
															 // earliest and latest are used to 

// this goes through all the spike time waves and re-creates copies using the indicies in the 
// idxMatrix. These waves are only needed to create boostrap STAs so they
// are just temporary spike time waves and get overwritten and finally deleted

// also the corresponding STAs are created calling STA from analogue. 
// this will be run in a mode such that the out-name does not end in "_STA" but in "_STB"
	DFREF dfSav= GetDataFolderDFR()

	// Create a free data folder to hold the results
	DFREF dfFree= NewFreeDataFolder()
	SetDataFolder dfFree

	VARIABLE		nWaves=DimSize(CountWave,0), ww, length, previous, range
	STRING		TempSTName
	previous=0		// number of entries in IndexWave previous to current ST wave
	MAKE/O/WAVE/N=(nWaves)	STB_References 	// will hold the references of all the individual STAs created 
											// will eventually be handed over to the AvgSTA_Threadsafe
	FOR (ww=0; ww< nWaves; ww+=1)
		//print 1,GetdataFolder(1)
		length=CountWave[ww]
		WAVE	currST = STWaveReferences[ww]
		WAVE	currI	 =	IWaveReferences[ww]
		
		TempSTName="W_temp_ST_"+num2istr(currentBSTround)
		
		MAKE/O/N=(length) $TempSTName
		WAVE	ST=$TempSTName
		ST[]	= IdxMatrix[previous+p][currentBSTround]	// filling with indicies
		ST[]	= currST[ST[p]]									// filling with the acctual spike times
		// At this point, it is decided whether the actual spike times according to the index are used,
		// or whether a random cyclic shift is added to those times in order to create the noise floor estimate
		IF (DoNoiseFloor)
			// shift the spike times by an random ammount that is sufficiently large to shift the spikes to a 
			// region of the input that is uncorrelated with the region that actually triggered the spike
			// to find out, what range of shift values is available, use earliest and latest
			IF (earliestST > latestST)				// meaningless values created to indicate that no values had been supplied by user
				range=wavemax(ST)-wavemin(ST) 
				// want at least 1 second shift and no more than range-1 s shift
				WAVE ST2=CyclicShift(ST, range/2+enoise(range/2-1,2))
			ELSE
				range=LatestST-earliestST 
				WAVE ST2=CyclicShift(ST, range/2+enoise(range/2-1,2), startTime=earliestST, endTime=LatestST)
			ENDIF
			ST[]=ST2[p]
		ENDIF
		Sort	ST,ST													// Sort to make sensible ST wave

		// now take the spike time wave with the boostrap sample of spike times
		// and the corresponding original current
		// and create a STA 
		STB_References[ww]=STAfromAnalogueWaveRef(currI, ST, 1,0,Suffix="B"+num2istr(ww)+num2istr(currentBSTround)+"_STB")

		previous+=length	
		KillWaves/Z ST, ST2

	ENDFOR
			
	// make average STA from the bootstrap STAs
	
	WAVE res=AvgSTA_ThreadSafe(STB_References,targetName="avgSTB") // no specific name needed, because we work in a 
																						// separate free data folder and 
																						// can use a generic name
																						//_round"+num2istr(currentBSTround))
	FOR (ww=0; ww< nWaves; ww+=1)
		WAVE toDel=STB_References[ww]
		KillWaves toDel
	ENDFOR
	WAVE staavg
	KillWaves STB_References, staavg
	return dfFree
END	//Bootstrap_MultiThread



Function FillBST_MultiThread(BSTResult, CountWave, IdxMatrix, STReferences, IReferences,  DoNoiseFloor,  earliestST, latestST)
	WAVE 		BSTResult, CountWave, IdxMatrix
	WAVE/WAVE	STReferences, IReferences
	VARIABLE	DoNoiseFloor,  earliestST, latestST

	
	Variable ncol= DimSize(BSTResult,1)
	Variable col
	Variable time2, ttime= stopMSTimer(-2)
	// Create a wave to hold data folder references returned by Worker.
	// /WAVE specifies the data type of the wave as "wave reference".
	Make/O/DF/N=(ncol) ddd
	VARIABLE numCores=ThreadProcessorCount
	// one processor will not be engaged to keep the computer responsive

	MultiThread/NT=(numCores-5)	ddd= Bootstrap_MultiThread(CountWave,IdxMatrix,p, STReferences, IReferences, DoNoiseFloor,  earliestST, latestST)
	
	

	
	for(col=0; col<ncol; col+=1)
		DFREF df= ddd[col]
		Duplicate/O df:avgSTB, DDummy
		BSTResult[][col] = DDummy[p]
		// now take over the info on avg sd , spikes and duration
		IF (col==0)
			STRING 	Notiz=Note(df:avgSTB)
			Note/K BSTResult
			Note/NOCR BSTResult, Notiz
		ENDIF
		//KillWaves/Z res
	endfor
	time2= stopMSTimer(-2)
	printf "%d rounds in %d sec;\r",ncol, (time2-ttime)*1e-6
	KillWaves ddd,DDummy

End

Threadsafe FUNCTION/WAVE CreateFilteredGain(WAVE M_STA,VARIABLE clmn, WAVE/C FFTofAC,VARIABLE DoNoiseFloor,VARIABLE globalRate,[WAVE/D FreqPoints, VARIABLE widthFac])
//WAVE		M_STA 		// matrix full of bootstrapped spike triggered averages 
//VARIABLE	clmn		// the clmn of interest in the present instance
//WAVE/C		FFTofAC	// fourier transform of autocorrelation, which was split at middle
//WAVE/D		FreqPoints	// vector with Frequency values at which filtered versions of 
//						// FT will be constructed
//VARIABLE	DoNoiseFloor,globalRate
//VARIABLE	widthFac		// width of the gaussian filter	; a value of 1 corresponds to 
						// Higgs and Spain's choice of f/2Pi
		// the optional parameters FreqPoints and widthFac  are supplöied as follows:
		// GaussFilter(FT,FreqPoints=FrequencyRange, widthFac=1)
		

		
			DFREF dfSav= GetDataFolderDFR()
	
			// Create a free data folder and set it as the current data folder
			SetDataFolder NewFreeDataFolder()

			Duplicate/O/R=[*][clmn] M_STA STA
			Redimension/N=(-1,0) STA
			VARIABLE SplitTime
			
			IF (DoNoiseFloor >0)
				SplitTime=0
			ELSE
				Wavestats/M=1/Q STA
				SplitTime=V_MaxLoc
			ENDIF
			IF ( SplitTime < DimOffset(STA,0)  || SplitTime > pnt2x(STA,numpnts(STA)-1) )
				Print "split Time is outside x-range" 
			ENDIF
			
			VARIABLE	dT=DimDelta(STA,0)
		
			VARIABLE SplitP=x2pnt(STA,SplitTime+dT/1000)
		
			VARIABLE truncLenght=2*floor(DimSize(STA,0)/2)
			MAKE/O/D/N=(truncLenght) out
			WAVE out

			out[0,DimSize(out,0)-1-SplitP] 						= STA[p+SplitP]	
			out[DimSize(out,0)-SplitP,DimSize(out,0)-1]		= STA[p+SplitP-DimSize(out,0)]	

			FFT/OUT=1 /DEST=W_G out
			KillWaves out
			W_G/=FFTofAC
			W_G*=cmplx(globalRate,0)			
			W_G=conj(W_G)			// added 03/08/2020 to obtain correct sign of the phase
	
			VARIABLE	df_in=DimDelta(FFTofAC,0)		// frequency step width in input
			VARIABLE	minf=DimOffset(FFTofAC,0)
			VARIABLE	maxf=minf+ (DimSize(FFTofAC,0) -1)*df_in
			VARIABLE	nPoints	
			IF (ParamIsDefault(FreqPoints))
			          	nPoints = min(floor((DimSize(FFTofAC,0) -1)/2), 50)		// max 50 Freq points, but no more than 
																		// half the number that is there already
				// linear distance in log requires a frequency factor 
				// to describe construction of the FreqList
				VARIABLE	freqFac =10^( log((maxf-minf)/df_in)/(nPoints-1))	// equivalent to the (n-1)th root of the quotient (maxf-minf)/df_in
				Make/D/O/N=(nPoints) FreqPoints
				Wave FreqPoints
				FreqPoints[]=minf+df_in*freqFac^p
				
			ELSE
				 nPoints= DimSize(FreqPoints,0)
				
			ENDIF
			IF (ParamIsDefault(widthFac))
				widthFac = 1
			ENDIF
	     	VARIABLE k,  center, width, weightSum
	     	DUPLICATE/C/D/O FFTofAC , GWeight
			Make/C/D/O/N=(nPoints) CmplxResult				// the complex valued result of weighting the noisy complex gain				         
	      FOR (k=0; k< nPoints; k+=1)
				width=FreqPoints[k]/2/Pi	* widthFac		// SD of Gaussian is 1/freq
				GWeight[]=cmplx( gauss(pnt2x(FFTofAC, p ),FreqPoints[k], width),0)
				WeightSum = mean(GWeight)*Dimsize(GWeight,0)
				GWeight/=WeightSum			// making sure the integral of the weight is one
				GWeight*=W_G	
				CmplxResult[k]=sum(GWeight)
			ENDFOR
		KillWaves/Z GWeight, W_G,STA
		SetDataFolder dfSav
		Return  CmplxResult
END

Threadsafe FUNCTION/WAVE	CyclicShift(SpikeTimeWave, Delta, [startTime,endTime])
WAVE		SpikeTimeWave	
VARIABLE	Delta, startTime, endTime

// spikes times (t0.....tn) are shifted forward by Delta d
// t0'=t0+d, t1'+d, .... tn'+d  those that are now after tn will be cycled by adding t0-tn
// e.g. tn' will end up at tn+d+t0-tn
// this is not too nice, as then there are two spikes at t0+d
// there is no good way around it, I could add the minimal ISI , one fundamental
// problem is that there is no way to know where to start placing the spike times
// unless the time at which the first spike can be expecte and the time at which the last spike can be expected are provided


// prepare the wave copy which will be manipulated and finally returned (by reference)
	Duplicate/O SpikeTimeWave, Temp_ST
	
	WAVE	ST=Temp_ST
	
	// just make sure the spike times are sorted:
	Sort ST,ST
	WaveStats/Q /M=1 ST
	IF (ParamIsDefault(StartTime))
		StartTime=V_min
	ENDIF
	
	IF (ParamIsDefault(EndTime))
		EndTime=V_max
	ENDIF
	
	IF (Delta> (EndTime-StartTime))	// shift would go cyclic more than one full round
		Delta -= floor(Delta/(EndTime-StartTime))*(EndTime-StartTime)
	ENDIF
	
	ST[]+=Delta
	// first check, if the min spike time is already above EndTime:
	IF (wavemin(ST)>EndTime)
		ST-=EndTime
	ELSEIF (wavemax(ST)<EndTime)
		// do nothing, apparently the shift has not moved a single spike beyond the allowed limits
	ELSE
		// some but not all spikes are over the allowed range
		FindLevel/EDGE=1 /P/Q ST, EndTime
		// all spikes with  indicies after V_LevelX have to be cycled to the begin
		IF (V_flag==0)		// level crossing was found
			ST[ceil(V_LevelX),DimSize(ST,0)-1]-=EndTime-StartTime
		ENDIF
		
		Sort ST, ST
	ENDIF
	
	Return ST
END // cyclic shift


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # #  	                         			 # # # # # # # # # # # # 
// # # # # # # #   VectorStrength-based Gain Bootstrap   # # # # # # # # # # # # 
// # # # # # # #  	                        			 # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

FUNCTION SineAmp_BalancedWrapBootstrap(WAVE FrequencyList)
// the conventional algorithm calculates a vector strength from 
// phases of spiketimes that were elicited in response to sinus components with 
// DIFFERENT amplitudes
// in a separate step, these amplitudes are averaged and then the gain is the 
// ratio of these two numbers (multiplied with 2 and the firing rate)

// here, this calculation is done on a spike basis: 
// in complex space, the individual contributions are not equal
// (complex number of magitude 1)
// but are scaled with 1/sine-Amplitude
// also, instead of dividing the sum of these components by the number of spikes
// and then multiplying with the firing rate (n_spikes/total time)
// we now only divide by total time

// this means, we need another data collection because we need a sine amplitude
// for every phase data point
// because tests suggest that FFT-based sine amplitude estimation is more reliable
// here only this fft-derived sine amp is used

// this is not compatible with the standard 
// analysis, it is processed here in one block


// Frequency List should be a column vector with frequencies that appear in the 
// notes of phase waves (spike phase, ending in "_SPhs")

// the Wave will be extended by further columns,
// containing the number of spikes, the mean VS, phase and the quantiles for 2.5%,50% and 97.5%)

	VARIABLE 	nFreqs=DimSize(FrequencyList,0)
	// makeing space for results
	Redimension/D/N=(nFreqs,13) FrequencyList
	// labelling results
	SetDimLabel 1,0,SineFreq,FrequencyList
	SetDimLabel 1,1,SpikeCount,FrequencyList
	SetDimLabel 1,2,meanPhs,FrequencyList
	SetDimLabel 1,3,meanVS,FrequencyList
	SetDimLabel 1,4,VS_2p5Perc,FrequencyList
	SetDimLabel 1,5,VS_median,FrequencyList
	SetDimLabel 1,6,VS_97p5Perc,FrequencyList
	SetDimLabel 1,7,meanVSGain,FrequencyList
	SetDimLabel 1,8,VSGain_2p5Perc,FrequencyList
	SetDimLabel 1,9,VSGain_median,FrequencyList
	SetDimLabel 1,10,VSGain_97p5Perc,FrequencyList
	SetDimLabel 1,11,avgSineAmp,FrequencyList
	SetDimLabel 1,12,totalDuration,FrequencyList
	
	STRING ListOfWaves
	VARIABLE	nWaves
	
	// looping through all frequencies
	VARIABLE 	ff, currfreq
	VARIABLE 	ww, t1,t2
	VARIABLE	rr, rmax=5// number repetitions for sine fit
	VARIABLE	Amp=0
	VARIABLE	totalDur=0
	VARIABLE 	nextIndex
	FOR (ff=0; ff< nFreqs; ff++)
		currfreq=FrequencyList[ff][0]
		WAVE Phases=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal= currfreq, up4NoteVal=currfreq)
		// here collect individual waves that went into the Phases and 
		// from them create the Amplitude wave 
		ListOfWaves=note(Phases)
		ListOfWaves=ReplaceString("\r", ListOfWaves, ";") 	// replace the carriage returns
												//with semicolon to have a standard wave list
		ListOfWaves=RemoveListItem(0, ListOfWaves,";") 		// the first line of the note is not a wave name
		WAVE/WAVE SPhsrefs=ListToWaveRefWave(ListOfWaves,1)
		// for total duration, also get corresponding current waves
		ListOfWaves=ReplaceString("_SPhs;", ListOfWaves, "_I;") 	// replace the carriage returns
		// now turn the list into References
		WAVE/WAVE Irefs=ListToWaveRefWave(ListOfWaves,1)
// now create amplitude wave by going through the individual waves
		Duplicate/O Phases, W_AmpWave
		Redimension/D W_AmpWave
		W_AmpWave=0
		totalDur=0
		nWaves=DimSize(SPhsrefs,0)
		nextIndex=0
		FOR (ww=0; ww<nWaves; ww++)
			WAVE CurrIWave=Irefs[ww]
			WAVE CurrSPhsWave=SPhsrefs[ww]
			Amp=numberbyKey( "Sine Amp (A)",note(CurrSPhsWave),":","\r")
						
			W_AmpWave[nextIndex,nextIndex+DimSize(CurrSPhsWave,0)-1]=Amp
			nextIndex+=DimSize(CurrSPhsWave,0)
			totalDur+=rightx(CurrIWave)
		ENDFOR
		
		IF (DimSize(Phases,0)!=DimSize(W_AmpWave,0))
			Abort "Amplitudes and Phases have to have the same number of points"
		ENDIF
		
		[WAVE QuantilesReturned,VARIABLE SpikeCount,VARIABLE VS,VARIABLE Phs]=SinAmp_WrapBootstrapVSKernel(Phases,W_AmpWave)
		FrequencyList[ff][1]= SpikeCount
		FrequencyList[ff][2]= Phs
		FrequencyList[ff][3]= VS*mean(W_AmpWave)
		FrequencyList[ff][4]= QuantilesReturned[0]
		FrequencyList[ff][5]= QuantilesReturned[1]
		FrequencyList[ff][6]= QuantilesReturned[2]
		FrequencyList[ff][11]= mean(W_AmpWave)

		// now turn vectorstrength into gain
		// for this, need average sine amplitude 
		// from CollectDataByNoteAndCriterium the note of the wave Phases
		// contains one line with a condition and then all other lines 
		// with the full path to the SPhs waves
		// The note of each of these waves contains the sine amplitude estimated from FFT of the current
		// and that seems to be the best way to estimate the sinus amplitude
		// 
		// instead: take corresponding current autocorrelation and fit (away from zero lag) with 
		// a sinusoid, then thake square root of the Amplitude
		FrequencyList[ff][12] =totalDur
		FrequencyList[ff][7] =VS*2/FrequencyList[ff][12]*FrequencyList[ff][1]
		
		FrequencyList[ff][8,10]= FrequencyList[ff][q-4]*2*FrequencyList[ff][1]/FrequencyList[ff][12]
		// gain is vector strength*2/sineamplitude*firing rate

	ENDFOR
	





END

FUNCTION [WAVE QuantilesReturned,VARIABLE SpikeCount,VARIABLE VS,VARIABLE Phs] SinAmp_WrapBootstrapVSKernel(WAVE Phases,WAVE Amps)

							// according to the indicies in the clmn-th column of indexwave 
							// this boostrap sample is then used to compute a vector strength and the average phase
							// those two are returned as real and imaginary part of a complex number

	MAKE/O/N=3 TargetQuantiles={0.025,0.5,0.975}
	WAVE QuantileWave=TargetQuantiles
	
	VARIABLE	nEntries=DimSize(Phases,0)
	VARIABLE	nboots=10000
	MAKE/C/O/N=(nBoots)	BootStatsC
			
	MultiThread  BootStatsC[]=Resample_AmpNormed_VSAndPhase(Phases,Amps) // returns in polar form magnitude and angle
	
	// get average vector from original sample
	Duplicate/FREE Phases, realPart, cmplxPart
	Redimension/D realPart, cmplxPart
	
	realPart		=	cos(Phases)/Amps
	cmplxPart		= 	sin(Phases)/Amps
	
	VARIABLE	rAVG	=	mean(realPart), cAVG=mean(cmplxPart)
	VARIABLE	meanVS		= 	sqrt(rAVG^2+cAVG^2), meanPhs = atan2(cAVG,rAVG)
	
	KillWaves/Z realPart, cmplxPart

	// now project every entry in the Bootstrap sample onto this vector
	BootStatsC[]=cmplx(real(BootStatsC[p])*cos((imag(BootStatsC[p])-(meanPhs))),0) // the phase does not matter any more
	Redimension/R BootStatsC
	// now go ahead with the same statistics as for case without phases

	STRING NoteInc, NoteString=""
	sprintf NoteInc, "Non-balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
	NoteString+= NoteInc+"\r"
	
	
	sprintf NoteInc, "Bootstrap distribution has the quantiles\r\t%1.3f\t%1.3f\t%1.3f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
	NoteString+= NoteInc+"\r"
				
	QuantilesFromSample(BootStatsC,QuantileWave,1) 	// before call QuantileWave contains the values between 0 & 1 at which  quantiles should be 
													// calculated, after the call, it  

	KillWaves/Z BootStatsC

	Return [QuantileWave, nEntries, meanVS, meanPhs]



END

FUNCTION WrapBootstrapVS(WAVE FrequencyList)
// Frequency List should be a column vector with frequencies that appear in the 
// notes of phase waves (spike phase, ending in "_SPhs")

// the Wave will be extended by further columns,
// containing the number of spikes, the mean VS, phase and the quantiles for 2.5%,50% and 97.5%)

	VARIABLE 	nFreqs=DimSize(FrequencyList,0)
	// makeing space for results
	Redimension/D/N=(nFreqs,13) FrequencyList
	// labelling results
	SetDimLabel 1,0,SineFreq,FrequencyList
	SetDimLabel 1,1,SpikeCount,FrequencyList
	SetDimLabel 1,2,meanPhs,FrequencyList
	SetDimLabel 1,3,meanVS,FrequencyList
	SetDimLabel 1,4,VS_2p5Perc,FrequencyList
	SetDimLabel 1,5,VS_median,FrequencyList
	SetDimLabel 1,6,VS_97p5Perc,FrequencyList
	SetDimLabel 1,7,meanVSGain,FrequencyList
	SetDimLabel 1,8,VSGain_2p5Perc,FrequencyList
	SetDimLabel 1,9,VSGain_median,FrequencyList
	SetDimLabel 1,10,VSGain_97p5Perc,FrequencyList
	SetDimLabel 1,11,avgSineAmp,FrequencyList
	SetDimLabel 1,12,totalDuration,FrequencyList
	
	STRING ListOfWaves
	VARIABLE	nWaves
	
	// looping through all frequencies
	VARIABLE 	ff, currfreq
	VARIABLE 	ww, t1,t2
	VARIABLE	rr, rmax=5// number repetitions for sine fit
	VARIABLE	squareAmp=0
	VARIABLE	totalDur=0
	FOR (ff=0; ff< nFreqs; ff++)
		currfreq=FrequencyList[ff][0]
		WAVE Phases=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal= currfreq, up4NoteVal=currfreq)
		[WAVE QuantilesReturned,VARIABLE SpikeCount,VARIABLE VS,VARIABLE Phs]=BootstrapVSKernel(Phases)
		FrequencyList[ff][1]= SpikeCount
		FrequencyList[ff][2]= Phs
		FrequencyList[ff][3]= VS
		FrequencyList[ff][4]= QuantilesReturned[0]
		FrequencyList[ff][5]= QuantilesReturned[1]
		FrequencyList[ff][6]= QuantilesReturned[2]
		FrequencyList[ff][11]= 0

		// now turn vectorstrength into gain
		// for this, need average sine amplitude 
		// from CollectDataByNoteAndCriterium the note of the wave Phases
		// contains one line with a condition and then all other lines 
		// with the full path to the SPhs waves
		// The note of each of these waves contains the sine amplitude estimated from FFT of the current
		// and that seems to be the best way to estimate the sinus amplitude
		// 
		// instead: take corresponding current autocorrelation and fit (away from zero lag) with 
		// a sinusoid, then thake square root of the Amplitude
		ListOfWaves=note(Phases)
		ListOfWaves=ReplaceString("\r", ListOfWaves, ";") 	// replace the carriage returns
												//with semicolon to have a standard wave list
		ListOfWaves=RemoveListItem(0, ListOfWaves,";") 		// the first line of the note is not a wave name
		WAVE/WAVE SPhsrefs=ListToWaveRefWave(ListOfWaves,1)
		ListOfWaves=ReplaceString("_SPhs;", ListOfWaves, "_I_AC;") 	// replace the carriage returns
		nWaves = ItemsInList(ListOfWaves,";")
		// now turn the list into References
		WAVE/WAVE ACrefs=ListToWaveRefWave(ListOfWaves,1)
		ListOfWaves=ReplaceString("_I_AC;", ListOfWaves, "_I;") 	// replace the carriage returns
		// now turn the list into References
		WAVE/WAVE Irefs=ListToWaveRefWave(ListOfWaves,1)
		
		// fit sine component to obtain sine amplitude, do a few times 
		totalDur=0
		FOR (ww=0; ww<nWaves; ww++)
//			WAVE CurrWave=ACrefs[ww]
			WAVE CurrSPhsWave=SPhsrefs[ww]
//			// now fit 2 Periods starting at 50 ms lag , do this a few times
//			// depending on how many periods fit in autocorrelation
//			rmax=min(5,floor(rightx(CurrWave)*currfreq/2)) 
//			squareAmp=0
//			FOR (rr=0; rr<rmax; rr++)
//				t1=50e-3+rr*2/currfreq
//				t2=t1+2/currfreq
//				MAKE/O/N=(4) fitCoefs={mean(CurrWave,t1,t2),0.5*(wavemax(CurrWave,t1,t2)-wavemin(CurrWave,t1,t2)),2*Pi*currfreq,0}
//				CurveFit/Q/G/H="0010"/N sin kwCWave=fitCoefs, currWave(t1,t2) /D 
//				squareAmp+=fitCoefs[1]
//
//			ENDFOR
//			squareAmp/=rmax
			// variant one:  use the results from the sine fit, 
			// but they have to be multiplied with sqrt(2)
		//	FrequencyList[ff][11]+=sqrt(squareAmp*2)//*DimSize(CurrSPhsWave,0)/FrequencyList[ff][1]
			// variant two: 		
			// take the sine amplitude from fft (multiplied with 2 divided by npnts)
			FrequencyList[ff][11]+=numberbyKey( "Sine Amp (A)",note(CurrSPhsWave),":","\r")
						// the Amp is weighed by spike contributions
			
			// now get duration from current wave
			WAVE CurrWave=Irefs[ww]
			totalDur+=rightx(CurrWave)
		ENDFOR
		FrequencyList[ff][11]/=nWaves
		FrequencyList[ff][12] =totalDur
		FrequencyList[ff][7,10]= FrequencyList[ff][q-4]*2/FrequencyList[ff][11]*FrequencyList[ff][1]/FrequencyList[ff][12]
		// gain is vector strength*2/sineamplitude*firing rate

	ENDFOR
	

KillWaves/Z fitCoefs
	
END

Function [Wave QuantileVS, Variable nSp, Variable meanVS, Variable meanPhs] BootstrapVSKernel(WAVE PhaseWave[,WAVE QuantileWave])
// here, the bootstrap is always performed with the complex vector strength vectors for each bootstrap sample being 
// projected onto the direction of the mean VS vector
// the option "do complex" (that existed before) is now set as standard and always implemented
// when the optional Amplitude wave is provided, averaging across all the phases is done differently:
// 
	IF (WaveDims(PhaseWave)>1)
		Abort "Bootstrap vector strength is not multi-dimensionality aware"
	ENDIF
	
	// bootstrap
	VARIABLE	isBalanced	= 0

	IF (paramisdefault(QuantileWave) )
		MAKE/O/N=3 TargetQuantiles={0.025,0.5,0.975}
		WAVE QuantileWave=TargetQuantiles
	ENDIF
	VARIABLE	nEntries=DimSize(PhaseWave,0)
	VARIABLE	nboots=100000
	MAKE/C/O/N=(nBoots)	BootStatsC
	
					
		IF (nEntries < 20)
			QuantileWave=NaN
			meanVS=NaN
			meanPhs=NaN
			Return [QuantileWave, nEntries, meanVS, meanPhs]
		ENDIF
			
		IF (isBalanced)
			IF (DimSize(PhaseWave,0) > 2^32 -1)
				Abort "Wave has too many entries for this balanced Bootstrap"
			ENDIF
			
			MAKE/U/I/N=(nBoots*nEntries) Indicies
			Indicies[]=mod(p,nEntries)
			Make/N=(NBoots*nEntries) rnd
			rnd=gnoise(1,2)
			Sort rnd, Indicies
			Redimension/N=(nEntries,nBoots) Indicies
			KillWaves/Z rnd
			
	
			MultiThread  	BootStatsC[]=ResampleAndGetVSandPhase(PhaseWave, indexWave=Indicies, clmn=p) // returns in polar form magnitude and angle
			KillWaves/Z Indicies
		ELSE		// not balanced
			MultiThread  	BootStatsC[]=ResampleAndGetVSandPhase(PhaseWave) // returns in polar form magnitude and angle
		ENDIF
		
	
		// get average vector from original sample
		VARIABLE/C avgVSandPhase= VSandPhasefromPhases(PhaseWave)
		meanVS = real(avgVSandPhase)
		meanPhs= imag(avgVSandPhase) // the average phase of the original sample
		// now project every entry in the Bootstrap sample onto this vector
		BootStatsC[]=cmplx(real(BootStatsC[p])*cos((imag(BootStatsC[p])-(meanPhs))),0) // the phase does not matter any more
		Redimension/R BootStatsC
		// now go ahead with the same statistics as for case without phases
		
		STRING NoteInc, NoteString=""
		IF (isBalanced)
			sprintf NoteInc,"Balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			NoteString+= NoteInc+"\r"
		ELSE
			sprintf NoteInc, "Non-balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			NoteString+= NoteInc+"\r"
		ENDIF
		
		sprintf NoteInc, "Bootstrap distribution has the quantiles\r\t%1.3f\t%1.3f\t%1.3f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
		NoteString+= NoteInc+"\r"
					
		QuantilesFromSample(BootStatsC,QuantileWave,1) 	// before call QuantileWave contains the values between 0 & 1 at which  quantiles should be 
														// calculated, after the call, it  

		KillWaves/Z BootStatsC

		Return [QuantileWave, nEntries, meanVS, meanPhs]
	
END

Function BootstrapVS(PhaseWave[,QuantileWave, quiet, doComplex])
WAVE		PhaseWave
WAVE		QuantileWave
VARIABLE	quiet, doComplex

// a typical use of this is with a accumulated Phase wave obtained through 
// CollectDataByNoteAndCriterium
// as in the following line (for 5 Hz sine wave)
// BootstrapVS(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal= 5, up4NoteVal=5),quiet=0,doComplex=1)
// the waves ending in "_SPhs" have been created previously by a call to 
// CreatePhaseWave



// use "docomplex=1" to get the correct confidence interval
	IF (WaveDims(PhaseWave)>1)
		Abort "Bootstrap vector strength is not multi-dimensionality aware"
	ENDIF
	
	// bootstrap
	VARIABLE	isBalanced	= 0
	IF (paramisdefault(quiet) )
		quiet=1
	ENDIF

	IF (paramisdefault(QuantileWave) )
		MAKE/O/N=3 TargetQuantiles={0.025,0.5,0.975}
		WAVE QuantileWave=TargetQuantiles
	ENDIF

	IF (paramisDefault(doComplex))
		doComplex = 1	
	ENDIF
	// one problem with the two-dimensional nature of the vector strength (magnitude, phase), is
	// that testing against zero amplitude is not done properly: 
	// especially for small numbers of entries, the average vector strength amplitude will always have 
	// a finite size (1/sqrt(N) for random phases)
	// and zero magnitude will typically be outside the MAGNITUDE distribution
	// In order to test properly, the bootstrap has to account for magnitude AND phase of the bootstrap sample
	// This is done in the following way: 
	// for the original sample, the average vector (magniutde AND phase) is determined
	// all bootstrap sample vectors are then PROJECTED onto this vectors axis
	// such that the resulting vector has the same direction. This means, if it points into the 
	// same direction (away from zero) as the average vector, it will have a positive length, otherwise a negative length
	// the boostrap statistics is the distribution of these PROJECTED lengths, not simply the magnitude

	VARIABLE	nEntries=DimSize(PhaseWave,0)
	VARIABLE	nboots=100000
	IF (!doComplex)
		MAKE/O/N=(nBoots)	BootStats
	ELSE
		MAKE/C/O/N=(nBoots)	BootStatsC
	ENDIF
					
		IF (nEntries < 2)
			QuantileWave=NaN
			return -1
		ENDIF
			
		IF (isBalanced)
			IF (DimSize(PhaseWave,0) > 2^32 -1)
				Abort "Wave has too many entries for this balanced Bootstrap"
			ENDIF
			
			MAKE/U/I/N=(nBoots*nEntries) Indicies
			Indicies[]=mod(p,nEntries)
			Make/N=(NBoots*nEntries) rnd
			rnd=gnoise(1,2)
			Sort rnd, Indicies
			Redimension/N=(nEntries,nBoots) Indicies
			KillWaves/Z rnd
			
	
			IF (!doComplex)
				MultiThread  	BootStats[]=ResampleAndGetVS(PhaseWave, indexWave=Indicies, clmn=p)
			ELSE
				MultiThread  	BootStatsC[]=ResampleAndGetVSandPhase(PhaseWave, indexWave=Indicies, clmn=p) // returns in polar form magnitude and angle
			ENDIF
			KillWaves/Z Indicies
		ELSE		// not balanced
			IF (!doComplex)
				MultiThread  	BootStats[]=ResampleAndGetVS(PhaseWave)
			ELSE
				MultiThread  	BootStatsC[]=ResampleAndGetVSandPhase(PhaseWave) // returns in polar form magnitude and angle
			ENDIF
		ENDIF
		
		IF (!doComplex)

			IF ((isBalanced) && !quiet)
				printf "Balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			ELSEIF (!quiet)
				printf "Non-balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			ENDIF
			
			IF (!quiet)
				Printf "Bootstrap distribution has the quantiles\r\t%1.3f\t%1.3f\t%1.3f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
			ENDIF
					
			QuantilesFromSample(Bootstats,QuantileWave,1)
			
			IF (!quiet)
				Printf "\r%f\r%f\r%f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
			ENDIF
			KillWaves/Z BootStats
		ELSE		// doComplex
			// get average vector from original sample
			VARIABLE avgPhase= imag(VSandPhasefromPhases(PhaseWave)) // the average phase of the original sample
			// now project every entry in the Bootstrap sample onto this vector
			BootStatsC[]=cmplx(real(BootStatsC[p])*cos((imag(BootStatsC[p])-(avgPhase))),0) // the phase does not matter any more
			Redimension/R BootStatsC
			// now go ahead with the same statistics as for case without phases
			
			IF ((isBalanced) && !quiet)
				printf "Balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			ELSEIF (!quiet)
				printf "Non-balanced bootstrap of %g data points with %g re-samples was done.\r",nEntries, nBoots
			ENDIF
			
			IF (!quiet)
				Printf "Bootstrap distribution has the quantiles\r\t%1.3f\t%1.3f\t%1.3f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
			ENDIF
					
			QuantilesFromSample(BootStatsC,QuantileWave,1)
			
			IF (!quiet)
				Printf "\r%f\r%f\r%f\r",QuantileWave[0],QuantileWave[1],QuantileWave[2]
			ENDIF
			KillWaves/Z BootStatsC

		ENDIF
END //BootstrapVS


ThreadSafe FUNCTION	ResampleAndGetVS(PhaseWave[, indexWave, clmn])
WAVE		PhaseWave
WAVE		IndexWave
VARIABLE	clmn			// from the  original wave (PhaseWave) a set of samples is drawn 
							// according to the indicies in the clmn-th column of indexwave 
							// this boostrap sample is then used to compute a vector strength
							// this scalar value is then returned

	Variable nEntries=DimSize(PhaseWave,0)
	IF (ParamIsDefault(clmn))
		Statsresample/N=(nEntries) PhaseWave
		WAVE	W_Resampled
	ELSE
		Duplicate/FREE PhaseWave, W_Resampled
		W_Resampled[]=PhaseWave[indexWave[p][clmn]]
		
	ENDIF
	
		Duplicate/FREE W_Resampled, realPart, cmplxPart
		Redimension/D realPart, cmplxPart
		
		realPart		=	cos(W_Resampled)
		cmplxPart	= 	sin(W_Resampled)
		
		VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
		VARIABLE	VS	= sqrt(rAVG^2+cAVG^2)
		
		KillWaves/Z realPart, cmplxPart, W_Resampled
	
		Return VS

	
END //ResampleAndGetVS

ThreadSafe FUNCTION/C	ResampleAndGetVSAndPhase(PhaseWave[, indexWave, clmn])
WAVE		PhaseWave
WAVE		IndexWave
VARIABLE	clmn			// from the  original wave (PhaseWave) a set of samples is drawn 
							// according to the indicies in the clmn-th column of indexwave 
							// this boostrap sample is then used to compute a vector strength and the average phase
							// those two are returned as real and imaginary part of a complex number

	Variable nEntries=DimSize(PhaseWave,0)
	IF (ParamIsDefault(clmn))
		Statsresample/N=(nEntries) PhaseWave
		WAVE	W_Resampled
	ELSE
		Duplicate/FREE PhaseWave, W_Resampled
		W_Resampled[]=PhaseWave[indexWave[p][clmn]]
		
	ENDIF
	
		Duplicate/FREE W_Resampled, realPart, cmplxPart
		Redimension/D realPart, cmplxPart
		
		realPart	=	cos(W_Resampled)
		cmplxPart	= 	sin(W_Resampled)
		
		VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
		VARIABLE/C	VSandPhase	= r2polar( cmplx(rAVG,cAVG) )
		
		KillWaves/Z realPart, cmplxPart, W_Resampled
	
		Return VSandPhase

	
END //ResampleAndGetVSAndPhase

ThreadSafe FUNCTION/C	Resample_AmpNormed_VSAndPhase(PhaseWave,AmpWave)
WAVE		PhaseWave
WAVE		AmpWave
		VARIABLE nEntries=DimSize(PhaseWave,0)
		Duplicate/O PhaseWave, Temp
		Redimension/N=(-1,2) Temp
		TEmp[][1]=AmpWave[p]
		
		Statsresample/MC/N=(nEntries) Temp
		WAVE	M_Resampled

		Duplicate/R=[*][0]/FREE M_Resampled, realPart, cmplxPart
		Redimension/D realPart, cmplxPart
		
		realPart[]	=	cos(M_Resampled[p][0])/M_Resampled[p][1]
		cmplxPart[]	= 	sin(M_Resampled[p][0])/M_Resampled[p][1]
		// this creates complex numbers that are scaled with 1/Amp
		
		VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
		VARIABLE/C	VSandPhase	= r2polar( cmplx(rAVG,cAVG) )
		
		KillWaves/Z realPart, cmplxPart, M_Resampled
	
		Return VSandPhase

	
END //Resample_AmpNormed_VSAndPhase

