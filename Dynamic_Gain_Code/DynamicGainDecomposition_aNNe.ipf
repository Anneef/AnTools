#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include "ANtools_extended"
#include "DynamicGainBootstrap_aNNe"
#include ":IFDL v4 Procedures:IFDL"

FUNCTION AssignThetaPhaseToAllSTInFolder()

	STRING	STListofFolder=WaveList("*_ST",";", "DIMS:1,MAXCOLS:1" )	// all the spike time waves 
																			// in the current datafolder 
																			// simply picked out by ending
																			
	STRING	IListofFolder= ReplaceString("_ST;" , STListofFolder, "_I;")
																			// this should all be existing current waves relatd to the 
																			// spike time waves in the original list
	VARIABLE	cW, nW=ItemsInList(STListofFolder,";")
	
	WAVE/Wave STrefs=ListToWaveRefWave(STListofFolder,0)
	WAVE/Wave Irefs=ListToWaveRefWave(IListofFolder,0)
	
	MAKE/O/N=(nW)/DF DummyReturn
	IF (Exists("root:BP_Butter_3p5to10p5_6pole")==1)
		WAVE FiltWave=root:BP_Butter_3p5to10p5_6pole 		// has to exist in the root
	ELSE
		Abort "The required filter wave was not found. Aborting"
	ENDIF
	MultiThread DummyReturn[]=AssignThetaPhaseToST_WORKER(STrefs[p],Irefs[p], FiltWave)
	// collect the folders in which the theta amplitudes and phase waves are
		STRING	ThetaPhaseName, ThetaAmpName, IName, IFiltName
		VARIABLE	INameLength												
	
	// now go through the folders and collect the waves 
		for(cW=0; cW<nW; cW+=1)
			DFREF df= DummyReturn[cW]
			IName=StringFromList(cW,IListofFolder,";")
			INameLength=strlen(IName)
			//IFiltName=Iname+"Filt"
			ThetaPhaseName=IName[0,INameLength-2]+"ThetaPhs"														
			ThetaAmpName=IName[0,INameLength-2]+"ThetaAmp"	
			Duplicate/O df:ThetaPhs, $ThetaPhaseName
			Duplicate/O df:ThetaAmp, $ThetaAmpName
			
		endfor
	
	KillWaves/Z DummyReturn, STrefs,Irefs
END



Threadsafe FUNCTION/DF AssignThetaPhaseToST_WORKER(STWave,IWave, FiltWave)
WAVE STWave, IWave
WAVE	FiltWave
	
	DFREF dfSav= GetDataFolderDFR()

	// Create a free data folder to hold the extracted and filtered plane 
	DFREF dfFree= NewFreeDataFolder()
	SetDataFolder dfFree

	// create the result wave that will remain in the free data folder until read out

	Duplicate/O  STWave, ThetaPhs, ThetaAmp
	WAVE	ThetaPhs
	WAVE	ThetaAmp
	
	
	// Filter the current in the band 3.5 to 10.5 Hz 
	// using a 6 pole Bessel Filter (Cascade)
	// once from start to 
	// the filter wave
	// HAS TO BE CORRECTLY MADE:
	// 1. a cascade filter 
	// 2. made for the correct sample frequency, i.e. the sample frequency of the current waves
	
	// The following is essentially a copy of the function 
	// IIRApplyFilterIIRcascadeCoefs(cascadeCoefs,input,outputName) but made in a threadsafe function
	// in order to speed things up a bit
		
	Variable nSections= DimSize(FiltWave,0)
		
	// We use a FIFO for the delayed intermediate w[n-1] and w[n-2] values.
	// The intermediate states are stored for each section.
	// This allows saving of state for segmenting the filtering.
	// This wave is the optional state input and output of Igor 5's FilterIIR operation.
	Make/O/D/N=(nSections,2) delayedIntermediate=0
	WAVE/D delayedIntermediate
	// delayedIntermediate[k][0] is w[n-1]
	// delayedIntermediate[k][1] is w[n-2]
	// where w[n]= ( x[n] - b1k*w[n-1] - b2k*w[n-2] ) / b0k

	Duplicate/O/D IWave,IFilt
	Wave output= IFilt

	Variable n, inLength=numpnts(IWave)
	for( n= 0; n < inLength; n += 1 )
		Variable sectionIndex, yn, xn=IWave[n]
		for( sectionIndex= 0; sectionIndex<nSections; sectionIndex+=1 )
			Variable a0=FiltWave[sectionIndex][0]
			Variable a1=FiltWave[sectionIndex][1]
			Variable a2=FiltWave[sectionIndex][2]
			Variable b0=FiltWave[sectionIndex][3]
			Variable b1=FiltWave[sectionIndex][4]
			Variable b2=FiltWave[sectionIndex][5]
			// compute wk[n]= ( x[n] - b1k*w[n-1] - b2k*w[n-2] ) / b0k
			Variable wn1= delayedIntermediate[sectionIndex][0]
			Variable wn2= delayedIntermediate[sectionIndex][1]
			Variable wn= (xn - b1*wn1 - b2*wn2 ) / b0
			// prepare for next filtering operation by shifting w's
			delayedIntermediate[sectionIndex][0]=wn
			delayedIntermediate[sectionIndex][1]=wn1
			// compute this sections's output (next section's input)
			// yk[n]= a0k*w[n] + a1k*w[n-1] + a2k*w[n-2]
			yn= a0*wn + a1*wn1 + a2*wn2
			xn= yn	// cascade
		endfor	// end section loop
		output[n]= yn
	endfor		// end input loop
	
	// end of copy of 
	// IIRApplyFilterIIRCascadeCoefs(cascadeCoefs,input,outputName)

	VARIABLE	nPnts=DimSize(IFilt,0)
	Duplicate/FREE/O/D IFilt, tempAmp, tempPhs
	tempAmp[0,npnts-1]=IFilt[nPnts-1-p]	// reverse time and apply filter again
	
	inLength=numpnts(tempAmp)
	for( n= 0; n < inLength; n += 1 )
		xn=tempAmp[n]
		for( sectionIndex= 0; sectionIndex<nSections; sectionIndex+=1 )
			a0=FiltWave[sectionIndex][0]
			a1=FiltWave[sectionIndex][1]
			a2=FiltWave[sectionIndex][2]
			b0=FiltWave[sectionIndex][3]
			b1=FiltWave[sectionIndex][4]
			b2=FiltWave[sectionIndex][5]
			// compute wk[n]= ( x[n] - b1k*w[n-1] - b2k*w[n-2] ) / b0k
			wn1= delayedIntermediate[sectionIndex][0]
			wn2= delayedIntermediate[sectionIndex][1]
			wn= (xn - b1*wn1 - b2*wn2 ) / b0
			// prepare for next filtering operation by shifting w's
			delayedIntermediate[sectionIndex][0]=wn
			delayedIntermediate[sectionIndex][1]=wn1
			// compute this sections's output (next section's input)
			// yk[n]= a0k*w[n] + a1k*w[n-1] + a2k*w[n-2]
			yn= a0*wn + a1*wn1 + a2*wn2
			xn= yn	// cascade
		endfor	// end section loop
		output[n]= yn
	endfor		// end input loop
	KillWaves/Z delayedIntermediate	
	
	tempAmp=IFilt
	IFilt[0,npnts-1]=tempAmp[nPnts-1-p]	// back to causal time 

		
	// perform Hilbert transform, create "analytic signal" 
	// and derive instantaneous phase and amplitude
	
	HilbertTransform /DEST=DummyH IFilt
	tempAmp=sqrt(IFilt^2 + DummyH^2)
	tempPhs=atan2(-DummyH,IFilt)
	
	ThetaPhs[]=tempPhs[x2pnt(tempPhs, STWave[p])]
	ThetaAmp[]=tempAmp[x2pnt(tempAmp, STWave[p])]
		
	KillWaves/Z DummyH, tempAmp, tempPhs, IFilt
	return dfFree


END //AssignThetaPhaseToST_WORKER


Function AssignThetaPhaseToST(W_ST)
WAVE	W_ST		// spike time wave

// get the current wave belonging to this spike times
	STRING	STName=NameOfWave(W_ST), IName, IFiltName, ThetaPhaseName, ThetaAmpName
	VARIABLE	STNameLength=strlen(STName)
	
	IF 	(StringMatch(STName[STNameLength-3,STNameLength-1],"_ST" ))	// last to characters
																// in STName are "_ST"
		IName=STName[0,STNameLength-3]+"I"														
		IfiltName=STName[0,STNameLength-3]+"IFilt"														
		ThetaPhaseName=STName[0,STNameLength-3]+"ThetaPhs"														
		ThetaAmpName=STName[0,STNameLength-3]+"ThetaAmp"														

		WAVE I_Wave=$IName
		Duplicate/O/D  W_ST, $ThetaPhaseName, $ThetaAmpName
		WAVE	ThetaPhs=$ThetaPhaseName
		WAVE	ThetaAmp=$ThetaAmpName
		
		
		IF (!(Exists(IName)==1))
			// spike time wave does not yet exist
			Abort "Cannot find input for "+STName
		ENDIF
				
	ELSE
		Abort "Cannot decode names as spike time trace does not end in 'ST'"
	ENDIF
	
	// Filter the current in the band 3.5 to 10.5 Hz 
	// using a 6 pole Bessel Filter (Cascade)
	// once from start to 
	WAVE	FiltWave=root:BP_Butter_3p5to10p5_6pole 		// has to exist in the root
	// should be made for the correct sample frequency (here 100 kHz = 1e+5 Hz)
	IIRApplyFilterIIRCoefs(FiltWave,I_wave,IfiltName)
	WAVE IFilt=$IfiltName

	VARIABLE	nPnts=DimSize(IFilt,0)
	Duplicate/FREE/O/D IFilt, tempAmp, tempPhs
	tempAmp[0,npnts-1]=IFilt[nPnts-1-p]	// reverse time and apply filter again
	IIRApplyFilterIIRCoefs(FiltWave,tempAmp,IfiltName)
	tempAmp=IFilt
	IFilt[0,npnts-1]=tempAmp[nPnts-1-p]	// back to causal time 

		
	// perform Hilbert transform, create "analytic signal" 
	// and derive instantaneous phase and amplitude
	
	HilbertTransform /DEST=DummyH IFilt
	tempAmp=sqrt(IFilt^2 + DummyH^2)
	tempPhs=atan2(-DummyH,IFilt)
	
	ThetaPhs[]=tempPhs[x2pnt(tempPhs, W_ST[p])]
	ThetaAmp[]=tempAmp[x2pnt(tempAmp, W_ST[p])]

		
	// define 3 strata of data points: 
	// 1. time points during which the 4-10 Hz filtered input was in it highest third of values 
	// 	above its 66 percentile
	// 2. time points during which the 4-10 Hz filtered was in it lowest third of values
	// 	i.e. below its 33 percentile
	// 3. the middle values 
	
	
	// for every spike (in spike time wave)	
	// find out where the FILTERED input was at the spike time, in the highest 
	// , lowest or middle third of the data range
	
	
//	MAKE/O/N=2 TargetQuantiles={1/3,2/3}
//	QuantilesFromSample(I_wave,TargetQuantiles,1)	// written by AN, see ANTools.ipf (or Bootstrap.ipf)
	
	
	
//	ThetaPhs=0
//	ThetaPhs[]=( IFilt(W_ST[p])  < TargetQuantiles[0]) ? -1 : ThetaPhs[p]
//	ThetaPhs[]=( IFilt(W_ST[p])  > TargetQuantiles[1]) ?  1 : ThetaPhs[p]
	
	KillWaves/Z TargetQuantiles, DummyH, tempAmp, tempPhs, IFilt
	Return 1
END




FUNCTION StratByThetaPhase(n_strat, start_phase, identifier,SineFreq, doDisplay)
VARIABLE		n_strat, Start_Phase
STRING		Identifier
VARIABLE		SineFreq, doDisplay

		IF (n_strat < 2)
			Abort "At least 2 groups!"
		ENDIF
		// the function creates a wave that holds vector strength values (and CI and spike counts)
		// each of the n_strat column contains these data for a range of theta phase, together covering [-Pi, Pi]
		// the sineFreq identifies the waves by the frequency of the sine wave that was embedded in the stimulus
		
		// first: define the percentiles that determine the CI
		
		DFREF		rootFolderRf=GetDataFolderDFR()
		STRING		rootDFPath
		rootDFPath = GetDataFolder(1)		// OK
	
		
		VARIABLE	upperCI=0.975
		VARIABLE	lowerCI=0.025
		VARIABLE binWidth=2*Pi/n_strat
		MAKE/O/N=3 TrgtQuant={0.5,lowerCI,upperCI}
		
		// the first bin starts at start_phase
		// bring startphase to within [-Pi, Pi]
		start_phase= start_phase - floor((start_phase+Pi)/2/Pi)*2*Pi 
		MAKE/O/N=(n_strat+1) binBounds
		binBounds[]=start_phase+p*binWidth
		
		// now deal with possible conflicts, i.e. if one of the bins contains the boundary of the cyclic range 
		// (starts at leftBound>0, ends at rightBound<0) 
		
		// first bring all boundaries in the range
		binBounds = binBounds + Pi
		binBounds = binBounds - floor(binBounds/2/Pi)*2*Pi
		binBounds = binBounds - Pi
		
		// create Wave to hold results
		
			STRING	ResName= "BootStr_TPhase_"+identifier
			
			MAKE/O/N=(n_strat,5) $ResName
			WAVE	Result = $ResName
			
			SetDimLabel 1,0,Median,Result
			SetDimLabel 1,1,lowerCI,Result
			SetDimLabel 1,2,upperCI,Result
			SetDimLabel 1,3,nEntries,Result
			SetDimLabel 1,4,avgPhase,Result
			
			Note/K    Result, "Phase:"+num2str(binBounds[0])+"\r"
			Note/NOCR Result, "binWidth:"+num2str(binWidth) +"\r"
			
		// Fill wave
			VARIABLE ss // current stratum
			FOR (ss=0; ss<n_strat; ss+=1)
				// check whether stratum contains boundary of the cyclic range ( i.e. Pi)
				IF  ( (binBounds[ss+1] <= 0) && (binBounds[ss] > 0) )
				// need to join two datasets
				WAVE w1=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=-Pi,upCrit=binBounds[ss+1])
				WAVE w2=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=Pi)
				// do the joining in a free DF 
					SetDataFolder NewFreeDataFolder()
					Concatenate/NP {w1,w2}, joined
					
					TrgtQuant={0.5,lowerCI,upperCI}
					
					BootstrapVS(joined, QuantileWave=TrgtQuant, quiet=1, doComplex=1 )
					Result[ss][0,2]=TrgtQuant[q]
					Result[ss][3]	= DimSize(joined,0)
					Result[ss][4]	= imag(VSandPhasefromPhases(joined))
					
					KillWaves/Z joined
					SetDataFolder rootFolderRf
				ELSE
					TrgtQuant={0.5,lowerCI,upperCI}
					BootstrapVS(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1]), QuantileWave=TrgtQuant, quiet=1, doComplex=1 )
					Result[ss][0,2]=TrgtQuant[q]
					Result[ss][3]	= DimSize(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1]),0)
					Result[ss][4]	= imag(VSandPhasefromPhases(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1])))
				ENDIF
				
					
			ENDFOR
			
		// Display result if that was wanted	
			IF (DoDisplay)
				MAKE/O/T/N=(n_strat+1) $("root:ThetaGroupLabels"+Identifier)
				WAVE/T	labels=$("root:ThetaGroupLabels"+Identifier)
				MAKE/O/N=(n_strat+1) $("root:ThetaGroupTicks"+Identifier)
				WAVE	tickw=$("root:ThetaGroupTicks"+Identifier)
				tickw=binBounds
				STRING	S_Dummy
				FOR (ss=0; ss<=n_strat; ss+=1)
					sprintf S_Dummy,"%1.2f",abs(tickw[ss])
					labels[ss]=S_Dummy
					IF ( abs(abs(tickw[ss])-Pi) < 1e-3)
						labels[ss]="π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/2)<1e-3)
						labels[ss]="½π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/3)<1e-3)
						labels[ss]="⅓π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/4)<1e-3)
						labels[ss]="¼π"
					ELSEIF ( abs(abs(tickw[ss])-3*Pi/4)<1e-3)
						labels[ss]="¾π"
					ELSEIF ( abs(abs(tickw[ss])-3*Pi/4)<1e-3)
						labels[ss]="⅔π"
					ELSEIF ( abs(tickw[ss]) < 1e-3)
						labels[ss]="0"
					ENDIF
						
					IF (tickw[ss]<0)
					 labels[ss]="-"+labels[ss]
					ENDIF

				ENDFOR 
				
			
				STRING	Name =NameOfWave(Result)
				STRING	Name1=NameOfWave(Result)+"#1"
				STRING	Name2=NameOfWave(Result)+"#2"
				Display /W=(268.2,66.2,636.6,258.8) Result[0,n_strat-1][1] as "XFC "+Identifier+" "+GetDataFolder(0)
				AppendToGraph Result[0,n_strat-1][2], Result[0,n_strat-1][0]
				ModifyGraph userticks(bottom)={tickw,labels}
				ModifyGraph margin(top)=28,margin(right)=3,margin(left)=28,width=141.732,height=127.559
				ModifyGraph mode($Name)=5,mode($Name1)=5
				ModifyGraph mode($Name2)=3
				ModifyGraph marker($Name2)=9
				ModifyGraph lOptions($Name2)=2
				ModifyGraph rgb=(0,0,0)
				ModifyGraph msize($Name2)=10
				ModifyGraph mrkThick($Name2)=3
				ModifyGraph hbFill($Name)=2
				ModifyGraph usePlusRGB($Name)=1
				ModifyGraph useNegRGB($Name)=1
				ModifyGraph plusRGB($Name)=(0,0,0,13107)
				ModifyGraph negRGB($Name)=(0,0,0,13107)
				ModifyGraph toMode($Name)=1, toMode($Name1)=1
				ModifyGraph offset($Name)={wavemin(binBounds),0},offset($Name1)={wavemin(binBounds),0}
				ModifyGraph offset($Name2)={wavemin(binBounds)+binwidth/2,0}
				ModifyGraph muloffset={binwidth,0}
				ModifyGraph btLen=3,ftLen(left)=3, standoff(bottom)=0
				ModifyGraph rgb($Name)=(34952,34952,34952),useBarStrokeRGB($Name)=1,barStrokeRGB($Name)=(0,0,0,0)
				ModifyGraph useBarStrokeRGB($Name1)=1,barStrokeRGB($Name1)=(0,0,0,0)
				ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={9,5}
				Label left "Vector strength @ "+num2str(SineFreq)+" Hz"
				Label bottom "Theta Phase"
				SetAxis left 0,1
				STRING	textName
				STRING	readout
				FOR (ss=0; ss<n_strat; ss+=1)
					textName="scount"+num2istr(ss)
					TextBox/C/N=$textName/F=0/B=1/A=MC/X=(-50+100/(n_strat-1)/3+ss*100/(n_strat))/Y=53 "\\Zr080"+num2istr(Result[ss][%nEntries])
					textName="avgPhase"+num2istr(ss)
					sprintf readout,"%-1.2f",Result[ss][%avgPhase]
					TextBox/C/N=$textName/F=0/B=1/A=MC/X=(-50+100/(n_strat-1)/3+ss*100/(n_strat))/Y=-45 "\\Zr080"+readout

				ENDFOR
				TextBox/C/N=title/F=0/B=1/A=MC/X=-6/Y=60 GetDataFolder(1)

			ENDIF
END //StratByThetaPhase

FUNCTION StratByThetaPhaseAndAmp(n_strat, start_phase,CritSuffix, identifier,SineFreq, Min_ThetaAmp, Max_ThetaAmp, doDisplay)
VARIABLE		n_strat, Start_Phase
STRING		Identifier, CritSuffix
VARIABLE		SineFreq, doDisplay, Min_ThetaAmp, Max_ThetaAmp

		IF (n_strat < 2)
			Abort "At least 2 groups!"
		ENDIF
		// the function creates a wave that holds vector strength values (and CI and spike counts)
		// each of the n_strat column contains these data for a range of theta phase, together covering [-Pi, Pi]
		// the sineFreq identifies the waves by the frequency of the sine wave that was embedded in the stimulus
		// this is included in the wave note of the SPhs waves
		
		// the limits of ThetaAmp exclude spikes for which the amplitude of the theta component was not in the correct range
		
		// first: define the percentiles that determine the CI
		
		DFREF		rootFolderRf=GetDataFolderDFR()
		STRING		rootDFPath
		rootDFPath = GetDataFolder(1)		// OK
	
		
		VARIABLE	upperCI=0.975
		VARIABLE	lowerCI=0.025
		VARIABLE binWidth=2*Pi/n_strat
		MAKE/O/N=3 TrgtQuant={0.5,lowerCI,upperCI}
		
		// the first bin starts at start_phase
		// bring startphase to within [-Pi, Pi]
		start_phase= start_phase - floor((start_phase+Pi)/2/Pi)*2*Pi 
		MAKE/O/N=(n_strat+1) binBounds
		binBounds[]=start_phase+p*binWidth
		
		// now deal with possible conflicts, i.e. if one of the bins contains the boundary of the cyclic range 
		// (starts at leftBound>0, ends at rightBound<0) 
		
		// first bring all boundaries in the range
		binBounds = binBounds + Pi
		binBounds = binBounds - floor(binBounds/2/Pi)*2*Pi
		binBounds = binBounds - Pi
		
		// create Wave to hold results
		
			STRING	ResName= "BootStr_TPhase_"+identifier
			
			MAKE/O/N=(n_strat,5) $ResName
			WAVE	Result = $ResName
			
			SetDimLabel 1,0,Median,Result
			SetDimLabel 1,1,lowerCI,Result
			SetDimLabel 1,2,upperCI,Result
			SetDimLabel 1,3,nEntries,Result
			SetDimLabel 1,4,avgPhase,Result
			
			Note/K    Result, "Phase:"+num2str(binBounds[0])+"\r"
			Note/NOCR Result, "binWidth:"+num2str(binWidth) +"\r"
			
		// Fill wave
			VARIABLE ss // current stratus
			FOR (ss=0; ss<n_strat; ss+=1)
				// check whether status contains boundary of the cyclic range ( i.e. Pi)
				IF  ( (binBounds[ss+1] <= 0) && (binBounds[ss] > 0) )
				// need to join two datasets
				WAVE w1=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=-Pi,upCrit=binBounds[ss+1], CriteriumSuffix1=CritSuffix, lowCrit1=MIN_ThetaAmp, upCrit1=MAX_ThetaAmp)
				WAVE w2=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=Pi, CriteriumSuffix1=CritSuffix, lowCrit1=MIN_ThetaAmp, upCrit1=MAX_ThetaAmp)
				// do the joining in a free DF 
					SetDataFolder NewFreeDataFolder()
					Concatenate/NP {w1,w2}, joined
					
					TrgtQuant={0.5,lowerCI,upperCI}
					
					BootstrapVS(joined, QuantileWave=TrgtQuant, quiet=1 )
					Result[ss][0,2]=TrgtQuant[q]
					Result[ss][3]	= DimSize(joined,0)
					Result[ss][4]	= imag(VSandPhasefromPhases(joined))
					
					KillWaves/Z joined
					SetDataFolder rootFolderRf
				ELSE
					TrgtQuant={0.5,lowerCI,upperCI}
					BootstrapVS(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1], CriteriumSuffix1=CritSuffix, lowCrit1=MIN_ThetaAmp, upCrit1=MAX_ThetaAmp), QuantileWave=TrgtQuant, quiet=1 )
					Result[ss][0,2]=TrgtQuant[q]
					Result[ss][3]	= DimSize(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1], CriteriumSuffix1=CritSuffix, lowCrit1=MIN_ThetaAmp, upCrit1=MAX_ThetaAmp),0)
					Result[ss][4]	= imag(VSandPhasefromPhases(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1], CriteriumSuffix1=CritSuffix, lowCrit1=MIN_ThetaAmp, upCrit1=MAX_ThetaAmp)))
				ENDIF
				
					
			ENDFOR
			
		// Display result if that was wanted	
			IF (DoDisplay)
				MAKE/O/T/N=(n_strat+1) $("root:ThetaGroupLabels"+Identifier)
				WAVE/T	labels=$("root:ThetaGroupLabels"+Identifier)
				MAKE/O/N=(n_strat+1) $("root:ThetaGroupTicks"+Identifier)
				WAVE	tickw=$("root:ThetaGroupTicks"+Identifier)
				tickw=binBounds
				STRING	S_Dummy
				FOR (ss=0; ss<=n_strat; ss+=1)
					sprintf S_Dummy,"%1.2f",abs(tickw[ss])
					labels[ss]=S_Dummy
					IF ( abs(abs(tickw[ss])-Pi) < 1e-3)
						labels[ss]="π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/2)<1e-3)
						labels[ss]="½π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/3)<1e-3)
						labels[ss]="⅓π"
					ELSEIF ( abs(abs(tickw[ss])-Pi/4)<1e-3)
						labels[ss]="¼π"
					ELSEIF ( abs(abs(tickw[ss])-3*Pi/4)<1e-3)
						labels[ss]="¾π"
					ELSEIF ( abs(abs(tickw[ss])-3*Pi/4)<1e-3)
						labels[ss]="⅔π"
					ELSEIF ( abs(tickw[ss]) < 1e-3)
						labels[ss]="0"
					ENDIF
						
					IF (tickw[ss]<0)
					 labels[ss]="-"+labels[ss]
					ENDIF

				ENDFOR 
				
			
				STRING	Name =NameOfWave(Result)
				STRING	Name1=NameOfWave(Result)+"#1"
				STRING	Name2=NameOfWave(Result)+"#2"
				Display /W=(268.2,66.2,636.6,258.8) Result[0,n_strat-1][1] as "XFC "+Identifier+" "+GetDataFolder(0)
				AppendToGraph Result[0,n_strat-1][2], Result[0,n_strat-1][0]
				ModifyGraph userticks(bottom)={tickw,labels}
				ModifyGraph margin(top)=28,margin(right)=3,margin(left)=28,width=141.732,height=127.559
				ModifyGraph mode($Name)=5,mode($Name1)=5
				ModifyGraph mode($Name2)=3
				ModifyGraph marker($Name2)=9
				ModifyGraph lOptions($Name2)=2
				ModifyGraph rgb=(0,0,0)
				ModifyGraph msize($Name2)=10
				ModifyGraph mrkThick($Name2)=3
				ModifyGraph hbFill($Name)=2
				ModifyGraph usePlusRGB($Name)=1
				ModifyGraph useNegRGB($Name)=1
				ModifyGraph plusRGB($Name)=(0,0,0,13107)
				ModifyGraph negRGB($Name)=(0,0,0,13107)
				ModifyGraph toMode($Name)=1, toMode($Name1)=1
				ModifyGraph offset($Name)={wavemin(binBounds),0},offset($Name1)={wavemin(binBounds),0}
				ModifyGraph offset($Name2)={wavemin(binBounds)+binwidth/2,0}
				ModifyGraph muloffset={binwidth,0}
				ModifyGraph btLen=3,ftLen(left)=3, standoff(bottom)=0
				ModifyGraph rgb($Name)=(34952,34952,34952),useBarStrokeRGB($Name)=1,barStrokeRGB($Name)=(0,0,0,0)
				ModifyGraph useBarStrokeRGB($Name1)=1,barStrokeRGB($Name1)=(0,0,0,0)
				ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={9,5}
				Label left "Vector strength @ "+num2str(SineFreq)+" Hz"
				Label bottom "Theta Phase"
				SetAxis left 0,1
				STRING	textName
				STRING	readout
				FOR (ss=0; ss<n_strat; ss+=1)
					textName="scount"+num2istr(ss)
					TextBox/C/N=$textName/F=0/B=1/A=MC/X=(-50+100/(n_strat-1)/3+ss*100/(n_strat))/Y=53 "\\Zr080"+num2istr(Result[ss][%nEntries])
					textName="avgPhase"+num2istr(ss)
					sprintf readout,"%-1.2f",Result[ss][%avgPhase]
					TextBox/C/N=$textName/F=0/B=1/A=MC/X=(-50+100/(n_strat-1)/3+ss*100/(n_strat))/Y=-45 "\\Zr080"+readout

				ENDFOR
				TextBox/C/N=title/F=0/B=1/A=MC/X=-6/Y=60 GetDataFolder(1)
				TextBox/C/N=AmpCrit/F=0/B=1/A=MC/X=-6/Y=36 num2str(Min_ThetaAmp)+"<"+CritSuffix+"<"+num2str(Max_ThetaAmp)

			ENDIF
END



FUNCTION SubThreshGain(InWave_V, InWave_I, InWave_ST, STA_Dur,[V_thresh,deadtime])
WAVE		InWave_V, InWave_I, InWave_ST
VARIABLE	STA_Dur								// Duration of the STA used for the Gain calculation
												// needed only to adjust the window duration to this value
VARIABLE	V_thresh							// voltages above this are replaced by threshold (in V)
VARIABLE	deadtime							// if given, the measurement values are replaced with a 
												// straight line, from the point of crossing threshold, 
												// to the point "deadtime" seconds later
																										
// takes current and voltage traces, removes voltage values above the spike threshold (average of the 
// specific voltage trace as specified by the spike property assessment )
// current and voltage traces (after truncation) are cut into segments of the same length as the spike triggered average (here: 1s)
// standard Fourier transform approach to determine the transfer function between current and voltage
// in the present version the phase is inverted (missing conjugated complex)
														
														
	VARIABLE	WinLen=STA_Dur/DimDelta(InWave_V,0)
	// winlen has to be an even number (FFT)
	WinLen=round(1/DimDelta(InWave_I,0)/2-0.4)*2	// there were suprising round off errors when the simple floor(N/2)*2 was used
	VARIABLE	nWins= floor(DimSize(InWave_V,0)/WinLen)
	VARIABLE	nSpikes=DimSize(InWave_ST,0), sp
	// Duplicate inputs and output wave; they will be reshaped into a matricies and then overwritten with Fourierresults
	Duplicate/O  Inwave_V, Dummy_V
	Duplicate/O  Inwave_I, Dummy_I
	// Health check for equal length
	IF (DimSize(InWave_I,0) != DimSize(InWave_V,0))
		DoAlert 0,"Length In_V != Length In_I, Aborting"
		Return -1
	ENDIF
	
															
	// remove the spikes by replacing values > V_tresh with V_thresh
	IF (ParamIsDefault(V_thresh) )
		V_thresh=str2num(StringByKey("V/s)", StringFromList(1, note(InWave_ST),"\r" ),":"," "))
	ENDIF
	IF (ParamIsDefault(deadtime) ) // no deadtime given
		Dummy_V[]= (Dummy_V[p] > V_thresh) ? V_thresh : Dummy_V[p] // this replaces the high values
	ELSE
		VARIABLE deadpnts=round(deadtime/DimDelta(InWave_V,0))
		// have to go spike by spike and draw this line
		// first, get a list of times at which the threshold is crossed BEFORE each spike
		SpikesAtVdetectMultiThread(Inwave_V, InWave_ST,V_thresh,Suffix="_4Susz")
		// creates the wanted wave (see last comment), the name is the name of Inwave_ST, with additional suffix "_4subG"
		WAVE StartPoints=$(NameOfWave(Inwave_ST)+"_4Susz")
		// problems arise, if this wave contains NAN
		// this would happen, if the threshold was not crossed, i.e. only before a previous spike occured
		// if that happends, we have to start cutting out spikes from dummy_i and dummy_V
		
		VARIABLE p1,p2,y1,y2, slope
		FOR (sp=0; sp< nSpikes-1; sp++)
			If  (numtype(StartPoints[sp])==2) // is a NaN
				Abort "the threshold causes spike no. "+num2istr(sp)+" in wave"+nameofwave(InWave_V)+" to be skipped, that is not caught"
			ENDIF
			p1=x2pnt(Dummy_V,StartPoints[sp])
			p2=p1+deadpnts
			y1=Inwave_V[p1]
			slope=(Inwave_V[p2]-Inwave_V[p1])/(deadpnts)			
			Dummy_V[p1+1,p2-1]=y1+slope*(p-p1)
		ENDFOR
		// last round separately, because the last point of the wave might be further away than deadtime
		p1=x2pnt(Dummy_V,StartPoints[nspikes-1])
		IF (p1+deadpnts>DimSize(InWave_V,0)) // not enough points to find endpoint of line -> ise straight line
			Dummy_V[p1+1,DimSize(InWave_V,0)-1]=Inwave_V[p1]
		ELSE
			p2=p1+deadpnts
			y1=Inwave_V[p1]
			slope=(Inwave_V[p2]-Inwave_V[p1])/(deadpnts)			
			Dummy_V[p1+1,p2-1]=y1+slope*(p-p1)
		ENDIF

	ENDIF 
	// cut length to Matrix Nrows x Mcols  with N= winlen; M=nWins
	Redimension/D/N=(WinLen,nWins) Dummy_I, Dummy_V
	FFT/COLS /DEST=Dummy_V_FFT Dummy_V
	FFT/COLS /DEST=Dummy_I_FFT Dummy_I
	
	KillWaves Dummy_I,Dummy_V, StartPoints
	// by pure coincidence, the magnitude of the current component in any one of the higher frequencies
	// could turn out to be very very small. This would then cause a very large gain, purely due the variation in the 
	// magnitudes due to noise. 
	// this should be avoided
	// therefore: check the current FFT for unusually small magnitudes
	
	// find out whether there is such a small value

	MAKE/O/N=(round(2000*STA_Dur),DimSize(Dummy_I_FFT,1)) Dummy 	// only need to reduce the noise for frequencies below 2kHz, 
																				// those are in the first 2000*Sta_dur rows
	VARIABLE reps=0
	Dummy=cabs(Dummy_I_FFT)	// put the values of the (updated) I_FFT matrix into dummy (for all freqs up to 2000 Hz)
	DO
		wavestats/Q/M=1 Dummy
		Dummy_I_FFT[V_minRowLoc][V_minColLoc]+=cmplx(1e-9,1e-9)		// replace the entry with the smallest magnitude by 
																	// a magnitude of 1e-9. The phase does not matter much, 
																	// because the magnitude is so small
		Dummy[V_minRowLoc][V_minColLoc]+=1e-9	
		reps+=1
	WHILE ( (reps < 2000) && (V_min<1e-9))

	// calculate complex gain
	Dummy_V_FFT/= Dummy_I_FFT
	KillWaves Dummy_I_FFT
	KillWaves/Z Dummy

	STRING TargetName=NameOfWave(InWave_V)
	TargetName=TargetName[0,strlen(TargetName)-2]+"subG"
	// Average across columns
	MatrixOp/C/NTHR=3  /O  $TargetName = sumRows(Dummy_V_FFT)/numCols(Dummy_V_FFT) 	
	KillWaves Dummy_V_FFT

END // SubThreshGain

FUNCTION SpikesAtVdetect(InWave_V, InWave_ST,V_detect[,Suffix])
WAVE	InWave_V, InWave_ST		        // voltage and spike-time wave
VARIABLE	V_detect					// voltage threshold value to define modified spike time
										// in Volts, will be rounded to full millivolts
STRING Suffix							// if given, the new spike time wave carries this suffic after the name of InWave_ST

//	Variable timerRefNum, microsecs
//	timerRefNum = StartMSTimer
//	if (timerRefNum == -1)
//		Abort "All timers are in use"
//	endif

//V_detect = round(1000*V_detect)/1000 	// rounded to full millivolts

// travels back in time from the original spike time, 
// to find the latest time voltage crosses V_detect before the full blown spike

// complications: the time should be later than the previous full blown spike
// if the t_sp_detect[n] does fall before the previous full blown spike
// i.e. if t_sp_detect[n]<t_sp[n+1] that will be noted by putting a NaN in the result wave

// the function returns a wave with one entry for each spike in InWave_ST
// the result's name has a suffix derived from the threhsold. 
// therefore the threshold voltage is rounded to full millivolts

		VARIABLE		dt=DimDelta(InWave_V,0)
		
// health test: Spike times should define a narrow set of voltages in the voltage trace
// (assuming that a VOLTAGE threshold rather than a rate threshold was used as definition
// here I use that test to see whether the provided voltage wave and spike time wave are plausibly a valid pair
		Duplicate /O InWave_ST, Dummy
		WAVE 	dummy
		dummy[]=InWave_V[x2pnt(InWave_V,InWave_ST[p])]
		WAVESTATS/Q dummy
		VARIABLE	expectedVjitter=1500*dt/2*2 // dV/dt max * dt/2  * safety factor of 1.5
		IF ( (V_max-V_min) > expectedVjitter)
			DoAlert 0,"The range of threshold voltages related to the spike times\ris "+num2str((V_max-V_min)*1000)+" mV.\rThat seems too big, aborting..."
			KillWaves dummy
			Return -1
		
		ENDIF
		
		// create name for output
		STRING		SPdetectName
		IF (ParamIsDefault(Suffix))
			SPdetectName=NameOfWave(InWave_ST)+"_at_"+ReplaceString("-",  num2istr(V_detect*1000), "neg")
			SPdetectName=ReplaceString(".", SPdetectName, "p")
		ELSE
			SPdetectName=NameOfWave(InWave_ST)+Suffix
		ENDIF

		Duplicate /O InWave_ST, $SPdetectName
		WAVE	STdet=$SPdetectName
		SetScale d 0,0,"s", STdet
		SetScale/P x,0,1,"", STdet
		// has correct unit (seconds) and index scaling in x
		
		// now get a stretch of data before the original spike time
		// possibly upsample
		
		VARIABLE		nSpikes=DimSize(InWave_ST,0), sp, sPnt// Pnt number in Voltage wave where spike was detected
		VARIABLE	 	nPntsBefore, upsamplefac=round(100000*dt)// stretch of data before and after that is copied for analysis
																			// factor of how much the temporal sampling is increased
		VARIABLE		nPntsInDummy
		FOR (sp=0; sp< nSpikes; sp +=1)
			
						
			IF (sp >0)	// not the first spike)
				nPntsBefore=(InWave_ST[sp]-InWave_ST[sp-1])/dt
			ELSE
				nPntsBefore=InWave_ST[sp]/dt -1
			ENDIF
			
			sPnt=x2pnt(InWave_V,InWave_ST[sp])
			DUPLICATE/O/R=[sPnt-nPntsBefore,sPnt-1] InWave_V, dummy			
			// at this point upsample dummy and perform all further analysis on the upsampled versions
			Duplicate/O dummy, dummy_US
			Resample/up=(upsamplefac) dummy_US	// this does a sinc reconstruction
			WAVE	dummy_US					// upsampled version covering the exact same time
			
			nPntsInDummy=DimSize(dummy_US,0)
			FindLevel /EDGE=1 /P/Q/R=[nPntsInDummy-1,0] dummy_US,V_detect
			IF (V_Flag!=0)	// NO level crossing found
			   // => the wanted crossing lies even before the previous spike
			   // so put a NaN in the spike time
			   STdet[sp] = NaN
			ELSE
				STdet[sp] = pnt2x(dummy_US,floor(V_LevelX)) + dt/upsamplefac*(V_LevelX-floor(V_LevelX))
			ENDIF
			
	   ENDFOR	
	   
//	   	microsecs = StopMSTimer(timerRefNum)
//	Print microsecs/DimSize(InWave_ST,0), "microseconds per spike"


	   KillWaves dummy_US, dummy		
END // SpikesAtVdetect

Threadsafe FUNCTION WalkBack(voltageWav, V_detect, spikepoint, pointsbefore, dt)
WAVE			voltageWav	// wave in which to walk back to find
VARIABLE		V_detect
VARIABLE		spikepoint	// point n voltageWave where to start
VARIABLE		pointsbefore // max number of points to walk back
VARIABLE		dt
		FindLevel /EDGE=1 /P/Q/R=[spikepoint-1,spikepoint-pointsBefore] voltageWav,V_detect
		IF (V_Flag!=0)	// NO level crossing found
			   // => the wanted crossing lies even before the previous spike
			   // so put a NaN in the spike time
		   return NaN
		ELSE
			return pnt2x(voltageWav,floor(V_LevelX)) + dt*(V_LevelX-floor(V_LevelX))
		ENDIF   

END

FUNCTION SpikesAtVdetectMultiThread(InWave_V, InWave_ST,V_detect[,Suffix])
WAVE	InWave_V, InWave_ST		        // voltage and spike-time wave
VARIABLE	V_detect					// voltage threshold value to define modified spike time
										// in Volts, will be rounded to full millivolts
STRING Suffix							// if given, the new spike time wave carries this suffic after the name of InWave_ST

//	Variable timerRefNum, microsecs
//	timerRefNum = StartMSTimer
//	if (timerRefNum == -1)
//		Abort "All timers are in use"
//	endif
// in contrast to the standard SpikesAtVdetect, this procedure does not upsample the voltage waveform
// this results in differences on the order of half a sample period
// it is 600 times faster in a typical use case


//V_detect = round(1000*V_detect)/1000 	// rounded to full millivolts

// travels back in time from the original spike time, 
// to find the latest time voltage crosses V_detect before the full blown spike

// complications: the time should be later than the previous full blown spike
// if the t_sp_detect[n] does fall before the previous full blown spike
// i.e. if t_sp_detect[n]<t_sp[n+1] that will be noted by putting a NaN in the result wave

// the function returns a wave with one entry for each spike in InWave_ST
// the result's name has a suffix derived from the threhsold. 
// therefore the threshold voltage is rounded to full millivolts

		VARIABLE		dt=DimDelta(InWave_V,0)
		
// health test: Spike times should define a narrow set of voltages in the voltage trace
// (assuming that a VOLTAGE threshold rather than a rate threshold was used as definition
// here I use that test to see whether the provided voltage wave and spike time wave are plausibly a valid pair
		Duplicate /O InWave_ST,dummy, pointsBefore, spikepoints
		WAVE 	dummy
		dummy[]=InWave_V[x2pnt(InWave_V,InWave_ST[p])]
		WAVESTATS/Q dummy
		VARIABLE	expectedVjitter=1500*dt/2*2 // dV/dt max * dt/2  * safety factor of 1.5
		IF ( (V_max-V_min) > expectedVjitter)
			DoAlert 0,"The range of threshold voltages related to the spike times\ris "+num2str((V_max-V_min)*1000)+" mV.\rThat seems too big, aborting..."
			KillWaves dummy
			Return -1
		
		ENDIF
		
		// create name for output
		STRING		SPdetectName
		IF (ParamIsDefault(Suffix))
			SPdetectName=NameOfWave(InWave_ST)+"_at_"+ReplaceString("-",  num2istr(V_detect*1000), "neg")
			SPdetectName=ReplaceString(".", SPdetectName, "p")
		ELSE
			SPdetectName=NameOfWave(InWave_ST)+Suffix
		ENDIF

		Duplicate /O InWave_ST, $SPdetectName
		WAVE	STdet=$SPdetectName
		SetScale d 0,0,"s", STdet
		SetScale/P x,0,1,"", STdet
		// has correct unit (seconds) and index scaling in x
		
		// now get a stretch of data before the original spike time
		// possibly upsample
		
		VARIABLE		nSpikes=DimSize(InWave_ST,0), sp// Pnt number in Voltage wave where spike was detected
		VARIABLE	 	upsamplefac=round(100000*dt)// stretch of data before and after that is copied for analysis
																			// factor of how much the temporal sampling is increased
		VARIABLE		nPntsInDummy
		
		pointsBefore[1,nSpikes-1]=(InWave_ST[p]-InWave_ST[p-1])/dt
		pointsBefore[0]=InWave_ST[0]/dt -1
		
		pointsBefore*= upsamplefac
		
		spikepoints[]=x2pnt(InWave_V,InWave_ST[p])
		
		Multithread STdet[]=WalkBack(InWave_V, V_detect, spikepoints[p], pointsbefore[p], dt)
//	   	microsecs = StopMSTimer(timerRefNum)
//	Print microsecs/DimSize(InWave_ST,0), "microseconds per spike"

END // SpikesAtVdetect multithread

PROC WrapperForGainDecay_fixedVdetect(Abs_V_detect)
VARIABLE Abs_V_detect in mV

// done in a single folder at a time
// this is for the case that the gain is calculated based on spike triggered average input using spike times 
// obtained by using a fixed voltage to define the spike time
// the existence of the spikes was previously detected by using another criterium, e.g. crossing a rate threshold or a high voltage
// now the spike time is reestablished using a lower voltage

// this wrapper requires that AC_avg_scaled_splt_FFT and STA_avg_scaled already exist
		
		// Go through all voltage traces and the corresponding spike time traces and obtains new spike times
		XeqtInSubs("Xeqt4WList(\"*_V\",\"SpikesAtVdetect(~, §~rmvend_§_ST,Abs_V_detect)\")")
		// those are collected in waves ending in _ST_at_XYYY where YYY is the detection voltage in millivolts.
		// If it is negative, X is m for "minus", otherwise it is ""
		// this suffix is created here: 
		STRING Suff="_"+ReplaceString("-",  num2istr(Abs_V_detect*1000), "m")
		STRING STASuff="STA"+Suff
		STRING STSuff="ST_at"+Suff
		// Now recreate the spike-triggered average input
		// 
		STRING	CMDSTR=ReplaceString("??STSuff??","Xeqt4WList(\"*_I\",\"STAfromAnalogue( ~, §~rmvend_§_??STSuff?? , 1,0, Suffix=\\\"??STASuff??\\\")\")",STSuff)
		CMDSTR=ReplaceString("??STASuff??",CMDSTR,STASuff)
		Execute(CMDSTR)
		
		// create the average spike triggered average for a single folder
		STRING STAListe=WaveList("*_"+STASuff,";","")
		AvgSTAFromList(STAListe,STASuff+"_avg_scaled") 
		// using the Execute formalism to avoid the need to create a wavename etc. 
		Execute("Wavestats/Q "+STASuff+"_avg_scaled; SplitBeforeFFT("+STASuff+"_avg_scaled,V_maxloc)")
		// this create the split version with a name ending in "_splt"
		Execute("FFT/OUT=1/DEST="+STASuff+"_avg_scaled"+"_splt_FFT "+STASuff+"_avg_scaled"+"_splt")
		
		// calculate gain and filter in complex space
		Execute("Duplicate/O "+STASuff+"_avg_scaled"+"_splt"+"_FFT Gain_avg_scaled"+Suff)
		Execute("Gain_avg_scaled"+Suff+"/=AC_avg_scaled_splt_FFT;Gain_avg_scaled"+Suff+"*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0)")
		Execute("Gain_avg_scaled"+Suff+"= conj(Gain_avg_scaled"+Suff+")")
		Execute("GaussFilter(Gain_avg_scaled"+Suff+")")

END

PROC WrapperForGainDecay_Vdetect_re_Trial_VOnset(rel_trial_V_detect, rateThresh)
VARIABLE rel_trial_V_detect // in mV
VARIABLE rateThresh 		// in V/s

// done in a single folder at a time
// detection threshold varies trial by trial and is defined relative to the average onsetVoltage, 
// i.e. the average voltage at which the rate threshold  (rateThresh) is reached  (average within the given trial)

// the variable "rel_trial_V_detect" is interpreted relative to the trial-average onset voltage
// - 0.002 means 2 millivolts below the trial-average onset voltage 

// the existence of the spikes was previously detected by using another criterium, e.g. crossing a rate threshold or a high voltage
// now the spike time is reestablished using a lower voltage and with those modified spike times, the gain is calculated by spike triggerd average

// this wrapper requires that AC_avg_scaled_splt_FFT and STA_avg_scaled already exist
		
		// redo spike detection with criterium of crossing 0 and an onset detection as crossing rateThresh
		STRING Prefix="OU"
		STRING CMDSTR
		STRING WaveSuffix= "_at_avgVthr_"+ReplaceString(".",ReplaceString("-",  num2istr(rel_trial_V_detect*1000), "neg"),"p")
		
		// get spike times and thereby recompute also the onset times (times at which rate threshold is crossed) for the intended rateThresh
//		CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_V\",\" ReturnSpikeTimes(~,V_Threshold=0, MinISI=0.002, RateThresh4Onset = "+num2str(rateThresh)+")\")",Prefix)
//		Execute(CMDSTR)
			
		
		
		// Go through all voltage traces and the corresponding spike time traces and obtains new spike times
		// based on the average voltage at onset
		
		Xeqt4WList("*_V","VARIABLE/G avgVTHR=mean(§~rmvend_§_VTHR);SpikesAtVdetect(~, §~rmvend_§_ST,avgVTHR+rel_trial_V_detect,Suffix= WaveSuffix)")
		// those are collected in waves ending in _ST_at_XYYY where YYY is the detection voltage in millivolts.
		// If it is negative, X is m for "minus", otherwise it is ""
		// this suffix is created here: 

		STRING STASuff="STA"+WaveSuffix
		STRING STSuff="ST"+WaveSuffix
		// Now recreate the spike-triggered average input
		// 
		CMDSTR=ReplaceString("??STSuff??","Xeqt4WList(\"*_I\",\"STAfromAnalogue( ~, §~rmvend_§_??STSuff?? , 1,0, Suffix=\\\"??STASuff??\\\")\")",STSuff)
		CMDSTR=ReplaceString("??STASuff??",CMDSTR,STASuff)
		Execute(CMDSTR)
		
		// create the average spike triggered average for a single folder
		STRING STAListe=WaveList("*_"+STASuff,";","")
		AvgSTAFromList(STAListe,STASuff+"_avg_scaled") 
		// using the Execute formalism to avoid the need to create a wavename etc. 
		Execute("Wavestats/Q "+STASuff+"_avg_scaled; SplitBeforeFFT("+STASuff+"_avg_scaled,V_maxloc)")
		// this create the split version with a name ending in "_splt"
		Execute("FFT/OUT=1/DEST="+STASuff+"_avg_scaled"+"_splt_FFT "+STASuff+"_avg_scaled"+"_splt")
		
		// calculate gain and filter in complex space
		Execute("Duplicate/O "+STASuff+"_avg_scaled"+"_splt"+"_FFT Gain"+WaveSuffix)
		Execute("Gain"+WaveSuffix+"/=AC_avg_scaled_splt_FFT;Gain"+WaveSuffix+"*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0)")
		Execute("Gain"+WaveSuffix+"= conj(Gain"+WaveSuffix+")")
		Execute("GaussFilter(Gain"+WaveSuffix+")")

END

PROC GainDecayAcrossSubFolders(STA_Suffix)
STRING STA_Suffix // across all subfolders, STA waves with this suffix are collected and a joint STA is calcualted and from there the gain

	STRING Prefix="OU" // the prefix of all those STAs of individual trials

		
	// collect all spike triggered averages from all subfolders 
	STRING/G FolderSTAs=""
	STRING CMDSTR="Xeqt4WList(\"??Pref??*_"+STA_Suffix+"\",\"root:FolderSTAs+=\\\"§SUBFULL§~;\\\"\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	CMDSTR=ReplaceString("root:",CMDSTR,GetDataFolder(1))
	XeqtInSubs(CMDSTR)
	
	// create the AVG STA from across all folders. Result is in superordinate folder
	AvgSTAFromList(FolderSTAs,STA_Suffix)
	
	// calculate overall gain for data from all subfolders
	SplitBeforeFFT(AC_avg_scaled,0)
	FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt

	STRING STAName=STA_Suffix
	WAVESTATS/Q $STAName
	SplitBeforeFFT($STAName,V_maxloc)
	STAName=STA_Suffix+"_splt"
	STRING GainName="Gain"+STA_Suffix[3,strlen(STA_Suffix)-1]
	CMDSTR="FFT/OUT=1/DEST="+GainName+" "+STAName
	
	Execute(CMDSTR)
	STAName=STAName+"_splt_FFT"

	$GainName/=AC_avg_scaled_splt_FFT
	$GainName*=cmplx(str2num(StringByKey("total#spikes", note(STA_avg_scaled),":" ,"\r"))/str2num(StringByKey("totalduration", note(STA_avg_scaled),":" ,"\r")),0)
	$GainName= conj($GainName)	// fixes the sign of the phase
	GaussFilter($GainName)
	GainName=GainName+"_MgFlt"
	note $GainName,"Zero delay Gain\r"+STA_Suffix+"\r"+note(STA_avg_scaled)
	
END
