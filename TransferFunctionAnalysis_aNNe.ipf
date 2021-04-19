#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>
#include <Waves Average>
#include <Remove Points>
#include ":IFDL v4 Procedures:IFDL"
// Version Header
// 
// 20/10/2016	_  aNNe
// 	Changed STAfromAnalogue to take suffix (i.e. typically "STA") as optional parameter
// 	use has not changed for standard analysis, now it can be used from the new bootstrap routines
//
// 31/10/2016	_  aNNe
// 	Interpolation in Spike Statistics is adjusted according to sample rate. The interpolated rate
// 	aimed for is 100 kHz.
//
// 15/11/2016	_	aNNe
//		Code for the final version of bootstrap of confidence intervals is introduced into the main file
//		It relies on the same Functions as the simple analysis but performs it on the files with ending "STB"
//		Can be called with the command BST_Wrapper(#repetitions)
//
//		Removed spurious code from earlier realisations of the analysis (AvgAC, some old bootstrap version)
//		
//		Adjusted Code to be suitable for Omers data (where beginning and end of current waves is rubbish)
// 	For this call MakeAndAvgAc(FirstPoint=14532,LastPoint=915469)
//
//	23/08/2017	_	aNNe
//		Added NoiseFloor calculation to BST_Wrapper
//		BST_Wrapper(200,0)	// zero stands for NO noiseFloor
//		BST_Wrapper(200,1, earliestST=0.5, latestST=49.5)		// this is doing the noise floor and earliest and latest conceivable Spike times are provided
// 	that is not absolutely necessary but would help a bit
// 	if the STA is 1 s wide, the earliest possible ST is typically 0.5 and the latest duration - 0.5 seconds
//
//	14/05/2018 _ aNNe
//		Added Visualization to show vectorstrength
//
// 3/07/2018 _ aNNe
//  		Added in function ReturnSpikeTimes: Line: 822-831 :  quick health check: dV/dt should not be larger than say 2000 V/s
//
//	15/08/2018 _ aNNe
//			now all spike statistic waves (full width @half max, V threshold etc.) share the wave note with the spike time wave
//			placed functions for vector strength analysis and data statification here
// 23/08/2018
//			additional stratification functions included
//	
//	01/10/2018 _ aNNe
//			tools for stratification along theta oscillations included: AssignThetaPhaseToST(SpikeTimeWave)
//																							StratByThetaPhase(n_strat, start_phase, identifier,SineFreq doDisplay)
//			those should later be moved to a dedicated procedure file since they are not for TransferFunction analysis per se
//
//	11/03/2019 _ aNNe
//			Multithreading for Bootstrap
//
// 14/03/2019 _ aNNe 
//			Multitreading for AssignThetahase, now use "AssignThetaPhaseToAllSTInFolder()" to have all Theta Phase waves calcuated in parallel
// 
// 24/4/2019 _ aNNe 
//			Added CollectSpikePropertiesByFolder() to rapidly collect all spike properties
// 			it reads the wavenote of the spike time waves, those can be updated by calling spikestats
// 			from those it gets a value for each trial which are then averaged across all trials of the folder
//			for each feature XYZ, it creates a wave "ByFolder_XYZ" in the parent folder
// 
// 20/09/2019 _ aNNe
//			Bootstrap Confidence interval correction based on Bradley Efron & Trevor Hastie 2016 "Computer Age statistical inference" pg 192ff
//			going back to Efron 1982 
//
// 30/05/2020 _ aNNe 
//			Vector Strength confidence intervals fixed. The _complex_ value of the vector strength of the bootstrap 
//			sample is now taken into account, not only the magnitude. A 95% confidence _region_ in phase / magnitude
//		 	space is determined and the distance of this region from zero magnitude is taken as lower bound for the 
//			vector strength magnitude
//
// 03/08/2020 _ aNNe 
//			Fixed the sign of the Phase of the gain as calculated in GainCalculationWrapper and in BST_Wrapper
//			This required replacing gain with its conjugated complex.
//
// 01/04/2021 _ aNNe 
//			When using stratiying analysis, selecting spikes by one or two criteria (GainCalculationWrapperForSpikeSubset)
//			There is a new option to keep those ST waves (with suffix ST_wC) to run Bootstrap. The Bootstrap Wrapper how accepts 
//			a suffix to signal on which dataset botstrap should run (full or stratified)


MACRO GainCalculationWrapper()
// start in the folder above all cell folders
// goes through and creates the avg gain 

	STRING Prefix="OU"
	STRING CMDSTR
	STRING PathToFolder=GetDataFolder(1)
	
	CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_V\",\" ReturnSpikeTimes(~,V_Threshold=0, MinISI=0.005) \")",Prefix)
	XeqtInSubs(CMDSTR)
	CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G W_name= \\\"~\\\"; W_name=W_name[0,strlen(W_name)-2]+\\\"ST\\\"; STAfromAnalogue( ~, $W_Name, 1,0)\")",Prefix)
	XeqtInSubs(CMDSTR)
	
	// create autocorrelation traces for each trial in each subfolder
	XeqtInSubs("MakeAC()")
	
	// Collect all ACs from all subfolders
	STRING/G FolderACs= ""
	CMDSTR="Xeqt4WList(\"??Pref??*_AC\" ,\"root:FolderACs+= \\\"§SUBFULL§~;\\\"\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	CMDSTR=ReplaceString("root:",CMDSTR,GetDataFolder(1))
	XeqtInSubs(CMDSTR)
	
	// create the AVG AC from across all folders. Result is in superordinate folder
	AvgACFromList(FolderACs,"AC_avg_scaled")
	
	// collect all spike triggered averages from all subfolders 
	STRING/G FolderSTAs=""
	CMDSTR="Xeqt4WList(\"??Pref??*_STA\",\"root:FolderSTAs+=\\\"§SUBFULL§~;\\\"\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	CMDSTR=ReplaceString("root:",CMDSTR,GetDataFolder(1))
	XeqtInSubs(CMDSTR)
	
	// create the AVG STA from across all folders. Result is in superordinate folder
	AvgSTAFromList(FolderSTAs,"STA_avg_scaled")
	
	// calculate overall gain for data from all subfolders
	SplitBeforeFFT(AC_avg_scaled,0)
	WAVESTATS/Q STA_avg_scaled
	SplitBeforeFFT(STA_avg_scaled,V_maxloc)
	FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt
	FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt
	Duplicate /O STA_avg_scaled_splt_FFT, Gain_avg_scaled;
	Gain_avg_scaled/=AC_avg_scaled_splt_FFT
	Gain_avg_scaled*=cmplx(str2num(StringByKey("total#spikes", note(STA_avg_scaled),":" ,"\r"))/str2num(StringByKey("totalduration", note(STA_avg_scaled),":" ,"\r")),0)
	Gain_avg_scaled= conj(Gain_avg_scaled)	// fixes the sign of the phase
	GaussFilter(Gain_avg_scaled)
	note Gain_avg_scaled_MgFlt,note(STA_avg_scaled)
	
//	// now create avg AC and STA and gain in each subfolder
//	
//	CMDSTR="STRING/G FolderACs=\"\"; Xeqt4WList(\"??Pref??*_AC\",\"FolderACs+=\\\"~;\\\"\") ; AvgACFromList(FolderACs,\"AC_avg_scaled\")"
//	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
//	XeqtInSubs(CMDSTR)
//	CMDSTR="STRING/G FolderSTAs=\"\"; Xeqt4WList(\"??Pref??*_STA\",\"FolderSTAs+=\\\"~;\\\"\") ; AvgSTAFromList(FolderSTAs,\"STA_avg_scaled\")"
//	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
//	XeqtInSubs(CMDSTR)
//	
//	XeqtInSubs("SplitBeforeFFT(AC_avg_scaled,0);WAVESTATS/Q STA_avg_scaled;SplitBeforeFFT(STA_avg_scaled,V_maxloc);")
//	XeqtInSubs("FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt;FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt;Duplicate /O STA_avg_scaled_splt_FFT, Gain_avg_scaled;")
//	XeqtInSubs("Gain_avg_scaled/=AC_avg_scaled_splt_FFT; Gain_avg_scaled*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0) ")
//	XeqtInSubs("GaussFilter(Gain_avg_scaled);	note Gain_avg_scaled_MgFlt,note(STA_avg_scaled)")
END 
	
MACRO GainCalculationForSpikeSubset(Crit0_suff,Crit1_suff,Crit0_min, Crit0_max, Crit1_min, Crit1_max)
STRING 		Crit0_suff,Crit1_suff
VARIABLE		Crit0_min, Crit0_max, Crit1_min, Crit1_max
// selects spikes that fullfill a certain criterium
// i.e. a certain theta phase (and amplitude)
// start in the folder above all cell folders
// goes through and creates the avg gain 

// the criteria are in waves to be identified based on the base name of the spike time waves
// by replacing the suffix "ST" with the S_Crit0_suff string
// if the S_Crit1_suff string is empty, the second critrion is ignored


// use cases

// GainCalculationForSpikeSubset("ThetaPhs","ThetaAmpNorm",0, Pi/4, 0.5, 100)
// GainCalculationForSpikeSubset("ISI","",0.02, 2,0, 0)

	VARIABLE keepSpikeTimesForBootstrap=0
	
	
	STRING Prefix="OU"
	STRING CMDSTR
	STRING PathToFolder=GetDataFolder(1)
	
//	CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_V\",\" ReturnSpikeTimes(~,V_Threshold=0, MinISI=0.005) \")",Prefix)
//	XeqtInSubs(CMDSTR)

	// for every single Spike time wave, use the criteria to create a permanent wave containing a subset of the spikes and use this as the argument in the STAfromAnalogue call
	// replacing the full Spike time wave
// here edit criteria
	VARIABLE/G root:V_crit0_min=Crit0_min
	VARIABLE/G root:V_crit0_max=Crit0_max

	STRING/G		root:S_crit0_Suff = Crit0_suff
	
	// the note that will be applied to the resulting STA and gain
	STRING NoteAppendix
	NoteAppendix=Crit0_suff+"_min:="+num2str(Crit0_min)+"\r"
	NoteAppendix+=Crit0_suff+"_max:="+num2str(Crit0_max)+"\r"

	IF (strlen(Crit1_suff)>0)	// apply a second criteriun
		VARIABLE/G root:V_crit1_min=Crit1_min
		VARIABLE/G root:V_crit1_max=Crit1_max
	
		STRING/G		root:S_crit1_Suff = Crit1_suff
		STRING/G		root:S_crit0_Suff = Crit0_suff		
		IF (!keepSpikeTimesForBootstrap)
			// use spike times with criterium not as explicit wave but only as free waves handed over 
			CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G S_Source= \\\"~\\\"; String/G S_Crit0=S_Source[0,strlen(S_Source)-2]+root:S_crit0_Suff; String/G S_Crit1=S_Source[0,strlen(S_Source)-2]+root:S_crit1_Suff ;S_Source=S_Source[0,strlen(S_Source)-2]+\\\"ST\\\"; STAfromAnalogue( ~,	WaveSubsetByCriteria($S_Source, $S_Crit0,root:V_crit0_min,root:V_crit0_max, CriteriumWave1= $S_Crit1,lowCrit1=root:V_crit1_min, upCrit1=root:V_crit1_max,Logic=\\\"and\\\"), 1,0, Suffix=\\\"STAwC\\\")\")",Prefix)
		ELSE
			// first create STwC waves
			CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_ST\",\"STRING/G S_Source= \\\"~\\\"; String/G S_Crit0=S_Source[0,strlen(S_Source)-3]+root:S_crit0_Suff; String/G S_Crit1=S_Source[0,strlen(S_Source)-3]+root:S_crit1_Suff ; Duplicate/O 	WaveSubsetByCriteria(~, $S_Crit0,root:V_crit0_min,root:V_crit0_max, CriteriumWave1= $S_Crit1,lowCrit1=root:V_crit1_min, upCrit1=root:V_crit1_max,Logic=\\\"and\\\") , $(\\\"~\\\"+\\\"wC\\\")\")",Prefix)
			CMDSTR+=";"+ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G S_Source= \\\"~\\\";S_Source=S_Source[0,strlen(S_Source)-2]+\\\"STwC\\\"; STAfromAnalogue( ~,	$S_Source, 1,0, Suffix=\\\"STAwC\\\")\")",Prefix)
		ENDIF
		NoteAppendix+=Crit1_suff+"_min:="+num2str(Crit1_min)+"\r"
		NoteAppendix+=Crit1_suff+"_max:="+num2str(Crit1_max)+"\r"
	ELSE		// apply only one criterium
				// FOR instance for selecting by ISI 

		IF (!keepSpikeTimesForBootstrap)
			CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G S_Source= \\\"~\\\"; String/G S_Crit0=S_Source[0,strlen(S_Source)-2]+root:S_crit0_Suff; S_Source=S_Source[0,strlen(S_Source)-2]+\\\"ST\\\"; STAfromAnalogue( ~,	WaveSubsetByCriteria($S_Source, $S_Crit0,root:V_crit0_min,root:V_crit0_max), 1,0, Suffix=\\\"STAwC\\\")\")",Prefix)

		ELSE
		
			CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_ST\",\"STRING/G S_Source= \\\"~\\\"; String/G S_Crit0=S_Source[0,strlen(S_Source)-3]+root:S_crit0_Suff; Duplicate/O WaveSubsetByCriteria(~, $S_Crit0,root:V_crit0_min,root:V_crit0_max),$(\\\"~\\\"+\\\"wC\\\")\")",Prefix)
			CMDSTR+=";"+ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G S_Source= \\\"~\\\"; String/G S_Crit0=S_Source[0,strlen(S_Source)-2]+root:S_crit0_Suff; S_Source=S_Source[0,strlen(S_Source)-2]+\\\"STwC\\\"; STAfromAnalogue( ~,	$S_Source, 1,0, Suffix=\\\"STAwC\\\")\")",Prefix)
		ENDIF
	ENDIF
		
XeqtInSubs(CMDSTR)
	
//	// create autocorrelation traces for each trial in each subfolder
//	XeqtInSubs("MakeAC()")
//	
//	// Collect all ACs from all subfolders
//	STRING/G FolderACs= ""
//	CMDSTR="Xeqt4WList(\"??Pref??*_AC\" ,\"root:FolderACs+= \\\"§SUBFULL§~;\\\"\")"
//	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
//	CMDSTR=ReplaceString("root:",CMDSTR,GetDataFolder(1))
//	XeqtInSubs(CMDSTR)
//	
//	// create the AVG AC from across all folders. Result is in superordinate folder
//	AvgACFromList(FolderACs,"AC_avg_scaled")
//	
	// collect all spike triggered averages from all subfolders 
	STRING/G FolderSTAs=""
	CMDSTR="Xeqt4WList(\"??Pref??*_STAwC\",\"root:FolderSTAs+=\\\"§SUBFULL§~;\\\"\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	CMDSTR=ReplaceString("root:",CMDSTR,GetDataFolder(1))
	XeqtInSubs(CMDSTR)
	
//	// create the AVG STA from across all folders. Result is in superordinate folder
	AvgSTAFromList(FolderSTAs,"STA_avg_scaled_wC")
	
	// calculate overall gain for data from all subfolders
	SplitBeforeFFT(AC_avg_scaled,0)
	WAVESTATS/Q STA_avg_scaled_wC
	SplitBeforeFFT(STA_avg_scaled_wC,V_maxloc)
	FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt
	FFT/OUT=1/DEST=STA_avg_scaled_wC_splt_FFT STA_avg_scaled_wC_splt

	Duplicate /O STA_avg_scaled_wC_splt_FFT, Gain_avg_scaled_wC;

	Gain_avg_scaled_wC/=AC_avg_scaled_splt_FFT
	Gain_avg_scaled_wC*=cmplx(str2num(StringByKey("total#spikes", note(STA_avg_scaled_wC),":" ,"\r"))/str2num(StringByKey("totalduration", note(STA_avg_scaled_wC),":" ,"\r")),0)
	Gain_avg_scaled_wC= conj(Gain_avg_scaled_wC)
	GaussFilter(Gain_avg_scaled_wC)

	Note/K Gain_avg_scaled_wC_MgFlt
	note Gain_avg_scaled_wC_MgFlt,note(STA_avg_scaled_wC) + NoteAppendix
	note/K STA_avg_scaled_wC,note(Gain_avg_scaled_wC_MgFlt)
	
//	// now create avg AC and STA and gain in each subfolder
//	
//	CMDSTR="STRING/G FolderACs=\"\"; Xeqt4WList(\"??Pref??*_AC\",\"FolderACs+=\\\"~;\\\"\") ; AvgACFromList(FolderACs,\"AC_avg_scaled\")"
//	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
//	XeqtInSubs(CMDSTR)
//	CMDSTR="STRING/G FolderSTAs=\"\"; Xeqt4WList(\"??Pref??*_STA\",\"FolderSTAs+=\\\"~;\\\"\") ; AvgSTAFromList(FolderSTAs,\"STA_avg_scaled\")"
//	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
//	XeqtInSubs(CMDSTR)
//	
//	XeqtInSubs("SplitBeforeFFT(AC_avg_scaled,0);WAVESTATS/Q STA_avg_scaled;SplitBeforeFFT(STA_avg_scaled,V_maxloc);")
//	XeqtInSubs("FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt;FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt;Duplicate /O STA_avg_scaled_splt_FFT, Gain_avg_scaled_wC;")
//	XeqtInSubs("Gain_avg_scaled_wC/=AC_avg_scaled_splt_FFT; Gain_avg_scaled_wC*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0) ")
//	XeqtInSubs("GaussFilter(Gain_avg_scaled_wC);	note Gain_avg_scaled_wC_MgFlt,note(STA_avg_scaled)")
END 
	

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
			VARIABLE ss // current stratus
			FOR (ss=0; ss<n_strat; ss+=1)
				// check whether status contains boundary of the cyclic range ( i.e. Pi)
				IF  ( (binBounds[ss+1] <= 0) && (binBounds[ss] > 0) )
				// need to join two datasets
				WAVE w1=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=-Pi,upCrit=binBounds[ss+1])
				WAVE w2=CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=Pi)
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
					BootstrapVS(CollectDataByNoteAndCriterium("SPhs",NoteKey="Sine frequency (Hz)",low4NoteVal=SineFreq, up4NoteVal=SineFreq,CriteriumSuffix="ThetaPhs", lowCrit=binBounds[ss],upCrit=binBounds[ss+1]), QuantileWave=TrgtQuant, quiet=1 )
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




FUNCTION PhaseLocking(V_Wave, DoDisplay)
	WAVE			V_Wave
	Variable 	DoDisplay		// if >1 it gives a graph
	
	
	STRING 	IName, STName, VSName, ISIName, VoltName=NameOfWave(V_Wave)
	VARIABLE	VoltNameLength=strlen(VoltName)
	
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																// in VoltName are "_V"
		VSName=VoltName[0,VoltNameLength-3]+"_VS"														
		STName=VoltName[0,VoltNameLength-3]+"_ST"														
		WAVE ST=$STName
		IF (!(Exists(STName)==1))
			// spike time wave does not yet exist
			Abort "First get Spike time wave, then call phaselock visualization"
		ENDIF
		IName=VoltName[0,VoltNameLength-3]+"_I"														
		ISIName=VoltName[0,VoltNameLength-3]+"_ISI"														
		WAVE I_Wave=$IName
		WAVE ISI_Wave=$ISIName
		IF (!(Exists(IName)==1))
			// spike time wave does not yet exist
			Abort "Cannot find input for "+VoltName
		ENDIF
		IF (!(Exists(ISIName)==1))
			// spike time wave does not yet exist
			Abort "Cannot find ISI for "+VoltName
		ENDIF
				
	ELSE
		Abort "Cannot decode names as voltage trace does not end in 'V'"
	ENDIF
	
	Duplicate/O ST, $VSName
	WAVE	VecStr=$VSName
	Redimension/N=(-1,2) VecStr
	VecStr=0
	
	
	
	// get sine frequency that predominates in injected current --> max of FFT
	VARIABLE	L_Input=DimSize(I_wave,0)
	FFT/RP=[0,2*floor(L_Input/2)-1]/OUT=3/DEST=Temp_FFT  I_wave
	// input needs to have a correct x-scaling
	
	WAVESTATS/Q/M=1/R=[1,DimSize(Temp_FFT,0)] Temp_FFT
	VARIABLE	freq=round(V_maxLoc)
	
	KillWaves Temp_FFT
	
	VARIABLE		SpikeTime,RR, nSpikes=DimSize(ST,0)
	FOR (RR=0;RR<nSpikes; RR+=1)
		VecStr[RR][0]= sin(2*Pi*Freq*ST[RR])
		VecStr[RR][1]= cos(2*Pi*Freq*ST[RR])
		
	ENDFOR 
		
	// Visualization
	IF (DoDisplay)
		Display/L=VertCrossing/B=HorizCrossing VecStr[*][0] vs VecStr[*][1]
		ModifyGraph mode=3,marker=19,msize=2
		ModifyGraph height={Aspect,1}
		SetAxis VertCrossing -1,1
		SetAxis HorizCrossing -1,1
		TextBox/C/N=text0/F=0/A=MC "f\\Bsine\\M="+num2istr(freq)
		TextBox/C/N=text0/A=LT/X=-3.00/Y=-3.00
		ModifyGraph zColor($VSName)={ISI_Wave,*,*,YellowHot,1}
		ModifyGraph logZColor=1
		ModifyGraph useMrkStrokeRGB=1
		ModifyGraph mrkStrokeRGB=(43690,43690,43690)
		ModifyGraph freePos(VertCrossing)={0,HorizCrossing},freePos(HorizCrossing)={0,VertCrossing}
	ENDIF

	
END

FUNCTION/WAVE StratifyVS(VSWave,ISIWAVE, low,high)
WAVE		VSWave // contains x and y component of the vector that each spike represents in the Gauss-Plane
WAVE		ISIWave // contains ISIs
VARIABLE	low, high // limits for the ISI range we want to analyze // given in seconds
	
	DUPLICATE/O ISIWave, TempWeight
	Duplicate/O VSWave, TempVS
	
	TempWeight[]= ((ISIWave[p]>low) && (ISIWave[p]<high)) ? 1 : 0
	
	TempVS[][]*=TempWeight[p]	
	
	VARIABLE	nInclude=Sum(TempWeight)
	
	
	printf "%g spikes have ISI between %g and %g\r" ,nInclude, low, high
	
	MatrixOP/O TempDest=sumCols(TempVS)
	
	TempDest/=nInclude
	MAKE	Result= {sqrt((TempDest[0][0])^2+(TempDest[0][1])^2), atan2(TempDest[0][1],TempDest[0][0])}
	
	Killwaves  TempWeight, TempVS, TempDest
	
	Return Result
	
END

 FUNCTION/WAVE waveStratifyData(TargetData,CriteriumData, lowCriterium, highCriterium[,CriteriumData1, lowCriterium1, highCriterium1, logic])
	WAVE		TargetData 		// 1D wave with those data that will be returned after cases are filtered out
	WAVE		CriteriumData,CriteriumData1 	
									// 1D wave(s) that contain(s) the filtering criterium
									// the function returns the reference to a wave that contains only those entries of "TargetData"
									// for which the corresponding entries of "CriteriumData" are numerically above lowCriterium and below highCriterium
	VARIABLE	lowCriterium, highCriterium ,lowCriterium1, highCriterium1 
									// limits for the range we want to analyze 
									
	STRING	logic				// is "AND" or "OR"

	DFREF dfSav= GetDataFolderDFR()

	// Create a free data folder and set it as the current data folder
	SetDataFolder NewFreeDataFolder()

	Duplicate/O CriteriumData TempWeight
	VARIABLE	rr=0, nEntries=DimSize(TargetData,0), includeCount=0
		
	TempWeight[]= ((CriteriumData[p]>lowCriterium) && (CriteriumData[p]<highCriterium)) ? 1 : 0
	
	IF (!ParamIsDefault(CriteriumData1)) // a second criterium is handed over
	
		strswitch(logic)	// string switch
			case "AND":	// execute if case matches expression
				TempWeight[]= ((CriteriumData1[p]>lowCriterium1) && (CriteriumData1[p]<highCriterium1)) ? TempWeight[p] : 0
				break		// exit from switch
			case "OR":	// execute if case matches expression
				TempWeight[]= ((CriteriumData1[p]>lowCriterium1) && (CriteriumData1[p]<highCriterium1)) ? 1 : TempWeight[p]
				break
			default:			// optional default expression executed
				// DoAlert 0,"Logic has to be \"AND\" or \"OR\".\rThe second criterium is disregarded"
		endswitch
	ENDIF
	
	VARIABLE	nInclude=Sum(TempWeight)
	STRING	TargetName=NameOfWave(TargetData)+"_Strat"
	
	Make/O/N=(nInclude) $TargetName
	WAVE Target=$TargetName
		
	DO
		IF (TempWeight[rr])
				Target[inCludeCount]=TargetData[rr]
				includeCount+=1
		ENDIF	 
		
		rr+=1
	WHILE ( (rr<nEntries) && ( includeCount < nInclude) )
	
	Killwaves  TempWeight
	SetDataFolder dfSav

	
	RETURN Target
	

END	// waveStratifyData

FUNCTION StratifyData(TargetData,CriteriumData, lowCriterium, highCriterium, OutName)
	WAVE		TargetData 		// 1D wave with those data that will be returned after cases are filtered out
	WAVE		CriteriumData 	// 1D wave that contains the filtering criterium
									// the function returns the reference to a wave that contains only those entries of "TargetData"
									// for which the corresponding entries of "CriteriumData" are numerically above lowCriterium and below highCriterium
	VARIABLE	lowCriterium, highCriterium // limits for the ISI range we want to analyze // given in seconds
	STRING	OutName
	
	WAVE result=waveStratifyData(TargetData,CriteriumData, lowCriterium, highCriterium)
	Duplicate/O result,$OutName
	
END


FUNCTION	PairwiseDistance(XYLine,XYLineStack)
WAVE		XYLine,XYLineStack	// 

END


FUNCTION DistributeByDensity(Field, W_Pos)
WAVE		Field // field (x) gives a probability density (not normalized)
					// to find any of the // elements at the location x
WAVE		W_Pos	// handed over just to be filled by the procedure (and to signal the number of elements to be placed

VARIABLE	nParticle=DimSize(W_Pos,0)	


//get the positions at which the field is given
Duplicate/O field, W_xpos
WAVE W_xpos
W_xpos[]=pnt2x(Field,p)

// create the cumulative DF
Integrate Field /D=W_Integral
VARIABLE	Pmin=W_Integral[0], Pmax=W_Integral[DimSize(W_Integral,0)-1]
VARIABLE	Prange=Pmax-Pmin
// randomly draw nParticle numbers between the First and last value of the integral
// put these in ProbPos
Duplicate/O W_Pos, W_ProbPos
WAVE W_ProbPos

//ProbPos=enoise(Prange/2,2)+Prange/2+Pmin
W_ProbPos[]=Pmin+Prange/(NParticle-1)*p
// now use probpos as xvalue, find the xpos at which Integral crosses probpos
W_Pos[]=interp(W_ProbPos[p], W_Integral, W_xpos )


KillWaves W_xpos, W_Integral, W_ProbPos

END

FUNCTION EquiSample(XYSequence,nPoints, NNei)
WAVE	XYSequence 	// contains a LINEAR SEQUENCE of data point pairs in two columns
						// Will be OVERWRITTEN wich the equally 2D spaced representation
						// x-scaling will be represented by a third, new column, that holds the 
						// position of each point in units of the original index
						
VARIABLE	nPoints  // number of points (x-y pairs) to be returned
// the idea is to create a represenation of the Line that is given by the sequence of x-y pairs
// such that the linear density of points is constant along the line
VARIABLE	NNei		// number of neighbours on either side, using pairwise distance EXCLUDING central point
						// setting 0 regard the two straight lines from central point to  the two sequentially adjacent points
// Some health checks
	IF (DimSize(XYSequence,1)<2)
		DoAlert 0,"Requires a wave of x-y pairs with at least 2 columns!"
		Return -1
	ENDIF
	
	VARIABLE	SrcPnts=DimSize(XYSequence,0)
	IF (SrcPnts<10)	
		DoAlert 0,"Very few Points in Sequence. Aborting."
		Return -2
	ENDIF
	

// step: 1. calculate an local average distance of points by cumulating the cartesian x-y distances between
// points in the vicinity. 
// there could be variants: just local distance, distance from central point to n neihbours either side
// or even define the density at position i as the sum of distances 
// of pairs r_vec(i-k)-r_vec(i+k) for k=1 to n
// the later has an advantage for oscillatory noise on the line

// this steps also requires some scale in x and y and the concept used here is to match ranges
// so each distance in x is seen as a fraction f in [0,1] of the x-range in the data and accordingly for y

// step: 2. the use the local point distance as the required local density (which is the inverse of the actual distance)
//	then re-map the points to achieve an equal distance distribution along the line

// step 1: calculate local average distance of points

	MAKE/O/N=( DimSize(XYSequence,0) )  W_LocalDist
	CopyScales/P XYSequence , W_LocalDist		// thereby use the x scale for determining the positions of points
	// relative scales in x and y:
	// get extreme of first column
	WaveStats/Q/M=1 /R=[0,SrcPnts-1] XYSequence
	VARIABLE	xRange=V_Max-V_min
	WaveStats/Q/M=1 /R=[SrcPnts,2*SrcPnts-1] XYSequence
	VARIABLE	yRange=V_Max-V_min
	IF (min(abs(xRange),abs(yRange))<1e-19)	// at least one dimension is very flat
		DoAlert 0,"There is surprisingly small variation in the data"
		Return -3
	ENDIF
		W_LocalDist				   = 0 

	switch(NNei)
	
		case 0:		// Only the average distance from the central point to the two neighbours
			W_LocalDist[0,SrcPnts-2]  = sqrt( ((XYSequence[p] [1]-XYSequence[p+1][1])/yRange) ^2 + ((XYSequence[p][0]-XYSequence[p+1][0])/xRange) ^2 )
			W_LocalDist[1,SrcPnts-1] += sqrt( ((XYSequence[p] [1]-XYSequence[p-1][1])/yRange) ^2 + ((XYSequence[p][0]-XYSequence[p-1][0])/xRange) ^2 )
			W_LocalDist[0]			    *= 2
			W_LocalDist[SrcPnts-1]	 *= 2
	
			break		// exit from switch
		case 1:		// leaving out the central point, only look at distance between the two neighbours directly
			
			W_LocalDist[0] 		 	= 2*sqrt( ((XYSequence[ 0 ][1]-XYSequence[ 1 ][1])/yRange) ^2 + ((XYSequence[ 0 ][0]-XYSequence[ 1 ][0])/xRange) ^2 ) 
			W_LocalDist[SrcPnts-1]  = 2*sqrt( ((XYSequence[ p ][1]-XYSequence[p-1][1])/yRange) ^2 + ((XYSequence[ p ][0]-XYSequence[p-1][0])/xRange) ^2 ) 
			W_LocalDist[1,SrcPnts-2]=   sqrt( ((XYSequence[p-1][1]-XYSequence[p+1][1])/yRange) ^2 + ((XYSequence[p-1][0]-XYSequence[p+1][0])/xRange) ^2 )
			break
		default:			// optional default expression executed
			VARIABLE	runde
			FOR (runde=1; runde<=NNei; runde+=1)
				W_LocalDist[NNei,SrcPnts-1-NNei]  += sqrt( ((XYSequence[p-NNei][1]-XYSequence[p+NNei][1])/yRange) ^2 + ((XYSequence[p-NNei][0]-XYSequence[p+NNei][0])/xRange) ^2 )
			ENDFOR
			// now rescale those points where fewer points contributed to the cumulative distance
			FOR (runde=1; runde<NNei; runde+=1)
				W_LocalDist[runde]  			 		 *= NNei/runde
				W_LocalDist[SrcPnts-1-runde]      *= NNei/runde 
			ENDFOR

	endswitch 

// step2: place 'nPoints' Points in the xy plane along the line to obtain equal density
// create a wave that should hold the positions and hand it over to external function
	MAKE/O/N=(nPoints) Temp_Pos
	DistributeByDensity(W_LocalDist, Temp_Pos)
	// now Temp_Pos contains the positions (in the units of the index of XYSequence) that are required to 
	// represent XYSequence evenly in XY
	
	// next step is an interpolation of the data at those points and the construction of the resulting xysequence
	// now I quickly check that all positions are within the index range of the input (no extrapolations are required)
	
	Wavestats/Q/M=1 Temp_Pos
	IF ( (V_min < DimOffset(XYSequence,0) )|| (V_max > DimSize(XYSequence,0)-DimDelta(XYSequence,0)) )	// 
		DoAlert 0,"The algorithm requested extrapolation, strange!"
	ENDIF
	
	Duplicate/O/R=[][0] XYSequence, Temp_Seq, W_Xpos
	W_Xpos[]=pnt2x(Temp_Seq,p) // need explicit xpositions for the interp command
	Duplicate/O	Temp_Pos, TempXYT	// to hold interp results
	Redimension/N=(-1,3) TempXYT	// third column is for the implizit parameter (the x-scaling of the original input)
	TempXYT [ ][0]= interp(Temp_Pos[p],W_Xpos,Temp_Seq)
	Temp_Seq[ ]   = XYSequence[p][1]	// now for the y data	
	TempXYT [ ][1]= interp(Temp_Pos[p],W_Xpos,Temp_Seq)
	TempXYT [ ][2]= Temp_Pos[p]
	
	// now overwrite input sequence
	Duplicate/O TempXYT, XYSequence
	

	KillWaves W_LocalDist, Temp_Pos, TempXYT, Temp_Seq, W_Xpos
	Return 1
END

FUNCTION	EqualDistPhaseplotStack(VMatrix, dVMatrix, nPoints)
WAVE		VMatrix, dVMatrix		// matricies that hold a certain number of voltage or rate points sampled at
										// equidistant time intervals, indexed in the wave scale
										
// data is represented at equidistant interval sin the V - dV space using 
VARIABLE	nPoints
// points for each spike; resulting in 'nPoints' voltage and rate values per spike PLUS
// 'nPoints' time values. All stored in a 2D Matrix (nPoints x 3) where the columns hold 0-X; 1-Y and 2-Time
		// as a rule of thumb I would use 100o points for 20 ms of data for 1 spike


// the entire set for all spikes is returned in the 'Reduced' Array. Each layer is the 2D Matri of a spike
		
		VARIABLE	spike, nSpikes=DimSize(VMatrix,1)
		VARIABLE	nTimePnts= DimSize(VMatrix,0)
		
		MAKE/O/N=(nPoints,3,nSpikes) $(NameOfWave(VMatrix)+"_rdx")
		WAVE	Reduced = $(NameOfWave(VMatrix)+"_rdx")
		
		MAKE/O/N=(nTimePnts,2) TempXY // to hand over x and y sequences and retrieve xyt 
		CopyScales VMatrix, TempXY
		
		VARIABLE	resp
		FOR (spike=0; spike < nSpikes; spike +=1)
			TempXY[][0]	=  VMatrix[p][spike]
			TempXY[][1]	= dVMatrix[p][spike]
			if (EquiSample(TempXY,nPoints,5) < 0)
				spike=nSpikes
			endif
			// TempXY is now the reduced (nPoints x 3) data set
			Reduced[][][spike] = TempXY[p][q]
			SetScale/P x,DimOffset(VMatrix,0), DimDelta(VMatrix,0),"", TempXY
			Redimension/N=(nTimePnts,2) TempXY

		ENDFOR
		KillWaves TempXY

END // EqualDistPhaseplotStack






FUNCTION	ExtractSpikes(InWave_ST, InWave_V, PreThreshDur, PostThreshDur, [,Min_Vmx, InWave_Vmx, Min_mxRt, InWave_mxRt, upsamplefac])
WAVE		InWave_ST, InWave_V ,InWave_Vmx,  InWave_mxRt
VARIABLE Min_Vmx, Min_mxRt, upsamplefac, PreThreshDur, PostThreshDur
// the optional parameters allow selecting against certain minimal values of peak potential and peak rate of rise
// the upsampling is set to 3 by default. If no Upsamplin gis desired (e.g. in the case of a 100 kHz sampling rate)
// upsamplefac = 1 has to be used

IF (ParamIsDefault(upsamplefac))
	upsamplefac = 3
ELSEIF (upsamplefac <1 )
	upsamplefac =1
ENDIF


Duplicate/O InWave_ST, Dummy_Criterion
Dummy_Criterion=1

IF (ParamIsDefault(Min_Vmx ))
	Min_Vmx = -1 // -1 Volt - every signal should pass at least this
ELSEIF (ParamIsDefault(InWave_Vmx))	// if selection is intended, the corresponding wave holding the eak voltages for each spike has to be provided
	DoAlert 0, "Please povide the corresponding wave holding the peak potentials for each spike"
	Return -1
ELSE
	Dummy_Criterion *= (InWave_Vmx >= Min_Vmx)

ENDIF	


IF (ParamIsDefault(Min_mxRt ))
	Min_mxRt = -1000 // -1000 Volt/s - every signal should pass at least this
ELSEIF (ParamIsDefault(InWave_mxRt))	// if selection is intended, the corresponding wave holding the eak voltages for each spike has to be provided
	DoAlert 0, "Please povide the corresponding wave holding the peak rate of rise for each spike"
	Return -2
ELSE
	Dummy_Criterion *= (InWave_mxRt >= Min_mxRt)
ENDIF	


// calculate how many points and spikes are needed, create output wave and define the Variables to fill it 
VARIABLE	dt=DimDelta(InWave_V,0)
VARIABLE	npnts= round((PreThreshDur+PostThreshDur)/dt)	, newPnts=(npnts-1)*upsamplefac+1
npnts+=6// include more points before and after to have some lee-way for the interpolation boundary effects
	
VARIABLE	count=0, spike, nspikes = sum(Dummy_Criterion)
MAKE/O/N=(npnts) Dummy_Intermediate
MAKE/O/N=(newpnts, nspikes) $(NameOfWave(InWave_ST) +"_extr")
WAVE	output = $(NameOfWave(InWave_ST) +"_extr")
MAKE/O/N=(newpnts) $(NameOfWave(InWave_ST) +"_avg")
WAVE	avg = $(NameOfWave(InWave_ST) +"_avg")
SEtScale/P x,0,DimDelta(InWave_V,0),WaveUnits(InWave_V,0 ),output,avg

FOR (spike=0; spike<DimSize(InWave_ST,0); spike += 1)
	IF (Dummy_Criterion[spike] >0)
		Dummy_Intermediate[]= InWave_V[x2pnt(InWave_V, InWave_ST[spike]-PreThreshDur )-3 +p]
		SetScale/P x,(InWave_ST[spike]-PreThreshDur-3*dt) ,dt,"s", Dummy_Intermediate
		IF (upsamplefac >1)
			Resample/UP=3 Dummy_Intermediate
		ENDIF
		output[][count] = Dummy_Intermediate[x2pnt(Dummy_Intermediate,InWave_ST[spike]-PreThreshDur) + p]
		Redimension/N=(npnts)	Dummy_Intermediate	
		count+=1
	ENDIF
ENDFOR

MATRIXOP/O/NTHR=0 avg = sumRows(output)/numCols(output)

KillWaves/Z Dummy_Criterion, Dummy_Intermediate

END // ExtractSpikes
	
	
FUNCTION SubThreshGain(InWave_V, InWave_I, InWave_ST, STA_Dur)
WAVE		InWave_V, InWave_I, InWave_ST
VARIABLE	STA_Dur									// Duration of the STA used for the Gain calculation
														// needed only to adjust the window duration to this value
														
// takes current and voltage traces, removes voltage values above the spike threshold (average of the 
// specific voltage trace as specified by the spike property assessment 
// traces are cut into segments of the same length as the spike triggered average (here: 1s)
// standard Fourier transform approach to determine the transfer function between current and voltage
// in the present version the phase is inverted (missing conjugated complex)
														
														
	VARIABLE	WinLen=STA_Dur/DimDelta(InWave_V,0)
	// winlen has to be an even number (FFT)
	WinLen=round(1/DimDelta(InWave_I,0)/2-0.4)*2	// there were suprising round off errors when the simple floor(N/2)*2 was used
	VARIABLE	nWins= floor(DimSize(InWave_V,0)/WinLen)
	
	// Duplicate inputs and output wave; they will be reshaped into a matricies and then overwritten with Fourierresults
	Duplicate/O  Inwave_V, Dummy_V
	Duplicate/O  Inwave_I, Dummy_I
	// Health check for equal length
	IF (DimSize(InWave_I,0) != DimSize(InWave_V,0))
		DoAlert 0,"Length In_V != Length In_I, Aborting"
		Return -1
	ENDIF
	
															
	// remove the spikes by replacing values > V_tresh with V_thresh
	VARIABLE		V_thresh=str2num(StringByKey("V/s)", StringFromList(1, note(InWave_ST),"\r" ),":"," "))
	
	Dummy_V[]= (Dummy_V[p] > V_thresh) ? V_thresh : Dummy_V[p]
	
	// cut length to Matrix Nrows x Mcols  with N= winlen; M=nWins
	Redimension/D/N=(WinLen,nWins) Dummy_I, Dummy_V
	FFT/COLS /DEST=Dummy_V_FFT Dummy_V
	FFT/COLS /DEST=Dummy_I_FFT Dummy_I
	
	KillWaves Dummy_I,Dummy_V
	// by pure coincidence, the magnitude of the current component in any one of the higher frequencies
	// could turn out to be very very small. This would then cause a very large gain, purely due the variation in the 
	// magnitudes due to noise. 
	// this should be avoided
	// therefore: check the current FFT for unusually small magnitudes
	
	// find out whether there is such a small value

	MAKE/O/N=(round(1000*STA_Dur),DimSize(Dummy_I_FFT,1)) Dummy // only care about frequencies below 1kHz, 
																				// those are in the first 1000*Sta_dur rows
		Dummy[][]=cabs(Dummy_I_FFT[p][q])
	VARIABLE reps=0
	DO
		Dummy[][]=cabs(Dummy_I_FFT[p][q])
		wavestats/Q/M=1 Dummy
		Dummy_I_FFT[V_minRowLoc][V_minColLoc]+=cmplx(1e-9,1e-9)
		reps+=1
	WHILE ( (reps < 1000) && (V_min<1e-9))

	// calculate complex gain
	Dummy_V_FFT/= Dummy_I_FFT
	KillWaves Dummy_I_FFT

	STRING TargetName=NameOfWave(InWave_V)
	TargetName=TargetName[0,strlen(TargetName)-2]+"subG"
	// Average across columns
	MatrixOp/C/NTHR=3  /O  $TargetName = sumRows(Dummy_V_FFT)/numCols(Dummy_V_FFT) 	

END // SubThreshGain

FUNCTION SpikeStats(InWave_V, InWave_ST,[RateThresh])
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
			DoAlert 0,"The range of threshold voltages related to the spike times\ris "+num2str((V_max-V_min)*1000)+" mV.\rThat seems too big, aborting..."
			KillWaves dummy
			Return -1
		
		ENDIF

	// decompose name of inwaves to create new names: remove ending after last"_"
		VARIABLE		N_underScore=ItemsInList(NameOfWave(InWave_V),"_")
		STRING		SlopeName, ThreshName, VmaxName, mxRateName, VofMaxRateName, ISIName,FWHMName, AHPName, ThreshCrossTime
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
		FWHMName=SlopeName+"FWHM"
		ISIName=SlopeName+"ISI"
		AHPName=SlopeName+"AHP"
		ThreshCrossTime=SlopeName+"SonsetT"
		SlopeName=SlopeName+"RPD"

		Duplicate /O InWave_ST, $SlopeName, $ThreshName, $VmaxName, $mxRateName, $VofMaxRateName, $ISIName, $FWHMName, $AHPName, $ThreshCrossTime
		WAVE	slope=$SlopeName
		WAVE	thresh=$ThreshName
		WAVE	Vmx=$VmaxName
		WAVE	mxRt=$mxRateName
		WAVE	VmxRt=$VofMaxRateName
		WAVE	FWHM=$FWHMName
		WAVE	ISI=$ISIName  
		WAVE	AHP=$AHPName  
		WAVE	STX=$ThreshCrossTime  
		slope	=NaN
		thresh=NaN
		Vmx	=NaN
		mxRt	=NaN
		VmxRt	=NaN
		FWHM	=NaN
		AHP 	=NaN
		STX 	=NaN
		SetScale d 0,0,"1/s", slope
		SetScale d 0,0,"V", thresh, VmxRt, Vmx, AHP
		SetScale d 0,0,"V/s", mxRt
		SetScale d 0,0,"s", ISI
		SetScale d 0,0,"s", FWHM
		SetScale d 0,0,"s", STX
		// make sure all waves that were created have index scaling
		SetScale/P x,0,1, slope,thresh, VmxRt, Vmx,ISI, mxRt, FWHM, AHP, STX

// get local variation coefficient and append to wave note
		Note/K 	InWave_ST
		Note/NOCR InWave_ST, "local variation of spike times (LV):"
		Note/NOCR InWave_ST, num2str(ReturnLocalCV(InWave_ST))
// 10/2016 aNNe:  decided to append coefficient of variation as well.
		Note/NOCR InWave_ST, "; coef. of variation of spike times (CV):"
		Note/NOCR InWave_ST, num2str(ReturnCV(InWave_ST))

// find threshold as defined by rate threshold 
		IF (ParamIsDefault(RateThresh))
	           	RateThresh = 30 //(V/s)
		ENDIF

		VARIABLE		nSpikes=DimSize(InWave_ST,0), sp, sPnt// Pnt number in Voltage wave where spike was detected
		VARIABLE	 	nPntsBefore, nPntsAfter, upsamplefac=round(100000*dt)// stretch of data before and after that is copied for analysis
																			// factor of how much the temporal sampling is increased
		VARIABLE		nPntsInDummy, SlopeMaxWasApplied=0, SlopeMax=3000000
		FOR (sp=0; sp< nSpikes; sp +=1)
			
			// copy a stretch of voltages around the detection point
			// use a max of 2 ms before but no more than the distance to the previous spike
			// use a max of 1 ms after 
			
			IF (sp >0)	// not the first spike)
				nPntsBefore=min(0.002, InWave_ST[sp]-InWave_ST[sp-1])/dt
			ELSE
				nPntsBefore=0.002/dt
			ENDIF
			IF (sp+1 < nSpikes)	// not the last spike)
				nPntsAfter=min(0.001, InWave_ST[sp+1]-InWave_ST[sp])/dt
			ELSE
				nPntsAfter=0.001/dt
			ENDIF
			
			// calculate slope
			// check for drop below threshold slope

			
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
			
			IF (thresh[sp]!=thresh[sp]) // this is the case if threshold was could not be determined and is set to NaN
				FWHMval = NaN
				minAmp = -0.02		// -20 mV as a default
			ELSE
				HM=thresh[sp]+0.5*(Vmx[sp]-thresh[sp])
				FindLevel /EDGE=1 /Q/R=(InWave_ST[sp]-0.0015,InWave_ST[sp]+0.0015) InWave_V,HM
				IF (V_Flag==0)
					 FWHMval=0 -	V_LevelX
					 FindLevel /EDGE=2 /Q/R=(V_LevelX+0.0001,InWave_ST[sp]+0.0045) InWave_V,HM
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
			Variable endOfAHPSearch
			IF (sp< nSpikes-1)
				endOfAHPSearch=x2pnt(InWave_V,InWave_ST[sp+1])
			ELSE	
			   endOfAHPSearch=DimSize(Inwave_V,0)-1
			ENDIF
			FindPeak /B=(max(10,0.0001/dt))/M=(minAmp) /N /P/Q/R=[sPnt,endOfAHPSearch] InWave_V	
			IF (V_flag==0)	// minimum found
				AHP[sp]= V_PeakVal
			ENDIF
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
			RemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(thresh)
		ENDIF
		Note		 InWave_ST, "average threshold (rate crosses "+num2str(RateThresh)+" V/s):"
		Note/NOCR InWave_ST, num2str(Reporter) + " V"
		
		IF (mean(slope)!=mean(slope) ) // contains NaNs
			Duplicate/O slope, dummy
			RemoveNaNs(dummy)
			Reporter=Median(dummy)
		ELSE
			Reporter=Median(slope)
		ENDIF
		
		Note		 InWave_ST, "median phase slope at threshold:"
		Note/NOCR InWave_ST, num2str(Reporter) + " 1/s"
		
		Note		 InWave_ST, "average peak potential:"
		Note/NOCR InWave_ST, num2str(mean(Vmx)) + " V"
		Note		 InWave_ST, "average peak rate of rise (PRR):"
		Note/NOCR InWave_ST, num2str(mean(mxRt)) + " V/s"
		Note		 InWave_ST, "average voltage at PRR:"
		Note/NOCR InWave_ST, num2str(mean(VmxRt)) + " V"
		Note		 InWave_ST, "average ISI:"
		Note/NOCR InWave_ST, num2str(mean(ISI,1,nSpikes-1)) + " s"
		
		IF (mean(FWHM)!=mean(FWHM) ) // contains NaNs
			Duplicate/O FWHM, dummy
			RemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(FWHM)
		ENDIF
				
		Note		 InWave_ST, "average FWHM:"
		Note/NOCR InWave_ST, num2str(Reporter)  + " s"
		
		IF (mean(AHP)!=mean(AHP) ) // contains NaNs
			Duplicate/O AHP, dummy
			RemoveNaNs(dummy)
			Reporter=mean(dummy)
		ELSE
			Reporter=mean(AHP)
		ENDIF
			Note		 InWave_ST, "average AHP:"
		Note/NOCR InWave_ST, num2str(Reporter) + " V"
		
		
		// copy the wave note of the spike time wave to the notes of all derived waves
		note/K/NOCR slope, note(InWave_ST)
		note/K/NOCR thresh , note(InWave_ST)
		note/K/NOCR Vmx, note(InWave_ST)
		note/K/NOCR mxRt, note(InWave_ST)
		note/K/NOCR VmxRt, note(InWave_ST)
		note/K/NOCR FWHM, note(InWave_ST)
		note/K/NOCR ISI, note(InWave_ST)  
		note/K/NOCR AHP, note(InWave_ST)  
		note/K/NOCR STX, note(InWave_ST)
		
KillWaves/Z Dummy
END 

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
	STRING VofMaxRateName=BaseName+"VmxRt"
	STRING FWHMName=BaseName+"FWHM"
	STRING ISIName=BaseName+"ISI"
	STRING AHPName=BaseName+"AHP"
	STRING ThreshCrossTime=BaseName+"SonsetT"
	STRING SlopeName=BaseName+"RPD"
	STRING LVName= Basename +"LV"
		MAKE/N=(numFolders) /O $SlopeName, $ThreshName, $VmaxName, $mxRateName, $VofMaxRateName, $ISIName, $FWHMName, $AHPName, $ThreshCrossTime, $LVName
		WAVE	slope=$SlopeName
		WAVE	thresh=$ThreshName
		WAVE	Vmx=$VmaxName
		WAVE	mxRt=$mxRateName
		WAVE	VmxRt=$VofMaxRateName
		WAVE	FWHM=$FWHMName
		WAVE	ISI=$ISIName  
		WAVE	AHP=$AHPName  
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
		SetScale d 0,0,"V/s", mxRt
		SetScale d 0,0,"s", ISI
		SetScale d 0,0,"s", FWHM
		SetScale d 0,0,"", LV
		
		// make sure all waves that were created have index scaling
		SetScale/P x,0,1, slope,thresh, VmxRt, Vmx,ISI, mxRt, FWHM, AHP, LV

		
	VARIABLE	cST, nST
	VARIABLE	LVavg, THRavg, SLOPEavg, VMXavg, MXRTavg, VMXRTavg, ISIavg, FWHMavg, AHPavg
	STRING	STListe, Notiz	
	STRING	LVstring, THRstring, SLOPEstring, VMXstring, MXRTstring, VMXRTstring, ISIstring, FWHMstring, AHPstring
	STRING	currFolderName
	

	FOR (runde=0; runde< numFolders; runde+=1)
		currFolderName=StringFromList(runde,FolderListe,";")
		SetDataFolder $(currFolderName)
		LVavg=0
		THRavg=0
		SLOPEavg=0
		VMXavg=0		 
		MXRTavg=0
		VMXRTavg=0
		ISIavg=0
		FWHMavg=0
		AHPavg=0

		
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
			MXRTstring	= StringByKey("(PRR)",Notiz,":"," ")
			VMXRTstring	= StringByKey("PRR",Notiz,":"," ")
			ISIstring	= StringByKey("ISI",Notiz,":"," ")
			FWHMstring	= StringByKey("FWHM",Notiz,":"," ")
			AHPstring	= StringByKey("AHP",Notiz,":"," ")
			
			// cumulate numbers
			LVavg+=str2num(LVstring)
			THRavg+=str2num(THRstring)
			SLOPEavg+=str2num(SLOPEstring)
			VMXavg+=str2num(VMXstring)
			MXRTavg+=str2num(MXRTstring)
			VMXRTavg+=str2num(VMXRTstring)
			ISIavg+=str2num(ISIstring)
			FWHMavg+=str2num(FWHMstring)
			AHPavg+=str2num(AHPstring)
	
		ENDFOR // loop across waves in a given folder
		slope	[runde]=SLOPEavg/nST
		thresh[runde]=THRavg/nST
		Vmx	[runde]=VMXavg/nST
		mxRt	[runde]=MXRTavg/nST
		VmxRt	[runde]=VMXRTavg/nST
		FWHM	[runde]=FWHMavg/nST
		AHP 	[runde]=AHPavg/nST
		ISI 	[runde]=ISIavg/nST
		LV		[runde]=LVavg/nST
		
		SetDimLabel 0,runde,$(currFolderName),	slope,thresh,Vmx,mxRt,VmxRt,FWHM,AHP,LV,ISI	


		SetDataFolder topDf	
	ENDFOR // loop across folders
END

FUNCTION ReturnLocalCV(InWave)
// this returns a single number, the average local variation 
// according to equation 2 in 
// http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000433
// Shinomoto et al. 2009
WAVE		InWave // contains spike times!
		
		VARIABLE	LV, kk, nS=DimSize(InWave,0) // number of spikes
		LV=0
		FOR (kk=0; kk<nS-2; kk+=1)
			LV+= ( ((InWave[kk+1]-InWave[kk])-(InWave[kk+2]-InWave[kk+1]))/(InWave[kk+2]-InWave[kk]) )^2
			
		ENDFOR 
			LV*=3/(nS-2)		// ns is number of spikes, the original paper used number of ISI, hence here we divide by Nspikes-2
									// while Shinomoto divides by nIntervals-1
			Return LV
END		//ReturnLocalCV
FUNCTION ReturnCV(InWave)
// this returns a single number, the coefficient of variation of the inter-spike intervals

WAVE		InWave // contains spike times!
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

FUNCTION AnaIDProtocol()
// start in folder containing voltage (*_V) and current (_I) traces

	STRING	V_list=WaveList("*_V",";","")
	STRING	I_list=WaveList("*_I",";","")
	// healthcheck
	Variable Num_V=ItemsInList(V_list,";")
	Variable Num_I=ItemsInList(I_list,";")
	IF (Num_I != Num_V)
		DoAlert 0,"Number of voltage and current traces differs ("+num2istr(Num_V)+" vs "+num2istr(Num_I)+"."
		Num_V=min(num_V,num_I)	// continue with the minimal number
	ENDIF
	
	// attempt to obtain a few key parameters automatically
	// assume that two large transients in the stimulus wave betray the stimulus onset and offset
	// 
	VARIABLE BaselineCurrent
	MAKE/FREE/N=(num_V) W_OnsetT,W_OffsetT
	
	VARIABLE rr, thresh, min_dT
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
	
	VARIABLE onsetT=median(W_OnsetT)
	VARIABLE StimDuration=median(W_OffsetT)-median(W_OnsetT)
	KillWaves W_OnsetT,W_OffsetT
	
	printf "estiamted onset %g, duration %g \r",onsetT,StimDuration
	
	Make/O/N=(num_V) FirstSpikeLate, InitRate, Rate, StimAmp
	FirstSpikeLate=NaN
	InitRate = 0
	VARIABLE	BaseLineStim	
	FOR (rr=0; rr<num_V; rr++)
		WAVE stim = $(StringFromList(rr,I_list,";"))
		BaseLineStim+= mean(stim,DimOffset(stim,0),OnsetT-DimDelta(stim,0))
		StimAmp[rr] = mean(stim,OnsetT+DimDelta(stim,0),OnsetT+StimDuration-DimDelta(stim,0))
	ENDFOR
	BaseLineStim/=num_V
	StimAmp-=BaseLineStim
	STRING W_Name
	FOR (rr=0; rr<num_V; rr++)
		W_Name = StringFromList(rr,V_list,";")
		WAVE resp = $(W_Name)
		Rate[rr] = ReturnSpikeTimes(resp,V_Threshold=0, MinISI=0.002, Bool_RatePositive=1)
		IF (Rate[rr] > 0)
			WAVE ST=$(W_Name[0,strlen(W_Name)-2]+"ST")
			FirstSpikeLate[rr] = ST[0]-onsetT
			IF (Rate[rr] > 1)
				InitRate[rr] = 1/( ST[1]-ST[0])
			ENDIF
		ENDIF
	ENDFOR
	rate/=StimDuration
	IF ((stringmatch("s",WaveUnits(stim,0))) || (stringmatch("S",WaveUnits(stim,-1)))|| (stringmatch("SEC",WaveUnits(stim,-1)))|| (stringmatch("sec",WaveUnits(stim,-1))) )
		SetScale/P d,0,0,"Hz", InitRate, Rate
		SetScale/P d,0,0,"s",FirstSpikeLate
	ENDIF
	SetScale/P d,0,0,WaveUnits(stim,-1), StimAmp


END // AnaIDProtocol


FUNCTION ReturnSpikeTimes(InWave,[ V_Threshold, MinISI, Bool_RatePositive])
// this returns spike times in the time unit of the wave
// spike time is time of crossing Threshold, which defaults to 0 mV

// the diversity of spike shapes requires a detection criterion that is more sophisticated than a 
// simple crossing of a threshold from below
// main problem was noise in very slowly falling spikes

// the optional parameters work as follows:
// V_Threshold  - crossing this from below makes a sample a candidatre for spike time
// 					in the absence of other optional parameters thats the only criterion and the default for the
// 					threshold value is 0 mV
// MinISI  	    - a candidate spike time is rejected if it occurs sooner than MinISI (seconds) after the last accepted spike
// Bool_RatePositive  	    - a candidate spike time is rejected if it occurs
// 					in a region (2 pnts before + 2pnts after) in which the rate of rise is on average negative 

WAVE		InWave			// contains V_mem vs time
VARIABLE	V_Threshold, MinISI, Bool_RatePositive

IF (ParamIsDefault(V_Threshold ) )
	V_Threshold = 0
ENDIF
IF (ParamIsDefault(MinISI ) )
	MinISI = 0.0002 //(5 kHz instantaneous rate)
ENDIF
IF (ParamIsDefault(Bool_RatePositive ) )
	Bool_RatePositive = 0
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
	IF (wavemax(Dummy)>2000 ||  wavemin(Dummy)<-600 ) 
// there are probably artefacts
		KillWaves/Z Dummy
		DoAlert 0,"Presumed artefact in "+GetWavesDataFolder(InWave,2)
		print  GetWavesDataFolder(InWave,2)
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
		ST[]=(InWave[x2pnt(Inwave,ST[p])+4]> InWave[x2pnt(Inwave,ST[p])-4]) ? ST[p] : NaN
		RemoveNaNs(ST)
	ENDIF
	// now call spikestatistics
	// this is optional, but makes sense to do right now
	// if there are any spikes at all
	IF (DimSize(ST,0)>0)
		SpikeStats(InWave, ST,RateThresh=30)
	ENDIF
	//make a note about the threshold that was used:
	Note ST,"Threshold:"+num2str(V_Threshold)+" "+Waveunits(InWave,1)
	return DimSize(ST,0)
END

FUNCTION/C ReturnVectorStrengthVector(STWave,Freq)
WAVE		STWave			// contains spike times
VARIABLE	Freq	//Freq in Hz, ZeroPhase: for a sine starting at time ZERO, this is Pi/2
		
		
		VARIABLE	counter, NSpikes=DimSize(STWave,0)
		VARIABLE/C	vector=cmplx(0,0)
		
		FOR (counter=0; counter < NSpikes; counter+=1)
			vector+=exp(cmplx(0,2*Pi*STWave[counter]*Freq))
		ENDFOR
		
		RETURN		vector/NSpikes
END

Threadsafe FUNCTION/WAVE STAfromAnalogueWaveRef(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix, FirstPoint, LastPoint])
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp
												// here it will be only needed to determine with spikes 
												// can be used to calculate the STA
												// because there is sufficient current available before and after the spike
	VARIABLE	FP,LP	
	STRING	SX										

		IF (ParamIsDefault(FirstPoint))
			FP = 0
		ELSE
			FP=FirstPoint
		ENDIF
		IF (ParamIsDefault(LastPoint))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
			LP=DimSize(AnalogueTrace,0)-1
		ELSE
			LP = LastPoint
		ENDIF
		
		IF (ParamIsDefault(Suffix))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
			SX = "STA"
		
		ELSE
			
			SX=Suffix
		ENDIF
			
		STAfromAnalogue(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar, Suffix=SX, FirstPoint=FP, LastPoint=LP)

		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		STA_Name
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
		STA_Name=RemoveEnding(I_Name,Ending)+SX


		WAVE Dummy=$STA_Name
		
		Return Dummy
END



Threadsafe FUNCTION 	STAfromAnalogue(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix, FirstPoint, LastPoint])
// calculates Spike triggered average of the snipletts (length 'range') of the Analogue trace
// the snipletts are centered around the entries of spike Times
// range has to have the same unit as the analogue data 
// the setscale command assumes it is "s"
// if DoSpikeTrigVar is not ZERO, the spike triggered variance is calculated
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp
												// here it will be only needed to determine with spikes 
												// can be used to calculate the STA
												// because there is sufficient current available before and after the spike
												
		IF (ParamIsDefault(FirstPoint))
			FirstPoint = 0
		ENDIF
		// check for existance
		IF (!Waveexists(AnalogueTrace))
			Print  "Aborting STA from analogue, could not find AnalogueTrace"+ NameOfWave(AnalogueTrace)
				Return -1
		ELSEIF (!WaveExists(SpikeTimes))
			Print "Aborting STA from analogue, couldn't find SpikeTimes"+ NameOfWave(SpikeTimes)
				Return -1
		ENDIF
		// check for having more than zero spikes
		IF (DimSize(SpikeTimes,0)<1)
			Print "Aborting STA from analogue, SpikeTime wave was empty "+ NameOfWave(SpikeTimes) + NameOfWave(AnalogueTrace)
			Return -1
		ENDIF
		
		
		VARIABLE	AT_Zero= DimOffset(AnalogueTrace,0)
		VARIABLE	AT_deltaT=DimDelta(AnalogueTrace,0)
		VARIABLE	nInAna=DimSize(AnalogueTrace,0)-FirstPoint
		IF (!ParamIsDefault(LastPoint))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
													
			nInAna -= DimSize(AnalogueTrace,0)-1-min(DimSize(AnalogueTrace,0)-1,LastPoint)
		ELSE
			LastPoint=DimSize(AnalogueTrace,0)-1
		ENDIF
		
		IF (ParamIsDefault(Suffix))
			Suffix = "STA"
		ENDIF
		
		
//		IF ( abs(SpikeTimes[0]-AT_Zero) > 2 )
//			DoAlert 1,"Zero time of analogue trace is rather different from first spike time, ("+num2str(SpikeTimes[0]-AT_Zero)+"s). Continue anyway?"
//			IF (V_flag==2)	// Continue
//				Return -1
//			ENDIF
//		ENDIF
		
		VARIABLE		k, nSp=DimSize(Spiketimes,0)
		VARIABLE 	pRange=2*round(range/AT_deltaT/2)+1
		VARIABLE		p1, p2,n1, n2
		
		STRING		ST_Name=NameOfWave(SpikeTimes)
		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		STA_Name
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
//		IF (cmpstr("ST", StringFromList(ItemsInList(ST_Name,"_")-1,ST_Name,"_")) ==0)
			 STA_Name=RemoveEnding(I_Name,Ending)+Suffix
//		ELSE
//			STA_Name=ST_Name+"_STA"
//		ENDIF
		MAKE/D/O/N=(pRange) $STA_Name
		WAVE avg=$STA_Name
		avg=0
		SetScale/P x -(pRange-1)/2*AT_deltaT,AT_deltaT,"s", avg

		IF (DoSpikeTrigVar)
			MAKE/D/O/N=(pRange) $(NameOfWave(SpikeTimes)+"_STV")
			WAVE var=$(NameOfWave(SpikeTimes)+"_STV")
			var=0
			SetScale/P x -(pRange-1)/2*AT_deltaT,AT_deltaT,"s", avg, var
		ENDIF
		n1=0
		
		VARIABLE	firstST=SpikeTimes[nSp-1]			// first Spike time to be included in analysis
		VARIABLE	lastST=SpikeTimes[0]				// last Spike time to be included in analysis
			
		// calc avg
		FOR (k=0; k< nSp; k+=1)
			// not every spike can be used: the onces too close to begin or end of the stimulus
			// cannot be used for the STA as the input is not completely known over the neccessary range
			p1=x2pnt(AnalogueTrace, SpikeTimes[k] )-(pRange-1)/2
			p2=p1+pRange-1
			IF ( (p1>=FirstPoint) && (p2< FirstPoint+nInAna) ) // if input is available for the entire range needed for STA
				avg[0,pRange-1]+=AnalogueTrace[p+p1]
				firstST=min(firstST,SpikeTimes[k])
				lastST=max(lastST,SpikeTimes[k])
				n1+=1
			ENDIF
		ENDFOR
		avg/=n1  //divide by #spikes
		// if there are no points in ST, return -1, delete STA wave
		IF (n1<1)
			KillWaves avg
			Return -2
		ENDIF


// place some info in the Wave note 
		Wavestats/R=[FirstPoint, min(DimSize(AnalogueTrace,0)-1,LastPoint)]/Q AnalogueTrace
		String wavnote="STA from "+ NameOfWave(AnalogueTrace)
		wavnote +=" triggered on the "+num2istr(n1)+" spikes in "+NameOfWave(SpikeTimes)+".\r"
		wavnote +="spikerate:"+num2str(n1/nInAna/DimDelta(AnalogueTrace,0))+"\r"
		wavnote +="#spikes:"+num2str(n1)+"\r"
		wavnote +="duration:"+num2str(DimDelta(AnalogueTrace,0)*nInAna)+"\r" // sweep length
		wavnote +="AVG:"+ num2str(V_avg)+"\rSD:"+num2str(V_sdev)+"\r"
		
		Note/K avg
		Note  avg, wavnote
		
		// as the time intervall over which spikes occur could be much shorter than the total duration of the analogue trace,  only 
		// the time intervall between first and last included spike is used for rate calculation in the line above

		// calc stv
		IF (DoSpikeTrigVar)
			n2=0
			FOR (k=0; k< nSp; k+=1)
				p1=x2pnt(AnalogueTrace, SpikeTimes[k] )-(pRange-1)/2
				p2=p1+pRange-1
				IF ( (p1>=0) && (p2< nInAna) )
					var[0,pRange-1]+=(AnalogueTrace[p+p1]-avg[p])^2
					n2+=1
				ENDIF
			ENDFOR
			var/=n2-1
			IF (n2!=n1)
				// DoAlert 0,"Different number of spikes accounted for on the two iterations, abort"
				avg=0
				var=0
				return -1
			ENDIF
			Note/K var
			Note   var, wavnote

		ENDIF 	// if DoSpikeTrigVar
		


END

FUNCTION 	STAfromAnalogueNONThreadsafe(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix, FirstPoint, LastPoint])
// calculates Spike triggered average of the snipletts (length 'range') of the Analogue trace
// the snipletts are centered around the entries of spike Times
// range has to have the same unit as the analogue data 
// the setscale command assumes it is "s"
// if DoSpikeTrigVar is not ZERO, the spike triggered variance is calculated
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp
												// here it will be only needed to determine with spikes 
												// can be used to calculate the STA
												// because there is sufficient current available before and after the spike
												
		IF (ParamIsDefault(FirstPoint))
			FirstPoint = 0
		ENDIF
		// check for existance
		IF (!Waveexists(AnalogueTrace))
			Print  "Aborting STA from analogue, could not find AnalogueTrace"+ NameOfWave(AnalogueTrace)
			Return -1
		ELSEIF (!WaveExists(SpikeTimes))
			Print "Aborting STA from analogue, couldn't find SpikeTimes"+ NameOfWave(SpikeTimes)
			Return -1
		ENDIF
		// check for having more than zero spikes
		IF (DimSize(SpikeTimes,0)<1)
			Print "Aborting STA from analogue, SpikeTime wave was empty "+ NameOfWave(SpikeTimes) + NameOfWave(AnalogueTrace)
			Return -1
		ENDIF
		
		VARIABLE	AT_Zero= DimOffset(AnalogueTrace,0)
		VARIABLE	AT_deltaT=DimDelta(AnalogueTrace,0)
		VARIABLE	nInAna=DimSize(AnalogueTrace,0)-FirstPoint
		IF (!ParamIsDefault(LastPoint))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
													
			nInAna -= DimSize(AnalogueTrace,0)-1-min(DimSize(AnalogueTrace,0)-1,LastPoint)
		ELSE
			LastPoint=DimSize(AnalogueTrace,0)-1
		ENDIF
		
		IF (ParamIsDefault(Suffix))
			Suffix = "STA"
		ENDIF
		
		
//		IF ( abs(SpikeTimes[0]-AT_Zero) > 2 )
//			DoAlert 1,"Zero time of analogue trace is rather different from first spike time, ("+num2str(SpikeTimes[0]-AT_Zero)+"s). Continue anyway?"
//			IF (V_flag==2)	// Continue
//				Return -1
//			ENDIF
//		ENDIF
		
		VARIABLE		k, nSp=DimSize(Spiketimes,0)
		VARIABLE 	pRange=2*round(range/AT_deltaT/2)+1
		VARIABLE		p1, p2,n1, n2
		
		STRING		ST_Name=NameOfWave(SpikeTimes)
		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		STA_Name
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
//		IF (cmpstr("ST", StringFromList(ItemsInList(ST_Name,"_")-1,ST_Name,"_")) ==0)
			 STA_Name=RemoveEnding(I_Name,Ending)+Suffix
//		ELSE
//			STA_Name=ST_Name+"_STA"
//		ENDIF
		MAKE/D/O/N=(pRange) $STA_Name
		WAVE avg=$STA_Name
	
		avg=0
		SetScale/P x -(pRange-1)/2*AT_deltaT,AT_deltaT,"s", avg

		IF (DoSpikeTrigVar)
			MAKE/D/O/N=(pRange) $(NameOfWave(SpikeTimes)+"_STV")
			WAVE var=$(NameOfWave(SpikeTimes)+"_STV")
			var=0
			SetScale/P x -(pRange-1)/2*AT_deltaT,AT_deltaT,"s", avg, var
		ENDIF
		n1=0
		
		VARIABLE	firstST=SpikeTimes[nSp-1]			// first Spike time to be included in analysis
		VARIABLE	lastST=SpikeTimes[0]				// last Spike time to be included in analysis
			
		// calc avg
		FOR (k=0; k< nSp; k+=1)
			// not every spike can be used: the onces too close to begin or end of the stimulus
			// cannot be used for the STA as the input is not completely known over the neccessary range
			p1=x2pnt(AnalogueTrace, SpikeTimes[k] )-(pRange-1)/2
			p2=p1+pRange-1
			IF ( (p1>=FirstPoint) && (p2< FirstPoint+nInAna) ) // if input is available for the entire range needed for STA
				avg[0,pRange-1]+=AnalogueTrace[p+p1]
				firstST=min(firstST,SpikeTimes[k])
				lastST=max(lastST,SpikeTimes[k])
				n1+=1
			ENDIF
		ENDFOR
		avg/=n1  //divide by #spikes
		// if there are no points in ST, return -1, delete STA wave
		IF (n1<1)
			KillWaves avg
			Return -2
		ENDIF


// place some info in the Wave note 
		Wavestats/R=[FirstPoint, min(DimSize(AnalogueTrace,0)-1,LastPoint)]/Q AnalogueTrace
		String wavnote="STA from "+ NameOfWave(AnalogueTrace)
		wavnote +=" triggered on the "+num2istr(n1)+" spikes in "+NameOfWave(SpikeTimes)+".\r"
		wavnote +="spikerate:"+num2str(n1/nInAna/DimDelta(AnalogueTrace,0))+"\r"
		wavnote +="#spikes:"+num2str(n1)+"\r"
		wavnote +="duration:"+num2str(DimDelta(AnalogueTrace,0)*nInAna)+"\r" // sweep length
		wavnote +="AVG:"+ num2str(V_avg)+"\rSD:"+num2str(V_sdev)+"\r"
		
		Note/K avg
		Note  avg, wavnote
		
		// as the time intervall over which spikes occur could be much shorter than the total duration of the analogue trace,  only 
		// the time intervall between first and last included spike is used for rate calculation in the line above

		// calc stv
		IF (DoSpikeTrigVar)
			n2=0
			FOR (k=0; k< nSp; k+=1)
				p1=x2pnt(AnalogueTrace, SpikeTimes[k] )-(pRange-1)/2
				p2=p1+pRange-1
				IF ( (p1>=0) && (p2< nInAna) )
					var[0,pRange-1]+=(AnalogueTrace[p+p1]-avg[p])^2
					n2+=1
				ENDIF
			ENDFOR
			var/=n2-1
			IF (n2!=n1)
				// DoAlert 0,"Different number of spikes accounted for on the two iterations, abort"
				avg=0
				var=0
				return -1
			ENDIF
			Note/K var
			Note   var, wavnote

		ENDIF 	// if DoSpikeTrigVar
		


END

// Boostrap code for average across cells

FUNCTION/S 	ST_inputSectionsfromAnalogue(AnalogueTrace,SpikeTimes, Range)
WAVE		AnalogueTrace, SpikeTimes
VARIABLE	Range
// creates index of spike triggered sections (p1,p2, SD of entire input wave)


// the sections are centered around the entries of spike Times
// range has to have the same unit as the analogue data 
// the setscale command assumes it is "s"
		VARIABLE	AT_Zero= DimOffset(AnalogueTrace,0)
		VARIABLE	AT_deltaT=DimDelta(AnalogueTrace,0)
		VARIABLE	nInAna=DimSize(AnalogueTrace,0)
		
//		IF ( abs(SpikeTimes[0]-AT_Zero) > 10 )
//			DoAlert 1,"Zero time of analogue trace is rather different from first spike time, ("+num2str(SpikeTimes[0]-AT_Zero)+"s). Continue anyway?"
//			IF (V_flag==2)	// Continue
//				Return ""
//			ENDIF
//		ENDIF
		
		VARIABLE	k, nSp=DimSize(Spiketimes,0)
		VARIABLE 	pRange=2*round(range/AT_deltaT/2)+1
		VARIABLE	p1, p2,n1, n2
		VARIABLE	SD_I=sqrt(Variance(AnalogueTrace,0.75,44.25)) // this limitation in the time is only valid for Omers 45 s stimuli
		
		
		MAKE/D/O/N=(0,3) $(NameOfWave(SpikeTimes)+"_M")	// this allows me to add 
														// row by row
		WAVE BootIndx=$(NameOfWave(SpikeTimes)+"_M")
		
		n1=0
		
		VARIABLE	firstST=SpikeTimes[nSp-1]			// first Spike time to be included in analysis
		VARIABLE	lastST=SpikeTimes[0]				// last Spike time to be included in analysis
			
		// calc avg
		FOR (k=0; k< nSp; k+=1)
			p1=x2pnt(AnalogueTrace, SpikeTimes[k] )-(pRange-1)/2
			p2=p1+pRange-1
			IF ( (p1>=0) && (p2< nInAna) )
				BootIndx[DimSize(BootIndx,0)]={{p1},{p2},{SD_I}}    // add one row
				SetDimLabel 	0,(DimSize(BootIndx,0)-1), $(NameOfWave(AnalogueTrace)), BootIndx
				firstST=min(firstST,SpikeTimes[k])
				lastST=max(lastST,SpikeTimes[k])
				n1+=1
			ENDIF
		ENDFOR
// place some info in the Wave note 
		Wavestats/Q AnalogueTrace
		STRING wavnote="Spike triggered input indicies from "+ NameOfWave(AnalogueTrace)
		wavnote +=" triggered on the "+num2istr(n1)+" spikes in "+NameOfWave(SpikeTimes)+".\r"
		wavnote +="The avg and sdev of the analogue wave are "+ num2str(V_avg)+" and "+num2str(V_sdev) + " respectively.\r"
		wavnote +="The average spike rate is "+num2str((n1-1) / (lastST-firstST))
		
		Note/K BootIndx
		Note  BootIndx, wavnote
		Return NameOfWave(BootIndx)

END // ST_inputSectionsfromAnalogue 



//¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
//-------------    B O O T S T R A P   C O D E        --------------------
//________________________________________________________________________


FUNCTION BST_Wrapper(nRounds,DoNoiseFloor[, earliestST, latestST, FirstPoint,LastPoint, SpikeTimeSuffix])
VARIABLE		nRounds	//number of bootstrap rounds
// This calls all the necessary procedures for confidence interval bootstrap
// the waves AC_avg_scaled and FreqPoints have to exist before this can be run!
VARIABLE		DoNoiseFloor		// if DoNoiseFloor> 0, the spike times get shifted
										//and the outcome represents the noise floor rather than the confidence interval
VARIABLE		FirstPoint,LastPoint
// for Omers data this needs to be limited to the relevant part of the current
// points 14532 to 915469
// so for Omer this has to be called
// MakeAndAvgAC(FirstPoint=14532,LastPoint=915469)

VARIABLE		earliestST,latestST	// those are boundaries of possible spike times; needed in order to properly shift spike times, i.e. only needed for NoiseFloor 
	// if not provided they will be approximated by the time of the first and last spike

STRING		SpikeTimeSuffix // standard is "ST", but for special cases, it is usefull to define otherwise, for instance STwC ("with Criteria")


IF (ParamIsDefault(FirstPoint))
	FirstPoint = 0
ENDIF
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
		
		IF (ParamIsDefault(FirstPoint))
			PrepareST_List(SpikeTimeSuffix)
			// creates the IndexWave "BST_SpikeIndicies"
			// and the CountWave "BST_SpikeCount"
			// also creates 2 WAVES of WaveREFERENCES, which hold all STWaves and all IWaves
			// which are used, in the sequence in which they are listed in the 
		ELSE
			PrepareST_List(SpikeTimeSuffix,FirstPoint=FirstPoint,LastPoint=LastPoint) // awkward nomenclature, lhs are the names of the variables in the called procedure, 
																						 //	rhs is the local variables (the values) that are handed over
		ENDIF											

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
		IF (ParamIsDefault(FirstPoint))
			FillBST_MultiThread(BSTResult, CountWave, IdxMatrix, STReferences, IReferences, DoNoiseFloor, earliestST, latestST)
			//FillBST(BSTResult, CountWave, IdxMatrix, DoNoiseFloor, earliestST, latestST)
		ELSE
			FillBST_MultiThread(BSTResult, CountWave, IdxMatrix, STReferences, IReferences, DoNoiseFloor, earliestST, latestST, FirstPoint=FirstPoint,LastPoint=LastPoint) 
			//FillBST(BSTResult, CountWave, IdxMatrix, DoNoiseFloor, earliestST, latestST ,FirstPoint=FirstPoint,LastPoint=LastPoint) 
																						 // awkward nomenclature, lhs are the names of the variables in the called procedure, 
																						 //	rhs is the local variables (the values) that are handed over
		ENDIF											

		
		// collect gains from all those STAs, before some preparations
		
		VARIABLE		nFreqs=DimSize(Freqs,0)
		MAKE/O/D/N=(nFreqs,nRounds) BST_GainMtx
		WAVE	Gains= BST_GainMtx		// will hold the gains
		
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
		MAKE/O/WAVE/N=(nRounds) W_holder
		MUltiThread/NT = (ThreadProcessorCount -2) W_Holder[]=CreateFilteredGain(BSTResult,p, Local_AC_avg_scaled_splt_FFT,DoNoiseFloor,globalRate,FreqPoints=Freqs)
		
		
		FOR (rr=0; rr< nRounds; rr+=1)
			WAVE	W_G_MgFlt=W_Holder[rr]
			Gains[][rr]=W_G_MgFlt[p]
		ENDFOR
		KillWaves/Z W_D,W_G, W_D_splt, W_D_splt_FFT,  W_G_MgFlt, W_G_PhFlt
		STRING	NameOfResult
		// sort gains
		Variable upperLimit, lowerLimit
		IF (DoNoiseFloor >0)
			NameOfResult="Noise_Floor"
			upperLimit=0.90
			lowerLimit=0.95
		ELSE
			NameOfResult="Conf_Int"
			upperlimit=0.975
			lowerLimit=0.025
		ENDIF
		Duplicate/D/O Gains, $NameOfResult 
 		WAVE Conf_Int=$NameOfResult
		FOR (rr=0; rr< nFreqs; rr+=1)
			Duplicate/O/R=[rr][0,nRounds-1] Conf_Int Dum
			Redimension/N=(nRounds,0) Dum
			Sort Dum, Dum
			Conf_Int[rr][0,nRounds-1]= Dum[q]
		ENDFOR
		KillWaves W_Holder
		KillWaves/Z Dum, STReferences, IReferences
		
		// get the avg of the bootstrap distribution
		MatrixOP/O Boot_avg=averagecols(Conf_Int^t)^t
		
		// sort out the 5 and 95% and the median
		Conf_Int[][0]=Conf_Int[p][nRounds*lowerLimit]
		Conf_Int[][2]=Conf_Int[p][(nRounds-1)*upperLimit]
		Conf_Int[][1]=(Conf_Int[p][floor(nRounds/2)]+Conf_Int[p][ceil(nRounds/2-1)])/2
		Redimension/N=(-1,3) Conf_Int
		Note/K Conf_Int, Note_String
END


FUNCTION PrepareST_List(ST_Suff,[FirstPoint, LastPoint])
// this creates two waves to hold 
// BST_SpikeIndicies 	- of all spike time waves here are the indicies of usable spikes listed
// BST_SpikeCount  		- of all spike time waves this holds the wave names (in DimLabel) and it holds
//							  the number of usable spike times (the once too close to begin / end cannot be used

// this has to be called ONCE before ANY Bootstrap
STRING			ST_Suff					// suffix identifying spike times, without leading "_"
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp

	IF (ParamIsDefault(FirstPoint))
		FirstPoint = 0
	ENDIF
	IF (ParamIsDefault(LastPoint))
		LastPoint = inf
	ENDIF
	
	
	VARIABLE	range=1		// duration of STA (in seconds)
	
	STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
					cellFolderList = ReplaceString(",", cellFolderList, ";")
						cellFolderList= RemoveFromList("Packages;IGNORE;", cellFolderList, ";")

	VARIABLE		ii, nSubs=ItemsInList(cellFolderList,";")
	VARIABLE		kk, nSTs
	DFREF			rootFolderRf=GetDataFolderDFR()
	STRING		currentDFPath, STList, currIName, currSTName
	VARIABLE		minST, maxST, count
	currentDFPath = GetDataFolder(1)		
	
	
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
			
			// get DimOffset right, in case FirstPoint was > 0
			IF (!ParamIsDefault(FirstPoint))
				Duplicate/R=[FirstPoint, min(DimSize(CurrI,0)-1,LastPoint)]/O CurrI, W_CurDummy
				WAVE	currI=W_CurDummy
			ENDIF

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
			SetDimLabel 0,DimSize(BST_SpikeCount,0)-1,$(GetDataFolder(0)+":"+currSTName) BST_SpikeCount
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
			
FUNCTION	BootStrapSweeps(ResultMatrix, CountWave,IdxMatrix, currentBSTround, DoNoiseFloor, earliestST,latestST, [FirstPoint, LastPoint])
WAVE			ResultMatrix, CountWave, IdxMatrix
VARIABLE		currentBSTround			// the current bootstrap round
VARIABLE		DoNoiseFloor				// if > 0, the spike times get shifted and the outcome represents the noise floor rather than the confidence interval
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp
VARIABLE		earliestST,latestST		// boundaries of possible spike times; needed in order to properly shift spike times 
												// if not provided they will be approximated by the time of the first and last spike

// this goes through all the spike time waves and re-creates copies using the indicies in the 
// idxMatrix. These waves are only needed to create boostrap STAs so they
// are just temporary spike time waves and get overwritten and finally deleted

// also the corresponding STAs are created calling STA from analogue. 
// this will be run in a mode such that the out-name does not end in "_STA" but in "_STB"
	IF (ParamIsDefault(FirstPoint))
		FirstPoint = 0
	ENDIF
	IF (ParamIsDefault(LastPoint))
		LastPoint = inf
	ENDIF


	DFREF			rootFolderRf=GetDataFolderDFR()

	VARIABLE		nWaves=DimSize(CountWave,0), ww, length, previous
	STRING		currSTName, currIName, TempSTName, currFolderName
	previous=0		// number of entries in IndexWave previous to current ST wave
	VARIABLE		range	// will hold the range of spike values allowed
	
	FOR (ww=0; ww< nWaves; ww+=1)
		SetDataFolder rootFolderRf
		length=CountWave[ww]
		currSTName=GetDimLabel(CountWave, 0, ww )
		currFolderName=StringFromList(0,currSTName,":")
		IF (ItemsInList(currSTName,":")>1)		// if there are ANY subfolders
			SetDataFolder $(currFolderName)
		ENDIF
		currSTName=StringFromList(1,currSTName,":")

		WAVE	currST = $currSTName
		currIName=RemoveEnding(currSTName,"ST") + "I" 
		WAVE	currI=$currIName
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
		IF (ParamIsDefault(FirstPoint))
			STAfromAnalogue(currI, ST, 1,0,Suffix="B"+num2istr(currentBSTround)+"_STB")
		ELSE
			STAfromAnalogue(currI, ST, 1,0,Suffix="B"+num2istr(currentBSTround)+"_STB", FirstPoint=FirstPoint, LastPoint=LastPoint)
		ENDIF
		previous+=length	
		KillWaves/Z ST

	ENDFOR
		
	SetDataFolder rootFolderRf
	
	// make average STA from the bootstrap STAs
	AvgSTA(Suffix="B"+num2istr(currentBSTround)+"_STB")
	STRING resultName
	resultName="B"+num2istr(currentBSTround)+"_STB_avg"
	WAVE res=$resultName
	// not needed, only scaled version is relevant
	KillWaves/Z res
	// now get the scaled version
	resultName="B"+num2istr(currentBSTround)+"_STB_avg_scaled"
	WAVE res=$resultName
	ResultMatrix[][currentBSTround] = res[p]
	// now take over the info on avg sd , spikes and duration
	IF (currentBSTround==0)
		STRING 	Notiz=Note(res)
		Note/K ResultMatrix
		Note/NOCR ResultMatrix, Notiz
	ENDIF
	
	KillWaves/Z res
	
	// clean up (remove all STB waves from this round)
	FOR (ww=0; ww< nWaves; ww+=1)
		SetDataFolder rootFolderRf
		currSTName=GetDimLabel(CountWave, 0, ww )
		currFolderName=StringFromList(0,currSTName,":")
		IF (ItemsInList(currSTName,":")>1)		// if there are ANY subfolders
			SetDataFolder $(currFolderName)
		ENDIF

		//SetDataFolder $(currFolderName)
		currSTName=StringFromList(1,currSTName,":")
		currIName=RemoveEnding(currSTName,"ST") +"B"+  num2istr(currentBSTround)+"_STB"
		WAVE	currSTB= $currIName
		KillWaves/Z currSTB

	ENDFOR
	SetDataFolder rootFolderRf

END

Function FillBST(BSTResult, CountWave, IdxMatrix, DoNoiseFloor,  earliestST, latestST, [FirstPoint, LastPoint])
	WAVE 		BSTResult, CountWave, IdxMatrix
	VARIABLE	FirstPoint, LastPoint,DoNoiseFloor, earliestST, latestST
	
	// if DoNoiseFloor> 0, the spike times get shifted and the outcome represents the noise floor rather than the confidence interval
	// earliestST,latestST are boundaries of possible spike times; needed in order to properly shift spike times, i.e. only needed for NoiseFloor 
	// if not provided they will be approximated by the time of the first and last spike

	// this function will put the spike triggered average of all spike used durig a given BST round into the BSTResult
	Variable ncol= DimSize(BSTResult,1)
	Variable col
	Variable time2, ttime= stopMSTimer(-2)
	
	for(col=0; col<ncol; col+=1)
		IF (ParamIsDefault(FirstPoint))
			BootStrapSweeps(BSTResult, CountWave,IdxMatrix, col, DoNoiseFloor,earliestST, latestST)
		ELSE
			BootStrapSweeps(BSTResult, CountWave,IdxMatrix, col, DoNoiseFloor, earliestST, latestST, FirstPoint=FirstPoint,LastPoint=LastPoint) // awkward nomenclature, lhs are the names of the variables in the called procedure, 
																						 //	rhs is the local variables (the values) that are handed over
		ENDIF											
		time2= stopMSTimer(-2)
		printf " round %d; : %d sec;\r",col, (time2-ttime)*1e-6
		
	endfor


End

// - - - - -  MultiThreading versions - - - - - - - -

Threadsafe	FUNCTION/DF	Bootstrap_MultiThread(CountWave,IdxMatrix, currentBSTround, STWaveReferences, IWaveReferences, DoNoiseFloor,earliestST, latestST, [FirstPoint, LastPoint])
WAVE			CountWave, IdxMatrix
WAVE/WAVE	STWaveReferences, IWaveReferences
VARIABLE		currentBSTround			// the current bootstrap round
VARIABLE		DoNoiseFloor,earliestST, latestST // if DoNoiseFloor the spike times are shifted (by more than one orrelation time)
															 // earliest and latest are used to 
VARIABLE		FirstPoint, LastPoint	// to deal with the data from Omer, where only some part of the current can be used 
												// due to the limitations of pClamp

// this goes through all the spike time waves and re-creates copies using the indicies in the 
// idxMatrix. These waves are only needed to create boostrap STAs so they
// are just temporary spike time waves and get overwritten and finally deleted

// also the corresponding STAs are created calling STA from analogue. 
// this will be run in a mode such that the out-name does not end in "_STA" but in "_STB"
	IF (ParamIsDefault(FirstPoint))
		FirstPoint = 0
	ENDIF
	IF (ParamIsDefault(LastPoint))
		LastPoint = inf
	ENDIF
	
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
		STB_References[ww]=STAfromAnalogueWaveRef(currI, ST, 1,0,Suffix="B"+num2istr(ww)+num2istr(currentBSTround)+"_STB", FirstPoint=FirstPoint, LastPoint=LastPoint)

		previous+=length	
		KillWaves/Z ST, ST2

	ENDFOR
			
	// make average STA from the bootstrap STAs
	
	WAVE res=AvgSTA_ThreadSafe(STB_References,targetName="avgSTB") // no specific name needed, becuase we work in a 
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



Function FillBST_MultiThread(BSTResult, CountWave, IdxMatrix, STReferences, IReferences,  DoNoiseFloor,  earliestST, latestST, [FirstPoint, LastPoint])
	WAVE 		BSTResult, CountWave, IdxMatrix
	WAVE/WAVE	STReferences, IReferences
	VARIABLE	FirstPoint, LastPoint, DoNoiseFloor,  earliestST, latestST

	
	Variable ncol= DimSize(BSTResult,1)
	Variable col
	Variable time2, ttime= stopMSTimer(-2)
	// Create a wave to hold data folder references returned by Worker.
	// /WAVE specifies the data type of the wave as "wave reference".
	Make/O/DF/N=(ncol) ddd
	VARIABLE numCores=ThreadProcessorCount
	// one processor will not be engaged to keep the computer responsive
	IF (ParamIsDefault(FirstPoint))
		MultiThread/NT=(numCores-5)	ddd= Bootstrap_MultiThread(CountWave,IdxMatrix,p, STReferences, IReferences, DoNoiseFloor,  earliestST, latestST)
	ELSE
		MultiThread/NT=(numCores-5) ddd= Bootstrap_MultiThread(CountWave,IdxMatrix,p, STReferences, IReferences, DoNoiseFloor,  earliestST, latestST,FirstPoint=FirstPoint,LastPoint=LastPoint) 
																		// awkward nomenclature, lhs are the names of the variables in the called procedure, 
																		//	rhs is the local variables (the values) that are handed over
	ENDIF											
	

	
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

Threadsafe FUNCTION/WAVE	CreateFilteredGain(M_STA,clmn, FFTofAC,DoNoiseFloor,globalRate,[FreqPoints, widthFac])
WAVE		M_STA 		// matrix full of bootstrapped spike triggered averages 
VARIABLE	clmn		// the clmn of interest in the present instance
WAVE/C		FFTofAC	// fourier transform of autocorrelation, which was split at middle
WAVE/D		FreqPoints	// vector with Frequency values at which filtered versions of 
						// FT will be constructed
VARIABLE	DoNoiseFloor, globalRate
VARIABLE	widthFac		// width of the gaussian filter	; a value of 1 corresponds to 
						// Higgs and Spain's choice of f/2Pi
		// the optional parameters FreqPoints and widthFac  are supplöied as follows:
		// GaussFilter(FT,FreqPoints=FrequencyRange, widthFac=1)
		
			Duplicate/O/R=[*][clmn] M_STA STA

		
			DFREF dfSav= GetDataFolderDFR()
	
			// Create a free data folder and set it as the current data folder
			SetDataFolder NewFreeDataFolder()

			//Duplicate/O/R=[*][rr] BSTResult W_D
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
		
			STRING	OutName = "STA_splt"
			Duplicate/O/R=[0,2*floor(DimSize(STA,0)/2)-1] STA, $OutName
			WAVE out=$OutName
			out[0,DimSize(out,0)-1-SplitP] 						= STA[p+SplitP]	
			out[DimSize(out,0)-SplitP,DimSize(out,0)-1]		= STA[p+SplitP-DimSize(out,0)]	
			
			WAVE	STA_splt
			FFT/OUT=1 /DEST=W_G STA_splt
			KillWaves STA_splt
			W_G/=FFTofAC
			W_G*=cmplx(globalRate,0)			
			W_G=conj(W_G)			// added 03/08/2020 to obtain correct sign of the phase
	
			//  Will  produce a Magnite and Phase wave
			VARIABLE	df_in=DimDelta(FFTofAC,0)		// frequency step width in input
			VARIABLE	minf=DimOffset(FFTofAC,0)
			VARIABLE	maxf=minf+ (DimSize(FFTofAC,0) -1)*df_in
			VARIABLE	nPoints	
			IF (ParamIsDefault(FreqPoints))
			          	nPoints = min(floor((DimSize(FFTofAC,0) -1)/2), 50)		// max 50 Freq points, but no more than 
																		// half the number that is there already
				// linear distance in log requires a frequency factor 
				// to describe construction of the FreqList
				VARIABLE	freqFac =10^( log(maxf-minf)/(nPoints-1))
				Make/D/O/N=(nPoints) FreqPoints
				Wave FreqPoints
				FreqPoints[]=minf+freqFac^p
				
			ELSE
				 nPoints= DimSize(FreqPoints,0)
				
			ENDIF
			IF (ParamIsDefault(widthFac))
				widthFac = 1
			ENDIF
	     	VARIABLE k,  center, width, weightSum
	     	DUPLICATE/C/D/O FFTofAC , GWeight
			Make/D/O/N=(Dimsize(GWeight,0)) RealPart, ImagPart				// real valued wave to hold magnitude of weighted complex 
			DUPLICATE/O/D FreqPoints, Mag_out
	         

	      FOR (k=0; k< nPoints; k+=1)
				width=FreqPoints[k]/2/Pi	* widthFac		// SD of Gaussian is 1/freq
				GWeight[]=cmplx( gauss(pnt2x(FFTofAC, p ),FreqPoints[k], width),0)
				WeightSum = mean(GWeight)*Dimsize(GWeight,0)
				GWeight/=WeightSum			// making sure the integral of the weight is one
				GWeight*=W_G	
				RealPart[]=real(GWeight[p])				
				ImagPart[]=imag(GWeight[p])				
				Mag_out[k] =sqrt( sum(RealPart) ^2+  sum(ImagPart)^2)
			ENDFOR
		KillWaves/Z      GWeight , RealPart, ImagPart

		KillWaves/Z GWeight, RealPart, ImagPart, W_G
		SetDataFolder dfSav
		Return   Mag_out
END




// End of Bootstrap code
 	
		 





//¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
//-------------    F I L T E R  I N   C O M P L E X   S P A C E       --------------------
//________________________________________________________________________________________
FUNCTION	GaussFilter(FT,[FreqPoints, widthFac])
WAVE/C		FT 		// complex wave, fourier transform of wave
WAVE/D		FreqPoints	// vector with Frequency values at which filtered versions of 
						// FT will be constructed
VARIABLE	widthFac		// width of the gaussian filter	; a value of 1 corresponds to 
						// Higgs and Spain's choice of f/2Pi
		// the optional parameters FreqPoints and widthFac  are supplied as follows:
		// GaussFilter(FT, widthFac=1)

// check whether Input is actually a complex wave

	VARIABLE WavNumType=NumberByKey("NUMTYPE", WaveInfo(FT,0) ) 
	IF (mod(WavNumType,2)==0) 	// NOT complex
		DoAlert 0,"A complex wave should be supplied."
		RETURN -1
	ENDIF
	
//  Will  produce a Magnite and Phase wave
	VARIABLE	df_in=DimDelta(FT,0)		// frequency step width in input
	VARIABLE	minf=DimOffset(FT,0)
	VARIABLE	maxf=minf+ (DimSize(FT,0) -1)*df_in
	VARIABLE	nPoints	
		IF (ParamIsDefault(FreqPoints))
	           	nPoints = min(floor((DimSize(FT,0) -1)/2), 50)		// max 50 Freq points, but no more than 
																	// half the number that is there already
			// linear distance in log requires a frequency factor 
			// to describe construction of the FreqList
			VARIABLE	freqFac =10^( log(maxf-minf)/(nPoints-1))
			Make/D/O/N=(nPoints) FreqPoints
			Wave FreqPoints
			FreqPoints[]=minf+freqFac^p
			
		ELSE
			 nPoints= DimSize(FreqPoints,0)
			
		ENDIF
		IF (ParamIsDefault(widthFac))
			widthFac = 1
		ENDIF
	         VARIABLE k,  center, width, weightSum
	         DUPLICATE/C/D/O FT , GWeight
	         Make/D/O/N=(Dimsize(GWeight,0)) RealPart, ImagPart				// real valued wave to hold magnitude of weighted complex FT
	         DUPLICATE/O/D FreqPoints, Mag_out, Phase_out
	         

	         FOR (k=0; k< nPoints; k+=1)
// 	        IF ( FreqPoints[k] == 0 )			// catch the case when Frequency is zero and 
// 	        								// width would go to infinity
//	         	Mag_out = FT[x2pnt(FT,0)]
//	         	k+=1
//	         ENDIF		
		        width=FreqPoints[k]/2/Pi	* widthFac		// SD of Gaussian is 1/freq
			GWeight[]=cmplx( gauss(pnt2x(FT, p ),FreqPoints[k], width),0)
			WeightSum = mean(GWeight)*Dimsize(GWeight,0)
			GWeight/=WeightSum			// making sure the integral of the weight is one
			GWeight*=FT	
			RealPart[]=real(GWeight[p])				
			ImagPart[]=imag(GWeight[p])				
			Mag_out[k] =sqrt( sum(RealPart) ^2+  sum(ImagPart)^2)
			Phase_out[k]=imag(  r2polar(   cmplx(   sum(RealPart) , sum(ImagPart) )    )      )
	         ENDFOR
		KillWaves/Z      GWeight , RealPart, ImagPart
		STRING	OutName1 = NameOfWave(FT),OutName2 = NameOfWave(FT)
		OutName1=OutName1[0,24]+"_MgFlt"
		OutName2=OutName2[0,24]+"_PhFlt"
		IF ( (Exists(OutName1) ==1 ) || (Exists(OutName2) ==1 ) )	// outwaves already exist
//		DoAlert 1,"Overwrite existing filtered output waves?"
		VARIABLE V_overwrite=1
				IF 	(V_overwrite==1)		// DO overwrite	
					Duplicate/O   Mag_out, $OutName1
					Duplicate/O Phase_out, $OutName2				
				ENDIF
				KillWaves/Z Mag_out, Phase_out
		ELSE
				Rename Mag_out, $OutName1
				Rename Phase_out, $OutName2
		ENDIF
		KillWaves/Z GWeight, RealPart, ImagPart
END


FUNCTION PiecewiseFFT(InWave, winLen,sinSqrLen, overlap)
WAVE InWave
VARIABLE	winLen, SinSqrLen, Overlap

// creates pieces of length winLen from input 
// windows them with sine square onset (length sinsqrlen from zero to 1)
// calculates FFT and accumulates the FFT to calculate average FFT for the entire input
// lengths are in points

// overlap is the fraction of the entire window that is overlapping the next window

	VARIABLE		Ntot=DimSize(InWave,0)
	// attempts coverage for entire wave: 
	// suggest a window more or a window less 
	VARIABLE	stepForward = round((1-Overlap) * winLen)
	VARIABLE coverage= (Ntot-winlen)/ stepForward
	VARIABLE nFullWin=1+floor(coverage)
	coverage = coverage-floor(coverage) // fraction of a full window advance that does not fit in anymore 
	coverage = (1-coverage)*(1-Overlap) * winLen // length not covered by last full window
	coverage = coverage / Ntot // fraction not covered by last full window
	
	// i ftoo much is not covered try adjusting
	IF (coverage > 0.01) // more than 1 % not covered
		// check if enough full points are left over to add at least two points to every window
		IF ((coverage*Ntot) > 2*nFullWin)
			DoAlert 0, "Window length is increased from "+ num2istr(winlen) + " by "+num2istr(floor((coverage*Ntot) / nFullWin))+" points."
			winLen = winLen+ floor((coverage*Ntot) / nFullWin)
		ELSE
			DoAlert 1,"Coverage is only "+num2str(1 - coverage)+", proceed?"
			IF (V_Flag !=1)
				Return -1
			ENDIF
			
		ENDIF
	ENDIF
	
	VARIABLE		NWin
	MAKE/O/D/N=(winLen) Fenster
	WAVE	Fenster
	
	Fenster=1
	Fenster[0,SinSqrLen]=(sin(p/sinsqrLen/4*2*Pi))^2
	Fenster[winLen-1-SinSqrLen, winLen-1] = Fenster[SinSqrLen-(p-(winLen-1-SinSqrLen))]
	
	FOR (NWin=0; NWin<NFullWin; NWin+=1)
		Duplicate/O/D/R=[0+NWin*stepForward,winLen-1+NWin*stepForward] InWave, Dummy
		Wave Dummy
		// remove avg
		WAVESTATS/Q/M=1 Dummy
		Dummy-=V_avg
		Dummy*=Fenster
		FFT/OUT=1/DEST=Dummy_FFT  Dummy
		WAVE/C Dummy_FFT

		IF (NWin==0)
			Duplicate/O/D/C Dummy_FFT OUTWave
			WAVE/C OUTWave
		ELSE
			OUTWave+= Dummy_FFT
		ENDIF
	
	ENDFOR
	OUTWave/=NFullWin*winlen
	
	KillWaves/Z Fenster
END

FUNCTION PiecewiseGain(InWave, OutWave, OutName, winLen,sinSqrLen, overlap)
WAVE InWave, OUtwave
STRING OutName
VARIABLE	winLen, SinSqrLen, Overlap

// creates pieces of length winLen from input 
// windows them with sine square onset (length sinsqrlen from zero to 1)
// calculates FFT and accumulates the FFT to calculate average FFT for the entire input
// lengths are in points

// overlap is the fraction of the entire window that is overlapping the next window

	VARIABLE		Ntot=DimSize(InWave,0)
	// attempts coverage for entire wave: 
	// suggest a window more or a window less 
	VARIABLE	stepForward = round((1-Overlap) * winLen)
	VARIABLE coverage= (Ntot-winlen)/ stepForward
	VARIABLE nFullWin=1+floor(coverage)
	coverage = coverage-floor(coverage) // fraction of a full window advance that does not fit in anymore 
	coverage = (1-coverage)*(1-Overlap) * winLen // length not covered by last full window
	coverage = coverage / Ntot // fraction not covered by last full window
	
	// i ftoo much is not covered try adjusting
	IF (coverage > 0.01) // more than 1 % not covered
		// check if enough full points are left over to add at least two points to every window
		IF ((coverage*Ntot) > 2*nFullWin)
			DoAlert 0, "Window length is increased from "+ num2istr(winlen) + " by "+num2istr(floor((coverage*Ntot) / nFullWin/2)*2)+" points."
			winLen = winLen+ floor((coverage*Ntot) / nFullWin)
		ELSE
			DoAlert 1,"Coverage is only "+num2str(1 - coverage)+", proceed?"
			IF (V_Flag !=1)
				Return -1
			ENDIF
			
		ENDIF
	ENDIF
	
	VARIABLE		NWin
	MAKE/O/D/N=(winLen) Fenster
	WAVE	Fenster
	
	Fenster=1
	Fenster[0,SinSqrLen]=(sin(p/sinsqrLen/4*2*Pi))^2
	Fenster[winLen-1-SinSqrLen, winLen-1] = Fenster[SinSqrLen-(p-(winLen-1-SinSqrLen))]
	
	FOR (NWin=0; NWin<NFullWin; NWin+=1)
		Duplicate/O/D/R=[0+NWin*stepForward,winLen-1+NWin*stepForward] InWave, Dummy
		Wave Dummy
		// remove avg
		WAVESTATS/Q/M=1 Dummy
		Dummy-=V_avg
		Dummy*=Fenster
		FFT/OUT=1/DEST=Dummy_FFT  Dummy
		Duplicate/O/D/R=[0+NWin*stepForward,winLen-1+NWin*stepForward] Outwave, Dummy
		Wave Dummy
		// remove avg
		WAVESTATS/Q/M=1 Dummy
		Dummy-=V_avg
		Dummy*=Fenster
		FFT/OUT=1/DEST=Dummy_FFT2  Dummy
		
		WAVE/C Dummy_FFT2

		IF (NWin==0)
			Duplicate/O/D/C Dummy_FFT2 $OutName
			WAVE/C GainWAve=$OutName
			GainWave/=Dummy_FFT
		ELSE
			GainWAve+= Dummy_FFT2/Dummy_FFT
		ENDIF
	
	ENDFOR
	GainWAve/=NFullWin*winlen
	
	KillWaves/Z Fenster, Dummy, Dummy_FFT, Dummy_FFT2
END



FUNCTION/WAVE SplitBeforeFFT(InWave, SplitTime)
WAVE		InWave
VARIABLE	SplitTime
// Timepoint according to x-index
// tjis requires to wave to have appropriate x-scaling

// splits a function and puts it together such that the first point of the output corresponds
// to the point closest to the SplitTime in the input

// as the next step is application of a fast fourier transformation, the number
// of points in the output should be an even number
// if the input has an even number of points - everything is fine
// if the input has an odd number of points, one point has to be deleted. This is approached with the 
// assumption that the last and first points of the input are not neighbours (in a physical sense
// this assumption is in some cases wrong.
// starting from this assumption, it makes sense to delete one of the points (the last)

// check whether Split time is inside the x-index
	IF ( SplitTime < DimOffset(InWave,0)  || SplitTime > pnt2x(InWave,numpnts(InWave)-1) )
		Abort "The intended split time is outside the x-Range of in-wave."
	ENDIF
	VARIABLE	dT=DimDelta(Inwave,0)

	VARIABLE SplitP=x2pnt(InWave,SplitTime+dT/1000)

// split such that the entry closest to split time ends up at first point

// only using 'SplitTime' instead of 'SplitTime+dT/1000' in the case when 
// the split time lies exactly between two samples, would give a counter-intuitive
// result: the point b e f o r e  the split time would sometimes be the 
// first point after splitting
// shifting the SplitTime very slightly to the right this is avoided

	STRING	OutName = NameOfWave(InWave)
	OutName=OutName[0,25]+"_splt"
	Duplicate/O/R=[0,2*floor(DimSize(InWave,0)/2)-1] InWave, $OutName
	WAVE out=$OutName
	out[0,DimSize(out,0)-1-SplitP] 						= InWave[p+SplitP]	
	out[DimSize(out,0)-SplitP,DimSize(out,0)-1]		= InWave[p+SplitP-DimSize(out,0)]	
	// last line corrected 06/03/2016 before it read [p-splitP], which for the standard case of splitting 
	// at (nPnts-1)/2 gave only a 1 pnt offset. By some kind of weird coincidence this was exactly what is 
	// needed: it caused a duplication of the first point (the first point after the split time)
	// as the last point. That sounds bad until one remembers that the total number of points is - in my convention
	// an ODD number (central point plus a certain number of points left and right)
	// so we always discard a single point when the FFT is done. This is exactly the point that was duplicated
	// by accident. Perfect.
	// however, that is now replaced. (see intro to the function)
	Return out
END


FUNCTION CyclicSplit(InWave, SplitTime)
WAVE		InWave
VARIABLE	SplitTime		// Timepoint according to x-index

// splits a vector and puts it together resulting in a cyclic shift 

// check whether Split time is inside the x-index
	IF ( SplitTime < DimOffset(InWave,0)  || SplitTime > pnt2x(InWave,numpnts(InWave)-1) )
		Abort "The intended split time is outside the x-Range of in-wave."
		Return -1
	ENDIF
	VARIABLE SplitP=x2pnt(InWave,SplitTime)
// split such that the entry closest to split time end up at first point
	STRING	OutName = NameOfWave(InWave)
	OutName=GetWavesDataFolder(InWave,1)+OutName[0,25]+"_splt"
	
	Duplicate/O InWave, $OutName
	WAVE out=$OutName
	out[0,DimSize(InWave,0)-1-SplitP] 							= InWave[p+SplitP]	
	out[DimSize(InWave,0)-SplitP,DimSize(InWave,0)-1]		= InWave[p+SplitP-DimSize(InWave,0)]	
	
END

FUNCTION  STA2Gain(STAwave, t_corr_OU, sigma_OU, PeakOrZero, Frequencies)
WAVE		STAwave			// input STA
VARIABLE	t_corr_OU			// correlation time of the input
VARIABLE	sigma_OU			// standard deviation of the input
								// last two needed to calculate the power spectral density of the input
WAVE		Frequencies			// Frequencies (in Hz) at which the gauss smoothed response (mag and phase)
								// will be calculated 
													
STRING		PeakOrZero			// determines which time point in the STA should be treated as zero time
								// either "Zero" from time stamp or "peak" - time of max 
		STRING 			OutPhaseName, OutMagName
		OutPhaseName =	NameOfWave(STAwave)+"_Phs"
		OutMagName	    =	NameOfWave(STAwave)+"_Gain"
		
		VARIABLE	SplitTime
	// decide where to split (at zero time or at peak time	
	strswitch(PeakOrZero)	// string switch
		case "Peak":
		case "peak":
		case "PEAK":
			WaveStats/Q STAwave
			SplitTime=V_maxloc
			break						
		case "Zero":
		case "zero":
		case "ZERO":
			SplitTime=0
			break
		default:							
			DoAlert 0,"Please use 'zero' or 'peak' to determine where zero time is in the STA."
			RETURN -1						
	endswitch
	// if the STA has an uneven number of points, the last one is left out
	IF (mod(DimSize(STAwave,0),2)==1)
		Duplicate/O/R=[0,DimSize(STAwave,0)-2] STAwave, tempjunk
		WAVE	STAwave = tempjunk
	ENDIF
	SplitBeforeFFT(STAwave, SplitTime)		//OutName = NameOfWave(InWave)
										//OutName=OutName[0,25]+"_splt"
	STRING	ExpectedName=NameOfWave(STAwave)[0,25]+"_splt"
	IF (Exists(ExpectedName)!=1)		// wave does not exist
		DoAlert 0,"After splitting the STA the result '"+ExpectedName+"' could not be located"	
		RETURN -2
	ENDIF						
	WAVE	STA=$ExpectedName							
	// Fast Fourier Transform with NO Window, complex result
	
	FFT/OUT=1/DEST=STA_FFT STA
	GaussFilter(STA_FFT,FreqPoints=Frequencies, widthFac=1)
	// returns InName[0,24]+"_MagFilt" and  InName[0,24]+"_PhsFilt"
	ExpectedName="STA_FFT_MagFilt"
	WAVE	Out=$ExpectedName
	Rename Out,$OutMagName
	Out[]/=4*t_corr_OU*sigma_OU^2/ ( (2*Pi*Frequencies[p]*t_corr_OU)^2 +1)
	ExpectedName="STA_FFT_PhsFilt"
	WAVE	Out=$ExpectedName
	Rename Out,$OutPhaseName
	
	KillWaves/Z tempjunk, STA, STA_FFT
END



//---------------------------------------------------------------------
//=============   automation   ==============
//---------------------------------------------------------------------

FUNCTION Phaseplot(VoltWave)
// Creates a phaseplot for the entire recording
// color code signals time since last spike
// This function creates the following waves
// Differentiated voltage
// colorwave

WAVE		VoltWave		// contains membrane potential over time

// 1. get spike times
	ReturnSpikeTimes(VoltWave)
	// creates a wave that contains spike times and has a name equal to the 
	// name of the input with "_ST" added
	// if VoltWaves name ended on "_V" remove this part from the spike time name
	STRING		STName=NameOfWave(VoltWave)
	STRING		VoltName=STName
	VARIABLE	VoltNameLength=strlen(VoltName)
	
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																// in VoltName are "_V"
		STName=VoltName+"_ST"														
		WAVE ST=$STName
		Rename ST, $(STName[0,VoltNameLength-3]+"_ST")
		STName = NameOfWave(ST)
	ENDIF
	// even if the spike time wave was not renamed, it is know inside this  function as ST

//2. create voltage derivative
	STRING		DifName
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																// in VoltName are "_V"
		DifName=	VoltName[0,VoltNameLength-3]+"_dV"
	ELSE
		DifName= 	VoltName+"_dV"
	ENDIF
	WAVE	volt=$VoltName													
	Differentiate volt/D=$DifName
	WAVE	Dif=$DifName

//3. create color wave
	// make a wave that has as many points as the voltage trace
	STRING	CName
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																// in VoltName are "_V"
		CName=	VoltName[0,VoltNameLength-3]+"_Color"
	ELSE
		CName= 	VoltName+"_Color"
	ENDIF
	Duplicate/O volt,$CName
	WAVE	C1=$CName
	
	// it is 1 until first spike
	C1=1
	// for all spike times: reset C1 to 0 and then point by point count up the time since the last spike
	// until the next spike time
	VARIABLE	SCount, SNum=DimSize(ST,0)	// SCount is counter in loop
											// SNum is number of spikes 
	VARIABLE	SpikeTime1, SpikePoint1		// time at which a spike occurs
	VARIABLE	SpikeTime2, SpikePoint2		// and index of the time point, i.e. 
											// which entry in the voltage wave corresponds to the time
											
		SpikeTime1=ST[0]						// time of first spike
		SpikePoint1=x2pnt(volt, SpikeTime1 )		// based on the sampling interval of the voltage wave
											// that is known to IGOR, it is calculated
											// which point in the voltage wave is closest to this spike 	
	FOR (SCount=0; Scount < SNum-1; SCount+=1) // from first to second from last spike

		SpikeTime2=ST[Scount+1]				// time of next spike
		SpikePoint2=x2pnt(volt, SpikeTime2 )		
		C1[SpikePoint1,SpikePoint2]=(p-SpikePoint1)*DimDelta(volt,0)
											// count up time point by point in units of the sample time
											// of the voltage trace
		SpikePoint1=SpikePoint2				// prepare next round by handing over SpikePoints
	ENDFOR
	// now update until end of voltage trace
	C1[SpikePoint1,DimSize(C1,0)-1]=(p-SpikePoint1)*DimDelta(volt,0)

	// create a new color wave that lag 2 ms with respect to the current color wave
	STRING	C2Name=CName+"_plus2ms"
	Duplicate/O C1,$C2Name
	WAVE C2=$C2Name		// now inside the function the second color wave is called C2
	
	VARIABLE	pointOffset=round(0.002/DimDelta(volt,0))	// calculate how many points correspond to 2 milliseconds
	C2=1
	C2[pointOffset,DimSize(C2,0)-1]=C1[p-pointOffset]
	
	// now create graph
	Display /W=(34.8,42.2,559.2,354.8) Dif vs volt
	ModifyGraph margin(right)=99
	ModifyGraph rgb=(0,0,0)
	ModifyGraph zColor($DifName)={C2,*,0.1,YellowHot,1}
	ModifyGraph minor(bottom)=1
	Label left "dV\\BMem\\M/dt (\\U)"
	Label bottom "V\\Bmem\\M (\\U)"
	ColorScale/C/N=text0/F=0/M/H=13/A=MC/X=61.33/Y=-4.81 trace=$DifName
	AppendText "\\u#2Time since last spike (ms)"
	TextBox/C/N=text1/X=-36.00/Y=42.00/F=0/M/H=13/A=MC VoltName
	
	// kill the unnecessary C1 color wave (not delayed)
	KillWaves/Z C1	


END	// Phaseplot	

FUNCTION PhasePlot_SpikesOnly(VoltWave[,SpikeTimes, Trange])
WAVE			VoltWave, Spiketimes
VARIABLE		Trange
// Creates a phaseplot, uses only a stretch of data around each spike

// Trange is the time BEFORE AND AFTER the spike that is plotted
// If colouring is wanted, a wave is needed with as many entries as there are spikes
// then the x-th phaseplot can be  cloured according to the x-th value in coloursource
// using the function ColorTraceByWave() or the menu item color traces in Graphs
// both are enabled by loading ANTools.ipf

	IF (ParamIsDefault(Trange))
		Trange=2e-3		// 2ms
	ENDIF
	STRING		VoltName=NameOfWave(VoltWave)

	VARIABLE dt=DimDelta(VoltWave,0), Prange=ceil(Trange/dt)
	
	IF (ParamIsDefault(SpikeTimes))
		ReturnSpikeTimes(VoltWave)
		// creates a wave that contains spike times and has a name equal to the 
		// name of the input with "_ST" added
		// if VoltWaves name ended on "_V" remove this part from the spike time name
		STRING		STName=VoltName
		VARIABLE		VoltNameLength=strlen(VoltName)
		
		IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																	// in VoltName are "_V"
			STName=VoltName[0, VoltNameLength-2]+"ST"														
			WAVE ST=$STName
			STName = NameOfWave(ST)
		ENDIF
	// even if the spike time wave was not renamed, it is know inside this  function as ST
	ELSE
		WAVE	ST=SpikeTimes
	ENDIF
	
	

	// create voltage derivative
	STRING		DifName
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last two characters
																// in VoltName are "_V"
		DifName=	VoltName[0,VoltNameLength-3]+"_dV"
	ELSE
		DifName= 	VoltName+"_dV"
	ENDIF
	Differentiate voltwave/D=$DifName
	WAVE	Dif=$DifName

	
	// create a 3D Wave that holds the necessary pieces of voltage and derivative
	VARIABLE	SCount, SNum=DimSize(ST,0)	// SCount is counter in loop
											// SNum is number of spikes 
	MAKE/O/N=(2*PRange+1,SNum,2)	$(DifName[0,strlen(DifName)-3]+"PP")// for every spike PRange points before and after 
																							// the detection time and 2 signals (V and dV)
	WAVE	PP= 	$(DifName[0,strlen(DifName)-3]+"PP")																				
	VARIABLE	 SpikePoint1	
	
	Display /W=(34.8,42.2,559.2,354.8) as "Phaseplot "+VoltName

	FOR (SCount=0; Scount < SNum-1; SCount+=1) // from first to last spike

		SpikePoint1=x2pnt(VoltWave, ST[Scount] )-PRange/3*2		
		PP[][Scount][0]=voltWave[SpikePoint1+p]	
		PP[][Scount][1]=Dif[SpikePoint1+p]	
		AppendToGraph/L=VertCrossing/B=HorizCrossing PP[*][Scount][1] vs PP[*][Scount][0]
	ENDFOR

		
	// now create graph
	ModifyGraph margin(left)=10
	ModifyGraph rgb=(0,0,0)
	ModifyGraph freePos(HorizCrossing)={0,VertCrossing}
	ModifyGraph freePos(VertCrossing)={0,HorizCrossing}
	ModifyGraph minor(HorizCrossing)=1
	Label VertCrossing "dV\\BMem\\M/dt (\\U)"
	Label HorizCrossing "V\\Bmem\\M (\\U)"
	TextBox/C/N=text1/X=-36.00/Y=42.00/F=0/M/H=13/A=MC VoltName
	
	// kill the unnecessary differential
	KillWaves/Z Dif	
END // PhasePLot_spikesOnly

FUNCTION	InterSpikeInterval(STWave)
WAVE		STWave		// contains spike times
	
	STRING 		ISIName, STName = NameOfWave(STWave)
	VARIABLE	STNameLength=strlen(STName)
	
	IF 	(StringMatch(STName[STNameLength-3,STNameLength-1],"_ST" ))	// last to characters
																// in STName are "_ST"
		ISIName=	STName[0,STNameLength-4]+"_ISI"
	ELSE
		ISIName= 	STName+"_ISI"
	ENDIF

	Duplicate/O STWave,$ISIName
	WAVE	ISI=$ISIName
	ISI[1,DimSize(ISI,0)-1]-=STWave[p-1];		// subtract most recent spike time from current spike time
										// to get inter spike interval
												
	DeletePoints 0, 1, ISI					// remove the first point
										// as the waiting time before first spike is not knon
	// now ISI is calculated, print the coefficient of variation
	
	Wavestats/Q ISI
	printf ISIName+" coeff. of variation\t%f\r",V_sdev/V_avg
END		// InterSpikeInterval

//¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
//--------------------    S P I K E  C O R R E L A T I O N S       ----------------------
//_______________________________________________________________________________________


FUNCTION ConvolveEvents(Source,Sigma)
WAVE		Source
VARIABLE	Sigma

	IF (Sigma > 0)
		VARIABLE	x0,dt = DimDelta(Source,0)
			// let make Gaussian Template 6 Sigma wide
		MAKE/D/O/N=(round(4*sigma/dt)*2+1) Gaussian		// uneven number of points, central point has zero time
		WAVE	Gwave=$"Gaussian"
		x0=-round(4*sigma/dt)*dt
		SetScale/P x,x0, dt, "s", GWave
		GWave[]=1/sqrt(2)/Sigma
		GWave[]*=exp(-1/2*((x0+p*dt)/sigma)^2  )
		VARIABLE integr=sum(GWave)
		GWave/=integr
		Convolve/A	GWave, Source	// using acausal as the left wave is the "impulse response with zero delay equals central point
//		KillWaves/Z GWave
	ENDIF
END		//ConvolveEvents

	
FUNCTION BinEvents(inwav, dt, starttime,endtime)
WAVE		inwav
VARIABLE	dt, starttime, endtime
// takes a sequence of event times (symbolizing occurence of 1 spike)
// downsamples into time bins of width dt
// over the interval from start time to endtime
// if both are ZERO, the range of entries in inwave is used


	IF (exists(NameOfWave(inwav))==1)
		IF (DimSize(Inwav,0)<1)
			DoAlert 0, "No entries in wave"
			Return -2
		ENDIF
	ELSE
		DoAlert 0,"Wave does not exist"
		Return -1
	ENDIF
	
	VARIABLE	point, k, nBins, offset, npnts
	STRING		resName=NameOfWave(inwav)+"_bin"

	IF ((starttime==0) && (endtime == 0))
		WAVESTATS/Q/M=1 inwav
		starttime= V_min
		endtime = V_max
	ENDIF
	
	npnts=DimSize(inwav,0)
	offset=floor(  (starttime) /dt )*dt
	nBins= ceil((endtime-offset)/dt)
	
	MAKE/D/N=(nBins)/O $resName
	WAVE reswav=$resName
	SetScale/P x,offset, dt, "s", reswav

	reswav=0
	
	FOR (k=0; k<npnts; k+=1)
		point=floor( (inwav[k]-offset) / dt )
		reswav[point]+=1
	ENDFOR
	reswav/=dt
END 		// BinEvents

FUNCTION/S NormedCrossCorrResultName_AN(w1,w2)
WAVE		w1,w2
	// uses the built in Correlation
	// parameters are default i.e. linear correlation without any normalization or mean removal
	// afterwards the result is normalized such that every point is divided by the number of 
	// overlapping points that created that summed up to the result
	
	VARIABLE r, n1=max(DimSize(w1,0),DimSize(w2,0)), , n2=min(DimSize(w1,0),DimSize(w2,0))
	
	Redimension/D w1,w2
	STRING 		CCName="CC" + NameOfWave(w1)+"_"+NameOfWave(w2)+"_n"
	DUPLICATE/O w2,$CCName
	WAVE	CCn=$CCName
	
	Correlate w1, CCn
	
	CCn[0,n2-1]/=p+1		// first part of partial overlap (including 1 point with maximal overlap)
	IF (n1>n2)
		CCn[n2,n1-1]/=n2
	ENDIF
	CCn[n1,n1+n2-2]/=(n1+n2-1)-p
	return CCName
END

FUNCTION NormedCrossCorr_AN(w1,w2)
WAVE		w1,w2
	// uses the built in Correlation
	// parameters are default i.e. linear correlation without any normalization or mean removal
	// afterwards the result is normalized such that every point is divided by the number of 
	// overlapping points that created that summed up to the result
	// the result is returned in w2, which is overwritten!
	
	
	VARIABLE r, n1=max(DimSize(w1,0),DimSize(w2,0)), , n2=min(DimSize(w1,0),DimSize(w2,0))
	
	Redimension/D w1,w2
	
	Correlate w1, w2
	
	w2[0,n2-1]/=p+1		// first part of partial overlap (including 1 point with maximal overlap)
	IF (n1>n2)
		w2[n2,n1-1]/=n2
	ENDIF
	w2[n1,n1+n2-2]/=(n1+n2-1)-p
END


FUNCTION CondFiringWrapper(W_FolderNames,S_StimListName, nStims,duration, binning_dt, convolution_sigma, Initial_settling_time)
WAVE/T	W_FolderNames
STRING	S_StimListName
VARIABLE	nStims, duration
VARIABLE	binning_dt, convolution_sigma, Initial_settling_time		// all in seconds

	// the function wraps the application of "Correlate2Periods" for a whole set of cells and periods
	// to aid summary, the periods (illumination periods) are handed over in a matrix, where column corresponds to 
	// category (i.e. the illumination during period 1 and 4 was identical so 1 and 4 belong in the same column)
	// conditional rates are computed for all possible combinations so one has to be a bit careful in filling the matricies
	// combinatorial scaling applies!
	
	// for calculation, a simple loop over electrodes and periods is applied
	// for display the grouping is important
	
	VARIABLE	sorted=0	// set to 1 if the sequence always has to be rate1>rate2
	
	
	STRING		currentDFPath, S_ResultName
	currentDFPath = GetDataFolder(1)		// OK

	
	VARIABLE	nC1,nC2, Stim, nCells= DimSize(W_FolderNames,0)
	STRING	C1,C2, shortC1, shortC2
	
	STRING		AllResults=""
	
	// loop accross stimulus number
	FOR (Stim=0; Stim<nStims; Stim+=1)
						
		
		// loop accross first cell
		FOR (nC1=0; nC1< nCells; nC1+=1)
		//	C1=currentDFPath+W_FolderNames[nC1]+":"			// if folder list contains only partial path
			C1=W_FolderNames[nC1]+":"								// if folder list contains full path
			shortC1=StringFromList(ItemsInList(C1,":")-1, C1,":")
			WAVE/T	StimList=$(C1+S_StimListName)
			C1+=StimList[Stim]
			WAVE	w_C1=$C1
			
			// loop accross second cell
//			FOR (nC2=nC1+1; nC2<nCells; nC2+=1)
			FOR (nC2=0; nC2<nCells; nC2+=1)
				IF (nC2 == nC1)
					nC2+=1
				ENDIF
				
				IF (nC2 < nCells)
					//	C2=currentDFPath+W_FolderNames[nC2]+":"			// if folder list contains only partial path
					C2=W_FolderNames[nC2]+":"									// if folder list contains full path
					shortC2=StringFromList(ItemsInList(C2,":")-1, C2,":")
	
					WAVE/T	StimList=$(C2+S_StimListName)
					C2+=StimList[Stim]
					WAVE	w_C2=$C2
					
					IF (!(Sorted) || (DimSize(w_C1,0)>=DimSize(w_C2,0)) )	// only compute the conditional firing rate between two cells ONce and 
																		// always with the higher firing rate going first
					// the actual computation
						S_ResultName=RemoveEnding(NameOfWave(w_C1),"_1_ST")+"_"+RemoveEnding(NameOfWave(w_C2),"_1_ST")+"_CR"
						Correlate2Periods(w_C1,w_C2, duration, binning_dt,convolution_sigma,Initial_settling_time,0,1)
					ELSE				// i.e. Do sort, but the second sweep has more spikes
						S_ResultName=RemoveEnding(NameOfWave(w_C2),"_1_ST")+"_"+RemoveEnding(NameOfWave(w_C1),"_1_ST")+"_CR"
						Correlate2Periods(w_C2,w_C1, duration, binning_dt,convolution_sigma,Initial_settling_time,0,1)
					ENDIF
					// the zero as penultimate parameter assures that the conditional firing rate is computed, rather than the cross correlation
					// the resulting wave has the following name 
					// the zero as last parameter requests NON-symmetric normalization, i.e. division by the rate of the conditioning neuron 
					// NOT by the goemetric mean of the two rates from both neurons
					WAVE res=$S_ResultName
					IF (!(Sorted) || (DimSize(w_C1,0)>=DimSize(w_C2,0)) )
						S_ResultName="Per_"+num2istr(Stim)+"_"+shortC1+"_" +shortC2+"_CR"
					ELSE
						S_ResultName="Per_"+num2istr(Stim)+"_"+shortC2+"_" +shortC1+"_CR"
					ENDIF
	
					Duplicate/O res,$S_ResultName
					
					AllResults+=S_ResultName +";"
					KillWaves/Z res
				ENDIF			
			
			ENDFOR	// second cell
			
		ENDFOR // first cell
		
	ENDFOR	// Stim
				
END	





FUNCTION Correlate2Periods(ST1,ST2,duration, dt, ConvolSigma, transtime, NODC, SYMMETRIC)
											// dt is the time bin of spike train binning
VARIABLE	dt, ConvolSigma			// sigma (secs)  		- width of the Gaussian envelope, the binary spike train is convolved with SECONDS
VARIABLE	duration, transtime		// duration (secs)	- duration of the stimulus - only needed to calculate the average firing rate
											// transtime (secs) 	- time to wait for any stimulus onset effects in the spike rate to ebb off SECONDS

WAVE		ST1, ST2						// the spike times (secs) of the 2 electrodes from which spikes are analysed

VARIABLE	NODC, Symmetric			// if larger than 0, the DC component is removed from the waves prior to 
											// correlation, the result is NOT the conditional firing rate nü, but nü_cond-sqrt/(nü_1*nü_2)
											// for NODC=0, a different normalization for the result of the correlation has to be used:
											// the cross-correlated or Auto-Correlated signal is, target point by target point,
											// divided by the number of points that entered the calculation for that target point
											// this is achieved in the function NormedCrossCorr_AN(w1,w2)
											// Symmetric means division of the conditional rate by sqrt(rate1*rate2)
											// NON-Symmetric means division by rate1
	
	
	VARIABLE	s_first, s_last, scount1, scount2
	VARIABLE	earliestSpikeTime, latestSpikeTime
	
	// select which extension the results will have
	STRING		ext
	IF (NODC)
		ext="_CC"
	ELSE
		ext="_CR"
	ENDIF
	
	duration	=	duration - transtime


	scount1=DimSize(ST1,0)
	scount2=DimSize(ST2,0)
	
	// Decide which spikes to include based on the trastime, the time at the beginning that gets ignored
	// to account for stimulus onset effects (non-stationarity)
	
	IF (ST1[0] < transtime)
		FindLevel /EDGE=1 /P/Q ST1, transtime; 
		s_first=ceil(V_LevelX)
 
	ELSE
		s_first=0
	ENDIF
	s_last=scount1 - 1
	scount1=s_last-s_first+1
	Duplicate/O/R=[s_first,s_last] ST1 AN_TempSource1
	WAVE Source1= AN_TempSource1
	Source1-=transtime
	
	IF  (ST2[0] < transtime )
		FindLevel /EDGE=1 /P/Q ST2, transTime;  
		s_first=ceil(V_LevelX)
	ELSE
		s_first=0
	ENDIF
	s_last=scount2 -1
	scount2=s_last-s_first+1
	
	Duplicate/O/R=[s_first,s_last] ST2 AN_TempSource2
	WAVE Source2= AN_TempSource2
	Source2-=transtime
	
	
	earliestSpikeTime=min(Source1[0],Source2[0])
	latestSpikeTime=max(Source1[scount1-1],Source2[scount2-1])
	BinEvents(Source2,dt,earliestSpikeTime-4*ConvolSigma,latestSpikeTime+4*ConvolSigma)
	BinEvents(Source1,dt,earliestSpikeTime-4*ConvolSigma,latestSpikeTime+4*ConvolSigma)

	KillWaves Source1,Source2
	
	WAVE Source1=AN_TempSource1_bin
	WAVE Source2=AN_TempSource2_bin

	// convolves binary, binned spike trains with gaussian kernel
	// If the width of the Kernel is larger than zero , here I use "larger than dt"
	IF (ConvolSigma > dt)
		ConvolveEvents(Source1, ConvolSigma)
		ConvolveEvents(Source2, ConvolSigma)
	ENDIF
		
	STRING NoteStr="TransTime ="+num2str(transtime)+" s\r"
	NoteStr+="Duration= "+num2str(duration) + " s\r"
	NoteStr+="spikeCount= " + num2istr(scount1)+"\r"
	
	
	Duplicate/O Source1,$(NameOfWave(ST1)+"_AC");
	WAVE AC=$(NameOfWave(ST1)+"_AC")
	SetScale/P x,0,dt,"s", AC	// assures that the AC peak is at ZERO
	IF (NODC)
		Correlate/NODC Source1, AC;
	ELSE
		NormedCrossCorr_AN( Source1, AC)
	ENDIF
	DeletePoints 0, -(DimOffset(AC,0)+1)/DimDelta(AC,0), AC; 
	SetScale/P x,-1,DimDelta(AC,0),"s",AC;
	DeletePoints 2/DimDelta(AC,0),10000000,AC
	Note/K AC; Note AC,NoteStr

	Duplicate/O AC, $(NameOfWave(ST1)+"_AC_n");
	WAVE AC=$(NameOfWave(ST1)+"_AC_n")
	AC/=scount1/duration
	
	NoteStr="TransTime ="+num2str(transtime)+" s\r"
	NoteStr+="Duration= "+num2str(duration) + " s\r"
	NoteStr+="spikeCount= "+num2istr(scount2)+"\r"

	
	Duplicate/O Source2,$(NameOfWave(ST2)+"_AC");
	WAVE AC=$(NameOfWave(ST2)+"_AC")
	SetScale/P x,0,dt,"s", AC	// assures that the AC peak is at ZERO

	IF (NODC)
		Correlate/NODC Source2, AC;
	ELSE
		NormedCrossCorr_AN(Source2, AC);
	ENDIF
	DeletePoints 0, -( DimOffset(AC,0)+1)/DimDelta(AC,0), AC; 
	SetScale/P x,-1,DimDelta(AC,0),"s",AC;
	DeletePoints 2/DimDelta(AC,0),10000000,AC
	Note/K AC; Note AC,NoteStr

	Duplicate/O AC, $(NameOfWave(ST2)+"_AC_n");
	WAVE AC=$(NameOfWave(ST2)+"_AC_n")
	AC/=scount2/duration
	
	
	
	NoteStr="TransTime ="+num2str(transtime)+" s\r"
	NoteStr+="Duration= "+num2str(duration) + " s\r"
	NoteStr+="spikeCount= "  + num2istr(scount1)+" ; "+num2istr(scount2)+"\r"
	STRING	CCName=RemoveEnding(NameOfWave(ST1),"_1_ST")+"_"+RemoveEnding(NameOfWave(ST2),"_1_ST")+ext
	
	Duplicate/O Source2,$CCName;
	WAVE CC=$CCName;
	SetScale/P x,0,dt,"s",CC;	// make sure the correlated Function has no additional offset in time
	IF (NODC)
		Correlate/NODC Source1, CC;
	ELSE
		NormedCrossCorr_AN(Source1, CC)
	ENDIF
	VARIABLE	ResultDuration=0.2
	DeletePoints 0, -(DimOffset(CC,0)+ResultDuration)/DimDelta(CC,0), CC; 
	SetScale/P x,-ResultDuration,DimDelta(CC,0),"s",CC;
	DeletePoints 2*ResultDuration/DimDelta(CC,0),10000000,CC
	Note/K CC; Note CC,NoteStr

	IF (NODC)	// result is NOT a conditional rate but a cross-correlation 
		Duplicate/O CC,$(NameOfWave(ST1)+NameOfWave(ST2)+"_n");
		WAVE CC=$(NameOfWave(ST1)+NameOfWave(ST2)+ext+"_n");
	ENDIF
	IF (Symmetric)
		CC/=sqrt(scount1*scount2)/duration
	ELSE
		CC/=scount1/duration
	ENDIF
		
	// Cleaning up

	KillWaves Source1, Source2

END		// correlate2Periods





//============================================================================
//-----------------------------------------------------------------------------------------------------------------------------------------
//============================================================================
FUNCTION	VectorStrengthPlot(ListOfFrequencies, BasePrefix, BaseSuffix)
STRING		ListOfFrequencies, BasePrefix, BaseSuffix
// makes a plot of vector strength and phase for all Frequencies
// names of the voltage waves is constructed as
// prefix + Entry from a list + Suffix
// where Entry comes from the first Input to this function, a Semicolon separated list of frequencies
// and Suffix is the ending of a typical spike time wave name after the frequency
// current convention would require to call this function in this way:
// VectorStrengthPlot("50;100;200;500;750;","A","Hz1_ST")

// the analysis assumes that the deterministic sinusoidal component in the input starts at phase 0
// for time=0
// also it assumes that Frequencies are given with unit Hz 

	VARIABLE	NFreqs=ItemsInList(ListOfFrequencies,";")
	VARIABLE	Count, Freq
	STRING		STName, IName
	
	// names of the results
	STRING		VecSName=BasePrefix + "VecStr" + BaseSuffix
	STRING		PhName=BasePrefix + "Phase" + BaseSuffix
	STRING		FreqsName=BasePrefix + "Freqs" + BaseSuffix
	// create waves that hold the results
	MAKE/O/N=(NFreqs)  $VecSName, $PhName, $Freqsname
	WAVE	VecS =	$VecSName
	WAVE	Phase=	$PhName
	WAVE	Freqs=	$FreqsName
	// create a text wave to hold the labels for the frequency axis of the graph
	STRING	LabelName="F_Lbl_"+BaseSuffix
	MAKE/T/O/N=(NFreqs)  $LabelName
	WAVE/T	Labels=$LabelName
		
	VARIABLE/C	VecStrength
	FOR	(Count=0; Count<NFreqs; Count+=1)
		
		Freq=str2num(StringFromList(Count, ListOfFrequencies,";")) 	// turning the characters from the 
															// freqList into a number
															// representing a Frequency in Hz
															
		Labels[Count]=	StringFromList(Count, ListOfFrequencies,";")															
		Freqs[Count]  =	Freq
		
		STName= BasePrefix + StringFromList(Count, ListOfFrequencies,";") + BaseSuffix
		IF (exists(STName ) == 1)		// a wave with the name given in STName exists
			WAVE	ST=$STName
		ELSE
			DoAlert 0, "Could not find wave "+STName
			Return 0
		ENDIF
		// try to find a corresponding current wave (input to the cell)
		// this makes it possible to check whether indeed the input has a prominent component at
		// the specified frequency
		IName= RemoveEnding(STName,"_ST")
		IF (!stringMatch(STName,Iname ) ) // somethign was removed
			IName+="_I"
			IF (exists(Iname ) == 1)		// a wave with the name given in IName exists
				WAVE	Input=$IName
				FFT/OUT=3/DEST=Temp Input;				   // making a fourier transform of the input
				WAVESTATS/Q/R=[1,DimSize(Temp,0)-1] Temp // finding maximum in the FFT
				IF (abs(V_Maxloc-Freq)/Freq > 0.1)			    // checking distance between max in FFT
														    // and the freq, that is supposed to be in the input
					STRING	Alert
					Alert=  "On wave "+IName
					Alert+=" the largest spectral component\r is at"
					Alert+=num2str(V_maxloc)+" and not "
					Alert+=num2str(Freq)+"!"
					DoAlert 0,Alert
				ENDIF
				KillWaves /Z Temp
			ENDIF
		ENDIF
		// now calling the function that returns VectorStrength as a complex number 
		// in rectangular coordinates
		// put those into the prepared output waves (VecS and Phase)
		// by changing from rectangular coordinates into polar coordinate representation of complex numbers
		
	
	VecStrength= ReturnVectorStrengthVector(ST,Freq)	
	VecS[Count]=cabs(VecStrength)
	Phase[Count]= imag(r2polar(VecStrength))*360/2/Pi
	
	ENDFOR
	
	Display /W=(709.8,332,1105.2,539.6) VecS vs Freqs as "Dynamic Gain"+ BaseSuffix
	AppendToGraph/L=L_Phase Phase vs Freqs
	ModifyGraph userticks(bottom)={$FreqsName,$LabelName}
	ModifyGraph mode($VecSName)=4
	ModifyGraph marker($VecSName)=19
	ModifyGraph lStyle($PhName)=2
	ModifyGraph rgb($VecSName)=(0,0,0),rgb($PhName)=(52224,52224,52224)
	ModifyGraph log(bottom)=2
	ModifyGraph lblPos(left)=51,lblPos(L_Phase)=40
	ModifyGraph lblLatPos(L_Phase)=-1
	ModifyGraph btLen=3
	ModifyGraph freePos(L_Phase)={10,bottom}
	ModifyGraph axisEnab(left)={0,0.45}
	ModifyGraph axisEnab(L_Phase)={0.55,1}
	Label left "Vector strength"
	Label bottom "Sine frequency (Hz)"
	Label L_Phase "Phase (°)"
	SetAxis left 0,0.48
	SetAxis bottom 40,750
	ShowInfo
End			// vector strength

//============================================================================
//----------------------------------------------------------------------------
//============================================================================

FUNCTION	STA_Analysis(ListOfFrequencies, Prefix, STSuffix)
STRING		ListOfFrequencies, Prefix, STSuffix
// assumption is that input wave name ends on "_I", where Spike time name end on "_ST"
// call with cms analogous to "   STA_Analysis("50;100;200;500;750;", "A", "hz1_ST")   "


// this function creates a "continuous" dynamic gain curve from a spike triggered average 
// this is done for each input (50 Hz, 100 Hz)
// the duration of the STA is HARD-WIRED to be 1 s:
	VARIABLE	STAwidth=1		// one second. change if you whish

	VARIABLE	NFreqs=ItemsInList(ListOfFrequencies,";")
	VARIABLE	Count, Freq
	STRING		STName, IName
	
	STRING		FreqsName=	Prefix + "Freqs" + STSuffix
	WAVE		Freqs	    =	$FreqsName
	// create a text wave to hold the names of the individual smoothed gain curves
	// one for each frequency
	STRING	GainNameList=""
	STRING	xWaveList=""
	// the name of the wave for average gain and standard error of the gain 
	VARIABLE	STSuffixLength=strlen(STSuffix)
	STRING		Suffix		// general part of a name before the identifier (_ST or _STA or _V)
	
	IF 	(StringMatch(STSuffix[STSuffixLength-3,STSuffixLength-1],"_ST" ))	// last to characters
																// in STsuffix are "_ST"
		Suffix=STSuffix[0,STSuffixLength-4]
	ELSE
		Suffix=STSuffix
	ENDIF

	// Create Graph to which the individual gain curves and the average will be appended
	Display 		as	Prefix+"_GAIN_"+Suffix
	
	FOR	(Count=0; Count<NFreqs; Count+=1)
		Freq=str2num(StringFromList(Count, ListOfFrequencies,";")) 	// turning the characters from the 
															// freqList into a number
															// representing a Frequency in Hz
		Freqs[Count]=Freq
		STName= Prefix + StringFromList(Count, ListOfFrequencies,";") + STSuffix
		IF (exists(STName ) == 1)		// a wave with the name given in STName exists
			WAVE	ST=$STName
		ELSE
			DoAlert 0, "Could not find wave "+STName
			Return 0
		ENDIF
		// try to find a corresponding current wave (input to the cell)
		// this makes it possible to check whether indeed the input has a prominent component at
		// the specified frequency
		IName= RemoveEnding(STName,"_ST")
		IF (!stringMatch(STName,Iname ) ) // somethign was removed
			IName+="_I"
			IF (exists(IName ) == 1)		// a wave with the name given in IName exists
				WAVE	Input=$IName
				FFT/OUT=3/DEST=Temp Input;				   // making a fourier transform of the input
				WAVESTATS/Q/R=[1,DimSize(Temp,0)-1] Temp // finding maximum in the FFT
				IF (abs(V_Maxloc-Freq)/Freq > 0.1)			    // checking distance between max in FFT
														    // and the freq, that is supposed to be in the input
					STRING	Alert
					Alert=  "On wave "+IName
					Alert+=" the largest spectral component\r is at"
					Alert+=num2str(V_maxloc)+" and not "
					Alert+=num2str(Freq)+"!"
					DoAlert 0,Alert
				ENDIF
				KillWaves /Z Temp
			ENDIF
		ENDIF


// create a Spike triggered average input from the spike times and the input wave
		STAfromAnalogue( Input, ST, STAwidth,0)
		// this creates a wave with a name that corresponds to STName with added "_STA"
		// lets rename this according to convention, i.e. from "_ST_STA" make "_STA"
		VARIABLE	STNameLength=strlen(STName)
		STRING		STAName
		
			STAName=STName+"_STA"														
			WAVE STA=$STAName
		IF 	(StringMatch(STName[STNameLength-3,STNameLength-1],"_ST" ))	// last to characters
																	// in STName are "_ST"
			
																	
			Duplicate/O  STA, $(STName[0,STNameLength-4]+"_STA")
			KillWaves/Z STA
			WAVE	STA=$(STName[0,STNameLength-4]+"_STA")
			STAName = NameOfWave(STA)
		ENDIF
	
// create the input autocorrelation 
		STRING		ACName=IName+"_AC"
		Duplicate/O input,$ACName
		WAVE	AC=$ACName
		Correlate input, AC
// AC is a very long wave,  retain only the same region that is present in the STA (i.e. STAwidth)
		// find the point that corresponds to the same time as the first point in the STA
		VARIABLE	firstP=x2pnt(AC, DimOffset(STA,0) )
		IF (firstP >0)
			DeletePoints  0, firstP, AC
		ENDIF
		Redimension/N=(DimSize(STA,0)) AC
		SetScale/P x,DimOffset(STA,0),DimDelta(STA,0),WaveUnits(STA, 0 ),AC
	
// make circular shifted versions of the STA and the Autocorrelation (AC)
// in order to have zero time at point 0 (begin of wave)
		SplitBeforeFFT(AC, 0)
		WAVE AC=$(ACName+"_splt")
		SplitBeforeFFT(STA, 0)
		WAVE STA=$(STAName+"_splt")
		Redimension/D AC, STA				// use double precision to calculate the FFT
											// gives slightly higher safety against getting
										// a perfect zero in any of the entries

// calculate the complex wave gain, which is the quotient  FFT(STA) / FFT(AC)
	// can only do FFT for an even number of points, leave the last point out:
		VARIABLE	evenEnd=2*ceil((DimSize(STA,0)-1)/2)-1
		STRING		STAFFTName=RemoveEnding(STAName,"_STA")+"_Gain"
		FFT/OUT=1/RP=[0,evenEnd]/DEST=$STAFFTName STA
		WAVE/C	Gain=$STAFFTName
		FFT/OUT=1/RP=[0,evenEnd]/DEST=Temp AC
		Gain/=Temp
		Killwaves/Z Temp, AC,STA
		
// smooth the complex gain and return a wave for the magnitude and a wave for the phase

		IF (Exists("FrequencyRange")==1)
			WAVE FrequencyRange
			GaussFilter(Gain,FreqPoints=FrequencyRange, widthFac=1)
			xWaveList+="FrequencyRange;"
			AppendToGraph $(STAFFTName +"_MgFlt") vs FrequencyRange

		ELSE
			GaussFilter(Gain)
			WAVE	FreqPoints		// amde by the function GaussFilter
			xWaveList+="FreqPoints;"
			AppendToGraph $(STAFFTName +"_MgFlt") vs FreqPoints
		ENDIF
		GainNameList+=STAFFTName +"_MgFlt;"

	

	ENDFOR
	ModifyGraph rgb=(0,0,0)
	
// average the magnitude waves accross all the different inputs (50 Hz 100 Hz...)

	STRING 		avgName=Prefix+"_gain_avg_"+Suffix, SEMName=Prefix+"_gain_SEM_"+Suffix
	fWaveAverage(GainNameList, xWaveList, 3, 0, avgName, SEMName)
	WAVE	G_avg=$avgName
	WAVE	G_SEM=$SEMName
	IF (Exists("FrequencyRange")==1)
		AppendToGraph G_avg vs FrequencyRange
	ELSE
		AppendToGraph G_avg vs FreqPoints
	ENDIF
	ModifyGraph mode($avgName)=4,marker($avgName)=19,msize($avgName)=1.5;DelayUpdate
	ModifyGraph rgb($avgName)=(65280,43520,0);DelayUpdate
	ErrorBars $avgName Y,wave=($SEMName,$SEMName)
	ModifyGraph log=1
	Label left "Dynamic gain (Hz/A)";DelayUpdate
	Label bottom "Frequency (Hz)"
	ModifyGraph btLen=4

END		// STA Analysis



FUNCTION	NoiseFloor(ListOfFrequencies, Prefix, STSuffix, Tau_corr)
STRING		ListOfFrequencies, Prefix, STSuffix
VARIABLE	Tau_corr		// correlation time of the input
// this does an estimation of the region that we can trust the STA derived GAin
// it takes the spike times and shifts them forward by 5 correlation times plus  a random 
// phase of the sinusoid that is present
// this is done 20 times 
// at each frequency point the second largest gain value is taken to represent the limit to the top
// 5% of possible random STA gains

// the duration of the STA is HARD-WIRED to be 1 s:
	VARIABLE	STAwidth=1		// one second. change if you whish
	VARIABLE	repetitions=100	// number of shifted spike time wave to estimate noise floor 
								// cannot be choosen too high as with each reprtition the spike times 
								// are shifted by 5 tau corr
	IF (repetitions*tau_corr > 10 )
		DoAlert 0,"You are requesting many repetitions, please consider using completely random spike times\r rather	than shifting the original times"
	ENDIF							

	VARIABLE	NFreqs=ItemsInList(ListOfFrequencies,";")
	VARIABLE	Count, Freq
	STRING		STName, IName
	
	STRING		FreqsName=	Prefix + "Freqs" + STSuffix
	WAVE		Freqs	    =	$FreqsName
	// create a text wave to hold the names of the individual smoothed gain curves
	// one for each frequency
	STRING	GainNameList=""
	STRING	xWaveList=""
	// the name of the wave for average gain and standard error of the gain 
	VARIABLE	STSuffixLength=strlen(STSuffix)
	STRING		Suffix		// general part of a name before the identifier (_ST or _STA or _V)
	
	IF 	(StringMatch(STSuffix[STSuffixLength-3,STSuffixLength-1],"_ST" ))	// last to characters
																// in STsuffix are "_ST"
		Suffix=STSuffix[0,STSuffixLength-4]
	ELSE
		Suffix=STSuffix
	ENDIF

	STRING		MagName=""	// Name of the gauss Windowed Magnitude 
	STRING		NoiseFloorName=""
	// Create Graph to which the individual gain curves and the average will be appended
	
	FOR	(Count=0; Count<NFreqs; Count+=1)
		Freq=str2num(StringFromList(Count, ListOfFrequencies,";")) 	// turning the characters from the 
															// freqList into a number
															// representing a Frequency in Hz
		Freqs[Count]=Freq
		STName= Prefix + StringFromList(Count, ListOfFrequencies,";") + STSuffix
		NoiseFloorName=Prefix+"_noise_"+ StringFromList(Count, ListOfFrequencies,";")+Suffix

		IF (exists(STName ) == 1)		// a wave with the name given in STName exists
			WAVE	ST=$STName
		ELSE
			DoAlert 0, "Could not find wave "+STName
			Return 0
		ENDIF
		// try to find a corresponding current wave (input to the cell)
		// this makes it possible to check whether indeed the input has a prominent component at
		// the specified frequency
		IName= RemoveEnding(STName,"_ST")
		IF (!stringMatch(STName,Iname ) ) // somethign was removed
			IName+="_I"
			IF (exists(IName ) == 1)		// a wave with the name given in IName exists
				WAVE	Input=$IName
				FFT/OUT=3/DEST=Temp Input;				   // making a fourier transform of the input
				WAVESTATS/Q/R=[1,DimSize(Temp,0)-1] Temp // finding maximum in the FFT
				IF (abs(V_Maxloc-Freq)/Freq > 0.1)			    // checking distance between max in FFT
														    // and the freq, that is supposed to be in the input
					STRING	Alert
					Alert=  "On wave "+IName
					Alert+=" the largest spectral component\r is at"
					Alert+=num2str(V_maxloc)+" and not "
					Alert+=num2str(Freq)+"!"
					DoAlert 0,Alert
				ENDIF
				KillWaves /Z Temp
			ENDIF
		ENDIF

		VARIABLE	rndCount=0
		WAVESTATS/Q/M=1 ST		// get time  of first and last spike into V_min and V_max
		FOR (rndCount=0; rndcount<repetitions ; rndCount+=1)
			WAVE		ST_Temp=RandomInInterval(ST)
//			WAVE		ST_Temp=CyclicShift(ST,Tau_corr*5+abs(enoise(V_max-V_min-Tau_corr*5)) )
			// create a Spike triggered average input from the spike times and the input wave
			STAfromAnalogue( Input, ST_Temp, STAwidth,0)
			// this creates a wave with a name that corresponds to STName with added "_STA"
			// lets rename this according to convention, i.e. from "_ST_STA" make "_STA"
			STRING	STAName=NameOfWave(ST_Temp)+"_STA"
			WAVE STA=$STAName
		
	// create the input autocorrelation once
			IF (rndCount==0)
				STRING		ACName=IName+"_AC"
				Duplicate/O input,$ACName
				WAVE	AC=$ACName
				Correlate input, AC
	// AC is a very long wave,  retain only the same region that is present in the STA (i.e. STAwidth)
				// find the point that corresponds to the same time as the first point in the STA
				VARIABLE	firstP=x2pnt(AC, DimOffset(STA,0) )
				IF (firstP >0)
					DeletePoints  0, firstP, AC
				ENDIF
				Redimension/N=(DimSize(STA,0)) AC
				SetScale/P x,DimOffset(STA,0),DimDelta(STA,0),WaveUnits(STA, 0 ),AC
				SplitBeforeFFT(AC, 0)
				WAVE AC=$(ACName+"_splt")
				Redimension/D AC
			ENDIF
	// make circular shifted versions of the STA and the Autocorrelation (AC)
	// in order to have zero time at point 0 (begin of wave)
			SplitBeforeFFT(STA, 0)
			WAVE STA=$(STAName+"_splt")
			Redimension/D STA				// use double precision to calculate the FFT
												// gives slightly higher safety against getting
											// a perfect zero in any of the entries
	
	// calculate the complex wave gain, which is the quotient  FFT(STA) / FFT(AC)
		// can only do FFT for an even number of points, leave the last point out:
			VARIABLE	evenEnd=2*ceil((DimSize(STA,0)-1)/2)-1
			STRING		STAFFTName=RemoveEnding(STAName,"_STA")+"_Gain"
			FFT/OUT=1/RP=[0,evenEnd]/DEST=$STAFFTName STA
			WAVE/C	Gain=$STAFFTName
			FFT/OUT=1/RP=[0,evenEnd]/DEST=Temp AC
			Gain/=Temp
			Killwaves/Z Temp,STA
			
	// smooth the complex gain and return a wave for the magnitude and a wave for the phase
	
			IF (Exists("FrequencyRange")==1)
				WAVE FrequencyRange
				GaussFilter(Gain,FreqPoints=FrequencyRange, widthFac=1)
			ELSE
				GaussFilter(Gain)
				WAVE	FreqPoints		// made by the function GaussFilter
			ENDIF
			MagName = STAFFTName +"_MgFlt"
			WAVE	Mag= $MagName
			IF (rndCount==0)
				Duplicate/O Mag , $NoiseFloorName
				WAVE	Noise=$NoiseFloorName
				Redimension/N=(-1,repetitions) Noise
			ELSE
				Noise[][rndCount]=Mag[p]
			ENDIF		
		ENDFOR // rndCount
		
		// now out of all different realisations find the largest value for each Frequency 
		VARIABLE	FreqEntry, NFreq=DimSize(Noise,0)
		FOR (FreqEntry=0; FreqEntry < NFreq; FreqEntry+=1)
			Duplicate /O /R=[FreqEntry][*] Noise, Temp_Row
			Sort/R Temp_Row, Temp_Row	// sort from largest to smallest
			Noise[FreqEntry][0]=Temp_row[0][round(repetitions/20)]
			Noise[FreqEntry][1]=Temp_row[0][round(repetitions/2)]
			Noise[FreqEntry][2]=Temp_row[0][repetitions-round(repetitions/20)]
		ENDFOR	
		Redimension/N=(NFreq,3) Noise
		
	ENDFOR // across Frequencies
		KillWaves/Z ST_Temp, AC, Temp_Row
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


FUNCTION/WAVE RandomInInterval(SpikeTimeWave)
WAVE		SpikeTimeWave	
// receives a wave of (spike) time points
// returns a wave containing time points 
// same number as input, but randomly placed in the interval between 
// first and last spike time from the input wave

	// prepare the wave copy which will be manipulated and finally returned (by reference)
	Duplicate/O SpikeTimeWave, Temp_ST
	
	WAVE	ST=Temp_ST
	
	// just make sure the spike times are sorted:
	Sort ST,ST
	WaveStats/Q /M=1 ST
	
	ST=enoise((V_max-V_min)/2,2)
	ST+=V_min+(V_max-V_min)/2
	// all spikes with  indicies after V_LevelX have to be cycled to the begin
	Sort ST, ST
	
	Return ST
END // random in Interval




//________________________________________________________________
//===========================  group analysis =========================
//________________________________________________________________

FUNCTION AvgSTA([Suffix])
STRING		Suffix

		IF (ParamIsDefault(Suffix))
			Suffix = "STA"
		ENDIF


STRING		STA_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
		cellFolderList = ReplaceString(",", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";"), flagLocal=0
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK
VARIABLE	k, nSTAs, SD, nspikes, totspikes, duration, totduration=0

STRING		Notiz, Note_String

VARIABLE	FirstTime=1, totalSTANumber=0, AvgSD=0
// go through all subfolders

// catch case that there are no subfolders and that only in the current folder the average is created
	IF (nSubs == 0)
		nSubs=1
		flagLocal=1
	ENDIF	
	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		IF (flagLocal==0)
			FolderName=StringFromList(i, cellFolderList,";")
			SetDataFolder $(FolderName) 
		ELSE // local, no subfolders
			FolderName= GetDataFolder(1)
		ENDIF
		STA_List=WaveList("*_"+Suffix, ";", "" )
		nSTAs=ItemsInList(STA_List,";")
//		IF(Exists( currentDFPath+"STA_avg")==1)
//			WAVE	Avg=$(currentDFPath+"STA_avg")
//			Note/NOCR Avg,Foldername+"\r"
//		ENDIF
		FOR (k=0; k< nSTAs; k+=1)
			WAVE	STA=$(StringFromList(k, STA_List,";") )
			Notiz=Note(STA)
			Note_String=StringByKey("SD", Notiz,":" ,"\r")
			SD=str2num(Note_String)
			// read out standard deviation of the individual STA will be normalized with below
			IF (SD!=SD)	// SD is NaN
				DoAlert 0,"Incorrect wave note in "+GetWavesDataFolder(STA,2)
			ENDIF
			// now read out the number of spikes that went into the STA
			Note_String=StringByKey("#spikes", Notiz,":" ,"\r")
			nspikes=str2num(Note_String)
			
			// now read ou the duration of the trial for that STA
			Note_String=StringByKey("duration", Notiz,":" ,"\r")
			duration=str2num(Note_String)

			IF (FirstTime)
				Duplicate/O STA, $(currentDFPath+Suffix+"_avg")
				WAVE	Avg=$(currentDFPath+Suffix+"_avg")
				Avg=0
				FirstTime=0
				Note/K Avg
				Note/NOCR Avg, currentDFPath+"\r"
				Note/NOCR Avg,FolderName+"\r"
			ENDIF
			Avg+=STA/SD*duration
			totalSTANumber+=1
			AvgSD+=SD*duration
			totspikes+=nspikes
			totduration=totduration+duration
			Note/NOCR Avg,"\t"+ StringFromList(k, STA_List,";") +"\r"
		ENDFOR		// loop across STAs in one folder
		Note/K Avg
		Note/NOCR Avg, currentDFPath+"\r"
		Note/NOCR Avg,FolderName+"\r"
	ENDFOR			// loop across folders
	Avg/=totduration
	AvgSD/=totduration

	IF (flagLocal==0)
		VARIABLE/G avgInputSD=AvgSD
		
	ENDIF

	Note/NOCR Avg,"total#trial:"+num2str(totalSTANumber) +"\r"
	Note/NOCR Avg,"total#spikes:"+num2str(totspikes) +"\r"
	Note/NOCR Avg,"totalduration:"+num2str(totduration) +"\r"
	Note/NOCR Avg,"AVG_SD:"+num2str(AvgSD) +"\r"
	Duplicate/O Avg,$(currentDFPath+Suffix+"_avg_scaled")
	WAVE scaled = $(currentDFPath+Suffix+"_avg_scaled")
	scaled*=AvgSD						// creating an average STA with units "Ampere"

	SetDataFolder rootFolderRf

END		// avg STAs


//__________________________

// - - - - - - - - - - - - ThreadSafe Version with hand-over of wave references of all STAs

FUNCTION AvgSTAFromList(ListOfWaves,resultname)
STRING		ListOfWaves,resultname
#if (IgorVersion() >= 7.0) 
	
	Wave/WAVE TempFree = ListToWaveRefWave(ListOfWaves, 1)
	
	WAVE result=AvgSTA_ThreadSafe(TempFree, targetName=resultname)
#endif
END

Threadsafe FUNCTION/WAVE AvgSTA_ThreadSafe(WRefWave[, targetName])
WAVE/WAVE		WRefWave
// contains references to all waves that have to be averaged
STRING	targetName
// this version of the AvgSTA function receives the wave references of all the waves to be averaged
// this can be used to 
// 1. pre-select by certain criteria (number of spikes, spike features etc.)
// 2. run the Function in an pre-emptive Thread where it does not access the folders in the standard way

// this Function also RETURNS a wave reference!
// it also creates the waves avg and scaled
// This version has been tested to return the same result as the 

	VARIABLE		k, nSTAs=DimSize(WRefWave,0)
	VARIABLE		SD, nspikes, totspikes, duration, totduration=0
	
	STRING		Notiz, Note_String
	
	VARIABLE		totalSTANumber=0, AvgSD=0

// create the avg wave by duplication from the first individual STA
 			WAVE	STA=WRefWave[k]
			Duplicate/O STA, staavg
			WAVE Avg=staavg
			Avg=0
			Note/K Avg
 
		
		FOR (k=0; k< nSTAs; k+=1)
			WAVE	STA=WRefWave[k]
			Notiz=Note(STA)
			Note_String=StringByKey("SD", Notiz,":" ,"\r")
			IF (strlen(Note_String)==0)
				print "not note on "+GetWavesDataFolder(STA,2)
			ENDIF
			SD=str2num(Note_String)
			// read out standard deviation of the individual STA will be normalized with below
			
			// now read out the number of spikes that went into the STA
			Note_String=StringByKey("#spikes", Notiz,":" ,"\r")
			nspikes=str2num(Note_String)
			
	  		// now read ou the duration of the trial for that STA
			Note_String=StringByKey("duration", Notiz,":" ,"\r")
			duration=str2num(Note_String)

			Avg+=STA/(SD^2)*duration // possibly SD^2
			totalSTANumber+=1
			AvgSD+=(SD^2)*duration		//  SD^2 if used above
			totspikes+=nspikes
			totduration=totduration+duration
			Note/NOCR Avg,"\t"+ GetWavesDataFolder(STA,2) +"\r"
		ENDFOR		// loop across STAs 
	Avg/=totduration
	AvgSD/=totduration


	Note/NOCR Avg,"total#trial:"+num2str(totalSTANumber) +"\r"
	Note/NOCR Avg,"total#spikes:"+num2str(totspikes) +"\r"
	Note/NOCR Avg,"totalduration:"+num2str(totduration) +"\r"
	Note/NOCR Avg,"AVG_SD:"+num2str(AvgSD) +"\r"
	IF (!ParamIsDefault(targetName))
		Duplicate/O Avg,$targetName
		WAVE scaled=$targetname
	ELSE
		Duplicate/O Avg,scaled
		WAVE scaled 
	ENDIF
	
	scaled*=AvgSD						// creating an average STA with units "Ampere"
		 
	Return scaled
END		// avg STAs Threadsafe

// - - - - - - - - - - - - ThreadSafe Version with hand-over of wave references of all STAs

FUNCTION AvgACFromList(ListOfWaves,resultname)
STRING		ListOfWaves	// semicolon separated list 
STRING		resultname
#if (IgorVersion() >= 7.0) 
	
	Wave/WAVE TempFree = ListToWaveRefWave(ListOfWaves, 1)
	
	WAVE result=AvgAC_ThreadSafe(TempFree)
	Duplicate /O result,$resultname
#endif
END

Threadsafe FUNCTION/WAVE AvgAC_ThreadSafe(WRefWave)
WAVE/WAVE		WRefWave
// contains references to all waves that have to be averaged

// this version of the AvgAC function receives the wave references of all the waves to be averaged
// this can be used to 
// 1. pre-select by certain criteria (number of spikes, spike features etc.)
// 2. run the Function in an pre-emptive Thread where it does not access the folders in the standard way

// this Function also RETURNS a wave reference!
// it also creates the waves avgac and scaledavgac
// This version has been tested to return the same result as the standard MakeAndAvgAC
// it relies on waves and wavenotes that exist only after MakeAndAvgAC has been run once

	VARIABLE		k, nSTAs=DimSize(WRefWave,0)
	VARIABLE		var, nspikes, duration, totduration=0
	
	STRING		Notiz, Note_String
	
	VARIABLE		totalACNumber=0, Variance_avg=0

// create the avg wave by duplication from the first individual STA
 			WAVE	AC=WRefWave[k]
			Duplicate/O AC, avgac
 			avgac=0
 			Note/K avgac
		
		FOR (k=0; k< nSTAs; k+=1)
			WAVE	ac=WRefWave[k]
			Notiz=Note(AC)
			Note_String=StringByKey("variance", Notiz,":" ,"\r")
			var=str2num(Note_String)
			// read out standard deviation of the individual current,
			// as AC will be normalized with variance - see below
			
			
	  		// now read ou the duration of the trial for that AC
			Note_String=StringByKey("duration", Notiz,":" ,"\r")
			duration=str2num(Note_String)

			Avgac+=AC/var*duration
			Variance_avg+=var*duration

			totalACNumber+=1
			totduration=totduration+duration
			Note/NOCR Avgac,"\t"+ GetWavesDataFolder(ac,2) +"\r"
		ENDFOR		// loop across ACs 
	Avgac/=totduration
	Variance_avg/=totduration

	Note/NOCR Avgac,"total#trial:"+num2str(totalACNumber) +"\r"
	Note/NOCR Avgac,"totalduration:"+num2str(totduration) +"\r"
	Note/NOCR Avgac,"AVG_VARIANCE:"+num2str(Variance_avg) +"\r"

	Duplicate/O Avgac,scaledavgac
	WAVE scaledavgac 
	scaledavgac*=Variance_avg						// creating an average AC with units "Ampere^2"
	SetScale d 0,0,"A^2", scaledavgac

	Return scaledavgac


END		// avg ACs Threadsafe

FUNCTION AvgACForSingleFolder(DoPlot)
VARIABLE DoPlot

STRING		AC_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
		cellFolderList = ReplaceString(",", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK
VARIABLE	k, nACs, VAR
STRING		Notiz, VAR_String

VARIABLE	FirstTime=1, totalACNumber=0, Variance_avg=0, flagLocal=0

// catch case that there are no subfolders and that onyl in the current folder the average is created
	IF (nSubs == 0)
		nSubs=1
		flagLocal=1
	ENDIF	

// go through all subfolders
	IF (DoPlot)
			Display as "All ACs"
	ENDIF
	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		IF (flagLocal==0)
			FolderName=StringFromList(i, cellFolderList,";")
			SetDataFolder $(FolderName) 
		ELSE // local, no subfolders
			FolderName= GetDataFolder(1)
		ENDIF

		AC_List=WaveList("*_AC", ";", "" )
		nACs=ItemsInList(AC_List,";")
		IF(Exists( currentDFPath+"AC_avg")==1)
			WAVE	Avg=$(currentDFPath+"AC_avg")
			Note/NOCR Avg,Foldername+"\r"
		ENDIF
		FOR (k=0; k< nACs; k+=1)
			WAVE	AC=$(StringFromList(k, AC_List,";") )
			Notiz=Note(AC)
			VAR_String=StringByKey("variance", Notiz,":" ,"\r")
			VAR=str2num(VAR_String)
			IF (FirstTime)
				Duplicate/O AC, $(currentDFPath+"AC_avg")
				WAVE	Avg=$(currentDFPath+"AC_avg")
				Avg=0
				FirstTime=0
				Variance_avg=0
				Note/K Avg
				Note/NOCR Avg, currentDFPath+"\r"
				Note/NOCR Avg,FolderName+"\r"
			ENDIF
			Avg+=AC/VAR
			IF (DoPlot)
				AppendToGraph AC
			ENDIF
			totalACNumber+=1
			Variance_avg+=VAR
			Note/NOCR Avg,"\t"+ StringFromList(k, AC_List,";") +"\r"
		ENDFOR		// loop across ACs in one folder
	ENDFOR			// loop across folders
	Avg/=totalACNumber
	Variance_avg/=totalACNumber
	IF (flagLocal==0)
		VARIABLE/G avgInputVariance=Variance_avg
		
	ENDIF
	
	Note/NOCR Avg,"AVG_VARIANCE:"+num2str(Variance_avg) +"\r"
	SetDataFolder rootFolderRf

END		// AvgACForSingleFolder

FUNCTION AvgVecStrength()

VARIABLE	NFreqs=5		// number of different sine Freqs used
STRING		VecStr_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
			cellFolderList = ReplaceString(",", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK
VARIABLE	k, nVecStrs
STRING		Notiz

VARIABLE	FirstTime=1, totalVecStrNumber=0
// go through all subfolders
	
	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(i, cellFolderList,";")
		SetDataFolder $(FolderName) 
		VecStr_List=WaveList("*VecStr_*", ";", "" )
		nVecStrs=ItemsInList(VecStr_List,";")
		IF(Exists( currentDFPath+"VecStr_avg")==1)
			WAVE	Avg=$(currentDFPath+"VecStr_avg")
			Note/NOCR Avg,Foldername+"\r"
		ENDIF

		FOR (k=0; k< nVecStrs; k+=1)
			WAVE	VecStr=$(StringFromList(k, VecStr_List,";") )
			IF (FirstTime)
				Duplicate/O VecStr, $(currentDFPath+"VecStr_avg")
				WAVE	Avg=$(currentDFPath+"VecStr_avg")
				Redimension/N=(NFreqs) Avg
				Note/NOCR Avg, currentDFPath+"\r"
				Note/NOCR Avg,FolderName+"\r"
				Avg=0
				FirstTime=0
			ENDIF
			IF (DimSize(VecStr,0) == DimSize(Avg,0) )
				Avg+=VecStr
				totalVecStrNumber+=1
				Note/NOCR Avg,"\t"+ StringFromList(k, VecStr_List,";") +"\r"	
			ENDIF
		ENDFOR		// loop across VecStrs in one folder
	ENDFOR			// loop across folders
	Avg/=totalVecStrNumber
	SetDataFolder rootFolderRf

END		// acg VecStrs


FUNCTION ScaleVecStrength()

MAKE/N=(5)/O W_Temp_Freq={50,100,200,500,750}	
VARIABLE	NFreqs=DimSize(W_temp_Freq,0)		// number of different sine Freqs used
STRING		VecStr_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
			cellFolderList = ReplaceString(",", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK
VARIABLE	fIndex, k, nVecStrs, scale
STRING		Notiz, RateString, AmpString, VecStrName, ACName, STAName, ScaledName

VARIABLE	FirstTime=1, totalVecStrNumber=0, Rate, Amp
// go through all subfolders

// go frequency by frequency as sometimes a frequency is missing in a folder

	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(i, cellFolderList,";")
		SetDataFolder $(FolderName) 
		VecStr_List=WaveList("*VecStr_*", ";", "" )
		nVecStrs=ItemsInList(VecStr_List,";")
		IF(Exists( currentDFPath+"VecStrSCALED_avg")==1)
			WAVE	Avg=$(currentDFPath+"VecStrSCALED_avg")
			Note/NOCR Avg,Foldername+"\r"
		ENDIF
		VARIABLE kk

		FOR (k=0; k< nVecStrs; k+=1)
			VecStrName	=	StringFromList(k, VecStr_List,";")
			WAVE	VecStr=	$VecStrName
			IF (FirstTime)
				Duplicate/O VecStr, $(currentDFPath+"VecStrSCALED_avg")
				WAVE	Avg=$(currentDFPath+"VecStrSCALED_avg")
				Redimension/N=(NFreqs) Avg
				Note/NOCR Avg, currentDFPath+"\r"
				Note/NOCR Avg,FolderName+"\r"
				Avg=0
				FirstTime=0
			ENDIF
			IF (DimSize(VecStr,0) == DimSize(Avg,0) )
				ScaledName=ReplaceString("VecStr", VecStrName,"VecStrGain")
				Duplicate/O VecStr, $ScaledName
				WAVE	Scaled=$ScaledName			// to hold Vector Strength derived gain
				FOR(kk=0; kk<Dimsize(Scaled,0); kk+=1)	
					IF (Scaled[kk]!=Scaled[kk])
						print GetDataFolder(1),NameOfWave(Scaled)
					ENDIF
				ENDFOR
				FOR( fIndex=0; fIndex<NFreqs; fIndex+=1)
				// get for each frequency the sine amplitude (out of wavenote of autocorr)
				// and the firing rate (wavenote of STA)
					STAName= ReplaceString("VecStr", VecStrName, num2istr(W_Temp_Freq[fIndex]) ) +"A"
					// the last "A" turns the suffix ST into STA
					ACName = ReplaceString("_STA",STAName,"_I_AC")
					WAVE	STA=$STAName
					WAVE	AC = $ACName
					Notiz=Note(STA)
					Notiz= ReplaceString(" average spike rate is ",Notiz," average spike rate is:")
					Notiz= ReplaceString("\r",Notiz," \r")
					RateString=StringByKey("is", Notiz,":" ," ")
					Rate=str2num(RateString)
					
					Notiz=Note(AC)
					AmpString=StringByKey("SineAmp", Notiz,":" ,"\r")
					Amp=str2num(AmpString)
					
					Scaled[fIndex]*=2*Rate/Amp
				
				ENDFOR
			
				Avg+=Scaled
				totalVecStrNumber+=1
				Note/NOCR Avg,"\t"+ VecStrName +"\r"	
			ENDIF
		ENDFOR		// loop across VecStrs in one folder
	ENDFOR			// loop across folders
	Avg/=totalVecStrNumber
	SetDataFolder rootFolderRf

END		// scale Vector strengths



			
//________________________________________________________________


FUNCTION MakeAndAvgAC([FirstPoint,LastPoint])
VARIABLE		FirstPoint,LastPoint
// for Omers data this needs to be limited to the relevant part of the current
// points 14532 to 915469
// so for Omer this has to be called
// MakeAndAvgAC(FirstPoint=14532,LastPoint=915469)

IF (ParamIsDefault(FirstPoint))
	FirstPoint = 0
ENDIF
IF (ParamIsDefault(LastPoint))
	LastPoint = inf
ENDIF

STRING		Cur_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
		cellFolderList = ReplaceString(",", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK

STRING		IName, STAName, ACName
VARIABLE	k, nCurs, SD, duration, totduration=0

VARIABLE	FirstTime=1, totalSTANumber=0, Variance_avg=0
// go through all subfolders
	
	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(i, cellFolderList,";")
		SetDataFolder $(FolderName) 
		Cur_List=WaveList("*_I", ";", "" )
		nCurs=ItemsInList(Cur_List,";")
		IF(Exists( currentDFPath+"AC_avg")==1)
			WAVE	Avg=$(currentDFPath+"AC_avg")
			Note/NOCR Avg,Foldername+"\r"
		ENDIF
		FOR (k=0; k< nCurs; k+=1)
			IName	= StringFromList(k, Cur_List,";") 
			STAname=ReplaceString("_I", IName, "_STA")
			ACName=IName+"_AC"

			WAVE	STA=$STAName
			WAVE	Cur=$IName
			// get standard deviation to scale AC for averaging
			// also gets written in wavenote
			WAVESTATS/R=[FirstPoint, min(DimSize(Cur,0)-1,LastPoint)]/Q Cur
			SD=V_sdev
			duration=V_npnts*Dimdelta(Cur,0)
			
			Duplicate/R=[FirstPoint, min(DimSize(Cur,0)-1,LastPoint)]/O Cur,$ACName, W_CurDummy
			// if FirstPOint >0, the Dimoffset is now not corresponding to the 
			// dimOffset of the current wave. 
			// (normally it should be zero, but now it is pnt2x(AC, FirstPoint )
			// This needs to be corrected 
			WAVE	AC=$ACName
			SetScale/P x,DimOffset(Cur,0),DimDelta(Cur,0),"", AC, W_CurDummy

			Correlate/NODC W_CurDummy, AC
			// normalize for number of points in the original current wave
			AC/=V_npnts
// AC is a very long wave,  retain only the same region that is present in the STA (i.e. STAwidth)
			// find the point that corresponds to the same time as the first point in the STA
			VARIABLE	firstP=x2pnt(AC, DimOffset(STA,0) )
			IF (firstP >0)
				DeletePoints  0, firstP, AC
			ENDIF
			Redimension/N=(DimSize(STA,0)) AC
			SetScale/P x,DimOffset(STA,0),DimDelta(STA,0),WaveUnits(STA, 0 )+"^2",AC
			Redimension/D AC			
			// add information into note on duration and variance
			Note/K/NOCR AC,"duration:"+num2str(duration) +"\r"
			Note/NOCR AC,"variance:"+num2str(SD^2) +"\r"
			// try to obtain firing rate from corresponding _ST wave
			IName	= StringFromList(k, Cur_List,";") 
			STRING STname=ReplaceString("_I", IName, "_ST")
			IF (Exists(STname)==1)
				WAVE st=$STname
				Note/NOCR AC,"spikerate:"+num2str(DimSize(st,0)/Duration) +"\r"
			ENDIF
			
			IF (FirstTime)
				Duplicate/O STA, $(currentDFPath+"AC_avg")
				WAVE	Avg=$(currentDFPath+"AC_avg")
				Avg=0
				FirstTime=0
				Note/K Avg
				Note/NOCR Avg, currentDFPath+"\r"
				Note/NOCR Avg,FolderName+"\r"
			ENDIF
			Avg+=AC/SD^2*duration
			Variance_avg+=SD^2*duration
			totalSTANumber+=1
			totduration+=duration

			
			Note/NOCR Avg,"\t"+ StringFromList(k, Cur_List,";") +"\r"
		ENDFOR		// loop across ACs in one folder
	ENDFOR			// loop across folders
	Avg/=totduration
	Variance_avg/=totduration

	Note/NOCR Avg,"total#trial:"+num2str(totalSTANumber) +"\r"
	Note/NOCR Avg,"totalduration:"+num2str(totduration) +"\r"
	Note/NOCR Avg,"AVG_VARIANCE:"+num2str(Variance_avg) +"\r"
	
	Duplicate/O avg $(currentDFPath+"AC_avg_scaled")
	WAVE	Scaled=$(currentDFPath+"AC_avg_scaled")
	Scaled*=Variance_avg
	SetScale d 0,0,"A^2", Scaled

	SetDataFolder rootFolderRf

END		// Make and avg AC

FUNCTION MakeAC([FirstPoint,LastPoint])
VARIABLE		FirstPoint,LastPoint
// for Omers data this needs to be limited to the relevant part of the current
// points 14532 to 915469
// so for Omer this has to be called
// MakeAndAvgAC(FirstPoint=14532,LastPoint=915469)

IF (ParamIsDefault(FirstPoint))
	FirstPoint = 0
ENDIF
IF (ParamIsDefault(LastPoint))
	LastPoint = inf
ENDIF

STRING		Cur_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
		cellFolderList = ReplaceString(",", cellFolderList, ";")
		cellFolderList = RemoveFromList("Packages;IGNORE;", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK

STRING		IName, STAName, ACName
VARIABLE		k, nCurs, SD, duration, totduration=0
	
	FOR (i=0;i<=nSubs;i+=1)
		SetDataFolder rootFolderRf
		IF (nSubs > 0)	// if there are any subfolders
			FolderName=StringFromList(i, cellFolderList,";")
			SetDataFolder $(FolderName) 
		ENDIF
		Cur_List=WaveList("*_I", ";", "" )
		nCurs=ItemsInList(Cur_List,";")
		FOR (k=0; k< nCurs; k+=1)
			IName	= StringFromList(k, Cur_List,";") 
			STAname=ReplaceString("_I", IName, "_STA")
			ACName=IName+"_AC"

			WAVE	STA=$STAName
			WAVE	Cur=$IName
			// get standard deviation to scale AC for averaging
			// also gets written in wavenote
			WAVESTATS/R=[FirstPoint, min(DimSize(Cur,0)-1,LastPoint)]/Q Cur
			SD=V_sdev
			duration=V_npnts*Dimdelta(Cur,0)
			
			Duplicate/R=[FirstPoint, min(DimSize(Cur,0)-1,LastPoint)]/O Cur,$ACName, W_CurDummy
			// if FirstPOint >0, the Dimoffset is now not corresponding to the 
			// dimOffset of the current wave. 
			// (normally it should be zero, but now it is pnt2x(AC, FirstPoint )
			// This needs to be corrected 
			WAVE	AC=$ACName
			SetScale/P x,DimOffset(Cur,0),DimDelta(Cur,0),"", AC, W_CurDummy

			Correlate/NODC W_CurDummy, AC
			// normalize for number of points in the original current wave
			AC/=V_npnts
// AC is a very long wave,  retain only the same region that is present in the STA (i.e. STAwidth)
			// find the point that corresponds to the same time as the first point in the STA
			VARIABLE	firstP=x2pnt(AC, DimOffset(STA,0) )
			IF (firstP >0)
				DeletePoints  0, firstP, AC
			ENDIF
			Redimension/N=(DimSize(STA,0)) AC
			SetScale/P x,DimOffset(STA,0),DimDelta(STA,0),WaveUnits(STA, 0 )+"^2",AC
			Redimension/D AC			
			// add information into note on duration and variance
			Note/K/NOCR AC,"duration:"+num2str(duration) +"\r"
			Note/NOCR AC,"variance:"+num2str(SD^2) +"\r"
			// try to obtain firing rate from corresponding _ST wave
			IName	= StringFromList(k, Cur_List,";") 
			STRING STname=ReplaceString("_I", IName, "_ST")
			IF (Exists(STname)==1)
				WAVE st=$STname
				Note/NOCR AC,"spikerate:"+num2str(DimSize(st,0)/Duration) +"\r"
			ENDIF
			
		ENDFOR		// loop across ACs in one folder
	ENDFOR			// loop across folders

	SetDataFolder rootFolderRf

END		// Make AC


//
//			SplitBeforeFFT(AC, 0)
//			WAVE AC=$(ACName+"_splt")

//FUNCTION AvgAC(DoPlot)
//VARIABLE DoPlot
//
//STRING		AC_List
//STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
//		cellFolderList = ReplaceString(",", cellFolderList, ";")
//VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
//DFREF		rootFolderRf=GetDataFolderDFR()
//STRING		currentDFPath
//currentDFPath = GetDataFolder(1)		// OK
//VARIABLE	k, nACs, var
//STRING		Notiz, var_String
//
//VARIABLE	FirstTime=1, totalACNumber=0, Variance_avg=0, flagLocal=0
//
//// catch case that there are no subfolders and that onyl in the current folder the average is created
//	IF (nSubs == 0)
//		nSubs=1
//		flagLocal=1
//	ENDIF	
//
//// go through all subfolders
//	IF (DoPlot)
//			Display as "All ACs"
//	ENDIF
//	FOR (i=0;i<nSubs;i+=1)
//		SetDataFolder rootFolderRf
//		IF (flagLocal==0)
//			FolderName=StringFromList(i, cellFolderList,";")
//			SetDataFolder $(FolderName) 
//		ELSE // local, no subfolders
//			FolderName= GetDataFolder(1)
//		ENDIF
//
//		AC_List=WaveList("*_AC", ";", "" )
//		nACs=ItemsInList(AC_List,";")
//		IF(Exists( currentDFPath+"AC_avg")==1)
//			WAVE	Avg=$(currentDFPath+"AC_avg")
//			Note/NOCR Avg,Foldername+"\r"
//		ENDIF
//		FOR (k=0; k< nACs; k+=1)
//			WAVE	AC=$(StringFromList(k, AC_List,";") )
//			Notiz=Note(AC)
//			var_String=StringByKey("variance", Notiz,":" ," ")
//			var=str2num(var_String)
//			IF (FirstTime)
//				Duplicate/O AC, $(currentDFPath+"AC_avg")
//				WAVE	Avg=$(currentDFPath+"AC_avg")
//				Avg=0
//				FirstTime=0
//				Variance_avg=0
//				Note/K Avg
//				Note/NOCR Avg, currentDFPath+"\r"
//				Note/NOCR Avg,FolderName+"\r"
//			ENDIF
//			Avg+=AC/var
//			IF (DoPlot)
//				AppendToGraph AC
//			ENDIF
//			totalACNumber+=1
//			Variance_avg+=var
//			Note/NOCR Avg,"\t"+ StringFromList(k, AC_List,";") +"\r"
//		ENDFOR		// loop across ACs in one folder
//	ENDFOR			// loop across folders
//	Avg/=totalACNumber
//	Variance_avg/=totalACNumber
//	IF (flagLocal==0)
//		VARIABLE/G avgInputVariance=Variance_avg
//		
//	ENDIF
//	
//	Note/NOCR Avg,"AVG_VARIANCE:"+num2str(Variance_avg) +"\r"
//	
//	Duplicate/O avg $(currentDFPath+"AC_avg_scaled")
//	WAVE	Scaled=$(currentDFPath+"AC_avg_scaled")
//	Scaled*=Variance_avg
//	SetScale d 0,0,"A^2", Scaled
//
//	
//	SetDataFolder rootFolderRf
//
//END		// avg ACs


//##########################################################################################
//##########################################################################################
//#############                                              ###############################
//#############    vector strength and data stratification   ###############################
//#############                                              ###############################
//##########################################################################################
//##########################################################################################

FUNCTION CreatePhaseWave(V_Wave[,FirstPoint,LastPoint])
//This function returns a wave with the phase which have the x and y values of the phase lock vectors
//The V_Wave should be produce from an Input with a sinudoidal wave of 0 phase. 
//The program assumes that for each V_Wave there are equivalent I_wave(input), ISI_wave and ST_Wave. 
//If the value of DoDisplay is 1, the function will graph the vectors with Magnitude 1  
//and phase dependent of the lock with the sin wave of the input. 
//The color of the points is set depending of the Interspike interval(ISI) in a yellow-to-black scale
//being yellow the largest ISI 
//The difference between this function and the original from Andreas is that it creates another 
//wave with the phase of every vector
 
	WAVE			V_Wave
	VARIABLE		FirstPoint,LastPoint
	
	
		IF (ParamIsDefault(FirstPoint))
			FirstPoint = 0
		ENDIF

		VARIABLE	nInAna=DimSize(V_Wave,0)-FirstPoint
		IF (!ParamIsDefault(LastPoint))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
													
			nInAna -= DimSize(V_Wave,0)-1-min(DimSize(V_Wave,0)-1,LastPoint)
		ELSE
			LastPoint=DimSize(V_Wave,0)-1
		ENDIF


	STRING 	IName, STName, VSName, SPhaseName, VoltName=NameOfWave(V_Wave)
	VARIABLE	VoltNameLength=strlen(VoltName)
	
	IF 	(StringMatch(VoltName[VoltNameLength-2,VoltNameLength-1],"_V" ))	// last to characters
																// in VoltName are "_V"
		STName=VoltName[0,VoltNameLength-3]+"_ST"		
		SPhaseName=VoltName[0,VoltNameLength-3]+"_SPhs"	// spike phase													
		
		WAVE ST=$STName
		IF (!(Exists(STName)==1))
			// spike time wave does not yet exist
			Abort "First get Spike time wave, then call phaselock visualization"
		ENDIF
		IName=VoltName[0,VoltNameLength-3]+"_I"														
		WAVE I_Wave=$IName
		IF (!(Exists(IName)==1))
			// spike time wave does not yet exist
			Abort "Cannot find input for "+VoltName
		ENDIF
				
	ELSE
		Abort "Cannot decode names as voltage trace does not end in 'V'"
	ENDIF
	
	Duplicate/O ST, $SPhaseName
	WAVE	SPhase=$SPhaseName
	SPhase=0
	
	VARIABLE customFreq
	// get sine frequency that predominates in injected current --> max of FFT
	VARIABLE	L_Input=DimSize(I_wave,0)
	FFT/RP=[0,2*floor(L_Input/2)-1]/OUT=3/DEST=Temp_FFT  I_wave
	// input needs to have a correct x-scaling
	
	REDIMENSION/N=(150000) Temp_FFT
	// only  look for frequencies larger than 0.25 Hz
	VARIABLE firstFreqPoint=round(0.25/DimDelta(Temp_FFT,0))
	
	WAVESTATS/Q/M=1/R=[firstFreqPoint,DimSize(Temp_FFT,0)-1] Temp_FFT
	VARIABLE	freq=round(V_maxLoc)
	VARIABLE	p_freq=V_maxRowLoc
	
	// estimate sine amplitude from FFT
	VARIABLE Est_sine_amp = Temp_FFT[p_freq]/(floor(L_Input/2)) // normalize by 
																					// dividing by number of input points 
																					// and multiply by 2 (positive/neg freq.
	
	IF (freq==0)
		// if there was no sine signal embedded, the output is f=0, which we should not 
		// use to compute phases 
			DoAlert 0, "Ignored the sine freq in wave "+IName+".\r It was detected to be "+num2str(freq)+" Hz."
			customFreq=500
			DoAlert 1, "can you provide a sine frequency?"
			IF (V_flag==1)
				Prompt customFreq,"Enter sine frequency : "
				DoPrompt "Custom choice", customFreq
				if (V_Flag)
					freq=NaN
				else
					freq=	customFreq
				endif
				
			ELSE
				freq=NaN
			ENDIF
		ELSE
		// check if that one frequency component does stand out from the surrounding
		
		VARIABLE nextHighest=max(WaveMax(Temp_FFT,pnt2x(Temp_FFT,max(0,p_freq-6)),pnt2x(Temp_FFT,max(0,p_freq-2))), WaveMax(Temp_FFT,pnt2x(Temp_FFT,p_freq+2),pnt2x(Temp_FFT,p_freq+6)))
		//Also compute the standard deviation of the part of the wave surrounding that high point 
		// Excluding the point itself
		
		Duplicate/O/R=[p_freq-6,p_freq+6] Temp_FFT, TEmp_part
		DeletePoints 6,1,TEmp_part
		VARIABLE SD=sqrt(Variance(TEmp_part))
		// Use a combination of the two measures to exclude 
		// results from noisy FFTs
		IF ( (nextHighest/Temp_FFT[p_freq] > 0.5 ) || (Temp_FFT[p_freq]/SD < 3 ) )	// might need some adjustments in the two 
																					// hard coded limits 
		
			DoAlert 0, "Ignored the sine freq in wave "+IName+".\r It was detected to be "+num2str(freq)+" Hz."
			customFreq=500
			DoAlert 1, "can you provide a sine frequency?"
			IF (V_flag==1)
				Prompt customFreq,"Enter sine frequency : "
				DoPrompt "Custom choice", customFreq
				if (V_Flag)
					freq=NaN
				else
					freq=	customFreq
				endif
				
			ELSE
				freq=NaN
			ENDIF
		ENDIF	
	ENDIF 
	printf "The detected sine frequency is %g\r", freq
		
	KillWaves Temp_FFT, Temp_Part

	SPhase[]=atan2(sin(2*Pi*Freq*ST[p]),cos(2*Pi*Freq*ST[p]))
	
	VARIABLE/C VSandPhase=VSandPhasefromPhases(SPhase)
	
	//Delete the note 
	// write the Frequency of the sinusoidal wave in wave note
	// write total vector strength in the note 
	VARIABLE  spikerate=DimSize(ST,0)/nInAna/DimDelta(V_Wave,0)
	Note /K SPhase,"spikerate:"+num2str(spikerate)
	Note SPhase, "Sine frequency (Hz): " +num2str(Freq)
	Note SPhase, "Sine Amp ("+WaveUnits(I_Wave,-1)+"): " +num2str(Est_sine_amp)
	Note SPhase, "Vector Strength: " +num2str(real(VSandPhase))
	Note SPhase, "Avg Phase: " +num2str(imag(VSandPhase))
	Note SPhase, "Gain (Hz/"+WaveUnits(I_Wave,-1)+"): " +num2str(2*real(VSandPhase)/Est_sine_amp*spikerate)

END  // CreatePhaseWave


Function BootstrapVS(PhaseWave[,QuantileWave, quiet, doComplex])
WAVE		PhaseWave
WAVE		QuantileWave
VARIABLE	quiet, doComplex
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
	// In order to test properly, the bootstrap has to accoutn for magnitude AND phase of the bootstrap sample
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
END

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
		
		realPart	=	cos(W_Resampled)
		cmplxPart	= 	sin(W_Resampled)
		
		VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
		VARIABLE/C	VSandPhase	= r2polar( cmplx(rAVG,cAVG) )
		
		KillWaves/Z realPart, cmplxPart, W_Resampled
	
		Return VSandPhase

	
END //ResampleAndGetVSAndPhase


Threadsafe FUNCTION/C VSandPhasefromPhases(PhaseWave)
WAVE	PhaseWave

	Duplicate/FREE PhaseWave, realPart, cmplxPart
	Redimension/D realPart, cmplxPart
	
	realPart		=	cos(PhaseWave)
	cmplxPart	= 	sin(PhaseWave)
	
	VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
	VARIABLE	VS	= sqrt(rAVG^2+cAVG^2), phsAVG = atan2(cAVG,rAVG)
	
	KillWaves/Z realPart, cmplxPart

	VARIABLE/C returner=cmplx(VS,phsAVG)
	Return returner
	
END // VSandPhasefromPhases


Threadsafe FUNCTION VSfromPhases(PhaseWave)
WAVE	PhaseWave

	Duplicate/FREE PhaseWave, realPart, cmplxPart
	Redimension/D realPart, cmplxPart

	
	realPart		=	cos(PhaseWave)
	cmplxPart	= 	sin(PhaseWave)
	
	VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
	VARIABLE	VS	= sqrt(rAVG^2+cAVG^2)
	
	KillWaves/Z realPart, cmplxPart

	Return VS
END // VSfromPhases

FUNCTION/WAVE WaveSubsetByCriteria(sourceWave, CriteriumWave0, lowCrit0, upCrit0[, CriteriumWave1, lowCrit1, upCrit1, Logic] )
WAVE		sourceWave, CriteriumWave0, CriteriumWave1
								// suffix(es) of the waves containing the criterium for stratification, e.g. "mxRt" or "RPD"
VARIABLE		lowCrit0, upCrit0,lowCrit1, upCrit1					
								// upper and lower bound(aries) (exclusive) for the value in the criterum wave
STRING			Logic 			// "and" or "or" 
// takes a single wave and returns a variant of the wave in which only the entries at certain indicies are kept
// those indicies are identified by entries in one or two criteria waves:
// for eligible indicies, the respective entries in the criteria waves fulfill the criteria

// the function is multidimension-aware (for 2D)
// BUT IT ALWAYS RETURNS A 1D WAVE!
// that is due to the internal logic here, where elements are selected individually, not entire rows or columns
	// Check if Waves exist 
	IF (!WaveExists(sourceWave))
		Abort "SourceWave " + NameOfWave(sourceWave)+" not found"
	ELSE
		Wavestats/M=1/Q SourceWave
		VARIABLE SourcePnts=V_npnts+V_numNans+V_numINFs	// it should be possible to deal with NaNs and INFs without excluding 
																	// ALL waves that have them
																	// therefore I include them here
	ENDIF
	IF (!WaveExists(CriteriumWave0))
		Abort "CriteriumWave " + NameOfWave(CriteriumWave0)+" not found"
	ELSE
		Wavestats/M=1/Q CriteriumWave0
		VARIABLE Crit0Pnts=V_npnts+V_numNans+V_numINFs
	ENDIF
	 
	// Quick health check about wave sizes
	IF (SourcePnts != Crit0Pnts)
		Abort "Sourcewave "+NameOfWave(SourceWave)+" and criterium wave "+NameOfWave(CriteriumWave0)+" have different size"
	ENDIF
	IF (!ParamIsDefault(CriteriumWave1))
		// Check if Waves exist 
		IF (!WaveExists(CriteriumWave1))
			Abort "CriteriumWave " + NameOfWave(CriteriumWave1)+" not found"
		ELSE
			Wavestats/M=1/Q CriteriumWave1
			VARIABLE Crit1Pnts=V_npnts+V_numNans+V_numINFs
		ENDIF
		// Quick health check about wave sizes
	
		IF (SourcePnts != Crit1Pnts)
			Abort "Sourcewave "+NameOfWave(SourceWave)+" and criterium wave "+NameOfWave(CriteriumWave1)+" have different size"
		ENDIF
		// Quick health check regarding consistent number of criterium values
		IF ( (ParamIsDefault(lowCrit1)) || (ParamIsDefault(upCrit1)) )
			Abort "lower or upper bound for criterium 1 was not defined"
		ENDIF
		IF (ParamIsDefault(Logic))
			Abort "You are required to define the logic of combining the two criteria"
		ENDIF
	ENDIF
	
	// prepare a result wave
	DFREF		rootFolderRf=GetDataFolderDFR()
	
	// create target wave in free data folder
	
	SetDataFolder NewFreeDataFolder()
	
	MAKE/D/O/N=0 Target
		
	SetDataFolder rootFolderRf
		
	// create waves that hold entries reflecting fullfilling of criteria
	
	MATRIXOP/O Indicies0=within(CriteriumWave0,lowCrit0,upCrit0)	// Returns an array of the same dimensions as w  
																				// with the value 1 where the corresponding element of w 
																				// is between low and high (low <= w[i][j] < high) 
																				// and zero otherwise.
	Redimension/N=(SourcePnts)	 Indicies0	// make a linear wave (in case it was not 1D)
	
																	
	IF (!ParamIsDefault(CriteriumWave1))
		MATRIXOP/O Indicies1=within(CriteriumWave1,lowCrit1,upCrit1)	
		Redimension/N=(SourcePnts)	 Indicies1	// make a linear wave (in case it was not 1D)
		strswitch(Logic)	// string switch
			case "and":	// execute if case matches expression
				Indicies0*=Indicies1							
				break		// exit from switch
			case "or":	// execute if case matches expression
				MATRIXOP/O Indicies0=(Indicies0 || Indicies1)
				break
			case "AND":	// execute if case matches expression
				Indicies0*=Indicies1							
				break		// exit from switch
			case "OR":	// execute if case matches expression
				MATRIXOP/O Indicies0=(Indicies0 || Indicies1)
				break
			default:			// optional default expression executed
				Abort "Give a Logic as 'and' or 'or'"
		endswitch
	ENDIF
	Redimension/L/U Indicies0	// change data type from INT8 to unsigned INT64 (indecies)
	KillWaves Indicies1
	
	// get number of non-zero entries
	VARIABLE	nHits=sum(Indicies0)
	IF (nHits<=0)
	
	// this version should work for most
	Return target  // still an empty wave
	
	// sometiems at least one entry is required
	// returning a NAN might be adequate
	Target[0]={NaN}

		Return Target // contains one "NaN"
	ENDIF
	// load index into indicies0
	Indicies0[]*=p+1
	
	// remove entries that are negative:

	Sort Indicies0 Indicies0
	DeletePoints 0,SourcePnts-nHits,Indicies0
	Indicies0-=1
	// now it cointains the indicies of all items for which all criteria are fullfilled

	Redimension/N=(nHits) target
	target[]=SourceWave[Indicies0[p]]	
	
	Return target

END // WaveSubsetByCriteria





FUNCTION/WAVE CollectDataByNoteAndCriterium(Suffix[, NoteKey,low4NoteVal, up4NoteVal, CriteriumSuffix, lowCrit, upCrit, CriteriumSuffix1, lowCrit1, upCrit1	] )
STRING		Suffix
STRING		NoteKey		// key string that allows extraction of value from wave note
								// examples: "Sine frequency (Hz)" for the SPhs waves
								// or "average peak rate of rise (PRR)" for the ST waves
VARIABLE		low4NoteVal, up4NoteVal	
								// upper and lower bound (inclusive) for the numerical value extracted from the note

STRING		CriteriumSuffix,CriteriumSuffix1
								// suffix(es) of the waves containing the criterium for stratification, e.g. "mxRt" or "RPD"
VARIABLE		lowCrit, upCrit,lowCrit1, upCrit1					
								// upper and lower bound(aries) (exclusive) for the value in the criterum wave
								
// version 03/10/2018 _aNNe_ added the option to have a _second_ criterium, i.e. a second list of waves applying a second,
// separate criterium to select a subset of the data
								
	VARIABLE	DoStratify = 0		// binary switch to indicate whether criterium 0 is set (DoStratify+=2^0),criterium 1 is set (DoStratify+=2^1), etc.
	VARIABLE	DoNote = 0

	IF (!ParamIsDefault(NoteKey))
		DoNote = 1
		IF ( (ParamIsDefault(low4NoteVal)) || (ParamIsDefault(up4NoteVal)) )
			Abort "lower or upper bound for value from wave note was not defined"
		ENDIF
	ELSE
		NoteKey=""
		low4NoteVal= NaN
		up4NoteVal = NaN	
	ENDIF

	IF (!ParamIsDefault(CriteriumSuffix))
		DoStratify += 2^0
		IF ( (ParamIsDefault(lowCrit)) || (ParamIsDefault(upCrit)) )
			Abort "lower or upper bound for criterium was not defined"
		ENDIF
	ENDIF
	
	IF (!ParamIsDefault(CriteriumSuffix1))
		DoStratify += 2^1
		IF ( (ParamIsDefault(lowCrit1)) || (ParamIsDefault(upCrit1)) )
			Abort "lower or upper bound for criterium 2 was not defined"
		ENDIF
	ENDIF
	
	

	
	STRING		Wave_List
	STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
	cellFolderList = ReplaceString(",", cellFolderList, ";")
	
	VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
	DFREF		rootFolderRf=GetDataFolderDFR()
	STRING		rootDFPath
	rootDFPath = GetDataFolder(1)		// OK
	VARIABLE	k, nWaves, noteVal
	
	STRING		Notiz, Note_String, CriteriumWaveName, CriteriumWaveName1
	
	VARIABLE	FirstTime=1, totalVSNumber=0
	
	// create target wave in free data folder
	
	SetDataFolder NewFreeDataFolder()
	
	MAKE/D/O/N=0 Target
		
	SetDataFolder rootFolderRf
		
	
	// go through all subfolders
	


	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(i, cellFolderList,";")
		SetDataFolder $(FolderName) 
		Wave_List=WaveList("*_"+Suffix, ";", "" )
		nWaves=ItemsInList(Wave_List,";")
		// do not want to collect ALL vectorStrength waves, only those that belong to the 
		// targeted sine frequency 
		FOR (k=0; k< nWaves; k+=1)
			WAVE	Candidate=$(StringFromList(k, Wave_List,";") )

			Notiz=Note(Candidate)
			Note_String=StringByKey(NoteKey, Notiz,":" ,"\r")
			noteVal=str2num(Note_String)
			IF ( ( (noteVal>=low4NoteVal) && (noteVal<= up4NoteVal) ) || (!DoNote ) )// this wave we want
				IF (FirstTime)
					Note/K Target
					IF (!DoStratify)
						Concatenate/NP  {Candidate}, Target
						Note Target, "No Stratification\r"
					ELSEIF (DoStratify == 1)
						// finding corresponding criterium wave
						CriteriumWaveName=RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit)}, Target
						ENDIF
						Note/NOCR Target, "Stratification by "+ CriteriumSuffix+"\r"
						Note/NOCR Target, "Limits: "+num2str(lowCrit)+ ", "+ num2str(upCrit)+"\r\r"
					ELSEIF (DoStratify == 2)
						// finding corresponding criterium wave
						CriteriumWaveName=RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix1
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Second Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit1, upCrit1)}, Target
						ENDIF
						Note/NOCR Target, "Stratification by "+ CriteriumSuffix1+"\r"
						Note/NOCR Target, "Limits: "+num2str(lowCrit1)+ ", "+ num2str(upCrit1)+"\r\r"
					ELSEIF (DoStratify == 3)
						// finding two corresponding criterium waves
						CriteriumWaveName		= RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix
						CriteriumWaveName1	= RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix1
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"First Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSEIF (!Exists(CriteriumWaveName1)==1)
							DoAlert 0,"Second Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							WAVE Crit1=$CriteriumWaveName1
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit,CriteriumData1=crit1, lowCriterium1=lowCrit1, highCriterium1=upCrit1, logic="AND")}, Target

						ENDIF
						

						Note/NOCR Target, "Stratification by "+ CriteriumSuffix+"\r"
						Note/NOCR Target, "Limits: "+num2str(lowCrit)+ ", "+ num2str(upCrit)+"\r\r"
						Note/NOCR Target, "Stratification by "+ CriteriumSuffix1+"\r"
						Note/NOCR Target, "Limits: "+num2str(lowCrit1)+ ", "+ num2str(upCrit1)+"\r\r"
					ENDIF	// Stratify?

					FirstTime=0
					
				ELSE	// not First Time
					IF (!DoStratify)
						Concatenate/NP  {Candidate}, Target
					ELSEIF (DoStratify == 1)
						// finding corresponding criterium wave
						CriteriumWaveName=RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit)}, Target
						ENDIF
					ELSEIF (DoStratify == 2)
						// finding corresponding criterium wave
						CriteriumWaveName=RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix1
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit1, upCrit1)}, Target
						ENDIF
					ELSEIF (DoStratify == 3)
						// finding two corresponding criterium waves
						CriteriumWaveName		= RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix
						CriteriumWaveName1	= RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix1
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSEIF (!Exists(CriteriumWaveName1)==1)
							DoAlert 0,"Second Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							WAVE Crit1=$CriteriumWaveName1
							Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit,CriteriumData1=crit1, lowCriterium1=lowCrit1, highCriterium1=upCrit1, logic="AND")}, Target
						ENDIF
ENDIF // stratify?

				ENDIF
				Note/NOCR Target, GetWavesDataFolder(Candidate,2)+"\r"				
			ENDIF
		ENDFOR		// loop across target waves in one folder
	ENDFOR			// loop across folders

	SetDataFolder rootFolderRf
	
	return target

END		// CollectDataByNoteAndCriterium


FUNCTION/S ListDataByNote(Suffix, NoteKey,low4NoteVal, up4NoteVal)
STRING		Suffix
STRING		NoteKey		// key string that allows extraction of value from wave note
								// examples: "Sine frequency (Hz)" for the SPhs waves
								// or "average peak rate of rise (PRR)" for the ST waves
VARIABLE		low4NoteVal, up4NoteVal	
								// upper and lower bound (inclusive) for the numerical value extracted from the note

// Similar to CollectData
// here, all those waves with a certain name suffix are returned
// if they have, in their wave note, a pair of key and value and the value is between low and up
// the return is a LIST of full wavenames (incl. folders)
								
	VARIABLE	DoNote = 1


	
	STRING		Wave_List
	STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
	cellFolderList = ReplaceString(",", cellFolderList, ";")
	
	VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
	DFREF		rootFolderRf=GetDataFolderDFR()
	STRING		rootDFPath
	rootDFPath = GetDataFolder(1)		// OK
	VARIABLE	k, nWaves, noteVal
	
	STRING		Notiz, Note_String
	STRING		OutString=""
	
//	VARIABLE	FirstTime=1, totalVSNumber=0
	
	// create target wave in free data folder
	SetDataFolder rootFolderRf
		
	
	// go through all subfolders
	


	FOR (i=0;i<nSubs;i+=1)
		SetDataFolder rootFolderRf
		FolderName=StringFromList(i, cellFolderList,";")
		SetDataFolder $(FolderName) 
		Wave_List=WaveList("*_"+Suffix, ";", "" )
		nWaves=ItemsInList(Wave_List,";")
		// do not want to collect ALL candidate waves, only those that belong to the 
		// targeted sine frequency 
		FOR (k=0; k< nWaves; k+=1)
			WAVE	Candidate=$(StringFromList(k, Wave_List,";") )

			Notiz=Note(Candidate)
			Note_String=StringByKey(NoteKey, Notiz,":" ,"\r")
			noteVal=str2num(Note_String)
			IF ( ( (noteVal>=low4NoteVal) && (noteVal<= up4NoteVal) ) || (!DoNote ) )// this wave we want
				OutString+= GetWavesDataFolder(Candidate,2)	+";"
			ENDIF

		ENDFOR		// loop across target waves in one folder
	ENDFOR			// loop across folders
	SetDataFolder rootFolderRf
	
	return OutString

END		// ListDataByNoteAndCriterium
