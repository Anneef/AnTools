#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include "ANtools_extended"
#include "SpikeAnalysis_aNNe"
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # # 	                       # # # # # # # # # # # # 
// # # # # # # # 	Dynamic gain WRAPPERS   # # # # # # # # # # # # 
// # # # # # # # 	                       # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
FUNCTION GainCalculationFunc()
	STRING Prefix="A"
	STRING CMDSTR
	DFREF  topDf=GetDataFolderDFR()
	
	STRING	rootStr=GetDataFolder(1,topDF)
	
	STRING		FolderListe=ReplaceString(",", StringByKey("FOLDERS", DataFolderDir(1,topDf)+",",  ":", ";") ,";") 
	FolderListe= RemoveFromList("Packages;IGNORE;", FolderListe, ";")
	FolderListe=SortList(FolderListe,";",16)
	
	VARIABLE	ww,runde,numFolders=ItemsInList(FolderListe)
	
	STRING 	ListOfVWaves,ListOfIWaves,ListOfSTWaves="",NewName
	VARIABLE	NumEntries=0
	
	MAKE/O/WAVE/N=(0) AllVWRs
	// Spike times in subfolders
	FOR (runde=0; runde< numFolders; runde+=1)
		SetDataFolder $(StringFromList(runde,FolderListe,";"))
		ListOfVWaves=WaveList(Prefix+"*_V",";","")
		Wave/WAVE Vwr = ListToWaveRefWave(ListOfVWaves, 0)
		Concatenate {Vwr},AllVWRs
		NumEntries+=ItemsInList(ListOfVWaves, ";")
		SetDataFolder topDf	
	ENDFOR
	Make/DF/N=(NumEntries) AllSTs

	STRING WAVESInDF
	VARIABLE nwaves

	Multithread AllSTs[]=ReturnSpikeTimeWave(AllVWRs[p],V_Threshold=0, MinISI=0.002)
	// using multithreading with folders to be returned
	// now write STs into existing datafolders (otherwise they are immediately deleted
	FOR (ww=0; ww<NumEntries; ww++)
		DFREF df= AllSTs[ww]
		// find out how many waves are there
		WAVESInDF=ReplaceString(",",StringByKey("WAVES",DataFolderDir(2,df),":",";"),";") // gives semicolon speparated list of wave names
		nwaves=ItemsInList(WAVESInDF)

		WAVE V=AllVWRs[ww]
		FOR(runde=0;runde<nwaves;runde++)
			NewName=StringFromList(runde,WAVESInDF)
			SetDataFolder	df
			WAVE XX=$(NewName)
			SetDataFolder GetWavesDataFolderDFR(V)
			Duplicate/O XX, $(NewName)
		ENDFOR
	ENDFOR
	NumEntries=0
	SetDataFolder topDf
	MAKE/O/WAVE/N=(0) AllIWRs,AllSTWRs // will keep all wave references to all current and ST waves

	// calculate spike triggered average inputs
	FOR (runde=0; runde< numFolders; runde+=1)
		SetDataFolder $(StringFromList(runde,FolderListe,";"))
		ListOfIWaves=WaveList(Prefix+"*_I",";","")
		// check whether corresponding Spike times exist
		ListOfSTWaves=""
		FOR (ww=0;ww < ItemsInList(ListOfIWaves, ";"); ww++)
			ListOfSTWaves+=RemoveEnding( StringFromList(ww, ListOfIWaves,";") )+"ST;"
		ENDFOR
		Wave/WAVE STwr = ListToWaveRefWave(ListOfSTWaves, 0)	// throws error, if they don't exist
		Wave/WAVE Iwr = ListToWaveRefWave(ListOfIWaves, 0)	// throws error, if they don't exist
		Concatenate {Iwr},AllIWRs
		Concatenate {STwr},AllSTWRs
		NumEntries+=ItemsInList(ListOfIWaves, ";")
		SetDataFolder topDf	
	ENDFOR
	MAKE/WAVE/N=(NumEntries) AllSTAs

	Multithread AllSTAs[]=STAfromAnalogueWORKER(AllIWRs[p], AllSTWRs[p], 1,0) // returns 0 if ok, otherwise something negative
	// using multithreading with waveReferences to be returned
	// now write STAs into existing datafolders (otherwise they are immediately deleted)
	FOR (ww=0;ww < NumEntries; ww++)
		WAVE XX=AllSTAs[ww]
		IF (DimSize(xx,0)>1)
			NewName=NameOfWave(XX)
			DFREF df= GetWavesDataFolderDFR(AllIWRs[ww])
			SetDataFolder df
			Duplicate/O XX,$NewName	
		ELSE 
			print "Error in calculating STA"
		ENDIF	
	ENDFOR
	
	SetDataFolder topDf	
	MakeAC()
	// creates autocorrelation traces for each trial in each subfolder
	// also creates waves of wave references AllSTAWRs and AllACWRs
	WAVE/WAVE AllSTAWRs
	WAVE/WAVE AllACWRs
	
	STRING FolderSTAs=WaveRefWaveToList(AllSTAWRs,0)
	STRING FolderACs=WaveRefWaveToList(AllACWRs,0)
	AvgACFromList(FolderACs,"AC_avg_scaled")
	AvgSTAFromList(FolderSTAs,"STA_avg_scaled")
	// last two result in waves
	WAVE AC_avg_scaled
	WAVE STA_avg_scaled
	// calculate overall gain for data from all subfolders
	SplitBeforeFFT(AC_avg_scaled,0)
	// results in a wave
	WAVE AC_avg_scaled_splt
	WAVESTATS/Q STA_avg_scaled
	SplitBeforeFFT(STA_avg_scaled,V_maxloc)
	// results in a wave
	WAVE STA_avg_scaled_splt
	FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt
	FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt
	Duplicate/C/O STA_avg_scaled_splt_FFT, Gain_avg_scaled;
	Gain_avg_scaled/=AC_avg_scaled_splt_FFT
	Gain_avg_scaled*=cmplx(str2num(StringByKey("total#spikes", note(STA_avg_scaled),":" ,"\r"))/str2num(StringByKey("totalduration", note(STA_avg_scaled),":" ,"\r")),0)
	Gain_avg_scaled= conj(Gain_avg_scaled)	// fixes the sign of the phase
	GaussFilter(Gain_avg_scaled)
	WAVE Gain_avg_scaled_MgFlt
	note Gain_avg_scaled_MgFlt,note(STA_avg_scaled)
	
	KillWaves AllVWRs,AllIWRs,AllSTWRs,AllSTAs,AllSTs, AllSTAWRs,AllACWRs

END 

MACRO GainCalculationWrapper()
// start in the folder above all cell folders
// goes through and creates the avg gain 

	STRING Prefix="A"
	STRING CMDSTR
	STRING PathToFolder=GetDataFolder(1)
//	
//	// get spike times
//	CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_V\",\" ReturnSpikeTimes(~,V_Threshold=0, MinISI=0.002) \")",Prefix)
//	XeqtInSubs(CMDSTR)
//	// calculate spike triggered average currents
//	CMDSTR=ReplaceString("??Pref??","Xeqt4WList(\"??Pref??*_I\",\"STRING/G W_name= \\\"~\\\"; W_name=W_name[0,strlen(W_name)-2]+\\\"ST\\\"; STAfromAnalogue( ~, $W_Name, 1,0)\")",Prefix)
//	XeqtInSubs(CMDSTR)
//	
//	// create autocorrelation traces for each trial in each subfolder
//	XeqtInSubs("MakeAC()")
	
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
	
	// now create avg AC and STA and gain in each subfolder
	
	CMDSTR="STRING/G FolderACs=\"\"; Xeqt4WList(\"??Pref??*_AC\",\"FolderACs+=\\\"~;\\\"\") ; AvgACFromList(FolderACs,\"AC_avg_scaled\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	XeqtInSubs(CMDSTR)
	CMDSTR="STRING/G FolderSTAs=\"\"; Xeqt4WList(\"??Pref??*_STA\",\"FolderSTAs+=\\\"~;\\\"\") ; AvgSTAFromList(FolderSTAs,\"STA_avg_scaled\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	XeqtInSubs(CMDSTR)
	
	XeqtInSubs("SplitBeforeFFT(AC_avg_scaled,0);WAVESTATS/Q STA_avg_scaled;SplitBeforeFFT(STA_avg_scaled,V_maxloc);")
	XeqtInSubs("FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt;FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt;Duplicate /O STA_avg_scaled_splt_FFT, Gain_avg_scaled;")
	XeqtInSubs("Gain_avg_scaled/=AC_avg_scaled_splt_FFT; Gain_avg_scaled= conj(Gain_avg_scaled);Gain_avg_scaled*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0) ")
	XeqtInSubs("GaussFilter(Gain_avg_scaled);	note Gain_avg_scaled_MgFlt,note(STA_avg_scaled)")
END //GainCalculationWrapper()
	
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

	VARIABLE keepSpikeTimesForBootstrap=1
	
	
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
	
	// now create avg AC and STA and gain in each subfolder
	
	CMDSTR="STRING/G FolderACs=\"\"; Xeqt4WList(\"??Pref??*_AC\",\"FolderACs+=\\\"~;\\\"\") ; AvgACFromList(FolderACs,\"AC_avg_scaled\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	XeqtInSubs(CMDSTR)
	CMDSTR="STRING/G FolderSTAs=\"\"; Xeqt4WList(\"??Pref??*_STA\",\"FolderSTAs+=\\\"~;\\\"\") ; AvgSTAFromList(FolderSTAs,\"STA_avg_scaled\")"
	CMDSTR=ReplaceString("??Pref??",CMDSTR,Prefix)
	XeqtInSubs(CMDSTR)
	
	XeqtInSubs("SplitBeforeFFT(AC_avg_scaled,0);WAVESTATS/Q STA_avg_scaled;SplitBeforeFFT(STA_avg_scaled,V_maxloc);")
	XeqtInSubs("FFT/OUT=1/DEST=AC_avg_scaled_splt_FFT AC_avg_scaled_splt;FFT/OUT=1/DEST=STA_avg_scaled_splt_FFT STA_avg_scaled_splt;Duplicate /O STA_avg_scaled_splt_FFT, Gain_avg_scaled_wC;")
	XeqtInSubs("Gain_avg_scaled_wC/=AC_avg_scaled_splt_FFT; Gain_avg_scaled_wC*=cmplx(str2num(StringByKey(\"total#spikes\", note(STA_avg_scaled),\":\" ,\"\\r\"))/str2num(StringByKey(\"totalduration\", note(STA_avg_scaled),\":\" ,\"\\r\")),0) ")
	XeqtInSubs("GaussFilter(Gain_avg_scaled_wC);	note Gain_avg_scaled_wC_MgFlt,note(STA_avg_scaled)")
END 

MACRO GainCalculationForTrialSubset(Crit0_string,Crit0_min, Crit0_max)
STRING 		Crit0_string
VARIABLE		Crit0_min, Crit0_max
// selects trials that fullfill a certain criterium
// i.e. rate during trial is between rate min and rate max
// This should work to calculate the gain within a folder 
// but just as well to calculate the gain across all subfolders
// depending from where it is started
// as of writing, the second case has not been tested

// the criteria are in spike time wave's NOTES (notes of waves ending in "_ST") 
// if the S_Crit1_suff string is empty, the second critrion is ignored

// can only be called AFTER Spike times, spike-triggered average currents and autocorrelations have been computed

// use cases
// GainCalculationForTrialSubset("average firing rate",2, 5)

	VARIABLE keepSpikeTimesForBootstrap=1
	
	
	STRING Prefix="OU"
	STRING CMDSTR
	STRING PathToFolder=GetDataFolder(1)
	

	// for every single Spike time wave, check inside note whether Criterium is met
// here edit criteria
	VARIABLE/G root:V_crit0_min=Crit0_min
	VARIABLE/G root:V_crit0_max=Crit0_max

	STRING/G	root:S_crit0_String = Crit0_String
	
	// the note that will be applied to the resulting STA and gain
	STRING NoteAppendix
	NoteAppendix=Crit0_String+" MIN:="+num2str(Crit0_min)+"\r"
	NoteAppendix+=Crit0_String+" MAX:="+num2str(Crit0_max)+"\r"


	// apply only one criterium
	// FOR instance for selecting by rate 


	// read out which spike time wave conform to criterium using CollectDataByNoteAndCriterium
	// however, don't use the returned wave directlz, just the info in the wave note
	// there is a list of full paths to each conforming spike time wave
	
	// this runs automatically through all subfolders, if there are any
	
	STRING EligibleSTs=Note(CollectDataByNoteAndCriterium("ST", NoteKey=Crit0_string,low4NoteVal=Crit0_min, up4NoteVal=Crit0_max ))
	//turn into semicolon separated listMatch(
	EligibleSTs = ReplaceString("\r", EligibleSTs,";") 
	// remove first entry (just sazing it is not further stratified (only entire waves are taken)
	EligibleSTs = RemoveListItem(0, EligibleSTs,";") 
	
	// turn into list of STA waves
	STRING 	EligibleSTAs = ReplaceString("_ST;", EligibleSTs,"_STA;") 
	
	// turn into list of STA waves
	STRING 	EligibleACs = ReplaceString("_ST;", EligibleSTs,"_I_AC;") 

	AvgACFromList(EligibleACs,"AC_avg_scaled_wC")

	AvgSTAFromList(EligibleSTAs,"STA_avg_scaled_wC")
	
	// calculate gain
	SplitBeforeFFT(AC_avg_scaled_wC,0)
	WAVESTATS/Q STA_avg_scaled_wC
	SplitBeforeFFT(STA_avg_scaled_wC,V_maxloc)
	FFT/OUT=1/DEST=AC_avg_scaled_wC_splt_FFT AC_avg_scaled_wC_splt
	FFT/OUT=1/DEST=STA_avg_scaled_wC_splt_FFT STA_avg_scaled_wC_splt

	Duplicate /O STA_avg_scaled_wC_splt_FFT, Gain_avg_scaled_wC;

	Gain_avg_scaled_wC/=AC_avg_scaled_wC_splt_FFT
	Gain_avg_scaled_wC*=cmplx(str2num(StringByKey("total#spikes", note(STA_avg_scaled_wC),":" ,"\r"))/str2num(StringByKey("totalduration", note(STA_avg_scaled_wC),":" ,"\r")),0)
	Gain_avg_scaled_wC= conj(Gain_avg_scaled_wC)
	GaussFilter(Gain_avg_scaled_wC)

	Note/K Gain_avg_scaled_wC_MgFlt
	note Gain_avg_scaled_wC_MgFlt,note(STA_avg_scaled_wC) + NoteAppendix
	note/K STA_avg_scaled_wC,note(Gain_avg_scaled_wC_MgFlt)
	IF (keepSpikeTimesForBootstrap)
		// need to make a copy of each eligible spike time wave, give it suffix "STwC"
		Xeqt4List(EligibleSTs,"",0,"Duplicate ~, ~wC")
	ENDIF

	
END //GainCalculationForTrialSubset(Crit0_string,Crit0_min, Crit0_max)


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # # 	                       # # # # # # # # # # # # 
// # # # # # # # 	STA and AC computations   # # # # # # # # # # # # 
// # # # # # # # 	                       # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


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
 
		VARIABLE	SD_exponent=1 // the exponent of SD in normalization
										// 1 means STAs are scaled by input SD
										// 2 means STA is scaled by input variance
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

			Avg+=STA/(SD^SD_exponent)*duration // possibly SD^2
			totalSTANumber+=1
			AvgSD+=(SD^SD_exponent)*duration		//  SD^2 if used above
			totspikes+=nspikes
			totduration=totduration+duration
			Note/NOCR Avg,"\t"+ GetWavesDataFolder(STA,2) +"\r"
		ENDFOR		// loop across STAs 
	Avg/=totduration
	AvgSD/=totduration


	Note/NOCR Avg,"total#trial:"+num2str(totalSTANumber) +"\r"
	Note/NOCR Avg,"total#spikes:"+num2str(totspikes) +"\r"
	Note/NOCR Avg,"totalduration:"+num2str(totduration) +"\r"
	switch (SD_exponent)
	case 1:
		Note/NOCR Avg,"AVG_SD:"+num2str(AvgSD) +"\r"
		break
	case 2:
		Note/NOCR Avg,"AVG_VAR:"+num2str(AvgSD) +"\r"
		break
	endswitch	
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


Threadsafe FUNCTION/WAVE STAfromAnalogueWaveRef(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix])
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave
	VARIABLE	FP,LP	
	STRING	SX										

		
		IF (ParamIsDefault(Suffix))	// if there is a limitation as to how many points of the
													// stimulus can be used for analysis (in case of pClamp limitations)
			SX = "STA"
		
		ELSE
			
			SX=Suffix
		ENDIF
			
		STAfromAnalogue(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar, Suffix=SX)

		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		STA_Name
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
		STA_Name=RemoveEnding(I_Name,Ending)+SX


		WAVE Dummy=$STA_Name
		
		Return Dummy
END



Threadsafe FUNCTION STAfromAnalogue(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix])
// calculates Spike triggered average of the snipletts (length 'range') of the Analogue trace
// the snipletts are centered around the entries of spike Times
// range has to have the same unit as the analogue data 
// the setscale command assumes it is "s"
// if DoSpikeTrigVar is not ZERO, the spike triggered variance is calculated
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave

		STRING		ST_Name=NameOfWave(SpikeTimes)
		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
											
		// check for existance
		IF (!Waveexists(AnalogueTrace))
			Print  "Aborting STA from analogue, could not find AnalogueTrace"+ NameOfWave(AnalogueTrace)
				Return -3
		ELSEIF (!WaveExists(SpikeTimes))
			Print "Aborting STA from analogue, couldn't find SpikeTimes"+ NameOfWave(SpikeTimes)
				Return -2
		ENDIF
		// check for having more than zero spikes
		IF (DimSize(SpikeTimes,0)<1)
			Print "Aborting STA from analogue, SpikeTime wave was empty "+ NameOfWave(SpikeTimes) + NameOfWave(AnalogueTrace)
			Return -2
		ENDIF
		
		
		VARIABLE	AT_Zero= DimOffset(AnalogueTrace,0)
		VARIABLE	AT_deltaT=DimDelta(AnalogueTrace,0)
		VARIABLE	nInAna=DimSize(AnalogueTrace,0)

		
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
		

		STRING		STA_Name
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
			IF ( (p1>=0) && (p2< 0+nInAna) ) // if input is available for the entire range needed for STA
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
			Return -1
		ENDIF


// place some info in the Wave note 
		Wavestats/Q AnalogueTrace
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
		return 0		
END // STAfromAnalogue

Threadsafe FUNCTION/WAVE STAfromAnalogueWORKER(AnalogueTrace,SpikeTimes, Range,DoSpikeTrigVar,[ Suffix])
// calculates Spike triggered average of the snipletts (length 'range') of the Analogue trace
// the snipletts are centered around the entries of spike Times
// range has to have the same unit as the analogue data 
// the setscale command assumes it is "s"
// if DoSpikeTrigVar is not ZERO, the spike triggered variance is calculated
WAVE			AnalogueTrace, SpikeTimes
VARIABLE		Range, DoSpikeTrigVar
STRING		Suffix						// the suffix placed at the end of the name of the resulting wave

		STRING		ST_Name=NameOfWave(SpikeTimes)
		STRING		I_Name=NameOfWave(AnalogueTrace)
		STRING		Ending = StringFromList(ItemsInList(I_Name,"_")-1,I_Name,"_")
		STRING 		NameOfWRONGReturn=RemoveEnding(I_Name,Ending)+"WRONGSTA"
		DFREF dfSav= GetDataFolderDFR()
	
		// Create a free data folder and set it as the current data folder
		// for multithreading to work
		SetDataFolder NewFreeDataFolder()
		
		
		
		MAKE/N=1 $NameOfWRONGReturn
		WAVE Dummy=$NameOfWRONGReturn
		Dummy=NaN												
		// check for existance
		IF (!Waveexists(AnalogueTrace))
			Print  "Aborting STA from analogue, could not find AnalogueTrace"+ NameOfWave(AnalogueTrace)
				Return Dummy
		ELSEIF (!WaveExists(SpikeTimes))
			Print "Aborting STA from analogue, couldn't find SpikeTimes"+ NameOfWave(SpikeTimes)
				Return Dummy
		ENDIF
		// check for having more than zero spikes
		IF (DimSize(SpikeTimes,0)<1)
			Print "Aborting STA from analogue, SpikeTime wave was empty "+ NameOfWave(SpikeTimes) + NameOfWave(AnalogueTrace)
			Return Dummy
		ENDIF
		
		
		VARIABLE	AT_Zero= DimOffset(AnalogueTrace,0)
		VARIABLE	AT_deltaT=DimDelta(AnalogueTrace,0)
		VARIABLE	nInAna=DimSize(AnalogueTrace,0)

		
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
		

		STRING		STA_Name
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
			IF ( (p1>=0) && (p2< 0+nInAna) ) // if input is available for the entire range needed for STA
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
			// Restore the current data folder
			SetDataFolder dfSav

			Return Dummy
		ENDIF


// place some info in the Wave note 
		Wavestats/Q AnalogueTrace
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
				// Restore the current data folder
				SetDataFolder dfSav

				return Dummy
			ENDIF
			Note/K var
			Note   var, wavnote

		ENDIF 	// if DoSpikeTrigVar
		// Restore the current data folder
		SetDataFolder dfSav
		return avg
		
END // STAfromAnalogueWORKER

ThreadSafe FUNCTION/WAVE GetACWORKER(WAVE cur,WAVE STA,STRING ACName)	
		// Worker for MakeAC			
			DFREF dfSav= GetDataFolderDFR()
	
			// Create a free data folder and set it as the current data folder
			// for multithreading to work
			SetDataFolder NewFreeDataFolder()
			
			// get standard deviation to scale AC for averaging
			VARIABLE 	var=Variance(Cur)
			VARIABLE 	npnts=DimSize(Cur,0)
			VARIABLE	duration=npnts*Dimdelta(Cur,0)
						
			Duplicate/O Cur,$ACName, W_CurDummy
			WAVE	AC=$ACName
			SetScale/P x,DimOffset(Cur,0),DimDelta(Cur,0),"", AC, W_CurDummy

			Correlate/NODC W_CurDummy, AC
			// normalize for number of points in the original current wave
			AC/=npnts
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
			Note/NOCR AC,"variance:"+num2str(var) +"\r"
			// Restore the current data folder
			SetDataFolder dfSav			
			Return AC
END //GetACWORKER


FUNCTION MakeAC()
// goes through subfolders finds current waves and calculates the AC


STRING		Cur_List
STRING		FolderName, cellFolderList=StringByKey("FOLDERS", DataFolderDir(1),":", ";")+","
		cellFolderList = ReplaceString(",", cellFolderList, ";")
		cellFolderList = RemoveFromList("Packages;IGNORE;", cellFolderList, ";")
VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
DFREF		rootFolderRf=GetDataFolderDFR()
STRING		currentDFPath
currentDFPath = GetDataFolder(1)		// OK

STRING		IName, STANameList, ACNameList
VARIABLE		k, nCurs
	
	Cur_List=""
	MAKE/O/WAVE/N=(0) AllIWRs,AllSTAWRs,AllACWRs // will keep all wave references to all current and STA waves
	ACNameList=""

	FOR (i=0;i<=nSubs;i+=1)
		SetDataFolder rootFolderRf
		IF (nSubs > 0)	// if there are any subfolders
			FolderName=StringFromList(i, cellFolderList,";")
			SetDataFolder $(FolderName) 
		ENDIF
		Cur_List=WaveList("*_I", ";", "" )
		Wave/WAVE Iwr = ListToWaveRefWave(Cur_List, 0)

		nCurs=ItemsInList(Cur_List,";")
		// create List of STA names and AC names
		STANameList=""	// STAs are only collected for each folder, then immediately converted into wave references, Those are accumulated
		FOR (k=0; k< nCurs; k+=1)
			IName	= StringFromList(k, Cur_List,";") 
			STANameList+=ReplaceString("_I", IName, "_STA;")
			ACNameList+=IName+"_AC;"		// AC names are collected for all folders
		ENDFOR // loop across ACs in one folder
		Wave/WAVE STAwr = ListToWaveRefWave(STANameList, 0)
		Concatenate {Iwr},AllIWRs
		Concatenate {STAwr},AllSTAWRs
			
	ENDFOR	// loop across folders
	VARIABLE NumEntries=DimSize(AllSTAWRs,0)
	MAKE/WAVE/N=(NumEntries) AllACs

	MultiThread AllACs[]=GetACWORKER(AllIWRs[p],AllSTAWRs[p],StringFromList(p, ACNameList,";") )	
	// using multithreading with waveReferences to be returned
	// now write ACs into existing datafolders (otherwise they are immediately deleted)
	VARIABLE ww,npnts, duration
	STRING NewName, STname
	FOR (ww=0;ww < NumEntries; ww++)
		WAVE XX=AllACs[ww]
		WAVE Cur=AllIWRs[ww]
		NewName=NameOfWave(XX)
		DFREF df= GetWavesDataFolderDFR(Cur)	// get correct datafolder to write into
		SetDataFolder df
		Duplicate/O XX,$NewName	
		WAVE AC=$NewName
		AllACWRs[DimSize(AllACWRs,0)]={AC}
		// try to obtain firing rate from corresponding _ST wave
		IName	= NameOfWave(cur)
		STname=ReplaceString("_I", IName, "_ST")
		npnts=DimSize(Cur,0)
		duration=npnts*Dimdelta(Cur,0)

		IF (Exists(STname)==1)
			WAVE st=$STname
			
			Note/NOCR AC,"spikerate:"+num2str(DimSize(st,0)/Duration) +"\r"
		ENDIF

	ENDFOR

	SetDataFolder rootFolderRf
	KillWaves AllACs
END		// Make AC

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


FUNCTION	GaussFilter(FT,[FreqPoints, widthFac, maxNumFreqPoint, WarnForRealInput])
WAVE/C		FT 		// complex wave, fourier transform of wave
WAVE/D		FreqPoints	// vector with Frequency values at which filtered versions of 
						// FT will be constructed
VARIABLE	widthFac		// width of the gaussian filter	; a value of 1 corresponds to 
						// Higgs and Spain's choice of f/2Pi
VARIABLE	maxNumFreqPoint	// the prefered number of entries in the outpur, unless those number is larger than half the input number
		// the optional parameters FreqPoints, maxNumFreqPoint, and widthFac  are supplied as follows:
		// GaussFilter(FT, widthFac=1, FreqPoins=SomeExisting Wave, maxNumFreqPoint=100)
VARIABLE	WarnForRealInput // when set to 0, no warning is issued for real input waves
// check whether Input is actually a complex wave
	STRING	OutName1 = NameOfWave(FT),OutName2 = NameOfWave(FT)
	OutName1=OutName1[0,24]+"_MgFlt"
	OutName2=OutName2[0,24]+"_PhFlt"

	VARIABLE WavNumType=NumberByKey("NUMTYPE", WaveInfo(FT,0) ) 
	IF (mod(WavNumType,2)==0) 	// NOT complex
		IF (WarnForRealInput)
			DoAlert 0,"A complex wave should be supplied. Check your reasoning, the output will have meaningless phase."
		ENDIF
		Duplicate FT, DummyFT
		Redimension/C DummyFT
		WAVE/C FT = DummyFT
	ENDIF
	
//  Will  produce a Magnitude and Phase wave
	VARIABLE	df_in=DimDelta(FT,0)		// frequency step width in input
	VARIABLE	minf=DimOffset(FT,0)
	VARIABLE	maxf=minf+ (DimSize(FT,0) -1)*df_in
	VARIABLE	nPoints	
	
	// new logic for auomated frequency point calculation
	// minimal step size (from minf to the next point= should be the same as in the input
	// until July 2021 it was 1 Hz, which was coincidentially identical to the 
	// input df for all cases treated so far
	IF(ParamIsDefault(maxNumFreqPoint))
		maxNumFreqPoint=50
	ENDIF
	IF(ParamIsDefault(WarnForRealInput))
		WarnForRealInput=1
	ENDIF
		
	IF (ParamIsDefault(FreqPoints))
           	nPoints = min(floor((DimSize(FT,0) -1)/2), maxNumFreqPoint)		// max number of frequency points should not be more than 
																// half the number that is there already
		// linear distance in log requires a frequency factor 
		// to describe construction of the FreqList
//			VARIABLE	freqFac =10^( log(maxf-minf)/(nPoints-1))
		VARIABLE	freqFac =10^( log((maxf-minf)/df_in)/(nPoints-1))	// equivalent to the (n-1)th root of the quotient (maxf-minf)/df_in
																					// I assume that this is a better way to calculate this
																					// would be interesting to know how the two versions would compile
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
	KillWaves/Z    DummyFT  GWeight , RealPart, ImagPart
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

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
// # # # # # # #  	                         			 # # # # # # # # # # # # 
// # # # # # # #        Vectorn Strength-based Gain      # # # # # # # # # # # # 
// # # # # # # #  	                        			 # # # # # # # # # # # # 
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

FUNCTION CreatePhaseWave(V_Wave)
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
	
	

	VARIABLE	nInAna=DimSize(V_Wave,0)


	STRING 		IName, STName, VSName, SPhaseName, VoltName=NameOfWave(V_Wave)
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
	
	// estimate sine amplitude 
	VARIABLE FFT_derived_sineAmp = Temp_FFT[p_freq]/(floor(L_Input/2)) 	// normalize by 
																		// dividing by number of input points 
																		// and multiply with 2 (positive/neg freq.)
																
	STRING SineNote

	sprintf SineNote,"Sine Amp (%s):%g",WaveUnits(I_Wave,-1),FFT_derived_sineAmp
	
	KillWaves Temp_FFT, Temp_Part

	SPhase[]=atan2(sin(2*Pi*Freq*ST[p]),cos(2*Pi*Freq*ST[p]))
	
	VARIABLE/C VSandPhase=VSandPhasefromPhases(SPhase)
	
	//Delete the note 
	// write the Frequency of the sinusoidal wave in wave note
	// write total vector strength in the note 
	VARIABLE  spikerate=DimSize(ST,0)/nInAna/DimDelta(V_Wave,0)
	Note /K SPhase,"spikerate:"+num2str(spikerate)
	Note SPhase, "Sine frequency (Hz): " +num2str(Freq)
	Note SPhase, SineNote
	Note SPhase, "Vector Strength: " +num2str(real(VSandPhase))
	Note SPhase, "Avg Phase: " +num2str(imag(VSandPhase))
	Note SPhase, "Gain (Hz/"+WaveUnits(I_Wave,-1)+"): " +num2str(2*real(VSandPhase)/FFT_derived_sineAmp*spikerate)

END  // CreatePhaseWave



Threadsafe FUNCTION/C VSandPhasefromPhases(PhaseWave)
WAVE	PhaseWave	// contains all the phases (rad) that a series of events has with respect to a certain reference frequency
					// the vector strength analysis assigns a  vector strength and an average phase a  to this list
					// those two values are returned as the real and imaginary part of one complex number
					

	Duplicate/FREE PhaseWave, realPart, cmplxPart
	Redimension/D realPart, cmplxPart
	
	realPart		=	cos(PhaseWave)
	cmplxPart		= 	sin(PhaseWave)
	
	VARIABLE	rAVG	=	mean(realPart), cAVG=mean(cmplxPart)
	VARIABLE	VS		= 	sqrt(rAVG^2+cAVG^2), phsAVG = atan2(cAVG,rAVG)
	
	KillWaves/Z realPart, cmplxPart

	VARIABLE/C returner=cmplx(VS,phsAVG)
	Return returner
	
END // VSandPhasefromPhases


Threadsafe FUNCTION VSfromPhases(PhaseWave)
WAVE	PhaseWave	// contains all the phases (rad) that a series of events has with respect to a certain reference frequency
					// the vector strength analysis assigns a  vector strength and an average phase a  to this list
					// the vector strength is returned 
					

	Duplicate/FREE PhaseWave, realPart, cmplxPart
	Redimension/D realPart, cmplxPart

	
	realPart		=	cos(PhaseWave)
	cmplxPart		= 	sin(PhaseWave)
	
	VARIABLE	rAVG=mean(realPart), cAVG=mean(cmplxPart)
	VARIABLE	VS	= sqrt(rAVG^2+cAVG^2)
	
	KillWaves/Z realPart, cmplxPart

	Return VS
END // VSfromPhases

