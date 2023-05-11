#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// 05/04/17 aNNe  -  replaced names of Execute* procedures for shorter versions
//							ExecuteInAllSubFolders_aNNe --> XeqtInSubs
//							ExecuteForAllInList         --> Xeqt4List
//							ExecuteForSeries				 --> Xeqt4Series
//
// 03/07/18 aNNe  -  In function Xeqt4List: added Line 57: Sort the list : Case-insensitive alphanumeric sort that sorts wave0 and wave9 before wave10.
//
// 16/07/20 aNNe  -  NormTraces is now aware of the point range plotted
//
// 12/4/21  aNNe	  -  In the Xeqt4List and Xeqt4WList it is now possible to programmatically modify the strings from the list before they are introduced to the command line at
//						instances of '~'. to this end, use §~rmvendX§ , where 'X' is the letter or number or underscore, with which the ENDING starts that is being removed. 
//						importantly, it is the last occurrence of this letter (number etc.) in the string that is considered

Menu "Macros"
		Submenu "Repetitive Tasks"
		"Execute for List...",		EvokeExecForList ()
		"Execute for Series...",			EvokeExecForSeries ()
	end
end
Menu "Graph"

		"Add trace browser",		TraceScanner()
		"Add trace browser transp",		TraceScannerTransparency()
		"Color Traces",			ColorTraceByWave()
		"Normalize Traces",			NormTraces()
		"Set Opacity",			GlobalOpacity()
		"Distribute vertical axes", DistributeAxes()
end
		
// ***************************************************
// aNNe tools for fast processing of repetitive tasks
// ***************************************************
FUNCTION	EvokeExecForList()	
	STRING		targetStr, Liste
	STRING		SkipListe, Cmd
	VARIABLE		fullPath

	Prompt		targetStr, "Waves containing...."
	Prompt		SkipListe, "Semicolon separated list of waves to skip..."
	Prompt		Cmd, "Command to execute, use ~ in place of wave name and § in place of index"
	Prompt		fullPath, "Cput full path before wave name? 0 or 1"
	
	DoPrompt "ExecuteForList Dialog", targetStr, SkipListe, fullpath, Cmd
	IF (V_flag==1)		// cancelled
		Return -1
	ENDIF
	Liste= WaveList(targetStr, ";", "" )
	print "Xeqt4List(WaveList(\""+targetStr+"\", \";\", \"\" ),"+"\""+SkipListe+"\""+" , \""+ReplaceString("\"",Cmd, "\\\"")+"\")"
	Xeqt4List(Liste,SkipListe,FullPath, Cmd)
END

FUNCTION/S ReplaceTokenRemovingEnd(regEx, SourceString, Replacement)
	STRING	regEx,SourceString, Replacement
	// performs a string replacement in SourceString
	// searches for a match with the regExprStr "regEx"
	// then replaces with Replacement, but possibly modifies it before 
	IF (GrepString(SourceString,regEx )!=1)
		Return SourceString
	ENDIF
	
	// prepare the replacement OP
	STRING pre=""
	STRING sep=""
	STRING post=""
	STRING SpecialSep
	STRING TruncatedReplacement, EndStr
	VARIABLE	finish=0
	DO
		pre=""
		sep=""
		post=""
		
		SplitString/E=regEx SourceString, pre,sep,post
		// now replace the center with the Replacement string after removing its end
		// first make sure, the separation character has no special meaning in regEx:
		IF (GrepString(sep,sep)==1)	// all fine
			SpecialSep=Sep
		ELSE	// the sep string ha a special meaning in regEx
			SpecialSep="\\"+sep // add a leading escape character
			IF (GrepString(sep,SpecialSep)!=1)	// could not be repaired
				print " problem with special punctuation character: "+sep
				Return ""
			ENDIF
		ENDIF
		IF (GrepString(Replacement,SpecialSep ))   // the replacement contains at least one 
															 // character that separates the suffix
			SplitString/E="(.*)"+SpecialSep+"(.*)" Replacement, TruncatedReplacement,EndStr  // the first and last pattern here are the characters that are mapped to the 
																												// two parentheses 
																												// the first expression being greedy, this splits the string at the last instance of specialsep
			SourceString= pre+TruncatedReplacement+post
		ELSE
			finish=1
		ENDIF
	
	WHILE ( (GrepString(SourceString,regEx )==1) & (!finish))
	
	Return SourceString
END // ReplaceTokenRemovingEnd

FUNCTION/S ReplaceTokenOnlyEnd(regEx, SourceString, Replacement)
	STRING	regEx,SourceString, Replacement
	// performs a string replacement in SourceString
	// searches for a match with the regExprStr "regEx"
	// then replaces with Replacement, but possibly modifies it before 
	IF (GrepString(SourceString,regEx )!=1)
		Return SourceString
	ENDIF
	
	// prepare the replacement OP
	STRING pre=""
	STRING sep=""
	STRING post=""
	STRING SpecialSep
	STRING TruncatedReplacement, EndStr
	VARIABLE	finish=0
	DO
		pre=""
		sep=""
		post=""
		
		SplitString/E=regEx SourceString, pre,sep,post
		// now replace the center with the Replacement string after removing its end
		// first make sure, the separation character has no special meaning in regEx:
		IF (GrepString(sep,sep)==1)	// all fine
			SpecialSep=Sep
		ELSE	// the sep string ha a special meaning in regEx
			SpecialSep="\\"+sep // add a leading escape character
			IF (GrepString(sep,SpecialSep)!=1)	// could not be repaired
				print " problem with special punctuation character: "+sep
				Return ""
			ENDIF
		ENDIF
		IF (GrepString(Replacement,SpecialSep ))   // the replacement contains at least one 
															 // character that separates the suffix
			SplitString/E="(.*)"+SpecialSep+"(.*)" Replacement, TruncatedReplacement,EndStr  // the first and last pattern here are the characters that are mapped to the 
																												// two parentheses 
																												// the first expression being greedy, this splits the string at the last instance of specialsep
			SourceString= pre+EndStr+post
		ELSE
			finish=1
		ENDIF
	
	WHILE ( (GrepString(SourceString,regEx )==1) & (!finish))
	
	Return SourceString
END // ReplaceTokenOnlyEnd



FUNCTION	Xeqt4List(Liste,SkipList, FullPath, Cmd)
STRING		Liste, SkipList
STRING		Cmd
VARIABLE	fullPath
	
	// generates a string (from "Cmd" )to be executed, where the following replacements are made
	//	\~ 	--> ~
	// ~  	--> wavename from the list
	// §	--> index of the wave in the list 
	// §count§		--> total number of waves in the list 
	// §~rmvend_§	--> WaveName with characters removed up to and including the last "_"
	// §~onlyend_§	--> WaveName with characters up to but excluding the last "_"
	// the ascii "_" can be also any letter or digit
	
	// option: add full path name before wave name
	STRING	FP=GetDataFolder(1)
	VARIABLE	k, NumEntries=ItemsInList(SkipList, ";")
	FOR (k=0; k<NumEntries;k+=1)
		//print StringFromList(k, SkipList,";") 
		Liste=RemoveFromList(StringFromList(k, SkipList,";") , Liste,";")
	ENDFOR
	Liste=SortList(Liste, ";",16)
	NumEntries=ItemsInList(Liste, ";")
	STRING		cCmd, subCmd
	VARIABLE	m
	cmd=ReplaceString("\~", cmd, "+\+")

	FOR (k=0; k<NumEntries;k+=1)
		IF (fullPath)
			cCmd=ReplaceTokenRemovingEnd("(.*)§~rmvend([[:word:]]){1,1}§(.*)", cmd, FP+StringFromList(k, Liste,";"))
			// this allows use of strings in the list, AFTER THEIR ENDING WAS REMOVED
			// what counts as ending is defined by the appearance of a single letter (number or underscore)
			// the ending starts with the last occurrence of this letter in the string from the list
			// the regex says: the exact string "§~rmvend" followed by exactly one character (underscore or any character that is a letter or digit)
			// and a "§"
			// the expression before and after this one word character (last occurence in wavename) are read out as a subpattern. Those CANNOT contain a newline character!
			cCmd=ReplaceTokenOnlyEnd("(.*)§~onlyend([[:word:]]){1,1}§(.*)", cCmd, FP+StringFromList(k, Liste,";"))
			cCmd=ReplaceString("~", cCmd, FP+StringFromList(k, Liste,";") )
			cCmd=ReplaceString("+\+", cCmd, "~")

//			cCmd=ReplaceString("~", Cmd, FP+StringFromList(k, Liste,";") )
			

		ELSE
			cCmd=ReplaceTokenRemovingEnd("(.*)§~rmvend([[:graph:]]){1,1}§(.*)", cmd, StringFromList(k, Liste,";")) // see last case
			cCmd=ReplaceTokenOnlyEnd("(.*)§~onlyend([[:graph:]]){1,1}§(.*)", cCmd, StringFromList(k, Liste,";"))
			cCmd=ReplaceString("~", cCmd, StringFromList(k, Liste,";") )
			cCmd=ReplaceString("+\+", cCmd, "~")
//			cCmd=ReplaceString("~", Cmd, StringFromList(k, Liste,";") )
		ENDIF
		cCmd=ReplaceString("§count§", cCmd, num2istr(NumEntries) )
		cCmd=ReplaceString("§", cCmd, num2istr(k) )
		IF (strlen(cCmd) > 1200)

			m=0
			subCmd=StringFromList(m, cCmd,";")
			
			DO
				Execute subCmd
				m+=1
				subCmd=StringFromList(m, cCmd,";")
	
			WHILE (strlen(subCmd) > 0)
			
		ELSE	
			Execute cCmd
		ENDIF
		
	ENDFOR

END

// -   -   -   -   -   -   -   -   -   -
FUNCTION	Xeqt4WList(TargetStr, Cmd)
STRING		TargetStr
STRING		Cmd

// based on Xeqt4List, except, that the list of wave names is not explicitely handed over, instead
// the list is created based on the string handed over
	STRING Liste
	Liste= WaveList(targetStr, ";", "" )
	Xeqt4List(Liste,"",0, Cmd)

END

FUNCTION	Xeqt4TList(Cmd)
STRING		Cmd

// based on Xeqt4List, except, that the list of wave names is not explicitely handed over, instead
// the list is a list of TRACES of the top graph - only normal graph traces are included

// full paths are used



	STRING TListe, WListe=""
	TListe= TraceNameList("", ";", 1 )
	VARIABLE k, nTraces=ItemsInList(TListe,";")
	FOR (k=0; k< nTraces; k+=1)
		WListe+=GetWavesDataFolder(TraceNameToWaveRef("", StringFromList(k,TListe,";") ),2)+";"
	ENDFOR
	//WListe=SortList(WListe, ";",16)
	STRING		cCmd, subCmd
	VARIABLE	m

	FOR (k=0; k<nTraces;k+=1)
			cCmd=ReplaceString("§T§", cmd, StringFromList(k, TListe,";") ) // §T§ substitutes to the current trace name 
			cCmd=ReplaceString("§W§", ccmd, StringFromList(k, WListe,";") ) // §W§ substitutes to the current wave name
			
			cCmd=ReplaceString("§", cCmd, num2istr(k) ) // § substitutes to the current index
		IF (strlen(Cmd) > 800)

			m=0
			subCmd=StringFromList(m, cCmd,";")
			
			DO
				Execute subCmd
				m+=1
				subCmd=StringFromList(m, cCmd,";")
	
			WHILE (strlen(subCmd) > 0)
			
		ELSE	
			Execute cCmd
		ENDIF
		
	ENDFOR


END

// -   -   -   -   -   -   -   -   -   -
FUNCTION	Xeqt4SList(TargetStr, Cmd)
STRING		TargetStr
STRING		Cmd

// similar to Xeqt4WList, but here we collect STRINGS matching the target string, NOT WAVES
	STRING StrListe = StringList(targetStr, ";")
	
	// generates a string (from "Cmd" )to be executed, where the following replacements are made
	//	\~ 	--> ~
	// ~  	--> contents of string from the list
	// §	--> index of the string in the list 
	// §count§		--> total number of strings in the list 
	// §~rmvend_§	--> contents of string with characters removed up to and including the last "_"
	// the ascii "_" can be also any letter or digit
	

	StrListe=SortList(StrListe, ";",16)
	
	VARIABLE	k, NumEntries=ItemsInList(StrListe, ";")

	//create list of contents of the strings
	STRING listItem, Liste=""
	FOR (k=0; k<NumEntries;k+=1)
		SVAR locStr=$(StringFromList(k, StrListe,";"))
		Liste+=locStr+";"
	ENDFOR
	STRING		cCmd, subCmd
	VARIABLE	m
	cmd=ReplaceString("\~", cmd, "+\+")

	FOR (k=0; k<NumEntries;k+=1)

		cCmd=ReplaceTokenRemovingEnd("(.*)§~rmvend([[:graph:]]){1,1}§(.*)", cmd, StringFromList(k, Liste,";")) // see last case
		cCmd=ReplaceTokenOnlyEnd("(.*)§~onlyend([[:graph:]]){1,1}§(.*)", cCmd, StringFromList(k, Liste,";"))
		cCmd=ReplaceString("~", cCmd, StringFromList(k, Liste,";") )
		cCmd=ReplaceString("+\+", cCmd, "~")
		cCmd=ReplaceString("§count§", cCmd, num2istr(NumEntries) )
		cCmd=ReplaceString("§", cCmd, num2istr(k) )
		IF (strlen(cCmd) > 400)

			m=0
			subCmd=StringFromList(m, cCmd,";")
			
			DO
				Execute subCmd
				m+=1
				subCmd=StringFromList(m, cCmd,";")
	
			WHILE (strlen(subCmd) > 0)
			
		ELSE	
			Execute cCmd
		ENDIF
		
	ENDFOR

END

// -   -   -   -   -   -   -   -   -   -

FUNCTION	EvokeExecForSeries()	
	VARIABLE	FirstX=0, LastX=1, Step=1
	STRING		Cmd

	Prompt		FirstX, "First number...."
	Prompt		LastX, "Last number..."
	Prompt		Step, "Step size..."
	Prompt		Cmd, "Command tobe executed, use ~ in place of number"
	
	DoPrompt "ExecuteForList Dialog", FirstX, LastX, Step, Cmd
	IF (!V_flag)
		print "ExecuteForSeries("+num2istr(FirstX)+" , " +num2istr(LastX)+" , " +num2istr(Step)+" , \"" +ReplaceString("\"",Cmd, "\\\"")+"\")"
		Xeqt4Series(FirstX,LastX, Step, Cmd)
	ENDIF
END


FUNCTION Xeqt4Series(firstX,LastX, Step, Command)
VARIABLE	FirstX, LastX, Step
STRING		Command

STRING		cmd

VARIABLE	r
IF (Step>0)
	FOR (r=firstX; r<=Lastx; r+=Step)
		cmd=ReplaceString("\~", Command, "+\+")
		cmd=ReplaceString("~", cmd, num2str(r))
		cmd=ReplaceString("+\+", cmd, "~")
		Execute /Q cmd
	ENDFOR
ELSE
	FOR (r=firstX; r>=Lastx; r+=Step)
		cmd=ReplaceString("\~", Command, "+\+")
		cmd=ReplaceString("~", cmd, num2str(r))
		cmd=ReplaceString("+\+", cmd, "~")
		Execute /Q cmd
	ENDFOR
ENDIF	
END


Function XeqtInSubs(Cmd)
STRING Cmd
// the string Cmd is executed in every subfolder

// three key strings are defined: 
// §ROOT§ is the folder from which the comand is called
// §SUBFULL§ and §SUB§	are the full path or partial path name
// of the folder in which the string is currently executed
	
	DFREF		topDf=GetDataFolderDFR()
	STRING	substrFull,substrPartial, rootStr=GetDataFolder(1,topDF)
	STRING	cmdStr
	
	STRING		FolderListe=ReplaceString(",", StringByKey("FOLDERS", DataFolderDir(1,topDf)+",",  ":", ";") ,";") 
	FolderListe= RemoveFromList("Packages;IGNORE;", FolderListe, ";")
	FolderListe=SortList(FolderListe,";",16)
	VARIABLE	runde,numFolders=ItemsInList(FolderListe)
	
	FOR (runde=0; runde< numFolders; runde+=1)
		SetDataFolder $(StringFromList(runde,FolderListe,";"))
		subStrFull=rootStr+StringFromList(runde,FolderListe,";")+":"
		subStrPartial=StringFromList(runde,FolderListe,";")//+":"
		cmdStr=ReplaceString("§ROOT§", cmd, rootStr, 1)
		cmdStr=ReplaceString("§SUB§", cmdStr, subStrPartial, 1)
		cmdStr=ReplaceString("§SUBFULL§", cmdStr, subStrFull, 1)
		cmdStr=ReplaceString("§§", cmdStr, num2istr(runde))

		Execute cmdStr
		SetDataFolder topDf	
	ENDFOR

END

FUNCTION BatchRenameDataFolder(FindString,ReplacementStr,Prefix,Suffix)
STRING		FindString,ReplacementStr,Prefix,Suffix

// this renames subfolders in the present DataFolder

// Get List of DataFolders
	STRING		currdirStr, dirList=StringByKey("FOLDERS",  DataFolderDir(1)  ,":",";")
	// gives a list of Folder Names (potentially liberal)
	// COMMA separated
	VARIABLE	currDir, nDirs=ItemsInList(dirList,",")
	FOR (currDir=0; currDir<nDirs; currDir+=1)
		currdirStr=StringFromList(currDir,dirList,",")
		RenameDataFolder $currDirStr, $(Prefix+ ReplaceString(FindString, currdirStr, ReplacementStr)+ Suffix)
	ENDFOR
END

FUNCTION/WAVE CrawlTree()
// returns a wave containing refrences to all the waves in all folders under currently active folder

	DfREF StartDF=getDataFolderDFR()
	
	
	// a wave to hold references to all folders
	MAKE/O/N=0/DF AllRefs
	AddFolders(AllRefs, GetDataFolderDFR())
	// now this contains references to all folders
	
	SetDataFolder StartDF
	
	VARIABLE	ww,nw,ff, nf=DimSize(AllRefs,0)
	// now visit all the folders and get the waves
	// and print results
	MAKE/O/WAVE/N=0 AllWaves
	
	
	FOR (ff=0; ff< nf; ff++)
		DFREF cf=AllRefs[ff]
		DFREF returnTo=GetDataFolderDFR()
		SetDataFolder cf
		//print GetDataFolder(1)
		SetDataFolder returnTo
		 
		nw=CountObjectsDFR(cf, 1)
		FOR (ww=0; ww < nw; ww++)
			WAVE cw=cf:$(GetIndexedObjNameDFR(cf,1,ww))
			AllWaves[DimSize(AllWaves,0)]={cw}
			//print NameOfWave(cw), StringByKey("SIZEINBYTES",WaveInfo(cw,0),":",";")
		ENDFOR
	ENDFOR	
	
	nw=DimSize(AllWaves, 0)
	MAKE/O/N=(nw) WaveSizes
	MAKE/O/T/N=(nW) WaveNames
	MAKE/O/T/N=(nw) WaveHashes
	WaveSizes[]= str2num(StringByKey("SIZEINBYTES",WaveInfo(AllWaves[p],0),":",";"))
	WaveNames[]= GetWavesDataFolder(AllWaves[p],2)
	FOR (ww=0; ww < nw; ww++)
		WAVE currWave=AllWaves[ww]
		IF (Wavetype(currWave,1)<3)
				WaveHashes[ww] = WaveHash(currWave,3)
		ENDIF
	ENDFOR
	// Display w annotations
	Display/K=0 WaveSizes
	Label left "Size (Byte)";DelayUpdate
	ModifyGraph log(left)=1
	ModifyGraph mode=3,marker=19,msize=1.5,rgb=(0,0,0)
	
	// label outliers
	WAVESTATS/Q WaveSizes
	VARIABLE	Thresh=V_avg+V_sdev*1
	STRING	TagName
	FOR (ww=0; ww < nw; ww++)
		IF (WaveSizes[ww] > Thresh)
			TagName="tag_"+num2istr(ww)
			Tag/C/N=$TagName/F=0/B=1/H=12/A=LC/L=1 WaveSizes, ww,WaveNames[ww]

		ENDIF
	ENDFOR

			
END

FUNCTION/WAVE AddFolders(FolderRefsWave, FolderRef)
// this function adds the one FolderRef that is handed over 
// to the wave containing all previous folder refs
// and it calls itself with each of the folders inside itself as reference	  
WAVE/DF	FolderRefsWave // contains Df references to all previously discovered folders
DFRef FolderRef

FolderRefsWave[DimSize(FolderRefsWave,0)]={FolderRef}

// get number of folders inside present folder
VARIABLE ff, nf=CountObjectsDFR(FolderRef, 4 )
STRING	fName
FOR (ff=0; ff< nf; ff++)
	fName=GetIndexedObjNameDFR(FolderRef, 4, ff )
	SetDataFolder FolderRef:$fname
	//SetDataFolder $(":"+fName)
	AddFolders(FolderRefsWave, GetDataFolderDFR())
ENDFOR

END
	
// *********************************************************************
// ** * * * * * * *                                       * * * * * * ** 
// *** * * * * * * *         Data stratification         * * * * * * ***
// ** * * * * * * *                                       * * * * * * ** 
// *********************************************************************

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
			case "OR":
		//	case "or":	// execute if case matches expression
				MATRIXOP/O Indicies0=(Indicies0 || Indicies1)
				break
			case "and":
			//case "AND":	// execute if case matches expression
				Indicies0*=Indicies1							
				break		// exit from switch
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
	
	// sometimes at least one entry is required
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
	cellFolderList = ReplaceString("IGNORE;", cellFolderList, "")
	cellFolderList = ReplaceString("Packages;", cellFolderList, "")
	
	VARIABLE	i, nSubs=ItemsInList(cellFolderList,";")
	// for the case of finding NO subfolders, operate in the folder, from where the function was called, 
	IF (nsubs<1)
		cellfolderlist=GetDataFolder(1)	+";"
		nsubs=1
	ENDIF
	DFREF		rootFolderRf=GetDataFolderDFR()
	STRING		rootDFPath
	rootDFPath = GetDataFolder(1)		// OK
	VARIABLE	k, nWaves, noteVal
	
	STRING		Notiz, Note_String, CriteriumWaveName, CriteriumWaveName1
	
	VARIABLE	FirstTime=1
	
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
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit)}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit, upCrit)}, Target
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
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit1, upCrit1)}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit1, upCrit1)}, Target
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
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit,CriteriumData1=crit1, lowCriterium1=lowCrit1, highCriterium1=upCrit1, logic="AND")}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit, upCrit,CriteriumWave1=crit1, lowCrit1=lowCrit1, upCrit1=upCrit1, Logic="AND")}, Target

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
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit)}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit, upCrit)}, Target
						ENDIF
					ELSEIF (DoStratify == 2)
						// finding corresponding criterium wave
						CriteriumWaveName=RemoveEnding(StringFromList(k, Wave_List,";"),Suffix)+CriteriumSuffix1
						IF (!Exists(CriteriumWaveName)==1)
							DoAlert 0,"Criterium wave for " + GetWavesDataFolder(Candidate,2) + " not found\rTake all points."
							Concatenate/NP=0  {Candidate}, Target
						ELSE
							WAVE Crit=$CriteriumWaveName
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit1, upCrit1)}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit1, upCrit1)}, Target
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
							//Concatenate/NP=0  {waveStratifyData(Candidate,Crit,lowCrit, upCrit,CriteriumData1=crit1, lowCriterium1=lowCrit1, highCriterium1=upCrit1, logic="AND")}, Target
							Concatenate/NP=0  {WaveSubsetByCriteria(Candidate,Crit,lowCrit, upCrit,CriteriumWave1=crit1, lowCrit1=lowCrit1, upCrit1=upCrit1, Logic="AND")}, Target
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

END		// ListDataByNote
	
// *********************************************************************
// ** * * * * * * *                                       * * * * * * ** 
// *** * * * * * * *          aNNe tools for graphs      * * * * * * ***
// ** * * * * * * *                                       * * * * * * ** 
// *********************************************************************

FUNCTION TagTraceWithWaveName()
// Scans through all traces
// identifies the point where the avg distance to all other traces is maximal
// places a tag here carrying the wave name

STRING	TraceNL=TraceNameList("", ";", 1 )
VARIABLE	NTraces=ItemsInLIst(TraceNL,";"), trace
STRING	C_trace// current Trace
STRING	TempNameX,TempNameY
	

FOR	(trace=0;trace<NTraces;trace+=1)
	C_Trace=StringFromList(trace, TraceNL,";")
	TempNameX="W_tmp_x_"+num2istr(trace)
	TempNameY="W_tmp_y_"+num2istr(trace)
	WAVE	YW=GetCopyYW(C_Trace,TempNameX)
	WAVE	XW=GetCopyYW(C_Trace,TempNameY)
	IF (!WaveExists(XW)) // not plotted against another wave
			
		MAKE/N=(DimSize(YW,0)) TmpXW
		//	TempXW[]=DimOffset
	ENDIF
ENDFOR
END
	
FUNCTION/WAVE 	GetCopyYW(TraceName,NewName)
STRING		TraceName,NewName

	// find out  what the Y-coordinates of the trace are (displayed part) and 
	// copy those over in a new free wave with the name NewName
	WAVE Source=TraceNameToWaveRef("",TraceName)
	STRING	SourceName=NameOfWave(Source)
	VARIABLE	Instance=0
	IF (ItemsInList(TraceName,"#")>1)
		Instance=Str2num(StringFromList(1,TraceName,"#"))
	ENDIF
	MAKE/FREE/N=(2,4)	bnd		// keeps first and last index of each dimension
	STRING	YRangeStr=StringByKey("YRANGE", TraceInfo("",SourceName,Instance))
	VARIABLE NDim=ItemsInList(YRangeStr,"]")	// MultiDim
	VARIABLE Dim=0
	
	STRING	SubString

	FOR (Dim=0;Dim<nDim; Dim+=1)
		SubString=StringFromList(Dim, YRangeStr , "]")
		SubString=SubString[1,strlen(Substring)-1]
		IF (CmpStr("*", SubString ,0)==0)	
			bnd[0][Dim]=0
			bnd[1][Dim]=DimSize(Source,Dim)-1
		ELSEIF (ItemsInList(SubString,",")==1 )
			bnd[0][Dim]=str2num(SubString)
			bnd[1][Dim]=str2num(SubString)		
		ELSEIF (ItemsInList(SubString,",")==2 )
			bnd[0][Dim]=str2num(StringFromList(0,SubString,","))
			bnd[1][Dim]=str2num(StringFromList(0,SubString,","))
		ELSE // Unexpected
			DoAlert 0,"Something not working in extracting the subranges plotted here"
			Return Source
		ENDIF
	ENDFOR
	
	switch(NDim)	// numeric switch
		case 1:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]] Source, FreeWDummy
			break		// exit from switch
		case 2:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]] Source, FreeWDummy
			break		// exit from switch
		case 3:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]][bnd[0][2],bnd[1][2]] Source, FreeWDummy
			break		// exit from switch
		case 4:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]][bnd[0][2],bnd[1][2]][bnd[0][3],bnd[1][3]] Source, FreeWDummy
			break		// exit from switch
		default:			// optional default expression executed
			DoAlert 0,"Something not working in extracting the subranges plotted here"
			Return Source
	endswitch
	Redimension/N=-1 FreeWDummy // make it 1D

//	Rename FreeWDummy, $NewName
//	WAVE	ReturnW=$NewName
	Return FreeWDummy
END

FUNCTION/WAVE 	GetCopyXW(TraceName,NewName)
STRING		TraceName,NewName

	WAVE Source=TraceNameToWaveRef("",TraceName)
	STRING	SourceName=NameOfWave(Source)
	VARIABLE	Instance=0
	IF (ItemsInList(TraceName,"#")>1)
		Instance=Str2num(StringFromList(1,TraceName,"#"))
	ENDIF

	// find out whether there is ANY XWAVE
	STRING XWName=StringByKey("XWAVE",  TraceInfo("",SourceName,Instance),":",";")
	IF (Strlen(XWName)==0)	// no Xwave
			// no time to solve this
			// the flow should be: find the ONE dimension that is plotted
			// look at YRANGE!!
			// either there is only a single bracket with [*] - then it is the zero-th Dimension
			// otherwise find out in which dimension there is either a "*" or a range (there is a comma) or, 
			// NOT SURE THE FOLLOWING HAS TO BE ACCOUNTED FOR
			// if there is no such dimension, then a single point is plotted, in this case, the entry of hte 
			// wave at this one single point will be plotted against WHAT? I would not be sure how that is possible
			
	ELSE
			WAVE Source=$(StringByKey("XWAVEDF", TraceInfo("",SourceName,Instance))+":"+XWName)
	ENDIF

	// find out how what the Y-coordinates of the trace are (displayed part) and 
	// copy those over in a new free wave with the name NewName


	MAKE/FREE/N=(2,4)	bnd		// keeps first and last index of each dimension
	STRING	XRangeStr=StringByKey("XRANGE", TraceInfo("",SourceName,Instance))
	VARIABLE NDim=ItemsInList(XRangeStr,"]")	// MultiDim
	VARIABLE Dim=0
	
	STRING	SubString

	FOR (Dim=0;Dim<nDim; Dim+=1)
		SubString=StringFromList(Dim, XRangeStr , "]")
		SubString=SubString[1,strlen(Substring)-1]
		IF (CmpStr("*", SubString ,0)==0)	
			bnd[0][Dim]=0
			bnd[1][Dim]=DimSize(Source,Dim)-1
		ELSEIF (ItemsInList(SubString,",")==1 )
			bnd[0][Dim]=str2num(SubString)
			bnd[1][Dim]=str2num(SubString)		
		ELSEIF (ItemsInList(SubString,",")==2 )
			bnd[0][Dim]=str2num(StringFromList(0,SubString,","))
			bnd[1][Dim]=str2num(StringFromList(0,SubString,","))
		ELSE // Unexpected
			DoAlert 0,"Something not working in extracting the subranges plotted here"
			Return Source
		ENDIF
	ENDFOR
	
	switch(NDim)	// numeric switch
		case 1:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]] Source, FreeWDummy
			break		// exit from switch
		case 2:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]] Source, FreeWDummy
			break		// exit from switch
		case 3:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]][bnd[0][2],bnd[1][2]] Source, FreeWDummy
			break		// exit from switch
		case 4:	// execute if case matches expression
			DUPLICATE/FREE/R=[bnd[0][0],bnd[1][0]][bnd[0][1],bnd[1][1]][bnd[0][2],bnd[1][2]][bnd[0][3],bnd[1][3]] Source, FreeWDummy
			break		// exit from switch
		default:			// optional default expression executed
			DoAlert 0,"Something not working in extracting the subranges plotted here"
			Return Source
	endswitch
	Redimension/N=-1 FreeWDummy // make it 1D

	Rename FreeWDummy, $NewName
	WAVE	ReturnW=$NewName
	Return ReturnW
END


Function SoftColors(RndInit) //: GraphStyle
	VARIABLE RndInit
	PauseUpdate; Silent 1		// modifying window...
	VARIABLE		nTraces=ItemsInList(TraceNameList("", ";", 1 ),";") 
	MAKE/O/N=(max(nTraces,24))list,rnd;		// by always taking at least the number of different colours as random list (the max() here), 
														// the entire list of colors is available for random assignment even if there are fewer traces
	MAKE/O/N=(25,3) W_softC
	list[]=p
	IF (abs(RndInit)>1e-12)	
		SetRandomSeed(RndInit)
		rnd[0,]=enoise(10,2)
		Sort rnd,list
	ENDIF
	VARIABLE	num=0
	W_softC={{39321,45746,52428,58853,65535,39321,45746,52428,58853,65535,27499,34181,41891,50115,65535,3855,11308,20817,32382,45746,9766,16962,25957,36751,63993},{3855,11308,20817,32382,45746,21588,28527,36494,45489,55512,39321,45746,52428,58853,56283,27499,34181,41891,50115,58853,3855,11308,20817,32382,63993},{3855,11308,20817,32382,45746,3855,11308,20817,32382,45746,3855,11308,20817,32382,9766,39321,45746,52428,58853,65535,39321,45746,52428,58853,0}}
	VARIABLE	nColors=DimSize(W_softC,0)
	DO
		ModifyGraph/Z rgb[num]=(W_softC[mod(list[num],nColors)][0],W_softC[mod(list[num],nColors)][1],W_softC[mod(list[num],nColors)][2])
		num+=1
	WHILE	(num < nTraces)
	ModifyGraph/Z btLen=3
	KillWaves/Z list, rnd,W_SoftC
EndMacro

FUNCTION GlobalOpacity([opacity])
	VARIABLE	opacity // (0-1)
	// function is setting the opacity of all traces on the top graph to the entered value
	IF (ParamIsDefault(opacity))
		opacity = 1
		Prompt	opacity,"opacity value"
		DoPrompt "your choice", opacity
	
		IF (V_flag) // cancel clicked
			Abort
		ENDIF
	ENDIF
	STRING		TraceNames=TraceNameList("", ";", 1 )
	VARIABLE		num=0, nTraces=ItemsInList(TraceNames,";")
	VARIABLE		InfoPointer

	// Loop over all traces
	STRING Info
	VARIABLE	rr,tot
	DO
		// get tracemodifier str 
		Info=TraceInfo("",StringFromList(num,TraceNames,";"),0)
		InfoPointer=strsearch(Info, "RECREATION:",0)+11
		Info=Info[InfoPointer, strlen(Info)-1]
		// now inside the recreation string 
		// find every instance of "RGB(x)=(....);" and replace the color by the same color with opacity "opacity*2^16-1"
		Info=ReplaceRGBOpacities(Info,opacity)	
		Info=ReplaceString("(x)", Info, "["+num2istr(num)+"]")
		// execute every modifying command in the recreation str
		tot=ItemsInList(Info,";")
		rr=0
		DO
			Execute "ModifyGraph "+StringfromList(rr,Info,";")
			rr+=1
		WHILE (rr<tot)
		
		num+=1
	WHILE	(num < nTraces)
	
END

FUNCTION/S ReplaceRGBOpacities(InStr, opacity)
	STRING	InStr
	VARIABLE	opacity
	// takes a string (e.g. a graph macro) and replaces all occurances of 
	// RGB values with RGBO values
	opacity = round( opacity * (2^16-1) )
	
	VARIABLE	nChars=strlen(InStr)
	VARIABLE	inst,nInstances=-1, startpos=0
	VARIABLE	p1,p2
	DO 
		nInstances+=1
		startpos=strsearch(InStr,"RGB(x)=(",startpos,2)+8
	WHILE (startpos>7)
	
	STRING RGB
	VARIABLE	nEntries
	startpos=0

	FOR (inst=0; inst < nInstances; inst += 1)
		p1=strsearch(InStr,"RGB(x)=(",startpos,2)+7
		p2=strsearch(InStr,")",p1)
		RGB=InStr[p1,p2]
		nEntries=ItemsInList(RGB+",",",")
		IF ( (nEntries<3) || (nEntries > 4) )
			// strange, but simply ignore
		ELSEIF (nEntries == 3)
			// no opacity given yet
			RGB=RemoveEnding(RGB,")")
			RGB=RGB+","+num2istr(opacity)+")"
			InStr=InStr[0,p1-1]+RGB+InStr[p2+1,strlen(InStr)-1]			
		ELSE // opacity already given
			RGB=RGB[0,strsearch(RGB,",",strlen(RGB)-1,1)-1]
			RGB=RGB+","+num2istr(opacity)+")"
			InStr=InStr[0,p1-1]+RGB+InStr[p2+1,strlen(InStr)-1]			
		ENDIF

		startpos=p1+strlen(RGB)
	ENDFOR
	RETURN InStr
END	
		
FUNCTION SoftColors2(RndInit,[fraction, cyclelength])
// colors traces sequentially going 3 times throught the rainbow with different saturation
// the "fraction" can be used, especially without randomization, i.e. with rndinit=0
// to use only a smaller part of the entire color scale and in this way get only soft colors

// the random init allows to quickly generate different colorings. those are reproducible, as long as the number of traces
// stays the same

// cycle length = CL is used to color every CL-th trace with the same color

	
	VARIABLE RndInit, Fraction, cycleLength

	
	IF (ParamIsDefault(fraction))
		fraction=1
	ENDIF
	
	PauseUpdate; Silent 1		// modifying window...
	VARIABLE		nTraces=ItemsInList(TraceNameList("", ";", 1 ),";")
	IF (ParamIsDefault(cycleLength))
		cycleLength=nTraces
	ENDIF
	VARIABLE		nColors=round(max(18/fraction,min(nTraces,cycleLength)/fraction)) 
	MAKE/O/N=(nColors)list,rnd	
	MAKE/O/N=(18,3) W_softC
	VARIABLE nFixPnt=DimSize(W_softC,0)
	MAKE/O/N=(nFixPnt) FixPnts
	WAVE	FixPnts=FixPnts
	
	FixPnts[]=round(p*(nColors-1)/(nFixPnt-1))	// this is just linear, later these positions can be handed over
	W_softC={{65535,65535,63932,56392,43348,43561,65535,65535,62591,58187,36872,0,6104,43519,32000,32000,18812,0},{43478,59078,65535,65535,65535,34627,19622,5866,44376,62719,62719,62719,0,0,187,22687,32000,3702},{43348,43348,43348,43348,59035,65535,48407,5515,0,0,0,44344,53887,27285,0,0,0,32128}}

	MAKE/O/N=(nColors,3) myColors
		VARIABLE 	white=65535						// black is zero
	VARIABLE 	colorStep
	VARIABLE	clmn, fp							// counters in FOR loops: column and fixpoint
	// the following is done for 3 color columns one after the other:
	IF (nColors>nFixPnt)
		FOR (clmn=0; clmn <3; clmn +=1)
			// linear changes between the fixpoints
			FOR (fp=0; fp < nFixPnt -1; fp+=1  )
				myColors[FixPnts[fp],FixPnts[fp+1]][clmn]= W_softC[fp][clmn]+(p-FixPnts[fp])/(FixPnts[fp+1]-FixPnts[fp])*( W_softC[fp+1][clmn]- W_softC[fp][clmn] )
			ENDFOR
		ENDFOR
	ELSE
		Duplicate/O W_softC, myColors
	ENDIF
	list[]=p
	IF (abs(RndInit)>1e-12)	
		SetRandomSeed(RndInit)
		rnd[0,]=enoise(10,2)
		Sort rnd,list
	ENDIF
	VARIABLE	num=0
	// the actual coloring part 
	DO
		ModifyGraph/Z rgb[num]=(myColors[mod(list[mod(num,cycleLength)],nColors)][0],myColors[mod(list[mod(num,cycleLength)],nColors)][1],myColors[mod(list[mod(num,cycleLength)],nColors)][2])
		num+=1
	WHILE	(num < nTraces)
	ModifyGraph/Z btLen=3
	KillWaves/Z list, rnd,W_SoftC, myColors,FixPnts
END // soft color 2

FUNCTION ColorByGeo([cyclelength, invert])
VARIABLE cyclelength, invert
	
	
	ColorTab2Wave Geo // creates M_colors with 256 entries 
						  // I only intend to use entries 4 to 189 (a span of 186)
	WAVE 			M_colors
	VARIABLE 		num=0
	STRING			AllTraces=TraceNameList("", ";", 1 ), CurrTrace
	VARIABLE		NTraces= ItemsInList(AllTraces,";"), trNum
	
	IF (ParamIsDefault(cycleLength))
		cycleLength=nTraces
	ENDIF
	IF (ParamIsDefault(invert))
		invert=0
	ENDIF

//	Variable		stepSize = 186/max(5,min(NTraces,cyclelength))   // the 5 assures that the colors are not spread unnecessarily wide 
//														 						// for small numbers of traces

	Variable		stepSize = 186/(cyclelength-1)   // the 5 assures that the colors are not spread unnecessarily wide 
														 						// for small numbers of traces
	IF (!invert)
											 	
		DO
			ModifyGraph/Z rgb[num]=(M_colors[4+stepsize*mod(num,cyclelength)][0],M_colors[4+stepsize*mod(num,cyclelength)][1],M_colors[4+stepsize*mod(num,cyclelength)][2])
			num+=1
		WHILE	(num < nTraces)
	ELSE
		DO
			ModifyGraph/Z rgb[nTraces-1-num]=(M_colors[4+stepsize*mod(num,cyclelength)][0],M_colors[4+stepsize*mod(num,cyclelength)][1],M_colors[4+stepsize*mod(num,cyclelength)][2])
			num+=1
		WHILE	(num < nTraces)
	ENDIF		
	ModifyGraph/Z btLen=3
END

FUNCTION ColorTraceByWave([CB])
WAVE	CB // color basis


// the application case for this function is a graph with N traces 
// where each trace should be coloured according to a single value
// these values are stored in the CBasis wave
// 

	VARIABLE	CBasisIndx
	IF (ParamIsDefault(CB))	
		Prompt 	CBasisIndx, "which wave holds the basis for coloring", popup, WaveList("*",";","MAXCOLS:1" )
		DoPrompt "Color scale source", CBasisIndx
		IF (V_flag)		// user cancelled
			Return -1
		ENDIF
	
		
		STRING 	CBasisName = StringFromList(CBasisIndx-1,WaveList("*",";","MAXCOLS:1" ),";")
		WAVE		CB=$CBasisName
	ENDIF
	DFREF			FolderRf=GetDataFolderDFR()
	
	// - - - - - - - - - - - -   p r e p a r a t i o n  - - - - - - - - - - - - - - - - - - - 
	// Prepare PackageFolder to hold colorscale wave


	NewDataFolder/O/S Packages
	NewDataFolder/O/S AN_GraphCS
	
	
	STRING		AllTraces=TraceNameList("", ";", 1 ), CurrTrace
	VARIABLE		NTraces= ItemsInList(AllTraces,";"), trNum
	WAVESTATS/Q  CB
	IF (V_numNaNs>0)
		SetDataFolder FolderRf
		DoAlert 0,"Your color basis wave " + NameOfWave(CB)+" contains "+num2istr(V_numNaNs)+" NaN values. Aborting..."
		Return -1
	ENDIF
	VARIABLE		maxBasisValue=V_max, minBasisValue=V_min	, opacity=1
	
	STRING		GraphName = WinName(0,1) // top (index zero), Graph (bit 1)
	STRING		CScaleName= GraphName[0,min(27,strlen(GraphName)-1)]+"_ClSc"
	
	STRING		OriCScaleName, CScaleList= CTabList()
	VARIABLE		OriCScale
	Prompt		OriCScale, "Color scale", popup,  CScaleList
	Prompt		minBasisValue, "lower value of colorscale"
	Prompt		maxBasisValue, "upper value of colorscale"
	Prompt		opacity, "how opaque are traces"
	
	VARIABLE		SupportOpacity=NumberByKey("IGORVERS", IgorInfo(0))>6
	
	IF (SupportOpacity)	// at least igor 7 (supports opacity)
		DoPrompt		"Color Scale Settings", OriCScale,minBasisValue,maxBasisValue,opacity 
	ELSE
			DoPrompt		"Color Scale Settings", OriCScale,minBasisValue,maxBasisValue
	ENDIF
	
	IF (V_flag)		// user cancelled
		Return -1
	ENDIF


	OriCScaleName=StringFromList(OriCScale-1,CScaleList,";")
	ColorTab2Wave $OriCScaleName
	// this created the wave 
	WAVE M_Colors
	DUPLICATE/O M_Colors, $CScaleName
	WAVE CS=$CScaleName
	Redimension/N=(NTraces,-1) CS // one row per trace

	VARIABLE		nCls=DimSize(M_Colors,0)
	
	CS[][]=M_colors[min(max(0,(CB[p]-minBasisValue)/(maxBasisValue-minBasisValue)*(nCls-1)),nCls-1)][q]
	// the min(max )) construction catches the cases where the maxBasisValue < maxValue and minBasisValue > minValue
	FOR (trNum=0; trNum< NTraces; trNum+=1)
		CurrTrace = StringFromList(trNum,AllTraces,";")
		IF (SupportOpacity)
			ModifyGraph rgb($CurrTrace)=(CS[trNum][0],CS[trNum][1],CS[trNum][2], 65535*opacity)
		ELSE
			ModifyGraph rgb($CurrTrace)=(CS[trNum][0],CS[trNum][1],CS[trNum][2])
		ENDIF
	ENDFOR
	// Add colorscale to graph 
	IF (strlen(WaveUnits(CB, 0))>0)		// the source wave for the coloring has a unit 
			ColorScale/B=1/C/N=CS NameOfWave(CB)+"\\U"
	ELSE
			ColorScale/B=1/C/N=CS NameOfWave(CB)

	ENDIF	
	ModifyGraph margin(right)=85
	ColorScale/C/N=CS ctab={minBasisValue,maxBasisValue,$OriCScaleName,0};DelayUpdate
	ColorScale/C/N=CS tickLen=3.00;DelayUpdate
	//ColorScale/C/N=CS CBasisName
	ColorScale/C/N=CS/F=0

	SetDataFolder FolderRf

END //ColorTraceByWave
	
//////////////////////////////////////////////

FUNCTION TraceScannerTransparency()
// adds a forward, backward button and a display of the current trace number to the top graph
// this functionality is provided by the Function "SwitchTraceProc"

// different modes are possible: show only one trace at a time or show more than one
// in the latter case, a numerical offset between the indicies of the two traces has to be provided



		

		STRING	TopGraphName= StringFromList(0,WinList("*", ";","WIN:1"),";")
		
		// prepare the Variable sthat dictate how many traces are shown and which are visible together
		
		VARIABLE 	nTraces=ItemsInList(TraceNameList(TopGraphName, ";", 1 ), ";")
		VARIABLE		nShown=1, interval, maxIndex
		// nShown is number of simultaneously visible traces
		// interval is index distance between two successive traces that are simultaneously visible
		// maxIndex is maximal index of the first  visible trace
	
	
		Prompt		nShown, "number of traces visible"
		DoPrompt "Display setting I", nShown
		IF (V_Flag)
			Return -1								// User canceled
		ENDIF	
		IF (nShown <1)
			Return -1
		ELSEIF (nShown ==1)
			interval=0
		ELSE
			interval= round((nTraces)/nShown)
			Prompt		interval, "period of trace numbers"
			DoPrompt "Display setting II", interval
			IF (V_Flag)
				Return -1								// User canceled
			ENDIF	

		ENDIF
				
		IF (nShown > 1)
				IF ((interval  != (nTraces)/nShown) && (interval >1))
					DoAlert 0,"With "+num2istr(nTraces)+" traces and an interval of "+num2istr(interval)+"\rit might not cycle through all traces"
				ENDIF
				maxIndex=nTraces-1 - interval*(nShown-1)
		ELSE 
				maxIndex=nTraces-1
		ENDIF  
		
		
		// get size of this Graph
		GetWindow $TopGraphName wsize // Reads window dimensions into 
												//	V_left, V_right, V_top, and V_bottom
												// in points from the top left of the screen.
		VARIABLE		winWidth= V_right-V_left
		VARIABLE		winHeight= V_bottom-V_top
		VARIABLE		transparency
		
		IF (strlen(TopGraphName) < 1)
			DoAlert 0,"No Graph found"
			Return -1
		ENDIF
		

		STRING		Disp=TopGraphName+"_D_"+num2istr(nShown)+"_"+num2istr(interval)

		DoWindow/F $TopGraphName
				
		
		// add the controls
		VARIABLE	FontSize=max(6,min(32,(6+winHeight/45)))
		VARIABLE	ButtonSize=FontSize*2
		
		
		ModifyGraph margin(top)=18;	
		SetVariable $Disp,pos={winWidth-ButtonSize*4,5},size={ButtonSize,ButtonSize},bodyWidth=ButtonSize+1.2*fontsize*floor(log(ntraces-1)),proc=SwitchTraceProcTransp,title="trace #"
		SetVariable $Disp,limits={0,maxIndex,1},value= _NUM:0,fSize=FontSize
		// store currently selected trace in userdata
		SetVariable $Disp,userdata=num2istr(0)
		SetVariable $Disp,userdata(low)="0,0,0,"+num2istr(65535*0.15)+","
		SetVariable $Disp,userdata(high)="0,0,0,65535,"
		Slider Trans, vert=0, pos={winWidth-ButtonSize*2.5,10},size={ButtonSize*2,ButtonSize*0.8}, value=0
		Slider Trans, userdata=Disp // store the name of the corresponding set variable control in oder to access it, when the transparency setting changes
		Slider Trans, proc=TransparencySliderProc, limits={0,1,0.05},ticks=0, side=0
		
END
FUNCTION TraceScanner([VARIABLE nShown])
// adds a forward, backward button and a display of the current trace number to the top graph
// this functionality is provided by the Function "SwitchTraceProc"

// different modes are possible: show only one trace at a time or show more than one
// in the latter case, a numerical offset between the indicies of the two traces has to be provided



		

		STRING	TopGraphName= StringFromList(0,WinList("*", ";","WIN:1"),";")
		
		// prepare the Variable sthat dictate how many traces are shown and which are visible together
		
		VARIABLE 	nTraces=ItemsInList(TraceNameList(TopGraphName, ";", 1 ), ";")
		VARIABLE	interval, maxIndex
		// nShown is number of simultaneously visible traces
		// interval is index distance between two successive traces that are simultaneously visible
		// maxIndex is maximal index of the first  visible trace
	
	IF (ParamIsDefault(nShown))
		nShown=1		
		Prompt		nShown, "number of traces visible"
		DoPrompt "Display setting I", nShown
		IF (V_Flag)
			Return -1								// User canceled
		ENDIF	
		IF (nShown <1)
			Return -1
		ELSEIF (nShown ==1)
			interval=0
		ELSE
			interval= round((nTraces)/nShown)
			Prompt		interval, "period of trace numbers"
			DoPrompt "Display setting II", interval
			IF (V_Flag)
				Return -1								// User canceled
			ENDIF	

		ENDIF
	ELSE
		interval= round((nTraces)/nShown)
	ENDIF
				
		IF (nShown > 1)
				IF ((interval  != (nTraces)/nShown) && (interval >1))
						DoAlert 0,"With "+num2istr(nTraces)+" traces and an interval of "+num2istr(interval)+"\rit might not cycle through all traces"
				ENDIF
				maxIndex=nTraces-1 - interval*(nShown-1)
		ELSE 
				maxIndex=nTraces-1
		ENDIF  
		
		
		// get size of this Graph
		GetWindow $TopGraphName wsize // Reads window dimensions into 
												//	V_left, V_right, V_top, and V_bottom
												// in points from the top left of the screen.
		VARIABLE		winWidth= V_right-V_left
		VARIABLE		winHeight= V_bottom-V_top
		
		IF (strlen(TopGraphName) < 1)
			DoAlert 0,"No Graph found"
			Return -1
		ENDIF
		

	STRING		Disp=TopGraphName+"_D_"+num2istr(nShown)+"_"+num2istr(interval)

	DoWindow/F $TopGraphName
			
	
	// add the controls
	VARIABLE	FontSize=max(6,min(32,(6+winHeight/45)))
	VARIABLE	ButtonSize=FontSize*2
	
	
	ModifyGraph margin(top)=18;	
	SetVariable $Disp,pos={winWidth-ButtonSize*4,5},size={ButtonSize,ButtonSize},bodyWidth=ButtonSize+0.5*fontsize*floor(log(ntraces-1)),proc=SwitchTraceProc,title="trace #"
	SetVariable $Disp,limits={0,maxIndex,1},value= _NUM:0,fSize=FontSize
	Button ShowAll, title="All", pos={winWidth-ButtonSize*2.8,5},size={ButtonSize*0.8,ButtonSize*0.8}
	Button ShowAll, proc=TraceScannerButtonProc
	
END

FUNCTION NormTraces()

// the application case for this function is a graph with N traces 

// the function goes through all traces of the top graph and 
// applies a y-normalization according to the criteria entered via
// the user prompt
// the function is sensitive to the plotted point range
// this function does NOT cange the underlying data


	STRING		NormBasis="",NormBasisList="max;min;avg;var;SD;"
	STRING		Range="Global;Cursors;"
	STRING		NormType="additive;multiplicative;"
	VARIABLE	NormChoice=0
	VARIABLE	UseCursors=0, multiplicative=0
	
	Prompt 	NormChoice, "Criterium for normalization", popup, NormBasisList
	Prompt 	UseCursors, "Evaluation range", popup, Range
	Prompt 	multiplicative, "Normalization type", popup, NormType
	DoPrompt "Normalize traces", NormChoice,UseCursors,multiplicative
	IF (V_flag)		// user cancelled
		Return -1
	ENDIF
	NormBasis = StringFromList(NormChoice-1,NormBasisList,";")
	UseCursors = UseCursors-1
	multiplicative = multiplicative-1
	
	DFREF			FolderRf=GetDataFolderDFR()
	
	STRING		AllTraces=TraceNameList("", ";", 1 ), CurrTrace
	VARIABLE		NTraces= ItemsInList(AllTraces,";"), trNum
	VARIABLE 	p1,p2, pOffset
	STRING		str_Pstart, strPend
	FOR (trNum=0; trNum< NTraces; trNum+=1)
		CurrTrace = StringFromList(trNum,AllTraces,";")
		WAVE	CurrWave=TraceNameToWaveRef("", CurrTrace )
		// if traces are 1D parts of multidimensional waves
		IF (WaveDims(CurrWave)>1) // MultiDim wave, only part of it could be a trace
											// hence that part needs to be extracted to find the ranges of normalization
											
			WAVE Currwave=GetCopyYW(CurrTrace,"W_FreeW")
		ENDIF
		// deal with the case that the trace contains only a subset of the wave,
		// i.e. the cursor point does not reveal the offset (points 10 to 20 are plotted
		// cursor A is on the 4th displayed point pcsr(A)=3, so it actually corresponds
		// to the 13th point of the wave
		SplitString/E="^\[([[:digit:]]+),([[:digit:]]+)\]$" StringByKey("YRANGE", TraceInfo("", CurrTrace, 0 ) ,":",";"), str_Pstart, strPend
		IF (strlen(str_Pstart)) // if the string is not empty
			pOffset = str2num(str_Pstart)
		ELSE
			pOffset = 0
		ENDIF
		
		IF (UseCursors)
			p1=pcsr(A)+pOffset
			p2=pcsr(B)+pOffset
		ELSE
			p1=pOffset
			IF (strlen(strPend)) // if the string is not empty
				p2 = str2num(strPend)
			ELSE
				p2=DimSize(CurrWave,0)-1

			ENDIF

		ENDIF
		WAVESTATS/R=[p1,p2]/Q CurrWave

		IF (multiplicative)
			strswitch(NormBasis)	// string switch
				case "max":	// execute if case matches expression
				ModifyGraph muloffset[trNum]={0,abs(1/V_max)}
					break		// exit from switch
				case "min":	// execute if case matches expression
				ModifyGraph muloffset[trNum]={0,abs(1/V_min)}
					break		// exit from switch
				case "avg":	// execute if case matches expression
				ModifyGraph muloffset[trNum]={0,abs(1/V_avg)}
					break		// exit from switch
				case "var":	// execute if case matches expression
				ModifyGraph muloffset[trNum]={0,1/V_sdev^2}
					break		// exit from switch
				case "SD":	// execute if case matches expression
				ModifyGraph muloffset[trNum]={0,1/V_sdev}
					break		// exit from switch
			endswitch
		ELSE	//	(additive)
				strswitch(NormBasis)	// string switch
				case "max":	// execute if case matches expression
				ModifyGraph offset[trNum]={0,-V_max}
					break		// exit from switch
				case "min":	// execute if case matches expression
				ModifyGraph offset[trNum]={0,-V_min}
					break		// exit from switch
				case "avg":	// execute if case matches expression
				ModifyGraph offset[trNum]={0,-V_avg}
					break		// exit from switch
				case "var":	// execute if case matches expression
				ModifyGraph offset[trNum]={0,-V_sdev^2}
					break		// exit from switch
				case "SD":	// execute if case matches expression
				ModifyGraph offset[trNum]={0,-V_sdev}
					break		// exit from switch
			endswitch
	
		ENDIF
		
	ENDFOR


	SetDataFolder FolderRf

END //ColorTraceByWave

FUNCTION ExpFitAllTraces([includeString,excludestring, x1,x2])
// fit all traces between cursors
VARIABLE		x1,x2
STRING		includeString,excludestring
// trace names containing these strings are excluded (included)

		VARIABLE noA=0, noB=0, useX=0
		IF (strlen(Csrwave(a)) <= 0)
			noA=1
		ENDIF
		
		IF (strlen(Csrwave(b)) <= 0)
			noB=1
		ENDIF
		
		VARIABLE lx,rx
		
		IF (!ParamisDefault(x1))
			// override any potential cursorA
			useX=1
			lx=x1
			noA=0
		ENDIF
		IF (!ParamisDefault(x2))
			// override any potential cursorB
			useX=1
			rx=x2
			noB=0
		ENDIF
			
		STRING 	curTraceN, TraceNL=TraceNameList("",";",1)
		VARIABLE	totTraces=ItemsInList(TraceNL,";")
		VARIABLE	curTrace
		
		// take care of the include and exclude rules.
		
		
	
		
		IF (!ParamisDefault(includestring))
				FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
					curTraceN=StringFromList(curTrace, TraceNL,";")
					IF (!stringmatch(curTraceN,"*"+includeString+"*"))
						TraceNL=RemoveFromList(curTraceN, TraceNL,";")
						totTraces-=1
						curTrace-=1
					ENDIF
				ENDFOR
				totTraces=ItemsInList(TraceNL,";")
		ENDIF
		
		IF (!ParamisDefault(excludestring))
				FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
					curTraceN=StringFromList(curTrace, TraceNL,";")
					IF (stringmatch(curTraceN,"*"+excludestring+"*"))
						TraceNL=RemoveFromList(curTraceN, TraceNL,";")
						totTraces-=1
						curTrace-=1
					ENDIF
				ENDFOR
				totTraces=ItemsInList(TraceNL,";")
		
		ENDIF
		
		// store results in a wave named after the graph
		STRING 	windowName = WinName(0, 3)
		
		Make/O/N=(totTraces,4) $(windowName+"_exp")
		WAVE	results=$(windowName+"_exp")
		SetDimLabel 1,0,y0,results
		SetDimLabel 1,1,A,results
		SetDimLabel 1,2,tau,results
		SetDimLabel 1,3,deltaV,results
		

		
		// convert the pruned tracelist into wavereferences
		MAKE/O/WAVE/N=(totTraces) Tr_refs
		
		Tr_refs[]=TraceNameToWaveRef("", StringFromList(p, TraceNL,";") )
		
		// go through and fit traces
		// this will currently NOT work, if there is an xWAVE
		STRING auxname
		
		FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
			WAVE cW=Tr_refs[curTrace]
			// only for special case that info has to be obtained form another (parallel)
			// wave with related name (e.g. voltage belonging to the current
			// otherwise disable
			auxName=ReplaceString("_I", GetWavesDataFolder(cW,2), "_V")
		
			WAVE auxW=$auxName
			
			IF (noA)
				lx=DimOffset(cW,0)
			ELSEIF (useX)
				lx=min(x1,x2)
			ELSE
				lx=min(xcsr(B),xcsr(A))

			ENDIF
			IF (noB)
				rx=rightx(cW)-Dimdelta(cW,0)
			ELSEIF (useX)
				rx=max(x1,x2)
			ELSE
				rx=max(xcsr(B),xcsr(A))
			ENDIF
			
			CurveFit/Q/M=2/W=0 exp_XOffset, cW(lx,rx) /D 
			WAVE W_Coef
			//printf "%s\ry0=%g\rA=%g\rtau=%g\r\r",NameOfWave(cW),W_Coef[0],W_Coef[1],W_Coef[2]
			results[curTrace][0,2]=W_Coef[q]
			results[curTrace][3]=mean(auxW,lx,rx)
			
		ENDFOR
		print nameOfWave(results)
		KillWaves/Z W_Coef,Tr_refs
END // ExpFitAllTraces


FUNCTION DblExpFitAllTraces([includeString,excludestring])
// fit all traces between cursors
STRING		includeString,excludestring
// trace names containing these strings are excluded (included)
		VARIABLE noA=0, noB=0
		IF (strlen(Csrwave(a)) <= 0)
			noA=1
		ENDIF
		
		IF (strlen(Csrwave(b)) <= 0)
			noB=1
		ENDIF
		
		VARIABLE lx,rx
		
		STRING 	curTraceN, TraceNL=TraceNameList("",";",1)
		VARIABLE	totTraces=ItemsInList(TraceNL,";")
		VARIABLE	curTrace
		
		// take care of the include and exclude rules.
		
		IF (!ParamisDefault(includestring))
				FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
					curTraceN=StringFromList(curTrace, TraceNL,";")
					IF (!stringmatch(curTraceN,"*"+includeString+"*"))
						TraceNL=RemoveFromList(curTraceN, TraceNL,";")
						totTraces-=1
						curTrace-=1
					ENDIF
				ENDFOR
				totTraces=ItemsInList(TraceNL,";")
		ENDIF
		
		IF (!ParamisDefault(excludestring))
				FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
					curTraceN=StringFromList(curTrace, TraceNL,";")
					IF (stringmatch(curTraceN,"*"+excludestring+"*"))
						TraceNL=RemoveFromList(curTraceN, TraceNL,";")
						totTraces-=1
						curTrace-=1
					ENDIF
				ENDFOR
				totTraces=ItemsInList(TraceNL,";")
		
		ENDIF
		
		// store results in a wave named after the graph
		STRING 	windowName = WinName(0, 3)
		
		Make/O/N=(totTraces,6) $(windowName+"_dblexp")
		WAVE	results=$(windowName+"_dblexp")
		SetDimLabel 1,0,y0,results
		SetDimLabel 1,1,A1,results
		SetDimLabel 1,2,tau1,results
		SetDimLabel 1,3,A2,results
		SetDimLabel 1,4,tau2,results
		SetDimLabel 1,5,deltaV,results
		

		
		// convert the pruned tracelist into wavereferences
		MAKE/O/WAVE/N=(totTraces) Tr_refs
		
		Tr_refs[]=TraceNameToWaveRef("", StringFromList(p, TraceNL,";") )
		
		// go through and fit traces
		// this will currently NOT work, if ther eis an xWAVE
		STRING auxname
		
		FOR (curTrace=0; curTrace<totTraces; curtrace+=1)
			WAVE cW=Tr_refs[curTrace]
			// only for special case that info has to be obtained form another (parallel)
			// wave with related name (e.g. voltage belonging to the current
			// otherwise disable
			auxName=ReplaceString("_I", GetWavesDataFolder(cW,2), "_V")
		
			WAVE auxW=$auxName
			
			IF (noA)
				lx=DimOffset(cW,0)
			ELSE
				lx=min(xcsr(B),xcsr(A))

			ENDIF
			IF (noB)
				rx=rightx(cW)-Dimdelta(cW,0)
			ELSE
				rx=max(xcsr(B),xcsr(A))
			ENDIF
			CurveFit/Q/M=2/W=0/K={0.2079} exp_XOffset, cW(lx,rx) /D 
			WAVE W_Coef
			Redimension/N=(5) W_Coef
			W_Coef[4]=W_Coef[2]*2
			W_Coef[2]=W_Coef[2]/2
			W_Coef[1]=W_Coef[1]/2
			W_Coef[3]=W_Coef[1]
			CurveFit/G/N/Q/M=2/W=0/K={0.2079} dblexp_XOffset kwCWave=W_coef, cW(lx,rx) /D 
			W_Coef[2]=0.00012
			CurveFit/G/H="00100"/N/Q/M=2/W=0/K={0.2079} dblexp_XOffset kwCWave=W_coef, cW(lx,rx) /D 

			WAVE W_Coef
			//printf "%s\ry0=%g\rA=%g\rtau=%g\r\r",NameOfWave(cW),W_Coef[0],W_Coef[1],W_Coef[2]
			results[curTrace][0,4]=W_Coef[q]
			results[curTrace][5]=mean(auxW,lx,rx)
		ENDFOR
		print nameOfWave(results)
		KillWaves/Z Tr_refs //, W_Coef
END // DblExpFitAllTraces

			

FUNCTION DistributeAxes([spacinginpercent])
VARIABLE spacinginpercent

IF (ParamisDefault(spacinginpercent))
	spacinginpercent = 1
ENDIF

STRING AList=   AxisList("")
VARIABLE nAxes = ItemsInList(AList,";")

// get axis that are vertical (left or right)
VARIABLE aa
	STRING vAxes="", hAxes="", info, cAxis
	FOR (aa=0; aa< nAxes; aa++)
		cAxis = StringFromList(aa, AList,";")
		info =  axisinfo("",cAxis)
		info = StringByKey("AXTYPE",info,":",";")
		IF ( (StringMatch(info, "left" )) || (StringMatch(info, "right" )))
			vAxes+=cAxis+";"
		ELSE
			hAxes += cAxis+";"
		ENDIF
	ENDFOR
	// sort inverse:
	vAxes=SortList(vAxes,";",17)
	// now go through vertical axes and distribute them evenly in z
	VARIABLE	nvAxes=ItemsInList(vAxes,";")
	VARIABLE FractPerAxis=(1-(nvAxes-1)*spacinginpercent/100)/(nVAxes) // accounting for the nvAxis -1 gaps
	STRING RefhAxis = StringFromList(0,hAxes,";")
	FOR (aa=0; aa< nvAxes; aa++)
		cAxis = StringFromList(aa, vAxes,";")
			ModifyGraph axisEnab($cAxis)={aa*(fractPerAxis+spacinginpercent/100),aa*(fractPerAxis+spacinginpercent/100)+fractPerAxis}
			ModifyGraph freePos($cAxis)={0,$RefhAxis}
	
	ENDFOR
	ModifyGraph btLen=3

END

FUNCTION PlotTraces()
	// just to speed up plotting many waves, avoiding the cmd line messages


	STRING 		BaseName
	VARIABLE	n1,n2
	
	STRING	SourceWaveList=WaveList("*",";","DIMS:1,DIMS:2,MINCOLS:2" )
	
	Prompt 	BaseName, "Base Name of Wave names", popup, SourceWaveList
	Prompt	n1, "First index"
	Prompt	n2, "Last Index"
	DoPrompt "Please enter:", BaseName, n1, n2
	IF (V_flag==1) // cancel
		Return -1
	ENDIF
	VARIABLE	n
//	STRING		Title=GetDataFolder(0)+"_Waveforms"
	STRING		Title=BaseName+"_Waveforms"
	// differentiate between 2 cases:
	// the wave with the name given by basename does itself exist
	// this is the case of all waveforms being contained in the matrix S_x
	// in this case plot columns of the matrix given by the indicies
	// the traditional case is that indeed the individual waves exist and get plotted
	Display as Title
	DoWindow/C $Title
	IF (Exists(Basename)==1)		// if the basename itself points to a wave, this wave contails ALL the waveforms (in a matrix)
		WAVE base=$Basename
		IF (DimSize(base,1) <= max(n1,n2) )			
			STRING message=""
			message+="Max index for wave "
			message+=basename
			message+=" is "+ num2istr( Dimsize(base,1)-1)
			message+="."			
			DoAlert 0, message
		ENDIF
		IF ((n1==0)&&(n2==0))
			n1=0
			n2=DimSize(base,1)-1
		ENDIF
		FOR (n=min(n1,n2); n<=min(max(n1,n2),Dimsize(base,1)-1) ; n+=1)
			AppendToGraph base[*][n]
		ENDFOR
	ELSE						// no wave has the basename as a name, many waves with numerical suffixes contain the waveforms
		STRING 		CurrName
		FOR (n=min(n1,n2); n<=max(n1,n2); n+=1)
			CurrName=BaseName+num2istr(n)
			AppendToGraph $CurrName
		ENDFOR
	ENDIF
END

// **************************************************************************************
// *****************       aNNe tools for bulk - loading data files      ****************
// **************************************************************************************
Function/S DoOpenMultiFileDialog_aNNe()
// allows to select multiple files from a single directory and returns 
// a string with the files names ("\r" separated)
	Variable refNum
	String message = "Select one or more files"
	String outputPaths
	String fileFilters = "Data Files (*.txt,*.dat,*.csv,*.mat,*.abf,*.tif):.txt,.dat,.csv,.mat,.abf,.tif;"
	fileFilters += "All Files:.*;"

	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	outputPaths = S_fileName
	
	if (strlen(outputPaths) == 0)
		Print "Cancelled"
	else
		Variable numFilesSelected = ItemsInList(outputPaths, "\r")
		Variable i
		for(i=0; i<numFilesSelected; i+=1)
			//String path = StringFromList(i, outputPaths, "\r")
			//Printf "%d: %s\r", i, path	
		endfor
	endif
	
	return outputPaths		// Will be empty if user canceled
End

FUNCTION Load([Liste,PathName,IdentPosition, Prefix])
STRING		Liste, PathName, IdentPosition, Prefix
// this is an optional parameter. It is possible to call "Load()" and it works to go throught the hard wired lists
// if called "Load(Liste=".....") it will open the individual files named in the List
// the optional Parameter IdentPosition specifies where in the file name the identifier is located
// that should be used to name the loaded data
// a path can be specified in order to keep the Entries in Liste shorter
// IdentPosition = "last" or "first" or "all" 
// causes the first entries before the first "_" of the filename
// the last entries after the last "_" and before the "." 
// or the entire filename before the "."
// to be used as identifier
// Prefix = "Pre_" gives all names a prefix 
//
// a possible call is :
// Load(Liste=DoOpenMultiFileDialog_aNNe(), IdentPosition="last", Prefix="WT_")

	VARIABLE	k, ff, num
	STRING		CurrItem
	IF (ParamIsDefault(IdentPosition) )
		IdentPosition="all"
	ENDIF
	IF (ParamIsDefault(Prefix) )
		Prefix=""
	ENDIF



IF( ParamIsDefault(Liste) ) // no Liste handed over
	WAVE/T		Subjects=root:examples:subject_IDs
	num=DimSize(Subjects,0)	 			//number of subjects
	STRING	ResultWaveName, CurrFile, CCResultWaveName, SlidingCC, CurrExt
	STRING	NoteString
	
	FOR (k=0; k<num; k+=1)
		CurrItem=Subjects[k]
		CurrFile="E:AnalysisBerlin:meanCurves_repABR/"+CurrItem+"_ffr_ba.mat"
		ResultWaveName="ffr_ba_"+CurrItem
		MLLoadWave/O/Q/N=$ResultWaveName/C/G CurrFile
		 
		CurrFile="E:AnalysisBerlin:meanCurves_repABR/"+CurrItem+"_ffr_da.mat"
		ResultWaveName="ffr_da_"+CurrItem
		MLLoadWave/O/Q/N=$ResultWaveName/C/G CurrFile 
	ENDFOR
ELSE
	num=ItemsInList(Liste,"\r")
		FOR (k=0; k<num; k+=1)
			CurrFile=StringFromList(k, Liste,"\r")
			// getting
			ff= ItemsInList(CurrFile, ":") // how many folders deep - last is file name
			CurrItem=StringFromList(ff-1,CurrFile,":") // that's the filename w/o any path
			CurrExt= StringFromList(ItemsInList(CurrItem,".")-1, CurrItem,".") // this is the ending after the last "."

			// now construct the identifier out of the filename
			
			ff= ItemsInList(CurrItem,"_")		// how many segments of the name
			strswitch(IdentPosition)						// string switch
				case "all":
					ResultWaveName = CurrItem[0,strlen(CurrItem)-1-strlen(CurrExt)-1] // first character to just before last point
					break
				case "first":
					ResultWaveName=StringFromList(0,CurrItem,"_")
					break
				case "last":
					ResultWaveName=StringFromList(ff-1,CurrItem,"_")
					ResultWaveName=ResultWaveName[0,strlen(ResultWaveName)-1-strlen(CurrExt)-1]
					break
			endswitch
 			ResultWaveName=Prefix+ResultWaveName
			strswitch(CurrExt)						// string switch
				case "mat":	// 
					MLLoadWave/O/Q/R CurrFile
					FOR (ff=0; FF<ItemsInList( S_waveNames,";"); ff+=1)
						CurrFile=StringFromList(ff,S_WaveNames,";")
						Rename $CurrFile, $(CurrFile+"_"+CurrItem)
					ENDFOR
					S_WaveNames=""
					break
				case "tif":	// 
						ImageLoad/T=tiff/Q CurrFile
						// S_WaveNemes now contains the loaded wave and a trailing ";"
						Rename $(REmoveEnding(S_WaveNames)), $ResultWaveName
					break

				case "abf":	// axon binary files
					// the following command requires the path separator to be "\\" rather than ":"
					// but the spearator after the initial drive letter has to be a ":\\"
					CurrFile=ReplaceString(":", CurrFile, "\\\\")  
					CurrFile=ReplaceString("\\\\", CurrFile, ":\\\\",0,1)  // single replacement
					STRING cmdStr= "HT_ImportAbfFile(\""
					cmdStr+=CurrFile+"\",\"\",1,2,2,2,1) "  // the numbers have some relevance to the kind of data imported (gap-fee, episodic etc.
					print  cmdStr
					Execute cmdStr 
					S_WaveNames=""
					break
				case "dat":
					
					LoadWave/Q/N=L_Dummy /G CurrFile 
					FOR (ff=0; ff<ItemsInList( S_waveNames,";"); ff+=1)
						WAVE CurrWave=$(StringFromList(ff,S_WaveNames,";"))
						note/K CurrWave,"Loaded from "+CurrFile+"\r"
						note/NOCR CurrWave,"Wave #"+num2istr(ff)
						IF ((ff>0) || (Exists(ResultWaveName)==1))
							Rename CurrWave, $(ResultWaveName+"_"+num2istr(FF))
						ELSE
							Rename CurrWave, $(ResultWaveName)
						ENDIF
					ENDFOR
					S_WaveNames=""
					break
				case "txt":
					
					LoadWave/Q/N=L_Dummy /G CurrFile 
					FOR (ff=0; ff<ItemsInList( S_waveNames,";"); ff+=1)
						WAVE CurrWave=$(StringFromList(ff,S_WaveNames,";"))
						note/K CurrWave,"Loaded from "+CurrFile+"\r"
						note/NOCR CurrWave,"Wave #"+num2istr(ff)
						IF ((ff>0) || (Exists(ResultWaveName)==1))
							Rename CurrWave, $(ResultWaveName+"_"+num2istr(FF))
						ELSE
							Rename CurrWave, $(ResultWaveName)
						ENDIF
					ENDFOR
					S_WaveNames=""
					break
				case "csv":// for Charlie's Spectrograms, Nici's StopSignalTask
					
					LoadWave/Q/N=L_Dummy/L={0,0,0,1,1} /J /G CurrFile 
					FOR (ff=0; ff<ItemsInList( S_waveNames,";"); ff+=1)
						CurrFile=StringFromList(ff,S_WaveNames,";")
						IF (ff>0)
							Rename $CurrFile, $(ResultWaveName+"_"+num2istr(FF))
						ELSE
							Rename $CurrFile, $(ResultWaveName)
						ENDIF
					ENDFOR
					S_WaveNames=""
					break
				case "last":
					ResultWaveName=StringFromList(ff-1,CurrItem,"_")
					break
			endswitch
			
			
		ENDFOR
ENDIF



END

//*********************************************************************
//***************  AN tools for handling imaging stacks  **************
//*********************************************************************

FUNCTION CreateRGB_Int16(RedLayer,GreenLayer,BlueLayer, [RedMin, RedMax, GreenMin, GreenMax, BlueMin, BlueMax])
WAVE		RedLayer,GreenLayer,BlueLayer
VARIABLE	RedMin, RedMax, GreenMin, GreenMax, BlueMin, BlueMax
// check dimensions
VARIABLE		test= DIMSIZE(RedLayer,0)==DimSize(GreenLayer,0)
test = test && (DimSize(RedLayer,0) ==  DimSize(BlueLayer,0))
test = test && ( DIMSIZE(RedLayer,1)==DimSize(GreenLayer,1))
test = test && ( DIMSIZE(RedLayer,1)==DimSize(BlueLayer,1))

IF ( !test )
	Abort "x-y-Dimensions of the three input layers do not agree!"
ENDIF

IF (ParamIsDefault(RedMin))
	WAVESTATS/Q/M=1 RedLayer
	RedMin = V_Min
ENDIF
IF (ParamIsDefault(GreenMin))
	WAVESTATS/Q/M=1 GreenLayer
	GreenMin = V_Min
ENDIF
IF (ParamIsDefault(BlueMin))
	WAVESTATS/Q/M=1 BLueLayer
	BlueMin = V_Min
ENDIF
IF (ParamIsDefault(RedMax))
	WAVESTATS/Q/M=1 RedLayer
	RedMax = V_Max
ENDIF
IF (ParamIsDefault(GreenMax))
	WAVESTATS/Q/M=1 GreenLayer
	GreenMax = V_Max
ENDIF
IF (ParamIsDefault(BlueMax))
	WAVESTATS/Q/M=1 BLueLayer
	BlueMax = V_Max
ENDIF


test = (WaveDims(RedLayer) == 2)
test = test && (WaveDims(GreenLayer) == 2)
test = test && (WaveDims(BlueLayer) == 2)

IF ( !test )
	Abort "The three input layers are not all flat!"
ENDIF
// create outputStack


MAKE/D/O/N=(DimSize(BlueLayer,0), DimSize(BlueLayer,1),3) RGB_Stack,out
WAVE/D out;
WAVE/D RGB_Stack


// map data onto 0-2^16-1
out[][][0]=round( (RedLayer[p][q]-Redmin)/(Redmax-Redmin)*(2^16-1) )

out[][][1]=round( (GreenLayer[p][q]-Greenmin)/(GreenMax-Greenmin)*(2^16-1) )

out[][][2]=round( (BlueLayer[p][q]-Bluemin)/(Bluemax-Bluemin)*(2^16-1) )



// when custom chosen min and max values are used, the results could lie outside [0,4095]
	MatrixOP /O RGB_stack=clip(out,0,65535)

Redimension/U/W RGB_stack
//Print "Result is in RGB_stack"
//Display
//AppendImage out
KillWaves/Z out

END


// ================== miscellaneous tools ==========
FUNCTION condmean(source, cond, [relationsWaveName])
WAVE		source, cond
STRING	relationsWaveName

// returns the mean of the source, but only for a subset of entries 
// which entries are taken into account depends on cond
// if relations is not present, cond is assumed to be a boolean mask
// all indicies with non-zero entries in cond are used
// if relations is present, then each entry is expected to contain a 
// relation, i.e. a logical comparison with a missing left hand side
// i.e. "!=10" or ">5"
// all indicies are taken into account for which the entry in cond fullfills 
// ALL relations

Duplicate/O cond, dummyMaskTemp
WAVE dummyMaskTemp
dummyMaskTemp=1
Variable scale
IF (ParamIsDefault(relationsWaveName))
	dummyMaskTemp[]=cond[p]>0
	
	scale= sum(dummyMaskTemp)
	dummyMaskTemp*=source
	scale=sum(dummyMaskTemp)/scale
	KillWaves dummyMaskTemp
	
	Return 	scale
ENDIF

IF (exists(relationsWaveName ))
	WAVE/T relations=$relationsWaveName
ELSE
	DoAlert 0,relationsWaveName+" does not exist"
	Return NaN
ENDIF


VARIABLE	kk, target
STRING	rel, rel1, targetstr, Expression="([\\W][^-+]+)([[:digit:].-+]+)"
// () encloses a string
// \\W is anything not (underscore or any character that is a letter or digit)
// ^- is not "-"
// the last () is any decimal digit and the positive or negative sign
	
	FOR (kk=0; kk<DimSize(relations,0); kk+=1)
	rel1=relations[kk]
	
	SplitString /E="(^[=!<>]+)([[:digit:].-]+)$" rel1, rel , targetstr
	target = str2num(targetstr)
		strswitch(rel)						// string switch
			case ">":
				dummyMaskTemp[]*=(cond[p] > target)	
				break
			case "<":
				dummyMaskTemp[]*=(cond[p] < target) 	
				break
			case "==":
				dummyMaskTemp[]*=(cond[p] == target) 	
				break
			case "!=":
				dummyMaskTemp[]*=(cond[p] != target) 		
				break
			case ">=":
				dummyMaskTemp[]*=(cond[p] >= target) 		
				break
			case "<=":
				dummyMaskTemp[]*=(cond[p] <= target) 	
				break
		
		endswitch
	
	ENDFOR
	scale= sum(dummyMaskTemp)
	dummyMaskTemp*=source
	scale=sum(dummyMaskTemp)/scale
	KillWaves dummyMaskTemp
	
	Return 	scale

END

FUNCTION NormedCrossCorr(w1,w2)
WAVE		w1,w2
	// uses the built in Correlation
	// parameters are default i.e. linear correlation without any normalization or mean removal
	// afterwards the result is normalized such that every point is divided by the number of 
	// overlapping points that created that summed up to the result
	// the result is returned in w2, which is overwritten!
	
	
	VARIABLE n1=max(DimSize(w1,0),DimSize(w2,0)), , n2=min(DimSize(w1,0),DimSize(w2,0))
	
	Redimension/D w1,w2
	
	Correlate w1, w2
	
	w2[0,n2-1]/=p+1		// first part of partial overlap (including 1 point with maximal overlap)
	IF (n1>n2)
		w2[n2,n1-1]/=n2
	ENDIF
	w2[n1,n1+n2-2]/=(n1+n2-1)-p
END

//*************************************** MiSC statistics ******************+++


Function TraceScannerButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			DoWindow /F $(ba.Win)
			ModifyGraph hideTrace=0
			Button ShowAll disable=2
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function SwitchTraceProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	// update limits of the control (in case traces were added or deleted)
	// this is done before the action to avoid calls to non-existent traces
	VARIABLE nTraces=ItemsInList(TraceNameList(sva.win,";", 1 ))
	STRING cntrName=sva.ctrlName
	// Just in case there are double underscores:
	cntrName= ReplaceString("__", cntrName, "_")
	// extract number of traces to be shown and interval between them
	// those are the last two strings in the name, separated by underscores
	VARIABLE	nShown, interval
	VARIABLE	nUnderscores=ItemsInList(cntrName,"_")
	interval = str2num(StringFromList(nUnderscores-1, cntrName  ,"_"))
	nShown	= str2num(StringFromList(nUnderscores-2, cntrName  ,"_"))
	VARIABLE	rr


	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			DoWindow /F $(sva.Win)
			ModifyGraph hideTrace=1
			ModifyGraph hideTrace($(StringFromList(min(dval,nTraces-1),TraceNameList(sva.Win,";", 1 ))))=0
			FOR (rr=1; rr<nShown; rr+=1)
				// the following display more than 1 trace at a time, 
				// this is controlled via an intervalbetween the selected traces
				ModifyGraph hideTrace($(StringFromList(dval+interval*rr,TraceNameList(sva.Win,";", 1 ))))=0
			ENDFOR
			Button ShowAll disable=0
			break
		case -1: // control being killed
			break
	endswitch


	return 0
End

Function SwitchTraceProcTransp(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			// update limits of the control (in case traces were added or deleted)
			// this is done before the action to avoid calls to non-existent traces
			VARIABLE nTraces=ItemsInList(TraceNameList(sva.win,";", 1 ))
			STRING cntrName=sva.ctrlName
			// Just in case there are double underscores:
			cntrName= ReplaceString("__", cntrName, "_")
			// extract number of traces to be shown and interval between them
			// those are the last two strings in the name, separated by underscores
			VARIABLE	nShown, interval
			VARIABLE	nUnderscores=ItemsInList(cntrName,"_")
			interval = str2num(StringFromList(nUnderscores-1, cntrName  ,"_"))
			nShown	= str2num(StringFromList(nUnderscores-2, cntrName  ,"_"))
			
			VARIABLE	rr
			VARIABLE  highR,highG,highB,highO	// red green blue opacity of the highlighted trace(s)
			VARIABLE  lowR,lowG,lowB,lowO	// red green blue opacity of the suppressed trace(s)
			
			// get currently chosen trace number from dval of the control
			Variable newOne = sva.dval
			// get previously chosen trace number from userdata of the control
			Variable oldOne = str2num(sva.userdata)
			// now place current selection in userdata
			sva.userdata=num2istr(newOne)


			STRING	rgbString= GetUserData(sva.win, sva.ctrlName, "low" )
			lowR=str2num(StringFromList(0,rgbstring,","))
			lowG=str2num(StringFromList(1,rgbstring,","))
			lowB=str2num(StringFromList(2,rgbstring,","))
			IF (ItemsInList(rgbstring,",") <4)
				lowO=65535
			ELSE
				lowO=str2num(StringFromList(3,rgbstring,","))
			ENDIF

			rgbString= GetUserData(sva.win, sva.ctrlName, "high" )
			highR=str2num(StringFromList(0,rgbstring,","))
			highG=str2num(StringFromList(1,rgbstring,","))
			highB=str2num(StringFromList(2,rgbstring,","))
			IF (ItemsInList(rgbstring,",") <4)
				highO=65535
			ELSE
				highO=str2num(StringFromList(3,rgbstring,","))
			ENDIF

			DoWindow /F $(sva.Win)
			IF (newOne == oldOne) // when the end of the list is reached: reset all the traces 
			// because sometimes there are unexplained glitches, and traces are not properly updated
				ModifyGraph rgb=(lowR,lowG,lowB,lowO)

			ENDIF

			FOR (rr=0; rr<nShown; rr+=1)
				// the following display more than 1 trace at a time, 
				// this is controlled via an intervalbetween the selected traces
				ModifyGraph rgb[oldOne+interval*rr]=(lowR,lowG,lowB,lowO)
				ModifyGraph rgb[newOne+interval*rr]=(highR,highG,highB,highO)
			ENDFOR
			break
		case -1: // control being killed
			break
		default: // nothing happened, just mouse over
			// just update the color choices stored in userdata
			// extract rgba values for active trace (which is still the highlighted one)
			
			rgbstring=StringByKey("rgb(x)", TraceInfo(sva.win, StringFromList(sva.dval,TraceNameList(sva.Win,";", 1 )), 0 ),  "=", ";")
			// remove leading and trailing parentheses
			rgbstring=rgbstring[1,strlen(rgbstring)-2]+","

	endswitch


	return 0
End


			// extract rgba values for newOne trace (which is NOT the highlighted one)
			STRING rgbstring=StringByKey("rgb(x)", TraceInfo(sva.win, StringFromList(newOne,TraceNameList(sva.Win,";", 1 )), 0 ),  "=", ";")
			// remove leading and trailing parentheses
			rgbstring=rgbstring[1,strlen(rgbstring)-2]+","
			
						// extract rgba values for oldOne trace (which is still the highlighted one)
			rgbstring=StringByKey("rgb(x)", TraceInfo(sva.win, StringFromList(oldOne,TraceNameList(sva.Win,";", 1 )), 0 ),  "=", ";")
			// remove leading and trailing parentheses
			rgbstring=rgbstring[1,strlen(rgbstring)-2]+","


Function TransparencySliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa

	STRING TraceVarContrName=GetUserData(sa.win, sa.ctrlName, "" )
	
	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			Variable curval = sa.curval
			
			ControlInfo $TraceVarContrName
			VARIABLE selectedTrace=V_Value
			// extract rgba values for selected trace (which is the highlighted one)
			STRING rgbstring=StringByKey("rgb(x)", TraceInfo(sa.win, StringFromList(selectedTrace,TraceNameList(sa.Win,";", 1 )), 0 ),  "=", ";")
			// remove leading and trailing parentheses
			rgbstring=rgbstring[1,strlen(rgbstring)-2]+","
			VARIABLE  highR,highG,highB,highO	// red green blue opacity of the highlighted trace(s)
			highR=str2num(StringFromList(0,rgbstring,","))
			highG=str2num(StringFromList(1,rgbstring,","))
			highB=str2num(StringFromList(2,rgbstring,","))
			IF (ItemsInList(rgbstring,",") <4)
				highO=65535
			ELSE
				highO=str2num(StringFromList(3,rgbstring,","))
			ENDIF


			VARIABLE  lowR,lowG,lowB,lowO	// red green blue opacity of the suppressed trace(s)
			
			Variable Unselected = abs(selectedTrace -1) // for sure a different value than the selected one
			rgbstring=StringByKey("rgb(x)", TraceInfo(sa.win, StringFromList(Unselected,TraceNameList(sa.Win,";", 1 )), 0 ),  "=", ";")
			// remove leading and trailing parentheses
			rgbstring=rgbstring[1,strlen(rgbstring)-2]+","

			lowR=str2num(StringFromList(0,rgbstring,","))
			lowG=str2num(StringFromList(1,rgbstring,","))
			lowB=str2num(StringFromList(2,rgbstring,","))
			lowO=round(65535*curval)
			
			
			SetVariable $TraceVarContrName,userdata(low) =num2istr(lowR)+","+num2istr(lowG)+","+num2istr(lowB)+","+num2istr(lowO)+","
			SetVariable $TraceVarContrName,userdata(high)=num2istr(highR)+","+num2istr(highG)+","+num2istr(highB)+","+num2istr(highO)+","
			ModifyGraph rgb=(lowR,lowG,lowB,lowO)

			// prepare coloring the selected traces
			// need to find out how many
			
			TraceVarContrName= ReplaceString("__", TraceVarContrName, "_")	// remove possible accidental double underscores
			VARIABLE	nShown, interval
			VARIABLE	nUnderscores=ItemsInList(TraceVarContrName,"_")
			interval = str2num(StringFromList(nUnderscores-1, TraceVarContrName  ,"_"))
			nShown	= str2num(StringFromList(nUnderscores-2, TraceVarContrName  ,"_"))
			
			VARIABLE	rr


			FOR (rr=0; rr<nShown; rr+=1)
				// the following display more than 1 trace at a time, 
				// this is controlled via an intervalbetween the selected traces
				ModifyGraph rgb[selectedTrace+interval*rr]=(highR,highG,highB,highO)
			ENDFOR

			break
	endswitch

	return 0
End


// ***********************************************************************
// ** * * * * * * * *                                  * * * * * * * * **
// *** * * * * * * *   Statistics & Signal processing   * * * * * * * * **
// ** * * * * * * * *                                  * * * * * * * * **
// ***********************************************************************


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

FUNCTION	QuantilesFromSample(W_Sample,W_Quantiles,CertainThatNormal)
WAVE		W_Sample, W_Quantiles
INT		CertainThatNormal		// set to 1 if the underlying distribution is normal

// see Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, American Statistician 50, 361–365. 
// implemented here is method 8 from this paper as standard
// for CertainThatNormal, method 9 is implemented
// The resulting quantile estimates are approximately median-unbiased regardless of the distribution of x. 

		VARIABLE	num=DimSize(W_Sample,0)
		
		//	Some safety tests
		
		IF (num<4)
			DoAlert 0,"Quantile calculation requires more than 4 samples"
			Return -1
		ENDIF
		
		VARIABLE	nQ=DimSize(W_Quantiles,0), qq
		IF (nQ<1)
			DoAlert 0,"Quantile wave was empty"
			Return -2
		ENDIF
		
		Duplicate/O W_Sample, Sorted, p_k
		IF (CertainThatNormal)
			p_k[]=(p+1-3/8)/(num+1/4)		
		ELSE
			p_k[]=(p+1-1/3)/(num+1/3)
		ENDIF

		
		FOR (qq=0; qq<nQ; qq+=1)
			IF ( (W_Quantiles[qq]< p_k[0] ) || (W_Quantiles[qq]> p_k[num-1]  ) )
				DoAlert 0,"Cannot compute quantile "+num2str(W_Quantiles[qq])+" - sample size is too small."
				Return -3
			ENDIF
		ENDFOR
		
		
		
		WAVE Sorted
		WAVE	p_k
		
		
		Sort Sorted,Sorted
		
		W_Quantiles[]=interp(W_Quantiles[p],p_k,Sorted)
		KillWaves/Z Sorted,p_k
		
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
//	OutName=OutName[0,25]+"_splt" // using the liberal name length in IP 7
	OutName=OutName+"_splt"
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

// next copied out of removpoints.ipf, but threadsafe

Threadsafe Function thdsfRemoveNaNs(theWave)
	Wave theWave

	Variable p, numPoints, numNaNs
	Variable val
	
	numNaNs = 0
	p = 0											// the loop index
	numPoints = numpnts(theWave)				// number of times to loop

	do
		val = theWave[p]
		if (numtype(val)==2)					// is this NaN?
			numNaNs += 1
		else										// if not NaN
			theWave[p - numNaNs] = val			// copy to input wave
		endif
		p += 1
	while (p < numPoints)
	
	// Truncate the wave
	DeletePoints numPoints-numNaNs, numNaNs, theWave
	
	return(numNaNs)
End