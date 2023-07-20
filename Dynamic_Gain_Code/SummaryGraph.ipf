#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
MACRO CreateAvgSpikeByOnset(alignByRateThresh)
// to be called from INSIDE the cells folder
// this can easily be achieved for all folders by calling 
// XEQTinSubs("CreateAvgSpikeByOnset(YYY)")
VARIABLE	alignByRateThresh
	// the rate threshold at which the spike shapes are aligned 
	// to create an average spike
	
	VARIABLE/G rThresh=alignByRateThresh // a global copy to facilitate running the spike stats	
	// the range around the threshold that is taken into account is set by the number in the following line:
	// 'STAfromAnalogue(~,$W_Name, 8e-3'
	// means total duration is 8ms
	
	// re-create spike onset times
	// and update all spike statistics
	Xeqt4WList("OU*_V", "STRING/G W_Name=\"~\";W_Name=W_Name[0,strlen(W_Name)-2]+\"ST\"; SpikeStats(~,$W_Name,RateThresh=rThresh) ")
	MAKE/O/N=0 RPD_byWave
	Xeqt4WList("OU*_V","STRING/G W_Name=\"~\";W_Name=W_Name[0,strlen(W_Name)-2]+\"SonsetT\";STAfromAnalogue(~,$W_Name, 12e-3,0, Suffix=\"avgSpikeOn\");W_Name=W_Name[0,strlen(W_Name)-8]+\"avgSpikeOn\";STRING/G D_Name=W_Name+\"_DIF\" ;Differentiate/METH=0 $W_Name /D=$D_Name;Resample/RATE=100000 $W_Name,$D_Name; FindLevel/Q/EDGE=1/P $D_Name,rThresh;RPD_byWave[DimSize(RPD_byWave,0)]= {($D_Name[ceil(V_LevelX)]-$D_Name[ceil(V_LevelX)-1])/($W_Name[ceil(V_LevelX)]-$W_Name[ceil(V_LevelX)-1])};SetDimLabel 0,DimSize(RPD_byWave,0)-1,~,RPD_byWave")
	Xeqt4WList("OU*_avgSpikeOn","SetScale d 0,0,\"V\", ~")
	Xeqt4WList("OU*_avgSpikeOnDIF","SetScale d 0,0,\"V/s\", ~")
END

MACRO CreateAvgSpikeByVoltThreshold()
// to be called from INSIDE the cells folder
// this can easily be achieved for all folders by calling 
// XEQTinSubs("CreateAvgSpikeByVoltThreshold(YYY)")

	// it uses the spike time as originally defined when creating the *_ST waves
	// this is usually around 0 mV
	// to align spikes and average them
	
	// the range around the threshold that is taken into account is set by the number in the following line:
	// 'STAfromAnalogue(~,$W_Name, 8e-3'
	// means total duration is 8ms
	
	// re-create spike onset times
	// and update all spike statistics
	Xeqt4WList("OU*_V","STAfromAnalogue(~,§~rmvend_§_ST, 12e-3,0, Suffix=\"avgSpike\"); Differentiate/METH=0 §~rmvend_§_avgSpike /D=§~rmvend_§_avgSpike_DIF;Resample/RATE=100000 §~rmvend_§_avgSpike,§~rmvend_§_avgSpike_DIF")
	Xeqt4WList("OU*_avgSpike","SetScale d 0,0,\"V\", ~")
	Xeqt4WList("OU*_avgSpikeDIF","SetScale d 0,0,\"V/s\", ~")
END


FUNCTION SummaryGraph(FolderName, GraphNamePrefix,[Vmin,Vmax, dVmin,dVmax,GainMaxMag])
// to be called from inside the folder, but also providing the folders name
// this can be achieved using the keyword §SUB§ inside the XeqtInSubs command:
// XEQTinSubs("SummaryGraph(\"§SUB§\",\"NAMEPREFIX\",Vmin=-60e-3,Vmax=60e-3,dVmin=-100,dVmax=600,GainMaxMag=4.1e+11)")

// the summary layout can be created with the following lines (assuming the Title of summaries assigned below (GraphName)
// is  "SumTitle")
//Layout/T as "SummaryGraphs NAMEPREFIX"
//XEQTinSubs("AppendToLayout/T NAMEPREFIX_§SUB§; LayoutPageAction appendPage")

// it does require the average spike waveforms, to exist and to end in "avgSpikeOn"
STRING		GraphNamePrefix
STRING 	FolderName
VARIABLE 	vmin,Vmax, dVmin,dVmax,GainMaxMag

VARIABLE	ShowGains=1

IF (ParamIsDefault(Vmin))
	Vmin=-65e-3 // -65mV
ENDIF
IF (ParamIsDefault(Vmax))
	Vmax=-40e-3 // 40mV
ENDIF
IF (ParamIsDefault(dVmin))
	dVmin=-300 // 300V/s
ENDIF
IF (ParamIsDefault(dVmax))
	dVmax=600 // 600V/s
ENDIF
IF (ParamIsDefault(GainMaxMag))
	GainMaxMag=2.5e+11 // 250 Hz/nA
ENDIF
	
	
	
	STRING GraphName=GraphNamePrefix+"_"+FolderName
	
	Display/W=(200,50,450,300) as GraphName; DoWindow/C $(GraphName);
	STRING rememberFolder=GetDataFolder(1)
	// go through all datafolders of the same level and plot the gains from those folders
	VARIABLE	grey=round(0.5*(2^16-1))
	IF (ShowGains)
		SetDataFolder ::
		XeqtInSubs("AppendToGraph/L=GainLeft/B=GainBottom Gain_avg_scaled_MgFlt vs FreqPoints;")
		ModifyGraph rgb=(grey,grey,grey,(2^16-1))	
	ENDIF
//	ModifyGraph rgb=(0,0,0,round(0.5*(2^16-1)))		// use transparency --> problems after pdf export
	VARIABLE tt,nGreyTraces=ItemsInList(TraceNameList("", ";",1 ),";")
	SetDataFolder $rememberFolder
	// now plot avg spike waveforms
	Xeqt4WList("*_avgSpikeOn" ," AppendToGraph/L=VertCrossing/B=HorizCrossing  ~_DIF vs ~;AppendToGraph/L=SpikeLeft/B=SpikeBottom ~[2e-3/DimDelta(~,0),DimSize(~,0)-1] ; MOdifyGraph rgb[ItemsInList(TraceNameList(\"\", \";\", 1 ) ,\";\")-1]=(0,0,0);MOdifyGraph rgb[ItemsInList(TraceNameList(\"\", \";\", 1 ) ,\";\")-2]=(0,0,0)") 
	IF (ShowGains)
		WAVE Magnitude=Gain_avg_scaled_MgFlt
		WAVE FreqPoints=FreqPoints
		IF (!WaveExists(Magnitude))
			DoAlert 0,"Missing gain wave "+GetWavesDataFolder(Magnitude,2)
		ENDIF 
		IF (!WaveExists(FreqPoints))
			DoAlert 0,"Missing FreqPoints wave "+GetWavesDataFolder(FreqPoints,2)
		ENDIF 
		AppendToGraph/L=GainLeft/B=GainBottom Magnitude vs FreqPoints
		SetAxis GainBottom *,1000
		SetAxis GainLeft 0,GainMaxMag
	  	DO
	  		ModifyGraph rgb[tt]=(grey,grey,grey,(2^16-1))
	  		tt+=1
	  	WHILE (tt<= nGreyTraces-1)
	  	// now the last gain trace - the one from this folder
	  	ModifyGraph rgb[tt]=(0,0,0,1*(2^16-1))
		ModifyGraph freePos(GainLeft)={0,GainBottom}
		ModifyGraph freePos(GainBottom)={0,GainLeft}
		ModifyGraph axisEnab(GainLeft)={0.55,1}
		ModifyGraph axisEnab(GainBottom)={0.35,1}
		ModifyGraph log(GainBottom)=1,lblPos(GainLeft)=35,prescaleExp(GainLeft)=-9
		Label GainLeft "Gain(Hz/nA)";DelayUpdate
		ModifyGraph lblPos(GainBottom)=35;DelayUpdate
		Label GainBottom "Input Freq (Hz)"
	ELSE
		ModifyGraph rgb=(0,0,0,0.6*(2^16-1))

	ENDIF
	SetAxis HorizCrossing Vmin,Vmax
	SetAxis VertCrossing dVmin,dVmax;
	SetAxis SpikeLeft Vmin,Vmax
	// color all except the first nGreyTraces in another color (black)
	VARIABLE nTraces=ItemsInList(TraceNameList("", ";",1 ),";")
	tt=nGreyTraces


	ModifyGraph btLen=3
	ModifyGraph freePos(VertCrossing)={0,HorizCrossing}
	ModifyGraph freePos(HorizCrossing)={0,VertCrossing}
	ModifyGraph freePos(SpikeLeft)={-0.003,SpikeBottom}
	ModifyGraph freePos(SpikeBottom)={0,SpikeLeft}
	ModifyGraph axisEnab(VertCrossing)={0,0.45}
	ModifyGraph axisEnab(HorizCrossing)={0,0.45}
	ModifyGraph axisEnab(SpikeLeft)={0.55,1}
	ModifyGraph axisEnab(SpikeBottom)={0,0.2}
	ModifyGraph axisEnab(VertCrossing)={0,0.45}

	ModifyGraph noLabel(SpikeBottom)=2
	ModifyGraph lblPos(SpikeLeft)=35,prescaleExp(SpikeLeft)=3
	Label SpikeLeft "\U";DelayUpdate
	STRING Annotation="\\Z10trial\trate\twidth\tLV"
			Annotation+="\r"
			Annotation+="\t (Hz)\t (ms)\t"
	STRING AllWaves= SortList( WaveList("*_ST", ";", "" ),";",16)
	String currentw
	VARIABLE	avgLV=0, avgRate=0, avgWidth=0, avgPRR=0
	VARIABLE	LV,Rate, Width, PRR
	
	VARIABLE	ww, nW=ItemsInList(AllWaves,";")
	ww=0
	DO
		currentw=StringFromList(ww,AllWaves,";")
		rate=1/str2num(StringByKey("ISI", note($currentw),":"," " ))
		LV=str2num(StringFromList(0,StringByKey("(LV)", note($currentw),":"," " ),";"))
		width=1000*str2num(StringFromList(0,StringByKey("FWHM", note($currentw),":"," " ),";"))
		PRR=str2num(StringByKey("(PRR)", note($currentw),":"," " ))
		Annotation+="\\Z08\r"
		Annotation+=currentw
			Annotation+="\t"
		Annotation+=num2str(round(100*rate)/100)
			Annotation+="\t"
		Annotation+=num2str(round(100*width)/100)
			Annotation+="\t"
		Annotation+=num2str(round(100*LV)/100)
		avgLV+=LV
		avgRate+=Rate
		avgWidth+=width
		avgPRR+=PRR
		ww+=1
	While (ww<nw)
		Annotation+="\\Z08\r—————————————————————\r"
		Annotation+="AVGs"
			Annotation+="\t"
		Annotation+=num2str(round(100*avgrate/ww)/100)
			Annotation+="\t"
		Annotation+=num2str(round(100*avgwidth/ww)/100)
			Annotation+="\t"
		Annotation+=num2str(round(100*avgLV/ww)/100)
	
	TextBox/C/N=AnnoText/X=23.75/Y=-19.38/F=0/A=MC/T={87,118,154,197,216,234,252,288,324,360}  Annotation
	TextBox/C/N=Title/B=1/X=32.00/Y=48.00/F=0/A=MC FolderName

	// in case ACL desktops are used, put the window to the desktop 2

	SetWindow kwTopWin,userdata(ACL_desktopNum)=  "2"
	ModifyGraph width=396,height=432

END
