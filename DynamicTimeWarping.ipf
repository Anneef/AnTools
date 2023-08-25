#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

FUNCTION GetDTW(xSeq,ySeq [, doPath])
WAVE	xSeq, ySeq
VARIABLE		doPath		// set to 1, if warping path is also required

NVAR globRange =V_Range

MAKE/O/N=(2) wp
IF (ParamIsDefault(doPath))
	doPath = 0
ENDIF

IF (doPath)
	wp[0]=1
ELSE
	wp[0]=0
ENDIF

VARIABLE/G dtw_dist=0
VARIABLE loc_dist

	[loc_dist,wp]=DTW(xSeq,ySeq,range=globRange)
	dtw_dist = loc_dist
	IF (!doPath)
		KilLWaves/Z wp
	ENDIF
	
	return dtw_dist
	
END


FUNCTION [VARIABLE dtw_dist,WAVE warpingpath] DTW(WAVE xSeq,WAVE ySeq[, VARIABLE range])
//%DTW dynamic time warping for multidimensional time series
//%
//% Input
//% xSeq:  [n × d]               d dimensional time series of length n
//% ySeq:  [m × d]               d dimensional time series of length m
//%
//% Output
//% d:  [1 × 1]               dtw(xSeq,ySeq) with local Euclidean distance
//% p:  [L × 2]   (optional)  warping path of length L
//%
//%
//% Authorof the matlab code David Schultz, DAI-Lab, TU Berlin, Germany, 2016 
// ported to IGOR PRO by Andreas Neef 2020
// _aNNe_ introduced the limited warping range r to save time and 
// 		  avoid spurious alignments
// 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	IF (DimSize(xSeq,1) != DimSize(ySeq,1) )
		Abort ("DTW distance cannot be calculated if entries have different dimensionality")
	ENDIF
	
	INT DoPath=0
	IF (warpingpath[0]>=0)
		DoPath = 1
	ENDIF
	
	VARIABLE N=DimSize(xSeq,0), d=max(1,DimSize(xSeq,1) ) 
	VARIABLE M=DimSize(ySeq,0)
	
	VARIABLE nn,mm
	IF (ParamIsDefault(range) ) 
		range = max (M,N)				// that max will be checked at every application, not to 
											// exceed the meaningful range
	ENDIF
	range = max(range, abs(N-M)) // range has to at least allow to warp from M to N in order to
										  // map the sequences onto each other
	MAKE/O/N=(N,M) Dist; 
	WAVE Dist
	IF (d==1)
		Dist=inf

		Dist[0][0] = (xSeq[0]-ySeq[0])^2;
		//    C[1,N-1][0] = C[p-1][0] + (xSeq[p]-ySeq[0])^2
		//    C[0][1,M-1] = C[1][q-1] + (xSeq[0]-ySeq[q])^2
		//    C[1,N-1][1,M-1] = (xSeq[p]-ySeq[q])^2 + min(min(C[p-1][q-1],C[p-1][q]),C[p][q-1])
		// Alternative:
		FOR (nn=1; nn<min(range,N); nn++)
			Dist[nn][0] = Dist[nn-1][0] + (xSeq[nn]-ySeq[0])^2;
		ENDFOR
		FOR (mm=1; mm<min(range,M); mm++)
			Dist[0][mm] = Dist[0][mm-1] + (xSeq[0]-ySeq[mm])^2;
		ENDFOR
		FOR (nn=1; nn<N; nn++)
			FOR (mm=max(1,nn-range); mm<min(M,nn+range-1); mm++)
				Dist[nn][mm] = (xSeq[nn]-ySeq[mm])^2 + min(min(Dist[nn-1][mm-1],Dist[nn-1][mm]),Dist[nn][mm-1]);
			ENDFOR
		ENDFOR
	
	ELSE
		// For multivariate sequences, the distance has to be taken over all entries
		// here, the warping will be the same for all dimensions
		// this is the case in which the different entries all depend on time in the same 
		// way and are not independently varying over time
	
		// original function uses pdist2
		//	 D = pdist2(xSeq,ySeq).^2;
		//    % for a n1 
		//    C(:,1) = cumsum(D(:,1));
		//    C(1,:) = cumsum(D(1,:));
		//    for n = 2:N
		//        for m = 2:M
		//            C(n,m) = D(n,m) + min(min(C(n-1,m-1),C(n-1,m)),C(n,m-1));
		//        end
		//    end 
	
	
		// suggested  replacement for repmat
		// tmpB(1,:,:)=ySeq';
		//  tmpx=repmat(tmpB,[size(xSeq,1),1,1]);
		//  tmpy=repmat(xSeq,[1,1,size(tmpx,3)]);
		//  D(:,:)=sum(abs(tmpx-tmpy),2);	
		MAKE/O/N=(N,M,d) Temp
		WAVE Temp
		Temp=inf
		FOR (nn=0; nn<N; nn++)
			FOR (mm=max(0,nn-range); mm<min(M,nn+range); mm++)
				Temp[nn][mm][0,d-1]=xSeq[nn][r]-ySeq[mm][r] 
			ENDFOR
		ENDFOR
		 
		// Dist [p][q][r] now contains the distance between point p in xSeq and point q in ySeq 
		// the r-dimension is not shifted between the sequences
		// 
		// now add the differences along the dimensions 
		// and square the difference
	
		MatrixOP/O 	EucDist=abs(Temp)
		MatrixOP/O 	Temp=sumBeams(EucDist))
		MatrixOP/O 	EucDist=magSqr(Temp)
		KillWaves/Z Temp
		

		// now D should have dimension N × M
		//      C[0][0] = D[0][0]
		//      
		//      C[1,N-1][0] = C[p-1][0] + D[p][0]
		//      C[0][1,M-1] = C[1][q-1] + D[0][q]
		//      C[1,N-1][1,M-1] = D[] + min(min(C[p-1][q-1],C[p-1][q]),C[p][q-1]);
		
		Dist[0][0] = EucDist[0][0]
		FOR (nn=1; nn<min(range,N); nn++)
			Dist[nn][0] = Dist[nn-1][0] + (EucDist[nn][0])^2;
		ENDFOR
		FOR (mm=1; mm<min(range,M); mm++)
			Dist[0][mm] = Dist[0][mm-1] + (EucDist[0][mm])^2;
		ENDFOR
		FOR (nn=1; nn<N; nn++)
			FOR (mm=max(1,nn-range); mm<min(M,nn+range-1); mm++)
				Dist[nn][mm] = (EucDist[nn][mm])^2 + min(min(Dist[nn-1][mm-1],Dist[nn-1][mm]),Dist[nn][mm-1]);
			ENDFOR
		ENDFOR
	
	
	ENDIF
	dtw_dist = sqrt(Dist[N-1][M-1])
	
	//################### calculate warping path ######################
	IF (DoPath)
		nn = N-1
		mm = M-1
		// initializing at the target points
		
		//p = zeros(N+M-1,2);
		VARIABLE npts=N+M-1 // the maximal conceivable length of the path
		MAKE/O/N=(npts,2) warpingpath
		warpingpath = 0
		// for every step, the target path has a target index in xSeq (nn) and a target index in ySeq
		// at the last point of the target path, the final indicies of the two input sequences
		// are reached
		
		// one problem is to deal with the different index offsets in MATLAB (1 to N) versus any
		// sensible language (IGOR, Python, ...)
		// the final indicies in IGOR have to be M-1 and N-1 
		warpingpath[npts-1] = {{nn},{mm}}
		
		// now step back until 
		VARIABLE C_diag, C_r, C_d
		VARIABLE k = 1  // k is used to count back from the last entry
		
		
		IF (nn+mm > 0) // that condition means "both indicies reached the zero-th entry
		DO 
			k = k + 1		// k=1 would target the last entry, but that is already filled
	    					// so k=2 is the first entry we have to care about
			IF (nn == 0)
				mm = mm-1;
			ELSEIF (mm == 0)
				nn = nn-1;
			ELSE
				C_diag = Dist[nn-1][mm-1]
				C_r = Dist[nn][mm-1]
				C_d = Dist[nn-1][mm]
				IF (C_diag <= C_r)
					IF (C_diag <= C_d)	// C_diag is smallest --> take diagonal step, i.e.
						nn=nn-1			// one step along xSeq-sequence
						mm=mm-1			// and one along ySeq-sequence
					ELSE					// i.e. C_r >= C_diag > C_d
						nn=nn-1			// only advance along xSeq
					ENDIF
				ELSEIF (C_r <= C_d)		// i.e. C_r < C_diag && C_r <= C_d
					mm=mm-1				// only advance along ySeq
				ELSE						// i.e. C_d < C_r < C_diag 
					nn=nn-1				// only advance along xSeq
				ENDIF
			ENDIF
		        warpingpath[npts-k]={{nn},{mm}} 
	    WHILE (nn+mm > 0)  // that condition means "both indicies reached the zero-th entry
	    ENDIF
	    // now the last step: 
	    // the length of the path, until it reached the first entry in both, xSeq and ySeq, is k
	    
	    // lets keep the last k points of the path, which is currently M+N-1 points long
	    DeletePoints/M=0 0,(npts-k), warpingpath
	ELSE
		Redimension/N=(0) warpingpath
	ENDIF
END  // DTW - dynamic time warping


FUNCTION/WAVE GetSSG(matchstr, dofrechet,[nEpochs, eta, initindex, InitWave])
STRING	matchstr
VARIABLE	nEpochs
WAVE		eta
VARIABLE	InitIndex, dofrechet
WAVE 		InitWave

		STRING 		AllMAtches=WaveList(matchstr,";","")
    	WAVE/WAVE	X_waveref = ListToWaveRefWave(AllMAtches)
    	
		VARIABLE		Nseqs = DimSize(X_waveref,0)	// number of samples
		
		
		IF (ParamIsDefault(nEpochs) )
		    nEpochs = ceil(1000/Nseqs);
		ENDIF
		
		IF (ParamIsDefault(eta))
			 Variable Neta=max(1000,Nseqs)
		    MAKE/O/N=(Neta) eta
		    eta[] = 0.1-(0.1-0.01)/(Neta-1)*p
		ENDIF
		
		IF (ParamIsDefault(InitWave))			// no initial proxy for central sequence provided
		 	IF (ParamisDefault(InitIndex))	// no index provided to identify which one of the sample sequences
		 												// should serve as initial proxy
		 	// have to select a random index
		 		WAVE	InitWave = X_Waveref[round(abs(enoise(Nseqs-0.51)))]
			ELSE	// no wave directly provided but index provided
			
				WAVE 	InitWave = X_Waveref[InitIndex]
			ENDIF
		ENDIF
		MAKE/O/N=(1) frechetvars
		
		Duplicate/O InitWave, zentral
		    
		[zentral, frechetvars] = SSG(X_waveref, nEpochs, eta, dofrechet)
		Note/K zentral, "central wave to the sample\r"+ReplaceString(";",WaveRefWaveToList(X_waveref, 0), "\r")
		Return zentral // this is now the zentral wave
END


FUNCTION [WAVE zentral,WAVE frechetvars] SSG(WAVE/WAVE X_waveref, VARIABLE nEpochs, WAVE eta, VARIABLE CalcFrechet[, VARIABLE WeightBySimilarity])
// SSG Computes an approximate sample mean under dynamic time warping
// 
// SSG aims at finding a sample mean under dynamic time warping of  
// sample X = (x_1,...,x_N) containing time series x_i, each of arbitraty 
// length and uniform dimension d. Similarity between the
// d-dimensional data points x_ij is measured by means of the Euclidean 
// distance. Other local distances are not supported.
//
// Mathematically this can be formulated as an optimization problem.
// A sample mean is a time series z that satisfies 
//
//    F(z) = min F(x),
//            x
//
// where F is the Frechet function, defined by
//
//                  N
//    F(x) = (1/N) sum dtw^2(x,x_i).
//                 i=1
//
// 
// SSG is a stochastic subgradient method that aims at minimizing F.
// For this, the algorithm requires to choose an initial (mean candidate) 
// sequence. Details are discussed in [1].
// 
//
//
// Input
//   X_waveref        Wave containing references to the waves representing 
//							the time series x_i of size [n_i,d], where
//                   n_i is the length of time series x_i and d is the 
//                   uniform dimension of all time series.
//                   example:   X{1} = x_1;  // size [n_1,d]
//                              X{2} = x_2;  // size [n_2,d]
//
//   nEpochs         Number of epochs.
//   (o)             default: nEpochs = 1
//                   Hint: One epoch is often sufficient for larger
//                   datasets (approximately N > 300). For smaller datasets
//                   nEpochs should be decreased
//
//   eta         		wave representiong a learning rate schedule
//   (o)             At update k the learning rate is set to eta(k).
//                   If k > length(eta), then eta(end) is used.
//                   default: eta = linspace(0.1, 0.005, N)
//                            i.e. linear decreasing from 0.1 to 0.005
//                            within the first epoch
//
//   initindex    	(1) default: If initSequence is not defined, then a  
//   (o)                 random sample of X is selected as the initial mean
//                       candidate
//                   (2) If initSequence is the integer i > 0, then x_i is  
//                       choosen as the initial sequence
//                   (3) If initSequence <= 0, then the medoid of
//                       X is selected. 
//                       CAUTION: This may take very long on larger
//                                datasets since O(N^2) dtw distances have 
//                                to be computed, each consisting of 
//                                O(n_i*n_j) distance  computations of
//                                d-dimensional vectors.
//                   (4) If initSequence is a d-dimensional time series of 
//                       length n, then it is used as initial sequence.
//
// Output
//   z               Approximate mean time series of size [n x d], 
//                   where n is the length of the initial sequence.
//
//   f               Vector containing Frechet variations f(k) := F(z^(k)),
//   (o)             where z^(k) is the mean at epoch k.
//                   f(1) is the Frechet variation of the initial sequence
//               !!! CAUTION: Returning f doubles the computation time!
//                            Do not ask for f if not necessary.
//							Because IGOR does currently not allow for optional
//							return objects, the entries in wave frechet, as it is handed 
//							to this function, are used to determine whether Frechet variations
//							will be computed
//		weight by similarity 
//							When the centroid is updated, the degree of update (the impact of the current sample sequence ) 
//							is weighted with the Frechet Distance of that sample sequence  within all the sample sequences
//							With this option, outliers get less influence on the centroid
// (o) = optional
//__________________
//
// Credits
// * The SSG mean algorithm was introduced in [1].
//
// [1] Schultz and Jain - "Nonsmooth Analysis and Subgradient Methods
//     for Averaging in Dynamic Time Warping Spaces"
//
//
// Author of the matlab code David Schultz, DAI-Lab, TU Berlin, Germany, 2016
// ported to IGOR PRO by Andreas Neef 2020
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		NVAR globRange =V_Range
		VARIABLE		Nseqs = DimSize(X_waveref,0)	// number of samples
		VARIABLE    N_of_eta= DimSize(eta,0)
		WAVE 			CurrX = X_waveref[0]				// reference the first wave to 
		VARIABLE		d = DimSize(CurrX,1);      	// find out the dimensionality of the data
		VARIABLE		npts = DimSize(zentral,0)		// length of the zentral series
		VARIABLE		optFrechet							// store minimal value of Frechet distance ever achieved
		
		IF (CalcFrechet)
			Redimension/N=(nEpochs+1) frechetvars
			frechetvars = 0
			frechetvars[0] = Frechet(zentral,X_waveref)
			optFrechet = frechetvars[0]
		ENDIF
		
		IF (weightBySimilarity)
			MAKE/N=(Nseqs) FrechetDist
		ENDIF 
		
		// stochastic subgradient optimization
		VARIABLE k, nn
		VARIABLE	c // cumulative count of iterations performed
		DUPLICATE/O zentral, subgradient, zentral_opt	// zentral_opt stored zentral that achieved minimal frechet
		VARIABLE Dist
		MAKE/O/N=(DimSize(zentral,0),2) Path
		MAKE/O/W/N=(2) ValenceM
		MAKE/O/B/N=(2,2)  WarpM
		
		FOR  (k = 0;  k < nEpochs; k++)
		
		//		create random permutation of sample indicies
				MAKE/O/N=(NSeqs) indc=p, randInc=enoise(1,2)
				SORT indc, randInc
				KillWaves/Z randInc
				FOR (nn=0; nn < Nseqs; nn++)
					WAVE currX = X_waveref[indc[nn]]
		//    a random sequence is selected
		
		
				[Dist,Path] = dtw(zentral,currX,range=globRange);
		
		      [WarpM,ValenceM] = WarpingAndValenceMatrix(Path)
		      Duplicate/O zentral, subgradient
		         
				subgradient[0,npts-1][0,d-1] *= 2*ValenceM[p]
				MatrixOP/O subgradient=subgradient - 2 * WarpM x currX
		//        
		//        % determine learning rate for update no. c
			        c = (k)*Nseqs+nn;
		//	        print c, eta[min(c,N_of_eta-1)]
		// learnrate = eta[min(c,N_of_eta-1)] i.e. constant learnrate after the first N_of_eta iterations
		      zentral = zentral - eta[min(c,N_of_eta-1)] * subgradient; 	
				ENDFOR    
				DoUpdate/W=	Graph5															
				DoUpdate/W=	Graph7															
 
				IF (CalcFrechet)
					frechetvars[k+1] = Frechet(zentral,X_waveref)
					IF (frechetvars[k+1] < optFrechet)
						optFrechet = frechetvars[k+1]
						zentral_opt = zentral
					ENDIF
				ENDIF
				//printf "runde %g of %g\r",k+1,nEpochs
		ENDFOR 
		KillWaves/Z WarpM,ValenceM

END  // SSG - stochastic subgradient method

//
//
FUNCTION Frechet(xseq, X_waveref)
WAVE			xseq
WAVE/WAVE 	X_waveref
	
	NVAR globRange =V_Range
	VARIABLE		dist, Nseqs = DimSize(X_waveref,0)	// number of samples
	MAKE/O		wp={-1} // singalling to the dtw function, that a path should not be computed
	VARIABLE 	f = 0, nn
	//  sum dtw distance over all waves in sample
	FOR (nn=0; nn < Nseqs; nn++)
		WAVE 	currW=X_waveref[nn]
			Redimension/N=1 wp
			wp[0]=-1
	    	[dist,wp] = dtw(xseq,currW,range=globRange)
	    	f = f + dist^2
	ENDFOR 
	
	f/=Nseqs
	Return f
	KillWaves/Z wp
END // Frechet distance end
//
FUNCTION CharacterizeSamples(matchstr)
STRING	matchstr	
		STRING 		AllMAtches=WaveList(matchstr,";","")
    	WAVE/WAVE	X_waveref = ListToWaveRefWave(AllMAtches)
    	MAKE/O/N=1 $("Medio_of_"+ReplaceString("*",matchstr,""))
		WAVE MedoidSeq = $("Medoid_of_"+ReplaceString("*",matchstr,""))
		VARIABLE		Nseqs = DimSize(X_waveref,0)	// number of samples
		MAKE/O/N=(NSeqs,Nseqs) $("CrossDist_of_"+ReplaceString("*",matchstr,""))
		WAVE CrossDist = $("CrossDist_of_"+ReplaceString("*",matchstr,""))
		[MedoidSeq, CrossDist]= DTWCharacteristic(X_waveref)
END	
FUNCTION [WAVE MedoidSeq, WAVE CrossDist] DTWCharacteristic(WAVE/WAVE X_waveref) 
// MEDOIDSEQUENCE returns medoid of X
//  A medoid is an element of X that minimizes the Frechet function
//  among all elements in X
// CrossDist Returns a N x N Matrix with pairwise DTW distances  
// the frechet Distance of a sequence relative to the distribution is the 
// sum of all squared distances

	NVAR globRange =V_Range
	VARIABLE		dist, Nseqs = DimSize(X_waveref,0)	// number of samples
	MAKE/O		wp={-1} // singalling to the dtw function, that a path should not be computed
	VARIABLE 	f, f_min = inf , n1,n2
	Redimension/S/N=(NSeqs,NSeqs) CrossDist
	CrossDist = 0
	// fill cross distance matrix
	 FOR (n1=0; n1 < NSeqs; n1++)
		 FOR (n2=n1; n2 < NSeqs ; n2++)
			WAVE 	currW1=X_waveref[n1]
			WAVE 	currW2=X_waveref[n2]
			wp[0]={0}
			[dist,wp] = dtw(currW1,currW2,range=globRange)
			CrossDist[n1][n2] = dist
		ENDFOR
	ENDFOR
// fill other triangle
	MatrixOP/O 	CrossDist = CrossDist + CrossDist^t
	MatrixOP/O FrechetFromCross = magSqr(sumRows(CrossDist))
	Redimension/N=(-1,0)  FrechetFromCross
	WAVEStats/Q/M=1 FrechetFromCross
	Duplicate/O X_waveref[V_minLoc], MedoidSeq 
	NOTE/K MedoidSeq,"Medoid sequence\rIndex:"+num2istr(V_minLoc)+";\r"+ ReplaceString(";",WaveRefWaveToList(X_waveref, 0), "\r")
	KillWaves/Z wp
END
//

FUNCTION GetWarpingAndValenceMatrix(WM, VM, wp)
WAVE	WM, VM, wp

[WM, VM]= WarpingAndValenceMatrix( wp)
END



FUNCTION [WAVE WarpM, WAVE ValenceM] WarpingAndValenceMatrix(WAVE warpPath)
//% W is the (sparse) warping matrix of p
//% V is a vector representing the diagonal of the valence matrix 

VARIABLE len = DimSize(warpPath,0)

VARIABLE N = warpPath[len-1][0] +1 // Length of the first wave warped (in the case of SSG, thats the central one) 
VARIABLE M = warpPath[len-1][1] +1 // Length of the second wave warped
    
// original command: 
// W = sparse(p(:,1),p(:,2),ones(L,1),N,M);
// Igor does not support special coding for sparse 

	 // warpM contains only values ZERO or One,
	// so BYTE (/B) is an apropriate format for now
	// there is nothing less memory consuming
	MatrixOP/O WarpM = zeroMat(N,M,8) // 8 = 2^3 setting bit 3 to 1 means this is a BYTE wave
	
	VARIABLE ll
	
	FOR (ll=0; ll < len; ll++)
		WarpM[warpPath[ll][0]][warpPath[ll][1]] = 1
	ENDFOR
	 
// S = sparse(i,j,v) generates a sparse matrix S from the triplets i, j, and v such that S(i(k),j(k)) = v(k). 
// The max(i)-by-max(j) output matrix has space allotted for length(v) nonzero elements.
// If the inputs i, j, and v are vectors or matrices, they must have the same number of elements. Alternatively, the argument v and/or one of the arguments i or j can be scalars.
// example
// S = sparse(i,j,v,m,n) specifies the size of S as m-by-n.
//
	REDIMENSION/U/W ValenceM
	MAtrixOP/O ValenceM=sumRows(WarpM) 	// need to check how the results compare to Matlab!
																		// because the 
END // WarpingAndValenceMatrix
