#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../include/defs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"


typedef unsigned char	BitSequence;

BitSequence	*epsilon;

double
Pr(int u, double eta)
{
	int		l;
	double	sum, p;
	
	if ( u == 0 )
		p = exp(-eta);
	else {
		sum = 0.0;
		for ( l=1; l<=u; l++ )
			sum += exp(-eta-u*log(2)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
		p = sum;
	}
	return p;
}

double psi2(int m, int n, BitSequence const *epsilon) {
	int				i, j, k, powLen;
	double			sum, numOfBlocks;
	unsigned int	*P;
	
	if ( (m == 0) || (m == -1) )
		return 0.0;
	numOfBlocks = n;
	powLen = (int)pow(2, m+1)-1;
	if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
		return -1;
	}
	for ( i=1; i<powLen-1; i++ )
		P[i] = 0;	  /* INITIALIZE NODES */
	for ( i=0; i<numOfBlocks; i++ ) {		 /* COMPUTE FREQUENCY */
		k = 1;
		for ( j=0; j<m; j++ ) {
			if ( epsilon[(i+j)%n] == 0 )
				k *= 2;
			else if ( epsilon[(i+j)%n] == 1 )
				k = 2*k+1;
		}
		P[k-1]++;
	}
	sum = 0.0;
	for ( i=(int)pow(2, m)-1; i<(int)pow(2, m+1)-1; i++ )
		sum += pow(P[i], 2);
	sum = (sum * pow(2, m)/(double)n) - (double)n;
	free(P);
	
	return sum;
}

double Frequency(int n) {
	int		i;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
	
	sum = 0.0;
	for ( i=0; i<n; i++ ) {
		sum += 2*(int)epsilon[i]-1;
    }
	s_obs = fabs(sum)/sqrt(n);
	f = s_obs/sqrt2;
	p_value = erfc(f);

	return  p_value;
}

double ApproximateEntropy(int m, int n) {
	int				i, j, k, r, blockSize, seqLength, powLen, index;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P;
	
	seqLength = n;
	r = 0;
	
	for ( blockSize=m; blockSize<=m+1; blockSize++ ) {
		if ( blockSize == 0 ) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			powLen = (int)pow(2, blockSize+1)-1;
			if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
				return -1;
			}
			for ( i=1; i<powLen-1; i++ )
				P[i] = 0;
			for ( i=0; i<numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
				k = 1;
				for ( j=0; j<blockSize; j++ ) {
					k <<= 1;
					if ( (int)epsilon[(i+j) % seqLength] == 1 )
						k++;
				}
				P[k-1]++;
			}
			/* DISPLAY FREQUENCY */
			sum = 0.0;
			index = (int)pow(2, blockSize)-1;
			for ( i=0; i<(int)pow(2, blockSize); i++ ) {
				if ( P[index] > 0 )
					sum += P[index]*log(P[index]/numOfBlocks);
				index++;
			}
			sum /= numOfBlocks;
			ApEn[r] = sum;
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];
	
	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-1), chi_squared/2.0);
	
    return p_value;
}

double BlockFrequency(int M, int n) {
	int		i, j, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;
	
	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	for ( i=0; i<N; i++ ) {
		blockSum = 0;
		for ( j=0; j<M; j++ )
			blockSum += epsilon[j+i*M];
		pi = (double)blockSum/(double)M;
		v = pi - 0.5;
		sum += v*v;
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);

    return p_value;
}

double DiscreteFourierTransform(int n) {
	double	p_value, upperBound, N_l, N_o, d, *m = NULL, *X = NULL, *wsave = NULL;
	int		i, count, ifac[15];

	if ( ((X = (double*) calloc(n,sizeof(double))) == NULL) ||
		 ((wsave = (double *)calloc(2*n,sizeof(double))) == NULL) ||
		 ((m = (double*)calloc(n/2+1, sizeof(double))) == NULL) ) {
			if( X != NULL )
				free(X);
			if( wsave != NULL )
				free(wsave);
			if( m != NULL )
				free(m);
			return -1;
	}
	for ( i=0; i<n; i++ )
		X[i] = 2*(int)epsilon[i] - 1;
	
	__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(n, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */
	
	for ( i=0; i<n/2; i++ )
		m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2)); 
	count = 0;				       /* CONFIDENCE INTERVAL */
	upperBound = sqrt(2.995732274*n);
	for ( i=0; i<n/2; i++ )
		if ( m[i] < upperBound)
			count++;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));

    free(X);
    free(wsave);
    free(m);

    return p_value;
}

double LongestRunOfOnes(int n) {
	double			pval, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M, V[7];
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

	if ( n < 128 ) {
		// n is too small
		return -1;
	}
	if ( n < 6272 ) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if ( n < 750000 ) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
			V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}
	
	N = n/M;
	for ( i=0; i<N; i++ ) {
		v_n_obs = 0;
		run = 0;
		for ( j=0; j<M; j++ ) {
			if ( epsilon[i*M+j] == 1 ) {
				run++;
				if ( run > v_n_obs )
					v_n_obs = run;
			}
			else
				run = 0;
		}
		if ( v_n_obs < V[0] )
			nu[0]++;
		for ( j=0; j<=K; j++ ) {
			if ( v_n_obs == V[j] )
				nu[j]++;
		}
		if ( v_n_obs > V[K] )
			nu[K]++;
	}

	chi2 = 0.0;
	for ( i=0; i<=K; i++ )
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	pval = cephes_igamc((double)(K/2.0), chi2 / 2.0);

	return pval;
}

double RandomExcursions(int n) {
	int		b, i, j, k, J, x;
	int		cycleStart, cycleStop, *cycle = NULL, *S_k = NULL;
	int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
	int		counter[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	double	p_value, sum, constraint, nu[6][8];
	double	pi[5][6] = { {0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000}, 
						 {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
						 {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
						 {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
						 {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051} };
	
	if ( ((S_k = (int *)calloc(n, sizeof(int))) == NULL) ||
		 ((cycle = (int *)calloc(MAX(1000, n/100), sizeof(int))) == NULL) ) {
		if ( S_k != NULL )
			free(S_k);
		if ( cycle != NULL )
			free(cycle);
		return -1;
	}
	
	J = 0; 					/* DETERMINE CYCLES */
	S_k[0] = 2*(int)epsilon[0] - 1;
	for( i=1; i<n; i++ ) {
		S_k[i] = S_k[i-1] + 2*epsilon[i] - 1;
		if ( S_k[i] == 0 ) {
			J++;
			if ( J > MAX(1000, n/100) ) {
				free(S_k);
				free(cycle);
				return -2;
			}
			cycle[J] = i;
		}
	}
	if ( S_k[n-1] != 0 )
		J++;
	cycle[J] = n;

	constraint = MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
        free(S_k);
        free(cycle);
        return -3;
	}
	else {
		cycleStart = 0;
		cycleStop  = cycle[1];
		for ( k=0; k<6; k++ )
			for ( i=0; i<8; i++ )
				nu[k][i] = 0.;
		for ( j=1; j<=J; j++ ) {                           /* FOR EACH CYCLE */
			for ( i=0; i<8; i++ )
				counter[i] = 0;
			for ( i=cycleStart; i<cycleStop; i++ ) {
				if ( (S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1) ) {
					if ( S_k[i] < 0 )
						b = 4;
					else
						b = 3;
					counter[S_k[i]+b]++;
				}
			}
			cycleStart = cycle[j]+1;
			if ( j < J )
				cycleStop = cycle[j+1];
			
			for ( i=0; i<8; i++ ) {
				if ( (counter[i] >= 0) && (counter[i] <= 4) )
					nu[counter[i]][i]++;
				else if ( counter[i] >= 5 )
					nu[5][i]++;
			}
		}

		p_value = 1;

		for ( i=0; i<8; i++ ) {
			x = stateX[i];
			sum = 0.;
			for ( k=0; k<6; k++ )
				sum += pow(nu[k][i] - J*pi[(int)fabs(x)][k], 2) / (J*pi[(int)fabs(x)][k]);
			p_value = fmin(p_value, cephes_igamc(2.5, sum/2.0));
		}
        free(S_k);
        free(cycle);
        return p_value;
	}
}

double RandomExcursionsVariant(int n) {
	int		i, p, J, x, constraint, count, *S_k;
	int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double	p_value;
	
	if ( (S_k = (int *)calloc(n, sizeof(int))) == NULL ) {
		return -1;
	}
	J = 0;
	S_k[0] = 2*(int)epsilon[0] - 1;
	for ( i=1; i<n; i++ ) {
		S_k[i] = S_k[i-1] + 2*epsilon[i] - 1;
		if ( S_k[i] == 0 )
			J++;
	}
	if ( S_k[n-1] != 0 )
		J++;

	constraint = (int)MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
        free(S_k);
	    return -3;
	}
	else {
	    p_value = 1;
		for ( p=0; p<=17; p++ ) {
			x = stateX[p];
			count = 0;
			for ( i=0; i<n; i++ )
				if ( S_k[i] == x )
					count++;
			p_value = fmin(p_value, erfc(fabs(count-J)/(sqrt(2.0*J*(4.0*fabs(x)-2)))));
        }
        free(S_k);
        return p_value;
	}
}

double CumulativeSums(int n)
{
	int		S, sup, inf, z, zrev, k;
	double	sum1, sum2, p_value1, p_value2;

	S = 0;
	sup = 0;
	inf = 0;
	for ( k=0; k<n; k++ ) {
		epsilon[k] ? S++ : S--;
		if ( S > sup )
			sup++;
		if ( S < inf )
			inf--;
		z = (sup > -inf) ? sup : -inf;
		zrev = (sup-S > S-inf) ? sup-S : S-inf;
	}
	
	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value1 = 1.0 - sum1 + sum2;

	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value2 = 1.0 - sum1 + sum2;

    return fmin(p_value1, p_value2);
}

double Rank(int n) {
	int			N, i, k, r;
	double		p_value, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
	BitSequence	**matrix = create_matrix(32, 32);
	
	N = n/(32*32);
	if ( isZero(N) ) {
		p_value = 0.00;
	}
	else {
		r = 32;					/* COMPUTE PROBABILITIES */
		product = 1;
		for ( i=0; i<=r-1; i++ )
			product *= ((1.e0-pow(2, i-32))*(1.e0-pow(2, i-32)))/(1.e0-pow(2, i-r));
		p_32 = pow(2, r*(32+32-r)-32*32) * product;
		
		r = 31;
		product = 1;
		for ( i=0; i<=r-1; i++ )
			product *= ((1.e0-pow(2, i-32))*(1.e0-pow(2, i-32)))/(1.e0-pow(2, i-r));
		p_31 = pow(2, r*(32+32-r)-32*32) * product;
		
		p_30 = 1 - (p_32+p_31);
		
		F_32 = 0;
		F_31 = 0;
		for ( k=0; k<N; k++ ) {			/* FOR EACH 32x32 MATRIX   */
			def_matrix(32, 32, matrix, k, epsilon);
			R = computeRank(32, 32, matrix);
			if ( R == 32 )
				F_32++;			/* DETERMINE FREQUENCIES */
			if ( R == 31 )
				F_31++;
		}
		F_30 = (double)N - (F_32+F_31);

		chi_squared =(pow(F_32 - N*p_32, 2)/(double)(N*p_32) +
					  pow(F_31 - N*p_31, 2)/(double)(N*p_31) +
					  pow(F_30 - N*p_30, 2)/(double)(N*p_30));
		
		arg1 = -chi_squared/2.e0;

		p_value = exp(arg1);

		for ( i=0; i<32; i++ )				/* DEALLOCATE MATRIX  */
			free(matrix[i]);
		free(matrix);
	}

    return p_value;
}

double Universal(int n) {
	int		i, j, p, L, Q, K;
	double	arg, sqrt2, sigma, phi, sum, p_value, c;
	long	*T, decRep;
	double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
				8.1764248, 9.1723243, 10.170032, 11.168765,
				12.168070, 13.167693, 14.167488, 15.167379 };
	double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
				3.401, 3.410, 3.416, 3.419, 3.421 };
	
	/* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
	 * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
	 * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	L = 5;
	if ( n >= 387840 )     L = 6;
	if ( n >= 904960 )     L = 7;
	if ( n >= 2068480 )    L = 8;
	if ( n >= 4654080 )    L = 9;
	if ( n >= 10342400 )   L = 10;
	if ( n >= 22753280 )   L = 11;
	if ( n >= 49643520 )   L = 12;
	if ( n >= 107560960 )  L = 13;
	if ( n >= 231669760 )  L = 14;
	if ( n >= 496435200 )  L = 15;
	if ( n >= 1059061760 ) L = 16;
	
	Q = 10*(int)pow(2, L);
	K = (int) (floor(n/L) - (double)Q);	 		    /* BLOCKS TO TEST */
	
	p = (int)pow(2, L);
	if ( (L < 6) || (L > 16) )
	    return n;
	if ( ((double)Q < 10*pow(2, L)) )
	    return -2;
	if ( ((T = (long *)calloc(p, sizeof(long))) == NULL) )
	    return -3;

	/* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
	c = 0.7 - 0.8/(double)L + (4 + 32/(double)L)*pow(K, -3/(double)L)/15;
	sigma = c * sqrt(variance[L]/(double)K);
	sqrt2 = sqrt(2);
	sum = 0.0;
	for ( i=0; i<p; i++ )
		T[i] = 0;
	for ( i=1; i<=Q; i++ ) {		/* INITIALIZE TABLE */
		decRep = 0;
		for ( j=0; j<L; j++ )
			decRep += epsilon[(i-1)*L+j] * (long)pow(2, L-1-j);
		T[decRep] = i;
	}
	for ( i=Q+1; i<=Q+K; i++ ) { 	/* PROCESS BLOCKS */
		decRep = 0;
		for ( j=0; j<L; j++ )
			decRep += epsilon[(i-1)*L+j] * (long)pow(2, L-1-j);
		sum += log(i - T[decRep])/log(2);
		T[decRep] = i;
	}
	phi = (double)(sum/(double)K);

	arg = fabs(phi-expected_value[L])/(sqrt2 * sigma);
	p_value = erfc(arg);

	free(T);

	return p_value;
}

double Runs(int n) {
	int		S, k;
	double	pi, V, erfc_arg, p_value;

	S = 0;
	for ( k=0; k<n; k++ )
		if ( epsilon[k] )
			S++;
	pi = (double)S / (double)n;

	if ( fabs(pi - 0.5) > (2.0 / sqrt(n)) ) {
        // PI ESTIMATOR CRITERIA NOT MET!
		return 0.0;
	}
	else {

		V = 1;
		for ( k=1; k<n; k++ )
			if ( epsilon[k] != epsilon[k-1] )
				V++;
	
		erfc_arg = fabs(V - 2.0 * n * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2*n));
		p_value = erfc(erfc_arg);

		return p_value;
	}
}

double Serial(int m, int n) {
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	
	psim0 = psi2(m, n, epsilon);
	psim1 = psi2(m-1, n, epsilon);
	psim2 = psi2(m-2, n, epsilon);

	if (psim0 == -1 || psim1 == -1 || psim1 == -1)
	    return -1;
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m-1)/2, del1/2.0);
	p_value2 = cephes_igamc(pow(2, m-2)/2, del2/2.0);

	return fmin(p_value1, p_value2);
}

double
LinearComplexity(int M, int n)
{
	int       i, ii, j, d, N, L, m, N_, parity, sign, K = 6;
	double    p_value, T_, mean, nu[7], chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;
	
	N = (int)floor(n/M);
	if ( ((B_ = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((C  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((P  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ||
		 ((T  = (BitSequence *) calloc(M, sizeof(BitSequence))) == NULL) ) {
		if ( B_ != NULL )
			free(B_);
		if ( C != NULL )
			free(C);
		if ( P != NULL )
			free(P);
		if ( T != NULL )
			free(T);
		return -1;
	}


	for ( i=0; i<K+1; i++ )
		nu[i] = 0.00;
	for ( ii=0; ii<N; ii++ ) {
		for ( i=0; i<M; i++ ) {
			B_[i] = 0;
			C[i] = 0;
			T[i] = 0;
			P[i] = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0] = 1;
		B_[0] = 1;
		
		/* DETERMINE LINEAR COMPLEXITY */
		N_ = 0;
		while ( N_ < M ) {
			d = (int)epsilon[ii*M+N_];
			for ( i=1; i<=L; i++ )
				d += C[i] * epsilon[ii*M+N_-i];
			d = d%2;
			if ( d == 1 ) {
				for ( i=0; i<M; i++ ) {
					T[i] = C[i];
					P[i] = 0;
				}
				for ( j=0; j<M; j++ )
					if ( B_[j] == 1 )
						P[j+N_-m] = 1;
				for ( i=0; i<M; i++ )
					C[i] = (C[i] + P[i])%2;
				if ( L <= N_/2 ) {
					L = N_ + 1 - L;
					m = N_;
					for ( i=0; i<M; i++ )
						B_[i] = T[i];
				}
			}
			N_++;
		}
		if ( (parity = (M+1)%2) == 0 ) 
			sign = -1;
		else 
			sign = 1;
		mean = M/2.0 + (9.0+sign)/36.0 - 1.0/pow(2, M) * (M/3.0 + 2.0/9.0);
		if ( (parity = M%2) == 0 )
			sign = 1;
		else 
			sign = -1;
		T_ = sign * (L - mean) + 2.0/9.0;
		
		if ( T_ <= -2.5 )
			nu[0]++;
		else if ( T_ > -2.5 && T_ <= -1.5 )
			nu[1]++;
		else if ( T_ > -1.5 && T_ <= -0.5 )
			nu[2]++;
		else if ( T_ > -0.5 && T_ <= 0.5 )
			nu[3]++;
		else if ( T_ > 0.5 && T_ <= 1.5 )
			nu[4]++;
		else if ( T_ > 1.5 && T_ <= 2.5 )
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
	for ( i=0; i<K+1; i++ )
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	free(B_);
	free(P);
	free(C);
	free(T);

	return p_value;
}

double NonOverlappingTemplateMatchings(int m, int n) {
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, nu[6], *Wj = NULL; 
	FILE			*fp = NULL;
	double			sum, chi2, p_value, min_p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	N = 8;
	M = n/N;

	min_p_value = 1;

	if ( (Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == NULL ) {
		return -1;
	}
	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ) {
	    return -1;
	}
	if ( ((fp = fopen(directory, "r")) == NULL) ) {
	    return -2;
	}
	if ( ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL) ) {
		if ( sequence != NULL )
			free(sequence);
		return -3;
	}
	else {
		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

		for( jj=0; jj<MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++ ) {

			for ( k=0; k<m; k++ ) {
				fscanf(fp, "%d", &bit);
				sequence[k] = bit;
			}
			for ( k=0; k<=K; k++ )
				nu[k] = 0;
			for ( i=0; i<N; i++ ) {
				W_obs = 0;
				for ( j=0; j<M-m+1; j++ ) {
					match = 1;
					for ( k=0; k<m; k++ ) {
						if ( (int)sequence[k] != (int)epsilon[i*M+j+k] ) {
							match = 0;
							break;
						}
					}
					if ( match == 1 ) {
						W_obs++;
                        j += m-1;
                    }
				}
				Wj[i] = W_obs;
			}
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
			}

            p_value = cephes_igamc(N/2.0, chi2/2.0);
			if (p_value < min_p_value)
			    min_p_value = p_value;

			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
		}
	}
	
	if ( sequence != NULL )
		free(sequence);

	free(Wj);
    if ( fp != NULL )
        fclose(fp);

    return min_p_value;
}

double OverlappingTemplateMatchings(int m, int n) {
	int				i, k, match;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, j, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };
	BitSequence		*sequence;

	M = 1032;
	N = n/M;
	
	if ( (sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL ) {
	    return -1;
	}
	else
		for ( i=0; i<m; i++ )
			sequence[i] = 1;
	
	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for ( i=0; i<N; i++ ) {
		W_obs = 0;
		for ( j=0; j<M-m+1; j++ ) {
			match = 1;
			for ( k=0; k<m; k++ ) {
				if ( sequence[k] != epsilon[i*M+j+k] )
					match = 0;
			}
			if ( match == 1 )
				W_obs++;
		}
		if ( W_obs <= 4 )
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}
	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);

	free(sequence);

	return p_value;
}

void generateRandomFile(const char *filename, size_t numBytes) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        printf("Error opening file for writing.\n");
        return;
    }

    srand(time(NULL));

    for (size_t i = 0; i < numBytes; ++i) {
        unsigned char byte = rand() % 256;  
        fwrite(&byte, sizeof(byte), 1, file);
    }

    fclose(file);
    printf("Random file generated with %zu bytes.\n", numBytes);
}

int convertToBit(BYTE *x, int xBitLength, int bitsNeeded, int *bitsRead) {
	int		i, j, count, bit;
	BYTE	mask;

	count = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
			}
			else {
				bit = 0;
			}
			mask >>= 1;
			epsilon[*bitsRead] = bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}
	
	return 0;
}
void readsHexDigitsInBinaryFormat(FILE *fp, unsigned int N) {  // N - количество бит
    int bitsRead = 0;
    BYTE buffer[4];
    
    if ((epsilon = (BitSequence *) calloc(N, sizeof(BitSequence))) == NULL) {
        printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
        return;
    }

    // Чтение данных и конвертация в биты
    while (bitsRead < N) {
        if (fread(buffer, sizeof(unsigned char), 4, fp) != 4) {
            printf("READ ERROR: Insufficient data in file.\n");
            free(epsilon);
            return;
        }

        // Конвертация 4 байтов (32 бита) в биты
        if (convertToBit(buffer, 32, N, &bitsRead) == 1) {
            break;  // Прочитаны все нужные биты
        }
    }
}

typedef struct {
    double p_value;
    int index;
} TestResult;

void start(int tests[15], int n) {
	double values[15];
	printf("15 tests \n");
	int n1 = n / 16;
	values[0] = Runs(n1);
	values[1] = ApproximateEntropy(10, n1);
	values[2] = Frequency(n1);
	values[3] = BlockFrequency(128, n1);
	values[4] = CumulativeSums(n1);
	values[5] = LongestRunOfOnes(n1);
	values[6] = 1;//Rank(n1);
	values[7] = DiscreteFourierTransform(n1);
	values[8] = NonOverlappingTemplateMatchings(9, n1);
	values[9] = OverlappingTemplateMatchings(9, n1);
	values[10] = Universal(n1);
	values[11] = RandomExcursions(n1);
	values[12] = RandomExcursionsVariant(n1);
	values[13] = LinearComplexity(500, n1);
	values[14] = Serial(16, n1);

	for (int i = 0; i < 15; i++) {
		if (i == 0) {
			printf("Runs                             ");
		} else if (i == 1) {
			printf("ApproximateEntropy               ");
		 } else if (i == 2) {
			printf("Frequency                        ");
		} else if (i == 3) {
			printf("BlockFrequency                   ");
		} else if (i == 4) {
			printf("CumulativeSums                   ");
		} else if (i == 5) {
			printf("LongestRunOfOnes                 ");
		} else if (i == 6) {
			printf("Rank                             ");
		} else if (i == 7) {
			printf("DiscreteFourierTransform         ");
		} else if (i == 8) {	
			printf("NonOverlappingTemplateMatchings  ");
		} else if (i == 9) {
			printf("OverlappingTemplateMatchings     ");
		} else if (i == 10) {
			printf("Universal                        ");
		} else if (i == 11) {
			printf("RandomExcursions                 ");
		} else if (i == 12) {
			printf("RandomExcursionsVariant          ");
		} else if (i == 13) {
			printf("LinearComplexity                 ");
		} else if (i == 14) {
			printf("Serial                           ");
		}
		printf("p-value      :%f\n", values[i]);
	}
	
	TestResult results[15];

	for (int i = 0; i < 15; i++) {
        results[i].p_value = values[i];
        results[i].index = i;
    }

	// Сортировка по p-value (по возрастанию)
    for (int i = 0; i < 14; i++) {
        for (int j = i + 1; j < 15; j++) {
            if (results[i].p_value > results[j].p_value) {
                // Обмен p-value
                double temp_value = results[i].p_value;
                results[i].p_value = results[j].p_value;
                results[j].p_value = temp_value;

                // Обмен индексов
                int temp_index = results[i].index;
                results[i].index = results[j].index;
                results[j].index = temp_index;
            }
        }
    }

	 // Для остальных тестов ставим tests[i] = -1
    for (int i = 5; i < 15; i++) {
        tests[results[i].index] = -1; // Устанавливаем -1 для тестов с большими p-value
    }
	
	int n2 = n / 8;
	for (int i = 0; i < 15; i++) {
		if (tests[i] == 1) {
			if (i == 0) {
				values[0] = Runs(n2);
			} else if (i == 1) {
				values[1] = ApproximateEntropy(10, n2);
			} else if (i == 2) {
				values[2] = Frequency(n2);
			} else if (i == 3) {
				values[3] = BlockFrequency(128, n2);
			} else if (i == 4) {
				values[4] = CumulativeSums(n2);
			} else if (i == 5) {
				values[5] = LongestRunOfOnes(n2);
			} else if (i == 6) {
				values[6] = 1;//Rank(n2);
			} else if (i == 7) {
				values[7] = DiscreteFourierTransform(n2);
			} else if (i == 8) {
				values[8] = NonOverlappingTemplateMatchings(9, n2);
			} else if (i == 9) {
				values[9] = OverlappingTemplateMatchings(9, n2);
			} else if (i == 10) {
				values[10] = Universal(n2);
			} else if (i == 11) {
				values[11] = RandomExcursions(n2);
			} else if (i == 12) {
				values[12] = RandomExcursionsVariant(n2);
			} else if (i == 13) {
				values[13] = LinearComplexity(500, n2);
			} else if (i == 14) {
				values[14] = Serial(16, n2);
			}
		}
	}

	printf("5 tests \n");
	for (int i = 0; i < 15; i++) {
		if (tests[i] == 1) {
			if (i == 0) {
			printf("Runs                             ");
		} else if (i == 1) {
			printf("ApproximateEntropy               ");
		 } else if (i == 2) {
			printf("Frequency                        ");
		} else if (i == 3) {
			printf("BlockFrequency                   ");
		} else if (i == 4) {
			printf("CumulativeSums                   ");
		} else if (i == 5) {
			printf("LongestRunOfOnes                 ");
		} else if (i == 6) {
			printf("Rank                             ");
		} else if (i == 7) {
			printf("DiscreteFourierTransform         ");
		} else if (i == 8) {	
			printf("NonOverlappingTemplateMatchings  ");
		} else if (i == 9) {
			printf("OverlappingTemplateMatchings     ");
		} else if (i == 10) {
			printf("Universal                        ");
		} else if (i == 11) {
			printf("RandomExcursions                 ");
		} else if (i == 12) {
			printf("RandomExcursionsVariant          ");
		} else if (i == 13) {
			printf("LinearComplexity                 ");
		} else if (i == 14) {
			printf("Serial                           ");
		}
		printf("p-value     :%f\n", values[i]);
		}
	}

	TestResult results8[5];

	int counter = 0;
	for (int i = 0; i < 15; i++) {
		if (tests[i] == 1) {
			results8[counter].p_value = values[i];
			results8[counter].index = i;
			counter++;
		}
    }

	// Сортировка по p-value (по возрастанию)
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 5; j++) {
            if (results8[i].p_value > results8[j].p_value) {
                // Обмен p-value
                double temp_value = results8[i].p_value;
                results8[i].p_value = results8[j].p_value;
                results8[j].p_value = temp_value;

                // Обмен индексов
                int temp_index = results8[i].index;
                results8[i].index = results8[j].index;
                results8[j].index = temp_index;
            }
        }
    }

	 // Для остальных тестов ставим tests[i] = -1
    for (int i = 1; i < 5; i++) {
        tests[results8[i].index] = -1; // Устанавливаем -1 для тестов с большими p-value
    }

	int n3 = n;
	for (int i = 0; i < 15; i++) {
		if (tests[i] == 1) {
			if (i == 0) {
				values[0] = Runs(n3);
			} else if (i == 1) {
				values[1] = ApproximateEntropy(10, n3);
			} else if (i == 2) {
				values[2] = Frequency(n3);
			} else if (i == 3) {
				values[3] = BlockFrequency(128, n3);
			} else if (i == 4) {
				values[4] = CumulativeSums(n3);
			} else if (i == 5) {
				values[5] = LongestRunOfOnes(n3);
			} else if (i == 6) {
				values[6] = 1;//Rank(n3);
			} else if (i == 7) {
				values[7] = DiscreteFourierTransform(n3);
			} else if (i == 8) {
				values[8] = NonOverlappingTemplateMatchings(9, n3);
			} else if (i == 9) {
				values[9] = OverlappingTemplateMatchings(9, n3);
			} else if (i == 10) {
				values[10] = Universal(n3);
			} else if (i == 11) {
				values[11] = RandomExcursions(n3);
			} else if (i == 12) {
				values[12] = RandomExcursionsVariant(n3);
			} else if (i == 13) {
				values[13] = LinearComplexity(500, n3);
			} else if (i == 14) {
				values[14] = Serial(16, n3);
			}
		}
	}

	printf("1 tests \n");
	for (int i = 0; i < 15; i++) {
		if (tests[i] == 1) {
			if (i == 0) {
			printf("Runs                             ");
		} else if (i == 1) {
			printf("ApproximateEntropy               ");
		 } else if (i == 2) {
			printf("Frequency                        ");
		} else if (i == 3) {
			printf("BlockFrequency                   ");
		} else if (i == 4) {
			printf("CumulativeSums                   ");
		} else if (i == 5) {
			printf("LongestRunOfOnes                 ");
		} else if (i == 6) {
			printf("Rank                             ");
		} else if (i == 7) {
			printf("DiscreteFourierTransform         ");
		} else if (i == 8) {	
			printf("NonOverlappingTemplateMatchings  ");
		} else if (i == 9) {
			printf("OverlappingTemplateMatchings     ");
		} else if (i == 10) {
			printf("Universal                        ");
		} else if (i == 11) {
			printf("RandomExcursions                 ");
		} else if (i == 12) {
			printf("RandomExcursionsVariant          ");
		} else if (i == 13) {
			printf("LinearComplexity                 ");
		} else if (i == 14) {
			printf("Serial                           ");
		}
		printf("p-value      :%f\n", values[i]);
		}
	}
}


int main() {
	int n = 1000000  * 8;
    //generateRandomFile("data2.bin", (size_t) (n / 8));

	FILE *fp = fopen("data1.bin", "rb");

	readsHexDigitsInBinaryFormat(fp, n);

	fclose(fp);

	int tests[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

	start(tests, n);

	free(epsilon);

    return 0;
}
