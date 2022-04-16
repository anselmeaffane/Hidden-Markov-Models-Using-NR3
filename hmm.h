//Ajouter le 09/08/2020
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdio.h>      /* printf */
#include "nr3.h"

struct HMM {
	MatDoub a, b;
	VecInt obs;
	Int fbdone;
	Int mstat, nobs, ksym;
	Int lrnrm;
	MatDoub alpha, beta, pstate;
	VecInt arnrm, brnrm;
	Doub BIG, BIGI, lhood;
	HMM(MatDoub_I &aa, MatDoub_I &bb, VecInt_I &obs);
	void forwardbackward();
	void baumwelch();
	Doub loglikelihood() {return log(lhood)+lrnrm*log(BIGI);}
};

HMM::HMM(MatDoub_I &aa, MatDoub_I &bb, VecInt_I &obss) :
	a(aa), b(bb), obs(obss), fbdone(0),
	mstat(a.nrows()), nobs(obs.size()), ksym(b.ncols()),
	alpha(nobs,mstat), beta(nobs,mstat), pstate(nobs,mstat),
	arnrm(nobs), brnrm(nobs), BIG(1.e20), BIGI(1./BIG)  {
	Int i,j,k;
	Doub sum;

	/*Verify initial parameters*/
	//Initial data control display output into screen and file
	ofstream outstream;
	outstream.open("hmm_init_values.log");
	outstream << "a.nrows() = " << a.nrows() << endl;
	outstream << "mstat(a.nrows()) = mstat(" << a.nrows() << ")" << endl;
	outstream << "obs.size() = " << obs.size() << endl;
	outstream << "nobs(obs.size()) = nobs(" << obs.size() << ")" << endl;
	outstream << "b.ncols() = " << b.ncols() << endl;
	outstream << "ksym(b.ncols()) = ksym(" << b.ncols() << ")" << endl;
	outstream << "alpha(nobs,mstat) = alpha(" << nobs << "," << mstat << ")" << endl;
	outstream << "beta(nobs,mstat) = beta(" << nobs << "," << mstat << ")" << endl;
	outstream << "pstate(nobs,mstat) = pstate(" << nobs << "," << mstat << ")" << endl;
	outstream << "arnrm(nobs) = arnrm(" << nobs << ")" << endl;
	outstream << "brnrm(nobs) = brnrm(" << nobs << ")" << endl;
	outstream << "BIG = " << BIG << endl;
	outstream << "BIGI = " << BIGI << endl;

	cout << endl;
	outstream << endl;

	if (a.ncols() != mstat) throw("transition matrix not square");
	if (b.nrows() != mstat) throw("symbol prob matrix wrong size");
	for (i=0; i<nobs; i++) {
		if (obs[i] < 0 || obs[i] >= ksym) throw("bad data in obs");
		//Display Obs
		cout << "Observe = " << obs[i] << "\t";
		outstream << "Observe = " << obs[i] << "\t";
	}
	cout << endl;
	outstream << endl;

	for (i=0; i<mstat; i++) {
		sum = 0.;
		for (j=0; j<mstat; j++) sum += a[i][j];
		if (abs(sum - 1.) > 0.01) throw("transition matrix not normalized");
		for (j=0; j<mstat; j++) a[i][j] /= sum;
	}

	//Display A
	cout << endl;
	outstream << endl;

	for (i=0; i<mstat; i++) {
		for (j=0; j<mstat; j++) {
			cout << "a[" << i << "][" << j << "] = \t" << a[i][j] << "\t";
			outstream << "a[" << i << "][" << j << "] = \t" << a[i][j] << "\t";
		}
		cout << endl;
		outstream << endl;
	}
	cout << endl; outstream << endl;

	for (i=0; i<mstat; i++) {
		sum = 0.;
		for (k=0; k<ksym; k++) sum += b[i][k];
		if (abs(sum - 1.) > 0.01) throw("symbol prob matrix not normalized");
		for (k=0; k<ksym; k++) b[i][k] /= sum;
	}
	//Display B
	for (i=0; i<mstat; i++) {
		for (k=0; k<ksym; k++) {
			cout << "b[" << i << "][" << k << "] = \t" << b[i][k] << "\t";
			outstream << "b[" << i << "][" << k << "] = \t" << b[i][k] << "\t";
		}
		cout << endl;
		outstream << endl;
	}

	outstream.close();
	//outstream2.close();
}
void HMM::forwardbackward() {
	Int i,j,t;
	Doub sum,asum,bsum;
	//Ajouté le 13/08/2020 par Anselme Affane
	VecDoub pi(mstat,0.);
	pi[0]=0.51316; pi[1]=0.51316;
	ofstream outstream2;
	outstream2.open("forward-backward.log");
	cout << "===== Forward Algorithm (Matrice Alpha) =====" << endl;
	outstream2 << "===== Forward Algorithm (Matrice Alpha) =====" << endl;

	for (i=0; i<mstat; i++) {
		alpha[0][i] = b[i][obs[0]]*pi[i];//Ajout
		//alpha[0][i] = b[i][obs[0]];
		cout << "alpha [0][" << i << "] = \t" << (alpha[0][i]) << "\t";
		outstream2 << "alpha [0][" << i << "] = \t" << (alpha[0][i]) << "\t";
	}
	cout << endl;
	outstream2 << endl;

	arnrm[0] = 0;
	for (t=1; t<nobs; t++) {
		asum = 0;
		for (j=0; j<mstat; j++) {
			sum = 0.;
			for (i=0; i<mstat; i++) sum += alpha[t-1][i]*a[i][j]*b[j][obs[t]];
			alpha[t][j] = sum;
			asum += sum;
			cout << "alpha [" << t << "][" << j << "] = \t" << (alpha[t][j]) << "\t";
			outstream2 << "alpha [" << t << "][" << j << "] = \t" << (alpha[t][j]) << "\t";
		}
		cout << endl;
		outstream2 << endl;

		arnrm[t] = arnrm[t-1];
		if (asum < BIGI) {
			++arnrm[t];
			for (j=0; j<mstat; j++) alpha[t][j] *= BIG;
		}
	}

	cout << "P(O | Lambda) = \t" << asum << endl;
	outstream2 << "P(Observe Seq | Lambda) = \t" << asum << endl;

	cout << endl; outstream2 << endl;

	cout << "===== Backward Algorithm (Matrice Beta) =====" << endl;
	outstream2 << "===== Backward Algorithm (Matrice Beta) =====" << endl;

	for (i=0; i<mstat; i++) {
		beta[nobs-1][i] = 1.;
		cout << "beta [" << nobs-1 << "][" << i << "] = \t" << (beta[nobs-1][i]) << "\t";
		outstream2 << "beta [" << nobs-1 << "][" << i << "] = \t" << (beta[nobs-1][i]) << "\t";
		//cout << "beta [0][" << i << "] = \t" << (beta[nobs-1][i]) << "\t";
		//outstream2 << "beta [0][" << i << "] = \t" << (beta[nobs-1][i]) << "\t";
	}
	cout << endl;
	outstream2 << endl;

	brnrm[nobs-1] = 0;
	for (t=nobs-2; t>=0; t--) {
		bsum = 0.;
		for (i=0; i<mstat; i++) {
			sum = 0.;
			for (j=0; j<mstat; j++) sum += a[i][j]*b[j][obs[t+1]]*beta[t+1][j];
			beta[t][i] = sum;
			bsum += sum;

			cout << "beta [" << t << "][" << i << "] = \t" << (beta[t][i]) << "\t";
			outstream2 << "beta [" << t << "][" << i << "] = \t" << (beta[t][i]) << "\t";
		}
		brnrm[t] = brnrm[t+1];
		if (bsum < BIGI) {
			++brnrm[t];
			for (j=0; j<mstat; j++) {
				beta[t][j] *= BIG;

				cout << "beta [" << t << "][" << i << "] = \t" << (beta[t][i]) << "\t";
				outstream2 << "beta [" << t << "][" << i << "] = \t" << (beta[t][i]) << "\t";
			}
		}
		cout << endl;
		outstream2 << endl;
	}
	cout << endl; outstream2 << endl;

	cout << "===== Forward-Backwar Algorithm (Matrice Gamma) =====" << endl;
	outstream2 << "===== Forward-Backwar Algorithm (Matrice Gamma) =====" << endl;
	lhood = 0.;
	for (i=0; i<mstat; i++) lhood += alpha[0][i]*beta[0][i];
	lrnrm = arnrm[0] + brnrm[0];
	while (lhood < BIGI) {lhood *= BIG; lrnrm++;}
	for (t=0; t<nobs; t++) {
		sum = 0.;
		for (i=0; i<mstat; i++) sum += (pstate[t][i] = alpha[t][i]*beta[t][i]);
		// sum = lhood*pow(BIGI, lrnrm - arnrm[t] - brnrm[t]);
		for (i=0; i<mstat; i++) {
			pstate[t][i] /= sum;
			cout << "pstate [" << t << "][" << i << "] = \t" << (pstate[t][i]) << "\t";
			outstream2 << "pstate [" << t << "][" << i << "] = \t" << (pstate[t][i]) << "\t";
		}
		outstream2 << endl;
		cout << endl;
		cout << "P(O | Lambda) = \t" << sum << endl;
		outstream2 << "P(Observe Seq | Lambda) = \t" << sum << endl;
		outstream2 << endl;
		cout << endl;

	}
	fbdone = 1;
	outstream2.close();
}
void HMM::baumwelch() {
	Int i,j,k,t;
	Doub num,denom,term;
	MatDoub bnew(mstat,ksym);
	Doub powtab[10];

	ofstream outstream3;
	outstream3.open("baumwelch.log");
	cout << "===== Parameters Reestimation =====" << endl;
	outstream3 << "===== Parameters Reestimation =====" << endl;

	for (i=0; i<10; i++) powtab[i] = pow(BIGI,i-6);
	if (fbdone != 1) throw("must do forwardbackward first");
	for (i=0; i<mstat; i++) {
		denom = 0.;
		for (k=0; k<ksym; k++) bnew[i][k] = 0.;
		for (t=0; t<nobs-1; t++) {
			term = (alpha[t][i]*beta[t][i]/lhood)
				* powtab[arnrm[t] + brnrm[t] - lrnrm + 6];
			denom += term;
			bnew[i][obs[t]] += term;
		}
		for (j=0; j<mstat; j++) {
			num = 0.;
			for (t=0; t<nobs-1; t++) {
				num += alpha[t][i]*b[j][obs[t+1]]*beta[t+1][j]
					* powtab[arnrm[t] + brnrm[t+1] - lrnrm + 6]/lhood;
			}
			a[i][j] *= (num/denom);
		}
		for (k=0; k<ksym; k++) bnew[i][k] /= denom;
	}
	/*Print New A*/
	for (i=0; i<mstat; i++) {
        for (j=0; j<mstat; j++){
			cout << "New A [" << i << "][" << j << "] = \t" << (a[i][j]) << "\t";
			outstream3 << "New A [" << i << "][" << j << "] = \t" << (a[i][j]) << "\t";
            //cout << "New A = " << (a[i][j]) << "\t";
            //outstream3 << "New A = " << (a[i][j]) << "\t";
        }
        cout << endl;
        outstream3 << endl;
	}
	cout << endl; outstream3 << endl;
	/*Print New B*/
	for (i=0; i<mstat; i++) {
        for (k=0; k<ksym; k++){
			cout << "New B [" << i << "][" << k << "] = \t" << (bnew[i][k]) << "\t";
			outstream3 << "New B [" << i << "][" << k << "] = \t" << (bnew[i][k]) << "\t";
            //cout << "New B = " << (bnew[i][k]) << "\t";
            //outstream3 << "New B = " << (bnew[i][k]) << "\t";
        }
        cout << endl;
        outstream3 << endl;
	}
	cout << endl; outstream3 << endl;
	
	/*Print New Pi*/
	for (i=0; i<mstat; i++) {
		cout << "New Pi [" << i << "] = \t" << (pstate[0][i]) << "\t";
		outstream3 << "New Pi [" << i << "] = \t" << (pstate[0][i]) << "\t";
     }
	 cout << endl;
     outstream3 << endl;
     
	b = bnew;
	fbdone = 0;

	outstream3.close();
}
