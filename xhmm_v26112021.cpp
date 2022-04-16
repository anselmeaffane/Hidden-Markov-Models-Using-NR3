#include <iostream>
#include <fstream>
#include "nr3.h"
#include "hmm.h"
#include <string>
#include <vector>

using namespace std;

int main(){
    
    MatDoub A(2,2,0.);
    MatDoub B(2,27,0.);
    Int n, nn(183);
    VecInt O(nn);
	    
    //Define Transition Matrice
	A[0][0] = 0.47468; A[0][1] = 0.52532; 
	A[1][0] = 0.51656; A[1][1] = 0.48344; 
    //Define Emission Matrice
	B[0][0] = 0.03909; B[0][1] = 0.03537;B[0][2] = 0.03537;B[0][3] = 0.03909; B[0][4] = 0.03583;B[0][5] = 0.03630 ;B[0][6] = 0.04048; B[0][7] = 0.03537;B[0][8] = 0.03816;B[0][9] = 0.03909;B[0][10] = 0.03490;B[0][11] = 0.03723; B[0][12] = 0.03537;B[0][13] = 0.03909;B[0][14] = 0.03397; B[0][15] = 0.03397;B[0][16] = 0.03816;B[0][17] = 0.03676; B[0][18] = 0.04048; B[0][19] = 0.03443; B[0][20] = 0.03537; B[0][21] = 0.03955;B[0][22] = 0.03816;B[0][23] = 0.03723;B[0][24] = 0.03769; B[0][25] = 0.03955 ;B[0][26] = 0.03397 ;
	B[1][0] = 0.03735; B[1][1] = 0.03408;B[1][2] = 0.03455;B[1][3] = 0.03828; B[1][4] = 0.03782;B[1][5] = 0.03922;B[1][6] = 0.03688; B[1][7] = 0.03408;B[1][8] = 0.03875;B[1][9] = 0.04062;B[1][10] = 0.03735;B[1][11] = 0.03968; B[1][12] = 0.03548;B[1][13] = 0.03735 ;B[1][14] = 0.04062; B[1][15] = 0.03595;B[1][16] = 0.03641;B[1][17] = 0.03408; B[1][18] = 0.04062;B[1][19] = 0.03548;B[1][20] = 0.03922; B[1][21] = 0.04062;B[1][22] = 0.03455;B[1][23] = 0.03595;B[1][24] = 0.03408; B[1][25] = 0.03408;  B[1][26] = 0.03688;  
    //Define symbols V_k={a=0,b=1,c=2,d=3,e=4,f=5,g=6,h=7,i=8,j=9,k=10,l=11,m=12,n=13,o=14,p=15,q=16,r=17,s=18,t=19,u=20,v=21,w=22,x=23,y=24,z=25,espace=26}, Observation sequence
    // ovservation={y0=a, y1=}
    // also suppose that current research indicates a correlation between the size of tree 
    // growth rings and temperature For simplicity we only consider three different tree ring 
    // sizes small
    
	Int first;//Variables to select data into the dataset object

    // Upload 2 dimensional Gaussian data point from a file

	ifstream instream;
    instream.open("file.txt");
    for (n=0;n<nn;n++) {
	    instream >> first ;//select 1st column with "first" variable and selest the 2nd colums with "second" variable
		//O.assign(n, first);
		O.addVec(n, first);
    }			
//Create an Numerical Recipe Gaussian object call myHmm
	struct HMM myHmm(A, B, O);
	myHmm.forwardbackward(); //call forwardbackward() function from hmm.h
	myHmm.baumwelch();//Call baumwelch() function from hmm.h

}
