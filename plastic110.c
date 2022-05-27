#include<math.h>
#include<stdio.h>
#include<stdlib.h>
//%%%%%%%%%--Randon numbers Parameters--%%%%%%%%%%%%%
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NR_END 1
#define FREE_ARG char*
//#define _MATH_H
#define A1 1.0
#define A2 0.5
#define tau1 1.8
#define tau2 6.0
float ran1(long *idum);
void nrerror(char error_text[]);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define N 1000
#define EQ 2*N
#define Dmax 500
#define Nmax 400
#define Stdp 1 // 0 regular; 1 plasticidade
#define fator3 0.001
int conex1[N],conex2[N],list1[N][Nmax],list2[N][Nmax];
float matrix[N][N];
void derivs(float y[],float ydot[],float s[]);
int main(){		
				long idum;
				int i,j,k,l,m,R,ij[N],contD[N];
				float r,delta[2],deltaw[2],b,VR,h,h0,t,tr,tf,x[EQ],y[EQ],c1[EQ],c2[EQ],c3[EQ],ydot[EQ],s[N],tdisp[N],wij,tq1,wij_max,vec[N][Dmax],isi[EQ],tq2;
				FILE *fw,*fs,*fw1,*fs1,*raster;
				fw=fopen("matrixWij_tq2CI3.dat","w"); // tq1 primeiro tempo; tq2 segundo tempo
				fs=fopen("perfil_tq2CI3.dat","w"); // tq1 primeiro tempo; tq2 segundo tempo
				fw1=fopen("matrixWij_tq1CI3.dat","w"); // tq1 primeiro tempo; tq2 segundo tempo
				fs1=fopen("perfil_tq1CI3.dat","w"); // tq1 primeiro tempo; tq2 segundo tempo
				raster=fopen("raster.dat","w");
				
				//%%%%%%%%%%%%%%%%%%--Parametros--%%%%%%%%%%%%%
				b=70.0; 
				VR=-58.0;
				h=0.01;
			 	h0=h/2.728;
				tf=41000.0; //tempo final
			 	tr=0.0; // transiente      
				tq1=10000.0;  //tempo do perfil 1 ou 2
				tq2=40000;
				wij_max=1.0; // peso máximo da matriz
				wij=0.04; // [0:0.1] // peso inicial da matriz
				r=0.1; // [0:0.1] // alcance do acoplamento
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				for(i=0;i<N;i++){
								for(j=0;j<Nmax;j++){
												list1[i][j]=0;
												list2[i][j]=0;
								}
								for(j=0;j<N;j++)
												matrix[i][j]=0.0;
								conex1[i]=0;
								conex2[i]=0;
				}
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				R=(int)(r*N);
				for(i=0;i<N;i++)
								for(j=0;j<N;j++){
												k=abs(i-j);
												if((i!=j)&&((k<=R)||(k>=(N-R)))){
																list1[i][conex1[i]]=j; // Fornece o neurônio pós
																list2[j][conex2[j]]=i; // Lista invertida, fornece o neurônio pre
																matrix[i][j]=wij;
																conex1[i]++; conex2[j]++;
												}
								}
				//%%%%%%%%%%%%%%%%%%--initial-conditions--%%%%%%%%%%%
				for(i=0;i<N;i++){
								x[i]=-70+20.0*ran1(&idum);  
								x[i+N]=200+200*ran1(&idum);   
								y[i]=x[i];
								y[i+N]=x[i+N];
								ij[i]=5000;
								tdisp[i]=-40.0;
								contD[i]=0;
								for(j=0;j<Dmax;j++)
												vec[i][j]=0.0;
								isi[i]=0.0;
								isi[i+N]=0.0;
				}
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				m=1;
				for(t=0.0;t<tf;t+=h){
								for(i=0;i<N;i++){	
												s[i]=exp(-ij[i]*h0);
												ij[i]++;
								}
								//%%%%%%%%%%%%%--RUNGE-KUTTA-4--%%%%%%%%%%%%%    
								derivs(y,ydot,s);
								for(i=0;i<EQ;i++){
												c1[i]=h*ydot[i];
												y[i]=x[i]+c1[i]*0.5;
								}
								derivs(y,ydot,s);
								for(i=0;i<EQ;i++){
												c2[i]=h*ydot[i];
												y[i]=x[i]+c2[i]*0.5;
								}
								derivs(y,ydot,s);
								for(i=0;i<EQ;i++){
												c3[i]=h*ydot[i];
												y[i]=x[i]+c3[i]; 
								}
								derivs(y,ydot,s);
								for(i=0;i<EQ;i++)
												y[i]=x[i]+(c1[i]+h*ydot[i])/6.0+(c2[i]+c3[i])/3.0;
								//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
								for(i=0;i<N;i++){
												if(y[i]>-35.0){ 
												
												fprintf(raster,"%f %d\n",t,i);
																tdisp[i]=t;
																if(Stdp==1){
																				//%%%%%%%%%%%%%%--Stdp--%%%%%%%%%%%%%%%%%%%
																				for(j=0;j<conex1[i];j++){
																								delta[0]=tdisp[i]-tdisp[list1[i][j]];  	// pre-to-pós																			
																								delta[1]=tdisp[list2[i][j]]-tdisp[i];   // pós-to-pre
																								for(l=0;l<=1;l++)  //// Calculando os Delta W
																												if(delta[l]>=0)
																																deltaw[l]=A1*exp(-delta[l]/tau1);   
																												else
																																deltaw[l]=-A2*exp(delta[l]/tau2);
																								//// Alterando o peso na matrix
																								matrix[i][list1[i][j]]=matrix[i][list1[i][j]]+fator3*deltaw[0];
																								matrix[list2[i][j]][i]=matrix[list2[i][j]][i]+fator3*deltaw[1];
																								//// Limitando valor mínimo dos pesos da matriz // sempre positivo
																								if(matrix[i][list1[i][j]]<0)
																												matrix[i][list1[i][j]]=0;
																								if(matrix[list2[i][j]][i]<0)
																												matrix[list2[i][j]][i]=0;
																								//// Limitando valor máximo e mínimo dos elementos da matriz
																								if(matrix[i][list1[i][j]]>wij_max)
																												matrix[i][list1[i][j]]=wij_max;
																								if(matrix[list2[i][j]][i]>wij_max)
																												matrix[list2[i][j]][i]=wij_max;
																				}
																}
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
																y[i]=VR;
															 	y[i+N]+=b;
															 	ij[i]=0;
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
																if((t>tr)&&(contD[i]<Dmax)){
																				vec[i][contD[i]]=t;
																				if(contD[i]>0){
																								isi[i]=isi[i]+(t-vec[i][contD[i]-1]);
																								isi[i+N]=isi[i+N]+(t-vec[i][contD[i]-1])*(t-vec[i][contD[i]-1]);
																				}
																				contD[i]++;
																}
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
												}
												x[i]=y[i];
												x[i+N]=y[i+N];
								}
								//%%%%%%%%%%%%%--Perfil-Espacial--%%%%%%%%%%%%%
								if((t>=tq1)&&(m==1)){
												for(i=0;i<N;i++){
																fprintf(fs1,"%d\t%f\n",i,y[i]);
																for(j=0;j<N;j++)
																				fprintf(fw1,"%d\t%d\t%f\n",i,j,matrix[i][j]);
																fprintf(fw1,"\n");
												}
												
														
												m=2;
								}
												
								if((t>=tq2)&&(m==2)){
												for(i=0;i<N;i++){
																fprintf(fs,"%d\t%f\n",i,y[i]);
																for(j=0;j<N;j++)
																				fprintf(fw,"%d\t%d\t%f\n",i,j,matrix[i][j]);
																fprintf(fw,"\n");
												}				
												
												
												
												
												
												
												m=3;
								}
				}
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				x[0]=0.0;
				x[1]=0.0;
				x[4]=0.0;
				for(i=0;i<N;i++){
								contD[i]=contD[i]-1;
								isi[i]=isi[i]/contD[i];
								isi[i+N]=isi[i+N]/contD[i];
								x[0]=x[0]+isi[i];
								x[1]=x[1]+isi[i+N];
								for(j=0;j<conex1[i];j++)
												x[4]=x[4]+matrix[i][list1[i][j]];
				}
				x[0]=x[0]/N; 
				x[1]=x[1]/N;
				x[2]=sqrt(fabs(x[1]-x[0]*x[0]))/x[0]; 
				x[3]=1000/x[0];
				x[4]=x[4]/(N*R);
				printf(" gexM=%f; R=%d; CvM=%f; freqM=%f\n",x[4],R,x[2],x[3]);
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				return(0);
}
//%%%%%%%%%%%%%%%%%%%%%--DERIVS--%%%%%%%%%%%%%%%%%%%%%%%
void derivs(float y[],float ydot[],float s[]){
				int i,j;
				float yrev,a,C,gL,EL,deltaT,VT,tauW,gq[N],I0;
				I0=500.0;
				a=2.0;
				C=200.0;
				gL=12.0; 
				EL=-70.0;
				VT=-50.0;
				tauW=300.0;
				deltaT=2.0;
				yrev=0.0;
				for(i=0;i<N;i++){
								gq[i]=0.0;
								for(j=0;j<conex1[i];j++)
												gq[i]=gq[i]+(yrev-y[i])*s[list1[i][j]]*matrix[i][list1[i][j]];
								ydot[i]=(-gL*(y[i]-EL)+gL*deltaT*exp((y[i]-VT)/deltaT)-
																y[i+N]+I0+gq[i])/C;
								ydot[i+N]=(a*(y[i]-EL)-y[i+N])/tauW;
				}
}
//%%%%%%%%%%%%%%%%%%%%%%%--GERADOR-RANDOMICO--%%%%%%%%%%%%%%%
float ran1(long *idum){
				int j; long k;
				static long iy=0;
				static long iv[NTAB];
				float temp; 
				if(*idum<=0 || !iy){
								if(-(*idum)<1)
												*idum=1;
								else
												*idum=-(*idum);
								for(j=NTAB+7;j>=0;j--){
												k=(*idum)/IQ;
												*idum=IA*(*idum-k*IQ)-IR*k;
												if(*idum<0) 
																*idum +=IM;
												if(j<NTAB) 
																iv[j]=*idum;
								}
								iy=iv[0];
				}
				k=(*idum)/IQ;
				*idum=IA*(*idum-k*IQ)-IR*k;
				if(*idum<0)
								*idum += IM;
				j=iy/NDIV; iy=iv[j];
				iv[j]=*idum;
				if((temp=AM*iy)>RNMX)
								return RNMX;
				else 
								return temp;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void nrerror(char error_text[]){
				fprintf(stderr,"Numerical Recipes run-time error...\n");
				fprintf(stderr,"%s\n",error_text);
				fprintf(stderr,"...now exiting to system...\n");
				exit(1);
}
