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
double ran1(long *idum);
void nrerror(char error_text[]);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define N 1001
#define EQ 2*N
#define Dmax 400
#define fator3 0.001
#define STDP 1 // 1 plasticidade; 0 sem plasticidade
int conex1[N],conex2[N],list1[N][N],list2[N][N];
double matrix[N][N];
void derivs(double y[],double ydot[],double z[]);
int main(){
				long idum;
				int i,ii,j,jj,k,l,m,R,ij[N],contD[N],variable,iii,aux_con;
				double r,delta[2],deltaw[2],b,VR,h,h0,t,tr,tf,x[EQ],y[EQ],xant[i],
				c1[EQ],c2[EQ],c3[EQ],ydot[EQ],s[N],tdisp[N],wij,z[N],
				wij_max,vec[N][Dmax],isi[2];
				char file1[100];
				FILE *ff,*fg,*fw;
				wij_max=1.0;
				wij=0.0;
				r=0.5; //// [0:0.1]
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				b=70.0;
				VR=-58.0;
				h=0.01;
				h0=h/2.728;
				tf=1000.0;
				tr=0.0;

				scanf("%d", &variable);
				jj=variable;

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				for(i=0;i<N;i++){
								for(j=0;j<N;j++){
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
				for(i=0;i<N;i++){
								for(j=0;j<N;j++){
												k=abs(i-j);
												if((i!=j)&&((k<=R)||(k>=(N-R)))){
													aux_con=0;
													for(iii=0;iii<=conex1[i]; iii++)
													if(list1[i][iii]==j) aux_con=1;

																if(aux_con==0)
																{
																list1[i][conex1[i]]=j; // Fornece o neurÃ´nio pÃ³s
																conex1[i]++;
																}

																aux_con=0;
																for(iii=0;iii<=conex2[j]; iii++)
																if(list2[j][iii]==i) aux_con=1;

																if(aux_con==0)
																{
																list2[j][conex2[j]]=i; // Lista invertida, fornece o neurÃ´nio pre
																conex2[j]++;
																}

																matrix[i][j]=wij;
												}
								}
				}
				ii=jj;
				sprintf(file1,"diags_%d.dat",ii);
				ff=fopen(file1,"w");
				sprintf(file1,"pesos_%d.dat",ii);
				fg=fopen(file1,"w");
				sprintf(file1,"raster_%d.dat",ii);
				fw=fopen(file1,"w");
				//%%%%%%%%%%%%%%%%%%--initial-conditions--%%%%%%%%%%%
				for(i=0;i<N;i++){
								x[i]=-80.0+40.0*ran1(&idum);
								x[i+N]=400*ran1(&idum);
								if(ii==1){
												//%%%--fake-chimera--%%%%%
												if((i<300)||(i>700)){
																x[i]=-60.0;
																x[i+N]=200.0;
												}
								}
								if(ii==2){
												//%%%%%--inclinado--%%%%%
												x[i]=-70.0;
												x[i+N]=200+i*80.0/(N-1); ////i*80.0/(N-1);
								}
								if(ii==3){
												//%%%%--sincronizado--%%%%
												x[i]=-70.0;
												x[i+N]=200.0;
								}
								//%%%%%%%%%%%%%%%%%%%%%%%
								y[i]=x[i];
								y[i+N]=x[i+N];
								ij[i]=5000;
								tdisp[i]=-100.0;
								contD[i]=0;
								for(j=0;j<Dmax;j++)
												vec[i][j]=0.0;
				}
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				m=0; isi[0]=isi[1]=0.0;
				for(t=0.0;t<tf;t+=h){

								for(i=0;i<N;i++) z[i]=0.0;

								for(i=0;i<N;i++){
												s[i]=exp(-ij[i]*h0);
												ij[i]++;

												for(j=0;j<conex1[i];j++)
												{
												z[list1[i][j]]= z[list1[i][j]]+ (0.0-y[list1[i][j]])*s[i]*matrix[i][list1[i][j]];
												}
								}

								//fprintf(fw,"%f %f\n",t, z[0]);

								//%%%%%%%%%%%%%--RUNGE-KUTTA-4--%%%%%%%%%%%%%
								derivs(y,ydot,z);
								for(i=0;i<EQ;i++){
												c1[i]=h*ydot[i];
												y[i]=x[i]+c1[i]*0.5;
								}
								derivs(y,ydot,z);
								for(i=0;i<EQ;i++){
												c2[i]=h*ydot[i];
												y[i]=x[i]+c2[i]*0.5;
								}
								derivs(y,ydot,z);
								for(i=0;i<EQ;i++){
												c3[i]=h*ydot[i];
												y[i]=x[i]+c3[i];
								}
								derivs(y,ydot,z);
								for(i=0;i<EQ;i++)
												y[i]=x[i]+(c1[i]+h*ydot[i])/6.0+(c2[i]+c3[i])/3.0;
								//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
								for(i=0;i<N;i++){
												if(y[i]>=-35.0){
																tdisp[i]=t-h-(xant[i]*h)/(x[i]-xant[i]);
																if(t>=tr){
																				fprintf(fw,"%d\t%f\n",i,t);
																				fflush(fw);
																}

																//%%%%%%%%%%%%%%--Stdp--%%%%%%%%%%%%%%%%%%%
																if(STDP==1){
																				for(j=0;j<conex1[i];j++){
																								delta[0]=tdisp[i]-tdisp[list2[i][j]];  	// pre-to-pÃ³s
																								delta[1]=tdisp[list1[i][j]]-tdisp[i];   // pÃ³s-to-pre
																								for(l=0;l<=1;l++)  //// Calculando os Delta W
																												if(delta[l]>=0)
																																deltaw[l]=A1*exp(-delta[l]/tau1);
																												else
																																deltaw[l]=-A2*exp(delta[l]/tau2);
																								//// Alterando o peso na matrix
																								matrix[i][list1[i][j]]=matrix[i][list1[i][j]]+fator3*deltaw[1];
																								matrix[list2[i][j]][i]=matrix[list2[i][j]][i]+fator3*deltaw[0];
																								//// Limitando valor mÃ­nimo dos pesos da matriz // sempre positivo
																								if(matrix[i][list1[i][j]]<0)
																												matrix[i][list1[i][j]]=0;
																								if(matrix[list2[i][j]][i]<0)
																												matrix[list2[i][j]][i]=0;
																								//// Limitando valor mÃ¡ximo e mÃ­nimo dos elementos da matriz
																								if(matrix[i][list1[i][j]]>wij_max)
																												matrix[i][list1[i][j]]=wij_max;
																								if(matrix[list2[i][j]][i]>wij_max)
																												matrix[list2[i][j]][i]=wij_max;
																				}
																}
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
																y[i]=VR; y[i+N]+=b; ij[i]=0;
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
																if((t>tr)&&(contD[i]<Dmax)){
																				vec[i][contD[i]]=t;
																				if(contD[i]>0){
																								m++;
																								isi[0]+=(t-vec[i][contD[i]-1]);
																								isi[1]+=(t-vec[i][contD[i]-1])*(t-vec[i][contD[i]-1]);
																				}
																				contD[i]++;
																}
																//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
												}
												x[i]=y[i]; x[i+N]=y[i+N];
												xant[i]=x[i];
								}
				}
				//%%%%%%%%%%%%%%%%%--PARAMETRO-ORDEM--%%%%%%%%%%%%%%%%%%%%%%%
				if(m>0){
								isi[0]=isi[0]/m;
								isi[1]=isi[1]/m;
				}
				for(i=0;i<N;i++){
								ij[i]=0;
								//	fprintf(fw,"%d\t%f\n",i,y[i]);
								for(j=0;j<N;j++)
												fprintf(fg,"%d\t%d\t%f\n",i,j,matrix[i][j]);
				}
				fclose(fg); fclose(fw);
				x[0]=10.0*h;
				x[1]=tr+100.0;
				x[2]=tf-100.0;
				x[3]=0.0; l=0;
				/*
					 for(t=x[1];t<x[2];t+=x[0]){
					 l++;
					 x[5]=0.0;
					 x[4]=0.0;
					 for(i=0;i<N;i++){
					 s[i]=6.28*(t-vec[i][ij[i]])/(vec[i][ij[i]+1]-vec[i][ij[i]]);
					 if((vec[i][ij[i]+1]>0.0)&&(t>vec[i][ij[i]+1]))
					 ij[i]++;
					 x[4]=x[4]+cos(s[i]);
					 x[5]=x[5]+sin(s[i]);
					 }
					 x[4]=x[4]/N;
					 x[5]=x[5]/N;
					 x[3]+=sqrt(x[4]*x[4]+x[5]*x[5]);
					 }
					 x[0]=x[3]/l; // parametro de ordem global
					 */
				if(isi[0]>0.0){
								x[1]=1000.0/isi[0]; //frequencia media
								x[2]=sqrt(fabs(isi[1]-isi[0]*isi[0]))/isi[0]; // Cv
				}
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				k=0; x[0]=x[3]=x[4]=0.0;
				for(i=0;i<N;i++)
								for(j=0;j<N;j++)
												if(matrix[i][j]>0.0){
																k++;
																x[0]+=matrix[i][j];
																x[3]+=matrix[i][j]*matrix[i][j];
												}
				if(k>0){
								x[0]=x[0]/k;
								x[3]=x[3]/k;
								x[4]=sqrt(fabs(x[3]-x[0]*x[0]));
				}
				fprintf(ff,"%d\t%f\t%f\n",ii+1,x[0],x[4]);
				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				fclose(ff);

				return(0);
}
//%%%%%%%%%%%%%%%%%%%%%--DERIVS--%%%%%%%%%%%%%%%%%%%%%%%
void derivs(double y[],double ydot[],double z[]){
				int i;
				double a,C,gL,EL,deltaT,VT,tauW,I0;
				a=2.0;
				C=200.0;
				gL=12.0;
				EL=-70.0;
				VT=-50.0;
				tauW=300.0;
				deltaT=2.0;
				I0=500.0;
				for(i=0;i<N;i++){
								ydot[i]=(-gL*(y[i]-EL)+gL*deltaT*exp((y[i]-VT)/deltaT)-y[i+N]+I0+z[i])/C;
								ydot[i+N]=(a*(y[i]-EL)-y[i+N])/tauW;
				}
}
//%%%%%%%%%%%%%%%%%%%%%%%--GERADOR-RANDOMICO--%%%%%%%%%%%%%%%
double ran1(long *idum){
				int j; long k;
				static long iy=0;
				static long iv[NTAB];
				double temp;
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
