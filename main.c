#include "header.h"




#define numd 1
#define numN 1
#define numg1 1
#define numg2 1


//gamma0 is indexed from 0 to 9, and takes on the values
//{.001,.002,.004,.006,.008,.01,.15,.02,0.025,.03}
//gamma1 is defined in terms of a=gamma1/gamma0, a={0,1,...,9,10}
//_dx_ indicates the value of damage, with x=2d.  The index runs from 0
//to 20, and d takes on the values {0,0.5,1,...,19.5,20}

//human n=20000 g0=0.0006 g1=0.0001 d=12
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//don't change anything above these lines

#define Nrun 10					//number of networks
#define Ntrial 3000				//number of trajectories per network



int Nl[numN]={700};	//N value goes here
double g1l[numg1]={0.015};		//gamma0 value goes here
double g2l[numg2]={0.14};		//gamma1 value goes here
double dl[numd]={2.4};			//damage value goes here

//800,0.0035,0.0035,4
//don't change anything below these lines


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------












extern double ran2(int *seed);
extern void initialize_arrays();		//set zeros
extern void get_BA(int k);					//BA network
extern void get_Ran(int k);					//BA network

extern void randomize_BA_directions();	//randomize edge directions
extern void determine_supply(double q);			//get the required supply of each node

extern void randomly_kill_one_node(int numAlive,int timestep);
extern int propagate_deaths(int timestep);
void run_trajectories();




int iter;
int numAlive;
int qIter,run,trialIter;
int kBarabasi;

int minNumber,minPos,oldNumAlive;
int NumTrials;
double q_scale;
int nRuns;
int kB;
int gammaIter2,damageIter;
double ResProb;
double damage;

double KillProb;
double RessurectProb;
int maxVal;

int gammaIter,runStart;


int DeathVector[Nnodes];
double tempdub;
int NumTrials,nRuns;








int main(){





















	nRuns=Nrun;
	NumTrials=Ntrial;

	runStart=1;
	int Niter,g1iter,g2iter,diter;
	system("mkdir output");


	/*
	f1=fopen("params.txt","w");
	fprintf(f1,"N={");
	for(i=0;i<numN;i++){
		fprintf(f1,"%d",Nl[i]);
		if(i<numN-1){
			fprintf(f1,",");
		}
	}
	fprintf(f1,"}\n");
	fprintf(f1,"g1={");
	for(i=0;i<numg1;i++){
		fprintf(f1,"%f",g1l[i]);
		if(i<numg1-1){
			fprintf(f1,",");
		}
	}
	fprintf(f1,"}\n");
	fprintf(f1,"g2={");
	for(i=0;i<numg2;i++){
		fprintf(f1,"%f",g2l[i]);
		if(i<numg2-1){
			fprintf(f1,",");
		}
	}
	fprintf(f1,"}\n");
	fprintf(f1,"d={");
	for(i=0;i<numd;i++){
		fprintf(f1,"%f",dl[i]);
		if(i<numd-1){
			fprintf(f1,",");
		}
	}
	fprintf(f1,"}\n");
	fclose(f1);
	*/




	printf("\nrunning.....\n");
	printf("# of networks = %d\n",nRuns);
	printf("# of trajectories per network = %d\n",NumTrials);

	for(Niter=0;Niter<numN;Niter++){
		for(run=runStart;run<=nRuns+runStart-1;run++){
			N=Nl[Niter];
			seed=-2-659*run-661*N;
			kBarabasi=2;
			M=(N-1)*kBarabasi;

//			sprintf(fn2,"output_SF/edges_N%d_run%d.txt",N,run);
//			f2=fopen(fn2,"w");
			initialize_arrays();
			get_BA(kBarabasi);

//			for(i=0;i<M;i++){
//				fprintf(f2,"%d %d\n",edges.start[i],edges.end[i]);
//			}
//			fclose(f2);

			for(g1iter=0;g1iter<numg1;g1iter++){
				for(g2iter=0;g2iter<numg2;g2iter++){
					for(diter=0;diter<numd;diter++){


						damage=dl[diter];
						KillProb=g1l[g1iter];
						ResProb=g2l[g2iter];

//						printf("%d %.1f %.1f %.1f %d\n",N,1000*KillProb,1000*ResProb,damage,run);

						if(run==1){
							printf("\nN = %d\n",N);
							printf("g0 = %.2f x 10^(-3)\n",1000*KillProb);
							printf("g1 = %.2f g0\n",ResProb/KillProb);
							printf("d = %.2f\n\n",damage);
						}



						sprintf(fn2,"output/Lifetimes_run%d.txt",run);
						f2=fopen(fn2,"w");
						run_trajectories();
						fclose(f2);



					}//end of loop over damage
				}//end of loop over gamma1
			}//end of loop over gamma0









			/*
			sprintf(fn2,"output_Ran/edges_N%d_run%d.txt",N,run);
			f2=fopen(fn2,"w");
			initialize_arrays();
			get_Ran(kBarabasi);

			for(i=0;i<M;i++){
				fprintf(f2,"%d %d\n",edges.start[i],edges.end[i]);
			}
			fclose(f2);


			for(g1iter=0;g1iter<numg1;g1iter++){
				for(g2iter=0;g2iter<numg2;g2iter++){
					for(diter=0;diter<numd;diter++){


						damage=dl[diter];
						KillProb=g1l[g1iter];
						ResProb=g2l[g2iter];

						printf("Ra:  %d %.1f %.1f %.1f %d\n",N,1000*KillProb,1000*ResProb,damage,run);
						sprintf(fn1,"output_Ran/Time_of_death_N%d_g%d_g%d_d%d_run%d.txt",Niter,g1iter,g2iter,diter,run);
						f1=fopen(fn1,"w");
						//sprintf(fn2,"output_Ran/Lifetime_N%d_k%d_run%d_q%.2f_g%d_g%d_d%d.txt",N,kBarabasi,run,q,gammaIter,gammaIter2,damage);
						//f2=fopen(fn2,"w");
						run_trajectories();
						fclose(f1);
						//fclose(f2);
					}//end of loop over damage
				}//end of loop over gamma1
			}//end of loop over gamma0

			*/




		}//end of loop over runs
	}//end of loop over N
    return 0;
}






void run_trajectories(){


	q=0.5-.000001;
	determine_supply(q);
	for(trialIter=1;trialIter<=NumTrials;trialIter++){
		for(i=0;i<N;i++){
			nodes.status[i]=1;
			nodes.timeOfDeath[i]=0;
			nodes.wasResurrected[i]=0;
		}
		numAlive=N;
		iter=0;
		while(iter<((int)(damage*((double)N)/100.0))){
			maxVal=((int)((double)N)*ran2(&seed));
			if(nodes.status[maxVal]==1){
				nodes.status[maxVal]=0;
				numAlive-=1;
				iter+=1;
			}
		}
		iter=0;
		while(numAlive>0){			//while anyone is still alive
			//			randomly_kill_one_node(numAlive,iter);
			if(iter>0){
				for(i=0;i<N;i++){
					DeathVector[i]=0;
					nodes.wasResurrected[i]=0;
					tempdub=ran2(&seed);
					if(tempdub<KillProb){
						DeathVector[i]=-1;
					}
					else{
						if(tempdub<KillProb+ResProb){
							DeathVector[i]=1;
						}
					}

					if(nodes.status[i]==1 & DeathVector[i]==-1){
						nodes.status[i]=0;
						nodes.timeOfDeath[i]=-iter;
					}
					if(nodes.status[i]==0 & DeathVector[i]==1){
						nodes.status[i]=1;
						nodes.ressurectTime[i]=nodes.timeOfDeath[i];
						nodes.wasResurrected[i]=1;
						nodes.timeOfDeath[i]=0;
					}

				}
			}



			numAlive=propagate_deaths(iter);
			for(i=0;i<N;i++){
				if(nodes.wasResurrected[i]==1 & nodes.status[i]==0){
					nodes.timeOfDeath[i]=nodes.ressurectTime[i];
				}
			}

			if(((double)numAlive)<=((double)N)/100.0){
				fprintf(f2,"%d ",iter);
				break;
			}

			if(((double)numAlive)/((double)N)<0.00){
				break;
			}
			iter++;



			if(iter==20000){
				trialIter=NumTrials+1;
				printf("abandoning\n");
				fprintf(f2,"abandoned\n");
				break;

			}



		}

		/*
		maxVal=0;
		for(i=0;i<N;i++){
			fprintf(f1,"%d ",nodes.timeOfDeath[i]);
			if(nodes.timeOfDeath[i]<0){
				nodes.timeOfDeath[i]=-nodes.timeOfDeath[i];
			}
			if(nodes.timeOfDeath[i]>maxVal){
				maxVal=nodes.timeOfDeath[i];
			}

		}
		fprintf(f1,"\n");
		*/


	//	fprintf(f2,"%d ",maxVal);
	//	maxVal=0;
	//	for(i=0;i<N;i++){
	//		maxVal+=nodes.timeOfDeath[i];
	//	}
		//fprintf(f2,"%f\n",((double)maxVal)/((double)(N)));

	}

}



