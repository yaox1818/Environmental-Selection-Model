 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import utils.random.MathUtil;
import utils.random.Multinomial2;


public class SelectiveIndividualE extends Individual {
    
    private double hostFitness;
    int[] microbeGenesFitnessInHost;
    int[] hostGenesFitness;
    private int numOfMGenes;
    double[] microbeFitnessInHostRecords;
    double[] microbeFitnessInHostRecordsK;
    double[] hostFitnessRecords;
    double[] hostFitnessRecordsK;
    double k;

    public SelectiveIndividualE(Multinomial2 MultDist, int nomph, int noes, SelectiveSpeciesRegistry ssr, boolean HMS_or_TMS) {
        super(MultDist, nomph, noes);
        numOfMGenes=ssr.numOfTolGenes;
        microbeGenesFitnessInHost=new int[numOfMGenes];//此个体中微生物对所有特征的fitness
        hostGenesFitness=new int[numOfMGenes];//此个体中宿主对所有特征的fitness
        microbeFitnessInHostRecords=new double[noes];
        microbeFitnessInHostRecordsK=new double[noes];
        hostFitnessRecords=new double[noes];
        hostFitnessRecordsK=new double[noes];
		if (HMS_or_TMS){//true：为HMS，每个个体有不同的microbeGenesFitness
			for(int i=0;i<numOfMGenes;i++){
				microbeGenesFitnessInHost[i]=MathUtil.getNextInt(2)-1;
			}
		}else microbeGenesFitnessInHost=ssr.mfrH;
		hostGenesFitness=ssr.hfr;
        //k= Math.abs(MathUtil.getNextGaussian(parentalk[i]));
        k= MathUtil.getNextFloat(1);
        //k=MathUtil.getNextInt(10);
        //k=0;
        //k/=10;
        //microbeFitnessRecords为微生物对所有微生物的mf，hostFitnessRecords为宿主对所有OTU的fitness，microbeGenesFitness为微生物对各种trait的fitness，hostGenesFitness为宿主对各种trait的fitness
        ssr.getFitness(microbeFitnessInHostRecords,hostFitnessRecords, microbeGenesFitnessInHost,hostGenesFitness);//计算此个体中每种OTU自身和对宿主的fitness
        hostProvideResources(microbeFitnessInHostRecords,microbeFitnessInHostRecordsK,hostFitnessRecords,hostFitnessRecordsK,k,noes);
        hostFitness=ssr.getTotalFitness(microbiome,hostFitnessRecordsK);//计算此宿主个体的总fitness
        /*ssr.getMFitnessSelection(microbiome,microbeFitnessRecordsk);//因为原本的代码第一代宿主是根据环境随机生成的，现在的代码是第一代宿主会根据微生物的适应性再次生成代码所以不一样
        MultDist.updateProb(microbiome);
        MultDist.multisample(microbiome, nomph);
        hostFitness=ssr.getTotalFitness(microbiome,hostFitnessRecordsk);**/
        //hostFitness=ssr.getTotalFitness(microbiome,hostFitnessRecordsk);//计算此宿主个体的总fitness
        //microbeFitnessRecords为微生物对所有微生物的mf，hostFitnessRecords为宿主对所有OTU的fitness，microbeGenesFitness为微生物对各种trait的fitness，hostGenesFitness为宿主对各种trait的fitness
    }
    
    public double getHostFitness(double hostSelectionCoef){
        if(hostSelectionCoef==0)
            return hostFitness;
        else
            return Math.pow(hostSelectionCoef,hostFitness);
    }

    public void hostProvideResources(double[] microbeFitnessInHostRecords, double[] microbeFitnessInHostRecordsK,
                                     double[] hostFitnessRecords, double[] hostFitnessRecordsK, double k, int noes){
        for (int i = 0; i < noes; i++) {
            if (hostFitnessRecords[i]>0){
                microbeFitnessInHostRecordsK[i]= microbeFitnessInHostRecords[i]+hostFitnessRecords[i]*k;
                hostFitnessRecordsK[i]=hostFitnessRecords[i]*(1-k);
            }
            else if (hostFitnessRecords[i]<0){
                microbeFitnessInHostRecordsK[i]= microbeFitnessInHostRecords[i]+hostFitnessRecords[i]*k;
                hostFitnessRecordsK[i]=hostFitnessRecords[i]*(1+k);
            }
            else {
                microbeFitnessInHostRecordsK[i] = microbeFitnessInHostRecords[i];
                hostFitnessRecordsK[i] = hostFitnessRecords[i];
            }
        }
    }
    public double getCosTheta(){
        double innerProduct=0;
        double Mlength=0;
        double Hlength=0;
        for(int i=0;i<numOfMGenes;i++){
            innerProduct+=microbeGenesFitnessInHost[i]*hostGenesFitness[i];
            Mlength+=microbeGenesFitnessInHost[i]*microbeGenesFitnessInHost[i];
            Hlength+=hostGenesFitness[i]*hostGenesFitness[i];
        }
        return innerProduct/Math.sqrt(Hlength*Mlength);
    }
	
    public double[]  goodVersusBad(){
        double[] bacteriaContents=new double[3];
        for(int i=0;i<hostFitnessRecords.length;i++){
            if (hostFitnessRecords[i]<0)
                 bacteriaContents[0]+=getMicrobiome()[i];
            else if (hostFitnessRecords[i]==0)
                 bacteriaContents[1]+=getMicrobiome()[i];
            else
                 bacteriaContents[2]+=getMicrobiome()[i];
        }
        return bacteriaContents;
    }
  
    public String printBacteriaContents(){
       StringBuilder sb=new StringBuilder();
       double[] bacteriaContents=goodVersusBad();
       for (int i=0;i<bacteriaContents.length;i++){
            sb.append(bacteriaContents[i]).append('\t');
       }
       return sb.toString().trim();
    }
    
}
