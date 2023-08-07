/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import java.io.IOException;
import java.util.Set;

import utils.BinarySearch;
import utils.DiversityIndex;
import utils.VectorAddition;
import utils.random.MathUtil;
/**
 *
 * @author qz28
 */
public class SelectivePopulationE extends Population{
    
    private SelectiveIndividualE[] compositionOfIndividuals;
    private SelectiveSpeciesRegistry selSpReg;
	
    double hsCoef;
    private double averageHF;
    private double varianceHF;
	
	double[] fitnessList;
	double[] hostFitnessTotals;
	double[] microFitnessInEnv;
	double[] microFitnessInHost;
	double[] hostFitness;

    int[][] parentalMGeneFitnessInHost;
    int[][] parentalHGeneFitness;
    double[][] parentalMicrobeFitnessInHost;
    double[][] parentalMicrobeFitnessInHostK;
    double[][] parentalHostFitness;
    double[][] parentalHostFitnessK;
    double[] parentalK;

    private double a_diversityInEnvironment;

	
    
    public SelectivePopulationE(int numberOfMicrobePerHost, double[] environment, int noOfIndividual, double environmentalFactor, double pooledOrFiexd, int numberOfSamples, int sampleReplicates0, SelectiveSpeciesRegistry ssr, double hostSelectionCoef, boolean HMS_or_TMS) throws IOException {
            super(numberOfMicrobePerHost, environment, noOfIndividual, environmentalFactor, pooledOrFiexd, numberOfSamples, sampleReplicates0);
            compositionOfIndividuals = new SelectiveIndividualE[numberOfIndividual];
            fitnessList=new double[numberOfIndividual];
            hostFitnessTotals=new double[numberOfIndividual];
            hostFitness=new double[numberOfEnvironmentalSpecies];
            microFitnessInHost=new double[numberOfEnvironmentalSpecies];
            microFitnessInEnv=new double[numberOfEnvironmentalSpecies];
            parentalMGeneFitnessInHost =new int[numberOfIndividual][ssr.numOfTolGenes];
            parentalHGeneFitness=new int[numberOfIndividual][ssr.numOfTolGenes];
            parentalMicrobeFitnessInHost =new double[numberOfIndividual][numberOfEnvironmentalSpecies];
            parentalMicrobeFitnessInHostK =new double[numberOfIndividual][numberOfEnvironmentalSpecies];
            parentalHostFitness=new double[numberOfIndividual][numberOfEnvironmentalSpecies];
            parentalHostFitnessK=new double[numberOfIndividual][numberOfEnvironmentalSpecies];
            parentalK =new double[numberOfIndividual];
            selSpReg=ssr;
            hsCoef=hostSelectionCoef;
            averageHF=0;
            varianceHF=0;
            ssr.getFitness(microFitnessInEnv,ssr.mfrE);
            ssr.getFitness(microFitnessInHost,hostFitness, ssr.hfr, ssr.mfrH);
            for (int i=0;i<numberOfIndividual;i++){
                compositionOfIndividuals[i]=new SelectiveIndividualE(multiNomialDist,numberOfMicrobePerHost,numberOfEnvironmentalSpecies,ssr,HMS_or_TMS);
                //生成一系列个体，包括个体的hostfitness，个体对各种out的额fitness，各种otu在宿主体内的fitness
                fitnessList[i]=compositionOfIndividuals[i].getHostFitness(0);//记录每一个宿主个体的average score,取值-1到1.
            }           
        }
    
    @Override
    public void parentalInheritanceAndEnvironmentalAcquisition() {
        double runningTotals=0;//根据宿主的average score计算最终fitness放到hostFitnessTotals中
        for (int i=0;i<numberOfIndividual;i++){
            runningTotals+=Math.pow(hsCoef,fitnessList[i]);
            hostFitnessTotals[i]=runningTotals;
        }
        int[] oldAncestryIndex= new int[numberOfIndividual];
        System.arraycopy(ancestryIndex,0,oldAncestryIndex, 0,numberOfIndividual);
        for(int i=0;i<numberOfIndividual;i++){
            int index=BinarySearch.binarySearch(hostFitnessTotals, MathUtil.getNextFloat(runningTotals));//根据fitness随机选取亲代
            parentalIndex[i]=index;
            ancestryIndex[i]=oldAncestryIndex[parentalIndex[i]];
            System.arraycopy(getIndividuals()[index].microbiome,
            0,parentalContribution[i], 0,numberOfEnvironmentalSpecies);//将抽取到的亲本体内微生物赋给parentalContribution
            System.arraycopy(getIndividuals()[index].microbeGenesFitnessInHost,
                                0, parentalMGeneFitnessInHost[i], 0,getSpeciesRegistry().numOfTolGenes);//将抽取到的亲本的微生物对各种trait的fitness赋给parentalMGeneFitness
            System.arraycopy(getIndividuals()[index].hostGenesFitness,
                                0,parentalHGeneFitness[i], 0,getSpeciesRegistry().numOfTolGenes);//将抽取到的亲本宿主对各种trait的fitness赋给parentalHGeneFitness
            System.arraycopy(getIndividuals()[index].microbeFitnessInHostRecords,
                                0, parentalMicrobeFitnessInHost[i], 0,numberOfEnvironmentalSpecies);//将抽取到的亲本各种微生物的mf赋给parentalMicrobeFitness
            System.arraycopy(getIndividuals()[index].microbeFitnessInHostRecordsK,
                    0, parentalMicrobeFitnessInHostK[i], 0,numberOfEnvironmentalSpecies);
            System.arraycopy(getIndividuals()[index].hostFitnessRecords,
                                0,parentalHostFitness[i], 0,numberOfEnvironmentalSpecies);//将抽取到的亲本宿主对所有OTU的fitness赋给parentalHostFitness
            System.arraycopy(getIndividuals()[index].hostFitnessRecordsK,
                    0,parentalHostFitnessK[i], 0,numberOfEnvironmentalSpecies);
            parentalK[i]=getIndividuals()[index].k;
        }
        
        
        VectorAddition.additionOfVectors(environmentalContribution, percentageOfpooledOrFixed, 1 - percentageOfpooledOrFixed, microbiomeSum,
                                         initialEnvironment);//计算环境
        //添加选择
        getSpeciesRegistry().getMFitnessSelectionInEnv(environmentalContribution,microFitnessInEnv);
        
        for (int i = 0; i < numberOfIndividual; i++) {
			VectorAddition.additionOfVectors(mixedContribution, coefficient,//(1 - environmentalFactor) / this.numberOfMicrobePerHost
					environmentalFactor, parentalContribution[i],
					environmentalContribution);//计算每个个体的微生物来源
                        getSpeciesRegistry().getMFitnessSelectionInHost(mixedContribution, parentalMicrobeFitnessInHostK[i]);
			multiNomialDist.updateProb(mixedContribution);//根据mixedContribution和parentalMicrobeFitness生成的概率
			multiNomialDist.multisample(getIndividuals()[i].microbiome, numberOfMicrobePerHost);//根据updateProb生成的概率将微生物随机分配给microbiome
                        //getIndividuals()[i].microbeFitnessRecords=parentalMicrobeFitness[i];
                        //getIndividuals()[i].hostFitnessRecords=parentalHostFitness[i];
                        //getIndividuals()[i].microbeFitnessInHostRecords= parentalMicrobeFitnessInHost[i];
                        //getIndividuals()[i].hostFitnessRecords=parentalHostFitness[i];
                        //getIndividuals()[i].microbeGenesFitnessInHost= parentalMGeneFitnessInHost[i];//继承亲代个体对每种trait对自身的fitness
                        //getIndividuals()[i].hostGenesFitness=parentalHGeneFitness[i];//继承亲代个体对每种trait对宿主的fitness
            System.arraycopy(parentalMGeneFitnessInHost[i], 0,getIndividuals()[i].microbeGenesFitnessInHost, 0,getSpeciesRegistry().numOfTolGenes);
            System.arraycopy(parentalHGeneFitness[i], 0,getIndividuals()[i].hostGenesFitness, 0,getSpeciesRegistry().numOfTolGenes);
            System.arraycopy(parentalMicrobeFitnessInHost[i], 0,getIndividuals()[i].microbeFitnessInHostRecords, 0,numberOfEnvironmentalSpecies);
            System.arraycopy(parentalMicrobeFitnessInHostK[i], 0,getIndividuals()[i].microbeFitnessInHostRecordsK, 0,numberOfEnvironmentalSpecies);
            System.arraycopy(parentalHostFitness[i], 0,getIndividuals()[i].hostFitnessRecords, 0,numberOfEnvironmentalSpecies);
            System.arraycopy(parentalHostFitnessK[i], 0,getIndividuals()[i].hostFitnessRecordsK, 0,numberOfEnvironmentalSpecies);
            getIndividuals()[i].k= parentalK[i];//继承亲代的k值
            // getSpeciesRegistry().getFitness(getIndividuals()[i].microbeFitnessRecords,getIndividuals()[i].hostFitnessRecords, getIndividuals()[i].microbeGenesFitness,getIndividuals()[i].hostGenesFitness);
            /**for (int j = 0; j < numberOfEnvironmentalSpecies; j++) {
                if (getIndividuals()[i].hostFitnessRecords[j] > 0) {
                    getIndividuals()[i].microbeFitnessInHostRecordsK[j] = getIndividuals()[i].microbeFitnessInHostRecords[j] + getIndividuals()[i].hostFitnessRecords[j] * getIndividuals()[i].k;
                    getIndividuals()[i].hostFitnessRecordsK[j] = getIndividuals()[i].hostFitnessRecords[j] * (1 - getIndividuals()[i].k);
                } else if (getIndividuals()[i].hostFitnessRecords[j] < 0) {
                    getIndividuals()[i].microbeFitnessInHostRecordsK[j] = getIndividuals()[i].microbeFitnessInHostRecords[j] + getIndividuals()[i].hostFitnessRecords[j] * getIndividuals()[i].k;
                    getIndividuals()[i].hostFitnessRecordsK[j] = getIndividuals()[i].hostFitnessRecords[j] * (1 + getIndividuals()[i].k);
                }
                else {
                    getIndividuals()[i].microbeFitnessInHostRecordsK[j] = getIndividuals()[i].microbeFitnessInHostRecords[j];
                    getIndividuals()[i].hostFitnessRecordsK[j] = getIndividuals()[i].hostFitnessRecords[j];
                }
            }*/
            fitnessList[i]=getSpeciesRegistry().getTotalFitness(getIndividuals()[i].getMicrobiome(), getIndividuals()[i].hostFitnessRecordsK);
        }
    }
    
    public double averageHostFitness(){
        double totalF=0;
        for(double fitness:fitnessList){
            totalF+=fitness;
        }
        averageHF=totalF/fitnessList.length;
        return averageHF;
    }
    
    public double varianceHostFitness(){
        double totalF2=0;
        for(double fitness:fitnessList){
            totalF2+=fitness*fitness;
        }
        varianceHF=totalF2/fitnessList.length-averageHF*averageHF;
        return varianceHF;
    }
    
    public double cosOfMH(){
        double total=0;
        for(SelectiveIndividualE host:getIndividuals()){
            total+=host.getCosTheta();
        }
        return total/numberOfIndividual;
    }
    
    public String printOutHFitness() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numberOfIndividual; i++) {
			sb.append(Math.pow(hsCoef,fitnessList[i]-1)).append("\t");
		}
		return sb.toString().trim();
	}
    
    public String printOutMFitnessInHost() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numberOfIndividual; i++) {
			sb.append(Math.pow(selSpReg.mscih,selSpReg.getTotalFitness(getIndividuals()[i].getMicrobiome(), getIndividuals()[i].microbeFitnessInHostRecords)-1)).append("\t");
		}
		return sb.toString().trim();
	}
    public String printOutMFitnessInEnvironment() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < numberOfIndividual; i++) {
            sb.append(Math.pow(selSpReg.mscih,selSpReg.getTotalFitness(getIndividuals()[i].getMicrobiome(), microFitnessInEnv)-1)).append("\t");
        }
        return sb.toString().trim();
    }
    public String printBacteriaContents() {
                double[] bacteriaContents=new double[3];
                for (SelectiveIndividualE host:getIndividuals()){
                       VectorAddition.additionOfVectors(bacteriaContents,1,1, host.goodVersusBad(),bacteriaContents);
                }
                StringBuilder sb= new StringBuilder();
                for (int i=0;i<3;i++){
                       sb.append(bacteriaContents[i]).append("\t");
                }
                return sb.toString().trim();
    }
    public String printk() {
        double[] kop=new double[numberOfIndividual];
        for (int i = 0; i < numberOfIndividual; i++) {
            kop[i]=getIndividuals()[i].k;
        }
        StringBuilder sb= new StringBuilder();
        for (int i=0;i<numberOfIndividual;i++){
            sb.append(kop[i]).append("\t");
        }
        return sb.toString().trim();
    }

    public double alphaDiversityInEnvironment() {
        a_diversityInEnvironment = 0;
        a_diversityInEnvironment+= DiversityIndex.shannonWienerIndex(environmentalContribution);
        return a_diversityInEnvironment;
    }

    public double[]  divideBacteriaForHost(){
        double[] bacteriaClass=new double[3];
        for(int i=0;i<hostFitness.length;i++){
            if (hostFitness[i]<0)
                bacteriaClass[0]+=environmentalContribution[i];
            else if (hostFitness[i]==0)
                bacteriaClass[1]+=environmentalContribution[i];
            else
                bacteriaClass[2]+=environmentalContribution[i];
        }
        return bacteriaClass;
    }
    public double[]  divideBacteria(){
        double[] bacteriaClass=new double[3];
        for(int i=0;i<microFitnessInEnv.length;i++){
            if (microFitnessInEnv[i]<0)
                bacteriaClass[0]+=environmentalContribution[i];
            else if (microFitnessInEnv[i]==0)
                bacteriaClass[1]+=environmentalContribution[i];
            else
                bacteriaClass[2]+=environmentalContribution[i];
        }
        return bacteriaClass;
    }
    public String printBacteriaContentsInEnvForHost() {
        double[] bacteriaContents=new double[3];
        bacteriaContents=divideBacteriaForHost();
        StringBuilder sb= new StringBuilder();
        for (int i=0;i<3;i++){
            sb.append(bacteriaContents[i]).append("\t");
        }
        return sb.toString().trim();
    }
    public String printBacteriaContentsInEnv() {
        double[] bacteriaContents=new double[3];
        bacteriaContents=divideBacteria();
        StringBuilder sb= new StringBuilder();
        for (int i=0;i<3;i++){
            sb.append(bacteriaContents[i]).append("\t");
        }
        return sb.toString().trim();
    }
    @Override
    public SelectiveIndividualE[] getIndividuals(){
        return compositionOfIndividuals;
    }
    
    public SelectiveSpeciesRegistry getSpeciesRegistry(){
        return selSpReg;
    }
    
}
    
    
