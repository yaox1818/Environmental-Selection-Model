/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import utils.RandomSample;


public class SelectiveSpeciesRegistry extends SpeciesRegistry {
    int numOfTolGenes;
    private double numOfGenesPerMicrobe;
    double mscih;
    double mscie;
	int[] hfr;
	int[] mfrH;
	int[] mfrE;
    private ArrayList<Integer> geneCompositionRecords=new ArrayList<>();


    //initialNumberOfSpecies：微生物总种类  notg：25种性状  nogpm：5种性状  microbeSelectionCoef：参数Sm  fitnessToHost：25种性状的分数  fitnessToMicrobe：25种性状的分数
    public SelectiveSpeciesRegistry(int initialNumberOfSpecies,int notg, double nogpm, double microbeSelectionCoefInHost,double microbeSelectionCoefInEnv,int[] fitnessToHost, int[] fitnessToMicrobeInHost, int[] fitnessToMicrobeInEnv) {
        super(initialNumberOfSpecies);//微生物总种类
        numOfTolGenes=notg;//总性状数量
        numOfGenesPerMicrobe=nogpm;//每种微生物的性状数
        mscih=microbeSelectionCoefInHost;//Sm参数
        mscie=microbeSelectionCoefInEnv;//Sm参数
		hfr=fitnessToHost;//宿主性状分数，长度25的数组
		mfrH=fitnessToMicrobeInHost;//微生物性状分数，长度25的数组
		mfrE=fitnessToMicrobeInEnv;//微生物性状分数，长度25的数组
        List<Integer> geneIndex = new ArrayList<>();//元素只能为int类型，性状列表
        for (int i=0;i<numOfTolGenes;i++){
            geneIndex.add(i);//性状列表
        }
        for (int i=0;i<initialNumberOfSpecies;i++){//给每种OTU随机赋予性状
            Set<Integer> tempGenes=RandomSample.randomSample(geneIndex, (int)numOfGenesPerMicrobe);//根据geneIndex给每种OTU抽取性状以二进制的形式放到tempGenes中
            int geneBinaryComp=0;//
            for (int gene : tempGenes){
                geneBinaryComp+=Math.pow(2,gene);
            }
            geneCompositionRecords.add(geneBinaryComp);//记录每种OTU
        }
        
    }        
    public double getTotalFitness(double[] microbiome, double[] fitnessRecords ){//microbiome为宿主体内微生物的分布，fitnessRecords为宿主对各种OTU的fitness
        double totalFitness=0;
        double abundance=0;
        for (int i=0;i<microbiome.length;i++){
            abundance+=microbiome[i];
            totalFitness+=microbiome[i]*fitnessRecords[i];
        }
        return totalFitness/abundance;
    }

    /*public void variation(double probability){
        double p = MathUtil.getNextFloat(1);
        if (p<probability)

    }**/

    public void getFitness(double[] Fitness, int[] geneFitness){
        for(int i=0;i<Fitness.length;i++){
            int genotype=geneCompositionRecords.get(i);
            double fitness=0;
            int index=0;
                while(genotype>0){
                    if(genotype%2==1)
                        fitness+=geneFitness[index];
                    index++;
                    genotype=genotype>>1;
                }
            Fitness[i]=fitness/numOfGenesPerMicrobe;
        }
    }


    //根据geneCompositionRecords中所记录的trait计算每种OTU的fitness
    //microbeFitnessRecords为所有微生物OTU的mf，hostFitnessRecords为宿主对所有OTU的fitness，microbeGenesFitness为微生物对各种trait的fitness，hostGenesFitness为宿主对各种trait的fitness
    public void getFitness(double[] Fitness1,double[]Fitness2, int[] geneFitness1, int[] geneFitness2){
        for(int i=0;i<Fitness1.length;i++){
            int genotype=geneCompositionRecords.get(i);
            double fitness1=0;
            double fitness2=0;
            int index=0;
                while(genotype>0){
                    if(genotype%2==1){
                        fitness1+=geneFitness1[index];
                        fitness2+=geneFitness2[index];}
                    index++;
                    genotype=genotype>>1;
                }
            Fitness1[i]=fitness1/numOfGenesPerMicrobe;//根据每种otu的trait计算每种OTU的总fitness
            Fitness2[i]=fitness2/numOfGenesPerMicrobe;
        }
    }
    
    public void getMFitnessSelectionInHost(double[] microbiome, double[] fitnessToMicrobeRecords){
        double totalFitness=0;
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]=microbiome[i]*Math.pow(mscih,fitnessToMicrobeRecords[i]);
            totalFitness+=microbiome[i];
        }
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]/=totalFitness;
        }
    }

    public void getMFitnessSelectionInEnv(double[] microbiome,double[] fitnessToMicrobeRecords){
        double totalFitness=0;
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]=microbiome[i]*Math.pow(mscie,fitnessToMicrobeRecords[i]);
            totalFitness+=microbiome[i];
        }
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]/=totalFitness;
        }
    }
    
    public ArrayList<Integer> getGeneComposition(){
        return geneCompositionRecords;
    }
    
    
    
}
