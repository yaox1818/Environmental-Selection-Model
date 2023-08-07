/*******************************************************************************
 *
 * Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/


package microbiosima;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import utils.Distance;
import utils.DiversityIndex;
import utils.RandomSample;
import utils.VectorAddition;
import utils.random.MathUtil;
import utils.random.Multinomial2;

/**
 *
 * @author John
 */
public class Population {
	
	final int numberOfMicrobePerHost;//每个宿主体内的微生物数量
	final int numberOfEnvironmentalSpecies;//环境中微生物的OTU数量
	final int numberOfIndividual;
	final double environmentalFactor;
	final double percentageOfpooledOrFixed;
	
	double coefficient;
    double[][] parentalContribution; //
	double[] environmentalContribution;
	double[] mixedContribution;
	double[] initialEnvironment;
	double[] microbiomeSum;
	private double[] previous_microbiomes;
	int[] parentalIndex;
    int[] ancestryIndex;
	
	
	private double g_diversity;
	private double a_diversity;
	private double b_diversity;
	private double beta_diversity_coef;
    private double beta_diversity_coef_2;
    private double beta_diversity_coef_3;
    private double beta_diversity_coef_4;	
    private double weighted_b_diversity;
	
	private int numberOfGeneration;
	private int sampleReplicates;
	private int numberOfSamples;	

	private Individual[] compositionOfIndividuals;//Individual根据微生物来源随机生成微生物群落

	private List<Set<Integer>> samples = new ArrayList<>();
	
	List<Integer> host_index = new ArrayList<>();//ArrayList.add()参数只能是对象，所以用Integer而不是Int
	
	Multinomial2 multiNomialDist;

	public Population(int numberOfMicrobePerHost, double[] environment, int noOfIndividual, 
			double environmentalFactor,
			double pooledOrFiexd, int numberOfSamples, int sampleReplicates0) throws IOException {
		//宿主体内微生物数量、环境中各otu占比、宿主种群数量、环境占宿主微生物来源比例、亲代对环境的贡献
		this.numberOfMicrobePerHost = numberOfMicrobePerHost;//宿主体内微生物数量
		this.numberOfEnvironmentalSpecies = environment.length;//环境中微生物的种类
		this.numberOfIndividual = noOfIndividual;//宿主种群数量
		this.environmentalFactor = environmentalFactor;//环境占宿主微生物来源比例
		this.percentageOfpooledOrFixed = pooledOrFiexd;//NOTE: "or" sounds like a boolean variable
		//亲代对环境的贡献
		this.numberOfGeneration = 0;//代数
		this.numberOfSamples = numberOfSamples;//默认是0
		this.sampleReplicates = sampleReplicates0;//默认是0
		coefficient = (1 - environmentalFactor) / this.numberOfMicrobePerHost; //NOTE: change to this.系数=亲代占宿主微生物来源比例/宿主体内微生物总数
		initialEnvironment = environment;//初始环境
		
		compositionOfIndividuals = new Individual[numberOfIndividual];//初始化整个宿主种群
		parentalContribution = new double[noOfIndividual][numberOfEnvironmentalSpecies];
		//NOTE: change to numberOfIndividual 【宿主个数】【微生物种类】
		microbiomeSum = new double[numberOfEnvironmentalSpecies];
		environmentalContribution = new double[numberOfEnvironmentalSpecies];
		mixedContribution = new double[numberOfEnvironmentalSpecies];//微生物来源
		previous_microbiomes = new double[numberOfEnvironmentalSpecies];
                parentalIndex=new int[numberOfIndividual];
                ancestryIndex=new int[numberOfIndividual];
		
		multiNomialDist = new Multinomial2(numberOfEnvironmentalSpecies);
		
		beta_diversity_coef = 2.0  / numberOfIndividual
				/ (numberOfIndividual - 1)/numberOfMicrobePerHost;
        beta_diversity_coef_2 = 2.0  / numberOfIndividual
				/ (numberOfIndividual - 1);
        beta_diversity_coef_3 = 2.0  / numberOfSamples
				/ (numberOfSamples - 1)
				/ sampleReplicates;
        beta_diversity_coef_4 = 2.0  / numberOfSamples
				/ (numberOfSamples - 1)
				/ sampleReplicates/numberOfMicrobePerHost;
		//形成一个宿主列表和一个祖先数组，长度为宿主种群数量
		for (int i = 0; i < numberOfIndividual; i++) {
			host_index.add(i);
                        ancestryIndex[i]=i;
		}
		//形成一个列表？
		for (int i = 0; i < sampleReplicates; i++) {
			samples.add(RandomSample.randomSample(host_index,
					this.numberOfSamples));
		}

                multiNomialDist.updateProb(initialEnvironment);

		//按宿主的种群数量初始化一系列宿主，宿主体内的微生物是根据微生物来源中微生物占比为概率来随机生成的
		for (int i = 0; i < noOfIndividual; i++) {
			//multiNomialDist.updateProb(initialEnvironment);
			compositionOfIndividuals[i] = new Individual(
					multiNomialDist, this.numberOfMicrobePerHost,
					numberOfEnvironmentalSpecies);

		}
	}
    //microbiomeSum是什么
	public void sumSpecies() {//一代宿主中的所有微生物
		Arrays.fill(microbiomeSum, 0);//Arrays.fill是一个数组填充的方法，将microbiomeSum数组填满0.
		for (Individual host : getIndividuals()) {
			VectorAddition.additionOfVectors(microbiomeSum, 1, 1,
					microbiomeSum, host.getMicrobiome());
		}//VectorAddition.additionOfVectors是累计各OTU微生物数量的方法
		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
			microbiomeSum[i] = microbiomeSum[i] / numberOfIndividual
					/ numberOfMicrobePerHost;
		}//microbiomeSum为每种OTU占所有宿主体内微生物总和的比例
	}

	public double interGenerationDistance() {
		if (getNumberOfGeneration() == 0)
			return 0;
		else {
			return Distance.getDistance(previous_microbiomes, microbiomeSum);
		}
	}
    //environmentPopulationDistance
	public double environmentPopulationDistance() {
		if (getNumberOfGeneration() == 0)
			//TODO: can we merge these two??
			return Distance.getDistance(initialEnvironment, microbiomeSum);
		else
			return Distance.getDistance(environmentalContribution,
					microbiomeSum);
	}



	public void parentalInheritanceAndEnvironmentalAcquisition() {

        int[] oldAncestryIndex= new int[numberOfIndividual];
        System.arraycopy(ancestryIndex,0,oldAncestryIndex, 0,numberOfIndividual);//将上一代的ancestryIndex保存到oldAncestryIndex中
		for (int i = 0; i < numberOfIndividual; i++) {
            parentalIndex[i]=MathUtil.getNextInt(numberOfIndividual-1);//将0到149之间随机生成的整数放到parentalIndex[i]中
            ancestryIndex[i]=oldAncestryIndex[parentalIndex[i]];//通过parentalIndex作为索引从oldAncestryIndex获取数到ancestryIndex中
			System.arraycopy(compositionOfIndividuals[parentalIndex[i]].getMicrobiome(),//通过抽取到的亲代宿主按顺序放到parentalContribution中构成亲本来源
                                0,parentalContribution[i], 0,numberOfEnvironmentalSpecies);
		}//每次随机抽取一个宿主作为亲代，抽取宿主种群个数次
		//NOTE: Maybe we can merge some of these environmentalContribution, mixedContribution, parentalContribution
		//These fields are not used at other places
		//根据亲代贡献和环境自身来构建环境来源，microbiomeSum为所有宿主的OTU微生物占比，initialEnvironment为初始环境OTU占比
		VectorAddition.additionOfVectors(environmentalContribution,
				percentageOfpooledOrFixed, 1 - percentageOfpooledOrFixed, microbiomeSum,
				initialEnvironment);
		//根据环境贡献和亲本贡献来确定每一个宿主的微生物来源，parentalContribution[i]为此宿主的亲代（每个宿主亲代可能不同），environmentalContribution为环境来源（每个宿主都一样）
		for (int i = 0; i < numberOfIndividual; i++) {
			VectorAddition.additionOfVectors(mixedContribution, coefficient,
					environmentalFactor, parentalContribution[i],
					environmentalContribution);
			//根据宿主的微生物环境随机生成宿主体内微生物数量
			multiNomialDist.updateProb(mixedContribution);
			multiNomialDist.multisample(
					compositionOfIndividuals[i].getMicrobiome(),
					numberOfMicrobePerHost);
		}
	}
	
	public void getNextGen() {
		parentalInheritanceAndEnvironmentalAcquisition();
		for (int i = 0; i < sampleReplicates; i++) {
			samples.set(i,
					RandomSample.randomSample(host_index, numberOfSamples));
		}
		
		System.arraycopy(microbiomeSum, 0, previous_microbiomes, 0, microbiomeSum.length);

		numberOfGeneration++;
	}
        
    public double[][] sample(int i){
            Set<Integer> sampleId=samples.get(i);
            double[][] sampleMicrobiomes=new double[numberOfSamples][];
            int m=0;
            for (Integer index:sampleId){
                 sampleMicrobiomes[m]=getIndividuals()[index].getMicrobiome();
                 m++;
            }
            return sampleMicrobiomes;
    }
        
    public double[] pooledSample(int i){
            Set<Integer> sampleId=samples.get(i);
            double[] temp_sum = new double[numberOfEnvironmentalSpecies];
            for (Integer index : sampleId) {
		VectorAddition.additionOfVectors(temp_sum, 1, 1.0/numberOfSamples,
							temp_sum,
							getIndividuals()[index].getMicrobiome());
				}
            return temp_sum;
    }

	public double alphaDiversity(boolean sampleOrNot) {
                a_diversity = 0;
		if (sampleOrNot) {
			for (Individual host : getIndividuals()) {
				a_diversity+=DiversityIndex.shannonWienerIndex(host.getMicrobiome(),numberOfMicrobePerHost);
			}
			a_diversity /=  numberOfIndividual;
		}
                else {
			for (Set<Integer> sample : samples) {
				for (Integer index : sample) {
					a_diversity+=DiversityIndex.shannonWienerIndex(getIndividuals()[index].getMicrobiome(),numberOfMicrobePerHost);
				}
			}
			a_diversity /= (numberOfSamples * sampleReplicates);
		}
		return a_diversity;

	}

	
	public double betaDiversity(boolean sampleOrNot) {
                b_diversity=0;
		if (sampleOrNot) {
                        if (numberOfIndividual>1){
			for (int i = 1; i < numberOfIndividual; i++) {
				for (int j = 0; j < i; j++) {
					b_diversity+=DiversityIndex.PiIndex(getIndividuals()[i].getMicrobiome(), getIndividuals()[j].getMicrobiome());
				}
			}
			
			b_diversity *= beta_diversity_coef_2;
		} }
                else {
                        if (numberOfSamples>1){
                        for (int index=0;index<sampleReplicates;index++){
                            double[][] temp_microbiomes=sample(index);
                            for (int i = 1; i < numberOfSamples; i++) {
                                for (int j = 0; j < i; j++) {
                                    b_diversity += DiversityIndex.PiIndex(temp_microbiomes[i], temp_microbiomes[j]);
                                }
                            }
                        }
			b_diversity *= beta_diversity_coef_3; 
		}}
                return b_diversity;
	}

	public double gammaDiversity(boolean sampleOrNot) {
                g_diversity=0;
		if (sampleOrNot) {
			g_diversity = DiversityIndex.shannonWienerIndex(microbiomeSum);
		} else {
                        for(int index=0;index<sampleReplicates;index++){
                            double[] temp_sum=pooledSample(index);
                            g_diversity+=DiversityIndex.shannonWienerIndex(temp_sum);
                        }
			g_diversity = g_diversity / sampleReplicates;	
		}
                return g_diversity;
	}
        
    public double BrayCurtis(boolean sampleOrNot){
            weighted_b_diversity=0;
            if (sampleOrNot){
                for (int i=1;i<numberOfIndividual;i++){
                    for (int j=0;j<i;j++){
                        weighted_b_diversity+= DiversityIndex.BrayCurtis(getIndividuals()[i].getMicrobiome(), getIndividuals()[j].getMicrobiome());
                    }
                }
                weighted_b_diversity*=beta_diversity_coef;
            }
            else{
                for (int index=0;index<sampleReplicates;index++){
                            double[][] temp_microbiomes=sample(index);
                            for (int i = 1; i < numberOfSamples; i++) {
                                for (int j = 0; j < i; j++) {
                                    weighted_b_diversity += DiversityIndex.BrayCurtis(temp_microbiomes[i], temp_microbiomes[j]);
                                }
                            }
                        }
			weighted_b_diversity *= beta_diversity_coef_4;
            }
            return weighted_b_diversity;
    }

	public String printOut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < numberOfEnvironmentalSpecies; i++) {
			sb.append(microbiomeSum[i]).append("\t");
		}
		return sb.toString().trim();
	}
    
    public Individual[] getIndividuals(){
        return compositionOfIndividuals;
    }

	public int getNumberOfGeneration() {
		return numberOfGeneration;
	}

	public void resetGeneration() {
		numberOfGeneration = 0;
		
	}
	
    public boolean checkAncestry() {
        for (int i=1; i<numberOfIndividual;i++){
            if (ancestryIndex[i-1]!=ancestryIndex[i])
                return false;
            }
            return true;
        }
        
    public String indexPrintOut(int[] indexArray) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < indexArray.length; i++) {
			sb.append(indexArray[i]).append("\t");
		}
		return sb.toString().trim();
	}	

        
     


}

