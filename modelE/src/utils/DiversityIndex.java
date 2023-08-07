/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;

/**
 *
 * @author qz28
 */
public class DiversityIndex {
    
    public static double speciesRichness(double[] microbiome){
        double x=0;
        for(double abundance:microbiome){
            if (abundance>0)
                    x++;
        }
        return x;
    }
   
    
    public static double shannonWienerIndex(double[] microbiome, double total){
        double ind=0;
        for(double abundance:microbiome){
            if (abundance>0&&abundance<total)
                ind-=(abundance/total)*Math.log(abundance/total);
        }
        return ind;
    }
    
    public static double shannonWienerIndex(double[] microbiome){
        return shannonWienerIndex(microbiome,1);
    }
    
    public static double shannonWienerevenness(double[] microbiome, double total){
        return shannonWienerIndex(microbiome, total)/speciesRichness(microbiome);
    }
    
    public static double shannonWienerevenness(double[] microbiome){
        return shannonWienerevenness(microbiome,1);
    }
   
    public static double TrueDiversity(double[] microbiome, double q ,double total){
        if (q==0){
            return speciesRichness(microbiome);
        }
        else if (q==1){
            return Math.exp(shannonWienerIndex(microbiome,total));
        }
        else{
            double ind=0;
            for(double abundance:microbiome){
                ind+=Math.pow(abundance/total,q);
            }
            return Math.pow(ind, 1/(1-q));
        }
        
    }
    
    public static double TrueDiversity(double[] microbiome, double q){
        return TrueDiversity(microbiome,q,1);
    }
    
    public static double simpsonIndex(double[] microbiome, double total){
        double lamda=0;
        for (double abundance : microbiome){
            lamda+=abundance/total*abundance/total;
        }
        return 1-lamda;
    }
    
    public static double simpsonIndex(double[] microbiome){
        return simpsonIndex(microbiome,1);
    }
    
    public static double BrayCurtis(double[] microbiome1, double[] microbiome2,double total1, double total2){
        double k=0;
        for (int i=0;i<microbiome1.length;i++){
            k+=Math.abs(microbiome1[i]/total1-microbiome2[i]/total2);
        }
        return k/2;
    }
    
    public static double BrayCurtis(double[] microbiome1, double[] microbiome2){
        double k=0;
            for (int i=0;i<microbiome1.length;i++){
                k+=Math.abs(microbiome1[i]-microbiome2[i]);
            }
        return k/2;
    }
    
    public static double PiIndex(double[] microbiome1, double[] microbiome2){
        double k=0;
        for (int i=0;i<microbiome1.length;i++){
            if((microbiome1[i]>0&&microbiome2[i]<=0) || (microbiome2[i]>0&&microbiome1[i]<=0)){
                k++;
            }
        }
        return k/microbiome1.length;
    }
    
    public static double JaccardIndex(double[] microbiome1, double[] microbiome2){
        double k=0;
        double n=0;
        for (int i=0;i<microbiome1.length;i++){
            if(microbiome1[i]>0 || microbiome2[i]>0){
                n++;
                if(microbiome1[i]>0 && microbiome2[i]>0){
                    k++;
                }
            }
        }
        return 1-k/n;
    }
    
    public static double DiceIndex(double[] microbiome1, double[] microbiome2){
        double x=1-JaccardIndex(microbiome1,microbiome2);
        return 1-2*x/(1+x);
    }
    
    
    
}
