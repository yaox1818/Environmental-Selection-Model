/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utils;
import java.util.Arrays;
import utils.random.MathUtil;

/**
 *
 * @author qz28
 */
//根据两项之间的差来随机抽取
public class BinarySearch {
    public static int binarySearch(double[] toBeSearched, double obj){
        //double obj=MathUtil.getNextFloat(toBeSearched[toBeSearched.length-1]);
        if (toBeSearched.length==1)
            return 0;
        else{
            int k=(toBeSearched.length-1)/2;
            if (obj<=toBeSearched[k]){
                return binarySearch(Arrays.copyOfRange(toBeSearched, 0, k+1),obj);
            }
            else{
                return k+1+binarySearch(Arrays.copyOfRange(toBeSearched, k+1, toBeSearched.length),obj);
            }
        }
    }
    
}
