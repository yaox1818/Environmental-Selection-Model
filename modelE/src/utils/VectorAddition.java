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
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package utils;

import java.util.stream.IntStream;

/**
 *
 * @author John
 */
public class VectorAddition {
	//additionOfVectors是将
	public static void additionOfVectors(double[] sumResults, double weight1,
			double weight2, double[] v1, double[] v2) {
		//microbiomeSum, 1, 1, microbiomeSum, host.getMicrobiome()
		int dim1 = v1.length;
		int dim2 = v2.length;
		if (dim1 >= dim2) {
			for (int i = 0; i < dim2; i++) {
				sumResults[i] = weight1 * v1[i] + weight2 * v2[i];
			}
			for (int i = dim2; i < dim1; i++) {
				sumResults[i] = weight1 * v1[i];
			}
		} else {
			for (int i = 0; i < dim1; i++) {
				sumResults[i] = weight1 * v1[i] + weight2 * v2[i];
			}
			for (int i = dim1; i < dim2; i++) {
				sumResults[i] = weight2 * v2[i];
			}
		}
	}
	public static double[] additionOfVectorsMap(double[] sumResults, double weight1,
			double weight2, double[] v1, double[] v2) {

		sumResults = IntStream.range(0, sumResults.length).mapToDouble( n ->  (weight1 *  v1[n] + weight2 * v2[n]) ).toArray();
	
//		double[] wv1 = Arrays.stream(v1).parallel().map(n -> n * weight1).toArray();
//		double[] wv2 = Arrays.stream(v2).parallel().map(n -> n * weight2).toArray();
//		
//		double[] wv1 = Arrays.stream(v1).map(n -> n * weight1).toArray();
//		double[] wv2 = Arrays.stream(v2).map(n -> n * weight2).toArray();
		
//		Arrays.stream(v1).forEach( n -> n=n*weight1);
//		Arrays.stream(v2).forEach( n -> n=n*weight2);
		
//		double[] wv1 = IntStream.range(0, v1.length).mapToDouble(n -> v1[n] * weight1).toArray();
//		double[] wv2 = IntStream.range(0, v2.length).mapToDouble(n -> v2[n] * weight2).toArray();
//		DoubleStream map = Arrays.stream(v1).map(n -> n * weight1);
		
//		Arrays.setAll(sumResults, i -> wv1[i] + wv2[i]);
//		Arrays.parallelSetAll(sumResults, i -> wv1[i] + wv2[i]);
//		IntStream.range(0, sumResults.length).forEach( id -> sumResults[id] = weight1 *  v1[id] + weight2 * v2[id]); 

//		IntStream.range(0, sumResults.length).forEach( id -> sumResults[id] = wv1[id] + wv2[id]);
//		sumResults = IntStream.range(0, sumResults.length).mapToDouble( id ->  (weight1 *  v1[id] + weight2 * v2[id]) ).toArray();
		
//		sumResults = IntStream.range(0, sumResults.length).mapToDouble( id ->  (weight1 *  v1[id] + weight2 * v2[id]) ).toArray();
//		sumResults = IntStream.range(0, sumResults.length).mapToDouble( id ->  v1[id] + v2[id] ).toArray();
//		for (int i = 0; i < v1.length; i++) {
//			sumResults[i] = weight1 * v1[i] + weight2 * v2[i];
//		}
		return sumResults;
	}
	public static void additionOfVectorsFL(double[] sumResults, double weight1,
			double weight2, double[] v1, double[] v2) {

		
		for (int i = 0; i < v1.length; i++) {
			sumResults[i] = weight1 * v1[i] + weight2 * v2[i];
		}


	}

	public static void additionOfVectorsSum(double[] sumResults, double[] v1, double[] v2) {
//		sumResults = IntStream.range(0, sumResults.length).mapToDouble( id -> v1[id] + v2[id]).toArray();
//		System.out.println(Arrays.toString(sumResults));
		for (int i = 0; i < v1.length; i++) {
			sumResults[i] = v1[i] + v2[i];
		}

	}
}
