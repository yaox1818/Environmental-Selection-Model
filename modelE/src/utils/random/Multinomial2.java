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

package utils.random;

import java.util.Arrays;

//本类用于
public class Multinomial2 {

	private double[] distribution;
	private int rangeOne;

	public Multinomial2(int range) {
		this.rangeOne = range - 1;
		distribution = new double[range];

	}

    //形成一个数组？
	public void updateProb(double[] probabilities) {
		double cum_prob = 1;
		for (int i = 0; i < distribution.length; i++) {
			if (probabilities[i] == 0)
				distribution[i] = 0;
			else {
				distribution[i] = probabilities[i] / cum_prob;
				cum_prob = cum_prob - probabilities[i];
			}
		}
	}
    //sampledNumbers为宿主体内微生物群落，sampleSize为宿主体内微生物数量
	//作用是按各比例生成微生物数量，并保存到multisample中
	public double[] multisample(double[] sampledNumbers, int sampleSize) {
        //按占比随机生成各OTU的微生物数量
		for (int i = 0; i < rangeOne; i++) {
			int rand = Binomial.binomial_sampling(sampleSize, distribution[i]);//具体没看懂，但应该是根据概率生成某一OTU的数量
			sampledNumbers[i] = rand;
			sampleSize -= rand;

			if (sampleSize == 0) {
				Arrays.fill(sampledNumbers, i + 1, sampledNumbers.length, 0);
				break;
			}
		}
		sampledNumbers[rangeOne] = sampleSize;
		return sampledNumbers;
	}
}
