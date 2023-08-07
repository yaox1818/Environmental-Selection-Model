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

import utils.random.Multinomial2;

/**
 *
 * @author John
 */
public class Individual  {

	double[] microbiome;
	private int numberEnvironmentalSpecies;
	
	
	public Individual(Multinomial2 MultDist, int nomph, int noes) {//nomph为个体体内微生物数量，noes为环境中微生物种类数量
		numberEnvironmentalSpecies = noes;
		microbiome = new double[numberEnvironmentalSpecies];
                MultDist.multisample(microbiome, nomph);//作用是按各比例生成微生物数量，并保存到microbiome中
		numberMicrobePerHost = nomph;
	}

	public String microbial_sequences() {
		char[] microbiome_sequence = new char[numberEnvironmentalSpecies];
		for (int i = 0; i < numberEnvironmentalSpecies; i++) {
			if (microbiome[i] >=0.5) {
				microbiome_sequence[i] = '1';
			} else {
				microbiome_sequence[i] = '0';
			}
		}
		return new String(microbiome_sequence);
	}


	public double[] getMicrobiome() {
		return microbiome;
	}
        
	public String printOut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 1; i < numberEnvironmentalSpecies; i++) {
			sb.append(microbiome[i]).append("\t");
		}
		return sb.toString().trim();
	}

	@Deprecated
	private int numberMicrobePerHost;
	
}


