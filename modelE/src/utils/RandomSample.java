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

package utils;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import utils.random.MathUtil;

/**
 *
 * @author qz28
 */



//根据性状列表随机抽取不重复性状给OTU
public class RandomSample {
	//NOTE: What sampling scheme/algorithm is this? what kind of randomness does it provide
	public static <T> Set<T> randomSample(List<T> items, int m) {//geneIndex, (int)numOfGenesPerMicrobe，所有性状列表，每种微生物的性状数
//		Random rnd = new Random();
		HashSet<T> res = new HashSet<>(m);
		int n = items.size();
		for (int i = n - m; i < n; i++) {//抽取5个性状放到res中，性状不重复
//			int pos = rnd.nextInt(i + 1);
			int pos = MathUtil.getNextInt(i);//从0到i之间抽取一个整数
			T item = items.get(pos);//将随机抽取的性状储存到item
			if (res.contains(item))//如果抽到重复的性状，添加i性状
				res.add(items.get(i));
			else
				res.add(item);
		}
		return res;
	}

}
