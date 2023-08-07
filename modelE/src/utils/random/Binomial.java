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

/**
 *
 * @author John
 */
public class Binomial {
	public static double f_function(int x) {
		if (x == 0)
			return 0.08106146679532726;
		else if (x == 1)
			return 0.04134069595540929;
		else if (x == 2)
			return 0.02767792568499834;
		else if (x == 3)
			return 0.02079067210376509;
		else if (x == 4)
			return 0.01664469118982119;
		else if (x == 5)
			return 0.01387612882307075;
		else if (x == 6)
			return 0.01189670994589177;
		else if (x == 7)
			return 0.01041126526197209;
		else if (x == 8)
			return 0.009255462182712733;
		else if (x == 9)
			return 0.008330563433362871;
		else if (x >= 10) {
			long y = (x + 1) * (x + 1);
			return (1.0 / 12 - (1.0 / 360 - 1.0 / 1260 / y) / y) / (x + 1);
		}
		return 0;
	}

	public static int BTRD(int N, double p) {
		int m = (int) ((N + 1) * p);
		double r = p / (1 - p);
		double nr = (N + 1) * r;
		double npq = N * p * (1 - p);
		double npq2 = Math.sqrt(npq);
		double b = 1.15 + 2.53 * npq2;
		double a = -0.0873 + 0.0248 * b + 0.01 * p;
		double c = N * p + 0.5;
		double alpha = (2.83 + 5.1 / b) * npq2;
		double vr = 0.92 - 4.2 / b;
		double ur = 0.86 * vr;
//		Random generator = new Random();
		double v;
		double u;
		double f;
		int k = 0;
		double us;
		double km;
		double rou;
		double t;
		do {
			do {
				do {
					v = MathUtil.nextFloat();
					if (v <= ur) {
						u = v / vr - 0.43;
						return (int) ((2 * a / (0.5 - Math.abs(u)) + b) * u + c);
					}
					if (v >= vr) {
						u = MathUtil.nextFloat() - 0.5;
					} else {
						u = v / vr - 0.93;
						u = Math.signum(u) * 0.5 - u;
						v = MathUtil.nextFloat() * vr;
					}
					us = 0.5 - Math.abs(u);
					k = (int) ((2 * a / us + b) * u + c);
				} while (k < 0 || k > N);
				v = v * alpha / (a / (us * us) + b);
				km = Math.abs(k - m);
				if (km > 15)
					break;
				f = 1;
				if (m < k) {
					int i = m;
					while (i < k) {
						i = i + 1;
						f = f * (nr / i - r);
					}
				} else if (m > k) {
					int i = k;
					while (i < m) {
						i = i + 1;
						v = v * (nr / i - r);
					}
				}
				if (v <= f)
					return k;
			} while (v > f);
			v = Math.log(v);
			rou = (km / npq) * (((km / 3 + 0.625) * km + 1.0 / 6) / npq + 0.5);
			t = -km * km / (2 * npq);
			if (v < t - rou)
				return k;
		} while (v > t + rou);
		int nm = N - m + 1;
		double h = (m + 0.5) * Math.log((m + 1) / (r * nm)) + f_function(m)
				+ f_function(N - m);
		double nk = N - k + 1;
		if (v <= h + (N + 1) * Math.log(nm / nk) + (k + 0.5)
				* Math.log(nk * r / (k + 1)) - f_function(k)
				- f_function(N - k))
			return k;
		else {
			do {
				do {
					do {
						do {
							v = MathUtil.nextFloat();
							if (v <= ur) {
								u = v / vr - 0.43;
								return (int) ((2 * a / (0.5 - Math.abs(u)) + b)
										* u + c);
							}
							if (v >= vr) {
								u = MathUtil.nextFloat() - 0.5;
							} else {
								u = v / vr - 0.93;
								u = Math.signum(u) * 0.5 - u;
								v = MathUtil.nextFloat() * vr;
							}
							us = 0.5 - Math.abs(u);
							k = (int) ((2 * a / us + b) * u + c);
						} while (k < 0 || k > N);
						v = v * alpha / (a / (us * us) + b);
						km = Math.abs(k - m);
						if (km > 15)
							break;
						f = 1;
						if (m < k) {
							int i = m;
							while (i < k) {
								i = i + 1;
								f = f * (nr / i - r);
							}
						} else if (m > k) {
							int i = k;
							while (i < m) {
								i = i + 1;
								v = v * (nr / i - r);
							}
						}
						if (v <= f)
							return k;
					} while (v > f);
					v = Math.log(v);
					rou = (km / npq)
							* (((km / 3 + 0.625) * km + 1.0 / 6) / npq + 0.5);
					t = -km * km / (2 * npq);
					if (v < t - rou)
						return k;
				} while (v > t + rou);
				nk = N - k + 1;
			} while (v > h + (N + 1) * Math.log(nm / nk) + (k + 0.5)
					* Math.log(nk * r / (k + 1)) - f_function(k)
					- f_function(N - k));
			return k;
		}

	}


	public static int inversion_cdf(int N, double p) {
//		Random generator = new Random();
		double key = MathUtil.nextFloat();//MathUtil.nextFloat生成一个0到1之间的随机数
		double q = 1 - p;
		double s = p / q;
		double a = (N + 1) * s;//（N+1）*p/（1-p）
		double r = Math.pow(q, N);//Math.pow求次方
		for (int i = 0; i <= N; i++) {
			if (key <= r)
				return i;
			key = key - r;
			r = (a / (i + 1) - s) * r;
		}
		return 0;

	}


	//N为宿主体内微生物数量，p为微生物的一个与占比有关的数组，r为生成的一个OTU的微生物个数
	public static int binomial_sampling(int N, double p) {
		double limit = 10;
		if (p <= 0)
			return 0;
		if (p >= 1)
			return N;
		int r;
		if (p <= 0.5) {
			if (N * p >= limit){
				r=BTRD(N, p);
			}
			else{
				r=inversion_cdf(N, p);
			}
		} else {
			if (N * (1 - p) >= limit){
				r= N - BTRD(N, 1 - p);
			}
			else{
				r= N - inversion_cdf(N, 1 - p);
			}
		}
		return r;
	}

}
