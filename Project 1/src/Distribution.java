// Joseph Shaffer
// Shaffer.567
// CSE 6461 Project 1




import java.util.Random;
import org.jfree.chart.*;
import org.jfree.chart.plot.PlotOrientation;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import java.lang.Math;
import java.io.*;

public class Distribution {
	
	public static int testSize = 1000;
	
	// Returns random variable between 0 and 1
	public static double URV() {
		Random randomGenerator = new Random();
		return randomGenerator.nextDouble();
	}
	
	// Array of Uniform Random Variables
	public static double[] randoms(int n, int j, double x) {
		double[] random = new double[n];
		for (int i = 0; i < n; i++) {
			switch(j) {
			case 0: random[i] = URV(); break;
			case 1: random[i] = exponentialRV(1 / x); break;
			case 2: random[i] = poissonRV(x); break;
			}
			
		}
		return random;
	}
	
	// Calculates P(U > x) where count  = 1 - x
	// P(X > x) where count = 
	public static double greaterThanX(double count, double range) {
		return (count / range);
	}
	
	public static double exponentialRV(double lambda) {
		double x;
		do {
			x = URV();
		} while ((x == 0) || (x == 1));
		return (- lambda * Math.log(x));
	} 
	
	public static double poissonRV(double lambda) {
		double L = Math.exp(-lambda);
		double p = 1.0;
		int k = 0;
		do {
			k++;
			p *= URV();
		} while (p > L);
		return k - 1;
	}
	
	public static double factorial(double k) {
		if (k <= 1) {
			return 1;
		} else {
			return k * factorial(k-1);
		}
	}
	
	public static double calculatePoissonProb(double lambda, double k) {
		if (k <= 1) {
			return (Math.pow(lambda, k) * Math.exp(-lambda) / factorial(k));
		} else {
			return ((Math.pow(lambda, k) * Math.exp(-lambda) / factorial(k)) + calculatePoissonProb(lambda, k - 1));
		}
	}
	
	public static double PnMM1(double lambda, double mu, int n) {
		double po = 1 - (lambda / mu);
		return Math.pow(lambda / mu, (double) n) * po;
	}
	
	public static double[] simMM1 (double lambda, double mu) {
		double[] timeInStateN = new double[100];
		double endTime = 1.0e6;
		double arrivalTime = 1 / lambda;
		double serviceTime = 1 / mu;
		double time = 0.0;
		double t1 = 0.0;
		double t2 = endTime;
		int n = 0;
		int completed = 0;
		double throughput = 0.0;
		double wait = 0.0;
		double areaOfCurve = 0.0;
		double tOfLast = time;
		double numSys;

		
		while (time < endTime) {
			if (t1 < t2) {
				time = t1;
				timeInStateN[n] += (time - tOfLast);
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n++;
				tOfLast = time;
				t1 = time + exponentialRV(arrivalTime);
				if (n == 1) {
					t2 = time + exponentialRV(serviceTime);
					
				}
			}
			else {
				time = t2;
				timeInStateN[n] += (time - tOfLast);
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n--;
				tOfLast = time;
				completed++;
				if (n > 0) {
					t2 = time + exponentialRV(serviceTime);
				} else {
					t2 = endTime;
				}
			}
		}
		throughput = completed / time;
		numSys = areaOfCurve / time;
		wait = numSys / throughput;
		System.out.println("Theoretical Expected Number in M/M/1 System E(N): " + expectedNumberMM1(lambda, mu) + " customers");
		System.out.println("Expected Number in M/M/1 System E(N): " + numSys + " customers");
		System.out.println("Theoretical Expected Delay in M/M/1 System E(T): " + expectedDelayMM1(lambda, mu) + " minutes");
		System.out.println("Expected Delay in M/M/1 System E(T): " + wait + " minutes");
		System.out.println("By Little's Law E(N) = lambda * E(T): " + numSys + " = " + lambda + " * " + wait + " = " + wait * lambda);
		double[] pn = new double[100];
		for (int i = 0; i < timeInStateN.length; i++) {
			pn[i] = (timeInStateN[i] / time);
		}
		return pn; 
	}
		
	
	
	public static double expectedNumberMM1(double lambda, double mu) {
		return /*(Math.pow(lambda, 2)) / (mu * (mu - lambda));*/((lambda / mu) / (1 - (lambda / mu)));
	}
	
	public static double expectedDelayMM1(double lambda, double mu) {
		return /*lambda / (mu * (mu - lambda)); */1 / (mu - lambda);
	} 
	
	public static double PnMEk1(double lambda, double mu, double k, int x) {
		if (lambda > 0) {
			return (Math.pow((1 / lambda), k) * Math.pow((double) x, k - 1) * Math.exp(-(1 / lambda) * x) / factorial(k-1));
		} else {
			return ((Math.pow((double) x, k - 1) * Math.exp(- x / mu)) / (Math.pow(mu, k) * factorial(k - 1)));
		}
	}
	
	public static double[] simMEk1(double lambda, double mu, double k) {
		double[] timeInStateN = new double[100];
		double endTime = 1e06;
		double arrivalTime = 1 / lambda;
		double serviceTime = 1 / mu;
		double time = 0.0;
		double t1 = 0.0;
		double t2 = endTime;
		int n = 0;
		int completed = 0;
		double throughput = 0.0;
		double wait = 0.0;
		double areaOfCurve = 0.0;
		double tOfLast = time;
		double numSys;

		
		while (time < endTime) {
			if (t1 < t2) {
				time = t1;
				timeInStateN[n] += (time - tOfLast);
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n++;
				tOfLast = time;
				t1 = time + exponentialRV(arrivalTime);
				if (n == 1) {
					t2 = time + erlangRV(serviceTime, k);
				}
			}
			else {
				time = t2;
				timeInStateN[n] += (time - tOfLast);
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n--;
				tOfLast = time;
				completed++;
				if (n > 0) {
					t2 = time + erlangRV(serviceTime, k);
				} else {
					t2 = endTime;
				}
			}
		}
		throughput = completed / time;
		numSys = areaOfCurve / time;
		wait = numSys / throughput;
		System.out.println("Theoretical Expected Number in M/Ek/1 System E(N): " + expectedNumberMEk1(k, lambda / mu) + " customers");
		System.out.println("Expected Number in M/Ek/1 System E(N): " + numSys + " customers");
		System.out.println("Theoretical Expected Delay in M/Ek/1 System E(T): " + expectedDelayMEk1(k, lambda, mu) + " minutes");
		System.out.println("Expected Delay in M/Ek/1 System E(T): " + wait + " minutes");
		System.out.println("By Little's Law E(N) = lambda * E(T): " + numSys + " = " + lambda + " * " + wait + " = " + wait * lambda);
		double[] pn = new double[100];
		for (int i = 0; i < timeInStateN.length; i++) {
			pn[i] = (timeInStateN[i] / time);
		}
		return pn;
	}
	
	public static double enMEk1(double lambda, double mu, double k) {
		double endTime = 1000;
		double arrivalTime = 1 / lambda;
		double serviceTime = 1 / mu;
		double time = 0.0;
		double t1 = 0.0;
		double t2 = endTime;
		int n = 0;
		double areaOfCurve = 0.0;
		double tOfLast = time;

		while (time < endTime) {
			if (t1 < t2) {
				time = t1;
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n++;
				tOfLast = time;
				t1 = time + exponentialRV(arrivalTime);
				if (n == 1) {
					t2 = time + erlangRV(serviceTime, k);
				}
			}
			else {
				time = t2;
				areaOfCurve = areaOfCurve + n * (time - tOfLast);
				n --;
				tOfLast = time;
				if (n > 0) {
					t2 = time + erlangRV(serviceTime, k);
				} else {
					t2 = endTime;
				}
			}
		}
		return areaOfCurve/time;
	} 
	
	

	public static double expectedNumberMEk1(double k, double rho) {
		return (rho / (1 - rho)) * (1- ((rho / 2) * (1 - (1 / k))));
	}
	
	public static double expectedDelayMEk1(double k, double lambda, double mu) {
		double delay = ((k + 1) / (2 * k)) * ((lambda) / (mu * (mu - lambda)));
		return delay;
	}
	
	
	public static double expectedNumberMD1(double rho) {
		return (rho + 0.5 * (Math.pow(rho,  2.0) / (1.0 - rho)));
	}
	
	public static double erlangRV(double lambda, double k) {
		
		double sum = 0.0;
		for (int i = 1; i <= k; i++) {
			sum += exponentialRV(lambda / k);
		}
		return sum;
	}
	
	public static void main(String[] args) throws IOException {
		
		// Plot uniform random variable P(U > x) with calculated probability
		final XYSeries urvGreaterThanX = new XYSeries("Theoretical P(U > x)");
		for (double i = 0.5; i <= 1.0; i += 0.05) {
			
			urvGreaterThanX.add(i, greaterThanX(1-i, 1));
		}
		final XYSeries estimateURVGreaterThanX = new XYSeries("Calculated P(U > x)");
		double[] random0 = randoms(testSize, 0, 0);
		for (double i = 0.5; i <= 1.0; i = i + 0.05) {
			int count = 0;
			for (int j = 0; j < random0.length; j++) {
				if (random0[j] > i) {
					count++;
				}
			}
			estimateURVGreaterThanX.add(i, greaterThanX(count, testSize));
		}
		
		//System.out.println("Array of Uniform Random Variables: " + Arrays.toString(random0));
		
		/* Plot a sample size of random variables and plot the 
		 * results of U > x from the sample size of Uniform
		 * random variables
		 */
		final XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(urvGreaterThanX);
		dataset.addSeries(estimateURVGreaterThanX);
		JFreeChart urv = ChartFactory.createXYLineChart("P(U > x) vs x", "x", "P(U > x)", dataset, PlotOrientation.VERTICAL, true, true, false);
		int width = 640;
		int height = 480;
		File URVChart = new File("URVChart.jpeg");
		ChartUtilities.saveChartAsJPEG(URVChart, urv, width, height);
		
		
		// Plot calculated probability of exponential random variables
		final XYSeries exprvGreaterThanX = new XYSeries("Theoretical P(X > x)");
		for (double i = 0.0; i <= 1.0; i += 0.1) {
			exprvGreaterThanX.add(i, Math.log(Math.exp(-2 * i)));
		}
		
		final XYSeries estimateExpRVGreaterThanX = new XYSeries("Calculated P(X > x)");
		double[] random1 = randoms(testSize, 1, 2);
		for (double i = 0.0; i <= 1.0; i += 0.1) {
			int count = 0;
			for (int j = 0; j < random1.length; j++) {
				if (random1[j] > i) {
					count++;
				}
			}
			estimateExpRVGreaterThanX.add(i, Math.log(greaterThanX(count, testSize)));
		}
		
		//System.out.println("Array of Exponential Random Variables: " +  Arrays.toString(random1));
		
		/* Sample size of exponential random variables and plot sample
		 * size where random variables X are greater than x
		 */
		final XYSeriesCollection dataset1 = new XYSeriesCollection();
		dataset1.addSeries(exprvGreaterThanX);
		dataset1.addSeries(estimateExpRVGreaterThanX);
		JFreeChart expRV = ChartFactory.createXYLineChart("log(P(X > x)) vs x", "x", "log(P(X > x))", dataset1, PlotOrientation.VERTICAL, true, true, false);
		File EXPRVChart = new File("ExpRVChart.jpeg");
		ChartUtilities.saveChartAsJPEG(EXPRVChart, expRV, width, height);
		
		
		// Plot of calculated poisson probability 
		final XYSeries poissonGreaterThanX = new XYSeries("Theoretical P(Y > x)");
		for (int i = 0; i <= 10; i++) {
			poissonGreaterThanX.add(i, 1 - calculatePoissonProb(2.0, i));
		}
		
		/* Sample size of random poisson variables 
		 * plot sample size where Y > x
		 */
		final XYSeries estimatePoissonGreaterThanX = new XYSeries("Calculated P(Y > x)");
		double[] random2 = randoms(testSize, 2, 2);
		for (int i = 0; i <= 10; i++) {
			int count = 0;
			for (int j =0; j < random2.length; j++) {
				if (random2[j] > i) {
					count++;
				}
			}
			estimatePoissonGreaterThanX.add(i, greaterThanX(count, testSize));
		}
		
		//System.out.println("Array of Poisson Random Variables: " + Arrays.toString(random2));
		
		final XYSeriesCollection dataset2 = new XYSeriesCollection();
		dataset2.addSeries(poissonGreaterThanX);
		dataset2.addSeries(estimatePoissonGreaterThanX);
		JFreeChart poissonRV = ChartFactory.createXYLineChart("P(Y > x) vs x", "x", "P(Y > x)", dataset2, PlotOrientation.VERTICAL, true, true, false);
		File PoissonRVChart = new File("PoissonRVChart.jpeg");
		ChartUtilities.saveChartAsJPEG(PoissonRVChart, poissonRV, width, height);
		
		final XYSeries expectedMM1 = new XYSeries("Theoretical Pn for M/M/1");
		for (int i = 0; i < 100; i++) {
			expectedMM1.add(i, PnMM1(5, 6, i));
		}
		final XYSeries mm1Queue = new XYSeries("Calculated Pn for M/M/1");
		double[] pnMM1 = simMM1(5, 6);
		for (int i = 0; i < pnMM1.length; i ++) {
			mm1Queue.add(i, pnMM1[i]);
		}
		final XYSeriesCollection dataset3 = new XYSeriesCollection();
		dataset3.addSeries(mm1Queue);
		dataset3.addSeries(expectedMM1);
		JFreeChart mm1 = ChartFactory.createXYLineChart("Pn vs n", "n", "Pn", dataset3, PlotOrientation.VERTICAL, true, true, false);
		File MM1Chart = new File("MM1Chart.jpeg");
		ChartUtilities.saveChartAsJPEG(MM1Chart, mm1, width, height);
		
		
		
		
		final XYSeries MEk1Queue = new XYSeries("Calculated Pn for MEk1");
		double[] pnMEk1 = simMEk1(5, 6, 4);
		for (int i = 0; i < pnMEk1.length; i ++) {
			MEk1Queue.add(i, pnMEk1[i]);
		}
		final XYSeries mEk1Queue = new XYSeries("Theoretical Pn for MEk1");
		for (int i = 0; i <= pnMEk1.length; i ++) {
			mEk1Queue.add(i, PnMEk1(5, 6, 4, i));
		}
		final XYSeriesCollection dataset4 = new XYSeriesCollection();
		dataset4.addSeries(mEk1Queue);
		dataset4.addSeries(MEk1Queue);
		JFreeChart mek1 = ChartFactory.createXYLineChart("Pn vs n", "n", "Pn", dataset4, PlotOrientation.VERTICAL, true, true, false);
		File MEk1Chart = new File("MEk1Chart.jpeg");
		ChartUtilities.saveChartAsJPEG(MEk1Chart, mek1, width, height); 
		
		
		final XYSeries expectedNMEk1 = new XYSeries("Theoretical E(n) for M/Ek/1");
		final XYSeries expectedMEk1 = new XYSeries("Calculated E(n) for M/Ek/1");
		final XYSeries expectedMD1 = new XYSeries("E(n) for M/D/1");
		for (double i = 0.1; i <= 0.9; i += 0.1) {
			expectedMD1.add(i, expectedNumberMD1(i));
			expectedMEk1.add(i, enMEk1(5, 5 / i, 40));
			expectedNMEk1.add(i, expectedNumberMEk1(40, i));
		} 
		
		final XYSeriesCollection dataset5 = new XYSeriesCollection();
		dataset5.addSeries(expectedMEk1);
		dataset5.addSeries(expectedMD1);
		dataset5.addSeries(expectedNMEk1);
		JFreeChart EnMEk1 = ChartFactory.createXYLineChart("E(n) vs. Rho for M/Ek/1 and M/D/1 and k = 40, lambda = 5", "rho", "E(n)", dataset5, PlotOrientation.VERTICAL, true, true, false);
		File EnMEk1Chart = new File("EnMEk1Chart.jpeg");
		ChartUtilities.saveChartAsJPEG(EnMEk1Chart, EnMEk1, width, height); 
		} 
		
	}

