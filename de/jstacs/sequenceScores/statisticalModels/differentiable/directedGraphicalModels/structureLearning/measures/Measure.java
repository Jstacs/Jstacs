/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures;

import java.util.ArrayList;
import java.util.Arrays;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.Storable;
import de.jstacs.algorithms.graphs.MST;
import de.jstacs.algorithms.graphs.tensor.Tensor;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.InstanceParameterSet;

/**
 * Class for structure measures that derive an optimal structure with respect to
 * some criterion within a class of possible structures from data.
 * 
 * @author Jan Grau
 */
public abstract class Measure implements Cloneable, Storable, InstantiableFromParameterSet {
	/**
	 * The parameters of this measure
	 */
	protected MeasureParameterSet parameters;
	
	/**
	 * Creates a new {@link Measure} from its XML-representation.
	 * @param xml the XML-representation
	 * @throws NonParsableException the the XML could not be parsed
	 */
	protected Measure( StringBuffer xml ) throws NonParsableException {
		fromXML(xml);
	}
	
	/**
	 * Parses this {@link Measure} from its XML representation. Used in
	 * constructor {@link Measure#Measure(StringBuffer)}.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	protected void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXMLTag() );
		parameters = XMLParser.extractObjectForTags( xml, "parameters", MeasureParameterSet.class );
	}
	
	/**
	 * Creates a new {@link Measure} from its {@link MeasureParameterSet}.
	 * @param parameters the parameters
	 * @throws CloneNotSupportedException if the parameters could not be cloned
	 */
	protected Measure( MeasureParameterSet parameters ) throws CloneNotSupportedException {
		this.parameters = (MeasureParameterSet) parameters.clone();
	}
	
	/**
	 * Default constructor.
	 */
	protected Measure(){
		
	}

	/**
	 * Returns the XML-tag for storing this measure
	 * @return the tag
	 */
	public abstract String getXMLTag();
	
	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, parameters, "parameters");
		XMLParser.addTags(buf, getXMLTag() );
		return buf;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public final InstanceParameterSet<Measure> getCurrentParameterSet() throws Exception {
		return parameters;
	}
	
	/**
	 * Returns the name of the {@link Measure} and possibly some additional
	 * information about the current instance.
	 * 
	 * @return the name of the {@link Measure}
	 */
	public abstract String getInstanceName();

	/**
	 * Returns the optimal parents for the given data and weights. The returned
	 * array of parents <code>p</code> at each position <code>i</code> is build
	 * as follows:
	 * <ul>
	 * <li> <code>p[i][p.length - 1]</code> contains the index <code>i</code>
	 * itself</li>
	 * <li> <code>p[i][p.length - 2]</code> contains the &quot;most
	 * important&quot; parent</li>
	 * <li>...</li>
	 * <li> <code>p[i][0]</code> contains the &quot;least important&quot; parent</li>
	 * </ul>
	 * 
	 * @param fg
	 *            the data of the current (foreground) class
	 * @param bg
	 *            the data of the negative (background) class
	 * @param weightsFg
	 *            the weights for the sequences of <code>fg</code>
	 * @param weightsBg
	 *            the weights for the sequences of <code>bg</code>
	 * @param length
	 *            the length of the model, must be equal to the length of the
	 *            sequences
	 * 
	 * @return the the array <code>p</code> with the optimal parents
	 * 
	 * @throws Exception
	 *             if the lengths do not match or other problems concerning the
	 *             data occur
	 */
	public abstract int[][] getParents(DataSet fg, DataSet bg,
			double[] weightsFg, double[] weightsBg, int length)
			throws Exception;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Measure clone() throws CloneNotSupportedException {
		return (Measure) super.clone();
	}

	
	/**
	 * Prepares a matrix of pairwise association measures for the implementation of Kruskal's algorithm.
	 * @param fullMatrix the full matrix
	 * @return the reduced matrix for Kruskal's algorithm
	 * @see MST#kruskal(double[][])
	 */
	protected double[][] getMatrixForKruskal(double[][] fullMatrix){
		double[][] triang = new double[fullMatrix.length][];
		
		for(int i=0;i<triang.length;i++){
			triang[i] = new double[triang.length-i-1];
			for(int j=0;j<triang[i].length;j++){
				triang[i][j] = fullMatrix[i][i+j+1];
			}
		}
		return triang;
	}
	
	/**
	 * Helper method that converts the structure returned by {@link MST#kruskal(double[][])}
	 * to that returned by the {@link Measure#getParents(DataSet, DataSet, double[], double[], int)} method
	 * @param structure the original structure
	 * @param length the length of sequences
	 * @return the modified structure
	 */
	protected int[][] reStructure(int[][] structure, int length){
		int[][] dep = new int[length][];
		ArrayList<int[]> edges = new ArrayList<int[]>( structure.length );

		for( int counter3 = 0; counter3 < structure.length; counter3++ ) {
			edges.add( structure[counter3] );
			//System.out.println( Arrays.toString( dep2[counter3] ) );
		}

		int[] help;
		boolean[] used = new boolean[length];
		Arrays.fill( used, false );
		dep[0] = new int[]{ (int)0 };
		used[0] = true;

		do {
			for( int counter3 = 0; counter3 < edges.size(); ) {
				help = edges.get( counter3 );
				if( used[help[0]] || used[help[1]] ) {
					if( used[help[1]] ) {
						int counter2 = help[1];
						help[1] = help[0];
						help[0] = counter2;
					}
					dep[help[1]] = edges.remove( counter3 );
					used[help[1]] = true;
				} else {
					counter3++;
				}
			}
		} while( edges.size() > 0 );
		return dep;
	}
	
	/**
	 * Creates a new parent structure as defined by
	 * {@link #getParents(DataSet, DataSet, double[], double[], int)} from an
	 * order and a topological ordering of positions.
	 * 
	 * @param o
	 *            the topological ordering
	 * @param order
	 *            the order
	 * 
	 * @return the parent structure in a two-dimensional array
	 */
	protected static int[][] toParents(int[] o, byte order) {
		int[][] parents = new int[o.length][];
		for (int i = 0; i < parents.length; i++) {
			parents[o[i]] = new int[(order < i ? order : i) + 1];
			for (int j = i; j >= i - order && j >= 0; j--) {
				parents[o[i]][parents[o[i]].length - (i - j) - 1] = o[j];
			}
		}
		return parents;
	}

	/**
	 * Fills a {@link Tensor} <code>t</code> with the weights defined in
	 * <code>weights</code>.
	 * 
	 * @param t
	 *            the {@link Tensor} to be filled
	 * @param weights
	 *            the weights
	 */
	protected static void fillTensor(Tensor t, double[][] weights) {
		for (int i = 0; i < weights.length; i++) {
			for (int j = 0; j < weights[i].length; j++) {
				if (i != j) {
					// t.setValue(i, new int[]{j}, (byte)1, weights[i][j]);
					t.setValue((byte) 1, weights[i][j], i, j);
				} else {
					t.setRootValue(i, weights[i][i]);
				}
			}
		}
	}

	/**
	 * Fills a {@link Tensor} <code>t</code> with the weights defined in
	 * <code>weights</code>.
	 * 
	 * @param t
	 *            the {@link Tensor} to be filled
	 * @param weights
	 *            the weights
	 */
	protected static void fillTensor(Tensor t, double[][][] weights) {

		for (int i = 0; i < weights.length; i++) {
			for (int j = 0; j < weights[i].length; j++) {
				if (i != j) {
					for (int k = j + 1; k < weights[j].length; k++) {
						if (i != k) {
							// t.setValue(i, new int[]{j,k}, (byte)2,
							// weights[i][j][k]);
							t.setValue((byte) 2, weights[i][j][k], i, j, k);
						}
					}
				}
			}
		}
	}

	/**
	 * Computes the mutual information from <code>counts</code> counted on
	 * sequences with a total weight of <code>n</code>.
	 * 
	 * @param counts
	 *            the counts as returned by
	 *            {@link #getStatisticsOrderTwo(DataSet, double[], int, double)}
	 * @param n
	 *            the total weight
	 * 
	 * @return the mutual information
	 */
	protected static double[][][] getMI(double[][][][][][] counts, double n) {
		double[][][] mi = new double[counts.length][counts.length][counts.length];

		for (int i = 0; i < counts.length; i++) {
			for (int j = 0; j < counts[i].length; j++) {
				for (int k = j + 1; k < counts[i][j].length; k++) {
					for (int a = 0; a < counts[i][j][k].length; a++) {
						for (int b = 0; b < counts[i][j][k][a].length; b++) {
							for (int c = 0; c < counts[i][j][k][a][b].length; c++) {
								mi[i][j][k] += (counts[i][j][k][a][b][c] / n)
										* Math
												.log(n
														* counts[i][j][k][a][b][c]
														/ (counts[i][i][i][a][a][a] * counts[j][j][k][b][b][c]));
								mi[i][k][j] += (counts[i][j][k][a][b][c] / n)
										* Math
												.log(n
														* counts[i][j][k][a][b][c]
														/ (counts[i][i][i][a][a][a] * counts[j][j][k][b][b][c]));
							}
						}
					}
				}
			}
		}
		return mi;
	}

	/**
	 * Computes the conditional mutual information from <code>fgStats</code> and
	 * <code>bgStats</code> counted on sequences with a total weight of
	 * <code>n</code>.
	 * 
	 * @param fgStats
	 *            the counts in the foreground sequences as returned by
	 *            {@link #getStatisticsOrderTwo(DataSet, double[], int, double)}
	 * @param bgStats
	 *            the counts in the foreground sequences as returned by
	 *            {@link #getStatisticsOrderTwo(DataSet, double[], int, double)}
	 * @param n
	 *            the total weight
	 * 
	 * @return the conditional mutual information
	 */
	protected static double[][][] getCMI(double[][][][][][] fgStats,
			double[][][][][][] bgStats, double n) {
		double[][][] cmi = new double[fgStats.length][fgStats.length][fgStats.length];

		for (int i = 0; i < fgStats.length; i++) {
			for (int j = 0; j < fgStats[i].length; j++) {
				for (int k = j + 1; k < fgStats[i][j].length; k++) {
					for (int a = 0; a < fgStats[i][j][k].length; a++) {
						for (int b = 0; b < fgStats[i][j][k][a].length; b++) {
							for (int c = 0; c < fgStats[i][j][k][a][b].length; c++) {
								cmi[i][j][k] += (fgStats[i][j][k][a][b][c] / n)
										* Math
												.log((fgStats[i][j][k][a][b][c] * (fgStats[j][j][k][b][b][c] + bgStats[j][j][k][b][b][c]))
														/ ((fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]) * (fgStats[j][j][k][b][b][c])));
								cmi[i][j][k] += (bgStats[i][j][k][a][b][c] / n)
										* Math
												.log((bgStats[i][j][k][a][b][c] * (fgStats[j][j][k][b][b][c] + bgStats[j][j][k][b][b][c]))
														/ ((fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]) * (bgStats[j][j][k][b][b][c])));
								cmi[i][k][j] += (fgStats[i][j][k][a][b][c] / n)
										* Math
												.log((fgStats[i][j][k][a][b][c] * (fgStats[j][j][k][b][b][c] + bgStats[j][j][k][b][b][c]))
														/ ((fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]) * (fgStats[j][j][k][b][b][c])));
								cmi[i][k][j] += (bgStats[i][j][k][a][b][c] / n)
										* Math
												.log((bgStats[i][j][k][a][b][c] * (fgStats[j][j][k][b][b][c] + bgStats[j][j][k][b][b][c]))
														/ ((fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]) * (bgStats[j][j][k][b][b][c])));
							}
						}
					}
				}
			}
		}

		return cmi;
	}

	/**
	 * Computes the explaining away residual from <code>fgStats</code> and
	 * <code>bgStats</code> counted on sequences with a total weight of
	 * <code>nFg</code> and <code>nBg</code>, respectively.
	 * 
	 * @param fgStats
	 *            the counts in the foreground sequences as returned by
	 *            {@link #getStatisticsOrderTwo(DataSet, double[], int, double)}
	 * @param bgStats
	 *            the counts in the foreground sequences as returned by
	 *            {@link #getStatisticsOrderTwo(DataSet, double[], int, double)}
	 * @param nFg
	 *            the total weight in the foreground
	 * @param nBg
	 *            the total weight in the background
	 * 
	 * @return the explaining away residual
	 */
	public static double[][][] getEAR(double[][][][][][] fgStats,
			double[][][][][][] bgStats, double nFg, double nBg) {
		double pcFg = nFg / (nFg + nBg);
		double pcBg = nBg / (nFg + nBg);

		double[][][] ear = new double[fgStats.length][fgStats.length][fgStats.length];

		for (int i = 0; i < fgStats.length; i++) {
			for (int j = 0; j < fgStats[i].length; j++) {
				if (i != j) {
					for (int k = 0; k < fgStats[i][j].length; k++) {
						if (k != j && k != i) {
							for (int a = 0; a < fgStats[i][j][k].length; a++) {
								for (int b = 0; b < fgStats[i][j][k][a].length; b++) {
									for (int c = 0; c < fgStats[i][j][k][a][b].length; c++) {
										ear[i][j][k] += pcFg
												* ((fgStats[i][j][k][a][b][c] / nFg) * Math
														.log(nFg
																* fgStats[i][j][k][a][b][c]
																/ (fgStats[i][i][i][a][a][a] * fgStats[j][j][k][b][b][c])));
										ear[i][j][k] += pcBg
												* ((bgStats[i][j][k][a][b][c] / nBg) * Math
														.log(nBg
																* bgStats[i][j][k][a][b][c]
																/ (bgStats[i][i][i][a][a][a] * bgStats[j][j][k][b][b][c])));
										ear[i][j][k] -= ((fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]) / (nFg + nBg))
												* Math
														.log(((nFg + nBg) * (fgStats[i][j][k][a][b][c] + bgStats[i][j][k][a][b][c]))
																/ ((fgStats[i][i][i][a][a][a] + bgStats[i][i][i][a][a][a]) * (fgStats[j][j][k][b][b][c] + bgStats[j][j][k][b][b][c])));
									}
								}
							}
						}
					}
				}
			}
		}
		return ear;
	}

	/**
	 * Counts the occurrences of symbols of the {@link AlphabetContainer} of
	 * {@link DataSet} <code>s</code> using <code>weights</code>. The array
	 * <code>counts</code> is indexed as follows:<br>
	 * <code>counts[first index][second index][third index][symbol at first index][symbol at second index][symbol at third index]</code>
	 * .
	 * 
	 * @param s
	 *            the data
	 * @param weights
	 *            the weights
	 * @param length
	 *            the length of the sequences
	 * @param ess
	 *            the equivalent sample size
	 * 
	 * @return the array <code>counts</code> with the symbol occurrences
	 * 
	 * @throws Exception
	 *             if the lengths do not match or other problems concerning the
	 *             data occur
	 */
	protected static double[][][][][][] getStatisticsOrderTwo(DataSet s,
			double[] weights, int length, double ess) throws Exception {
		if (s.getElementLength() != length) {
			throw new Exception("Lengths do not match.");
		}
		AlphabetContainer alph = s.getAlphabetContainer();
		double currEss = 0;
		double[][][][][][] counts = new double[length][length][length][][][];
		for (int i = 0; i < counts.length; i++) {
			for (int j = 0; j < counts[i].length; j++) {
				for (int k = 0; k < counts[i][j].length; k++) {
					counts[i][j][k] = new double[(int) alph
							.getAlphabetLengthAt(i)][(int) alph
							.getAlphabetLengthAt(j)][(int) alph
							.getAlphabetLengthAt(k)];
					currEss = ess
							/ (alph.getAlphabetLengthAt(i)
									* alph.getAlphabetLengthAt(j) * alph
									.getAlphabetLengthAt(k));
					if (i == j && i == k) {
						currEss = ess / alph.getAlphabetLengthAt(i);
					} else if (i == j) {
						currEss = ess
								/ (alph.getAlphabetLengthAt(i) * alph
										.getAlphabetLengthAt(k));
					} else if (i == k || j == k) {
						currEss = ess
								/ (alph.getAlphabetLengthAt(i) * alph
										.getAlphabetLengthAt(j));
					}
					for (int l = 0; l < counts[i][j][k].length; l++) {
						for (int m = 0; m < counts[i][j][k][l].length; m++) {
							java.util.Arrays.fill(counts[i][j][k][l][m],
									currEss);
						}
					}
				}
			}
		}

		for (int i = 0; i < s.getNumberOfElements(); i++) {
			Sequence seq = s.getElementAt(i);
			double w = weights[i];
			for (int j = 0; j < seq.getLength(); j++) {
				counts[j][j][j][seq.discreteVal(j)][seq.discreteVal(j)][seq
						.discreteVal(j)] += w;
				for (int k = j + 1; k < seq.getLength(); k++) {
					counts[j][j][k][seq.discreteVal(j)][seq.discreteVal(j)][seq
							.discreteVal(k)] += w;
					counts[j][k][j][seq.discreteVal(j)][seq.discreteVal(k)][seq
							.discreteVal(j)] += w;
					counts[k][j][j][seq.discreteVal(k)][seq.discreteVal(j)][seq
							.discreteVal(j)] += w;
					counts[k][k][j][seq.discreteVal(k)][seq.discreteVal(k)][seq
							.discreteVal(j)] += w;
					counts[k][j][k][seq.discreteVal(k)][seq.discreteVal(j)][seq
							.discreteVal(k)] += w;
					counts[j][k][k][seq.discreteVal(j)][seq.discreteVal(k)][seq
							.discreteVal(k)] += w;
					for (int l = k + 1; l < seq.getLength(); l++) {
						counts[j][k][l][seq.discreteVal(j)][seq.discreteVal(k)][seq
								.discreteVal(l)] += w;
						counts[k][j][l][seq.discreteVal(k)][seq.discreteVal(j)][seq
								.discreteVal(l)] += w;
						counts[l][j][k][seq.discreteVal(l)][seq.discreteVal(j)][seq
								.discreteVal(k)] += w;
						counts[j][l][k][seq.discreteVal(j)][seq.discreteVal(l)][seq
								.discreteVal(k)] += w;
						counts[l][k][j][seq.discreteVal(l)][seq.discreteVal(k)][seq
								.discreteVal(j)] += w;
						counts[k][l][j][seq.discreteVal(k)][seq.discreteVal(l)][seq
								.discreteVal(j)] += w;
					}
				}
			}
		}

		return counts;
	}

	/**
	 * Counts the occurrences of symbols of the {@link AlphabetContainer} of
	 * {@link DataSet} <code>s</code> using <code>weights</code>. The array
	 * <code>counts</code> is indexed as follows:<br>
	 * <code>counts[first index][second index][symbol at first index][symbol at second index]</code>
	 * .
	 * 
	 * @param s
	 *            the data
	 * @param weights
	 *            the weights
	 * @param length
	 *            the length of the sequences
	 * @param ess
	 *            the equivalent sample size
	 * 
	 * @return the array <code>counts</code> with the symbol occurrences
	 * 
	 * @throws Exception
	 *             if the lengths do not match or other problems concerning the
	 *             data occur
	 */
	protected static double[][][][] getStatistics(DataSet s, double[] weights,
			int length, double ess) throws Exception {
		if (s.getElementLength() != length) {
			throw new Exception("Lengths do not match.");
		}
		AlphabetContainer alph = s.getAlphabetContainer();
		double currEss = 0;
		double[][][][] counts = new double[length][length][][];
		for (int i = 0; i < counts.length; i++) {
			for (int j = 0; j < counts[i].length; j++) {
				counts[i][j] = new double[(int) alph.getAlphabetLengthAt(i)][(int) alph
						.getAlphabetLengthAt(j)];
				currEss = i == j ? ess / alph.getAlphabetLengthAt(i) : ess
						/ (alph.getAlphabetLengthAt(i) * alph
								.getAlphabetLengthAt(j));
				for (int k = 0; k < counts[i][j].length; k++) {
					java.util.Arrays.fill(counts[i][j][k], currEss);
				}
			}
		}

		for (int i = 0; i < s.getNumberOfElements(); i++) {
			Sequence seq = s.getElementAt(i);
			double w = weights[i];
			for (int j = 0; j < seq.getLength(); j++) {
				counts[j][j][seq.discreteVal(j)][seq.discreteVal(j)] += w;
				for (int k = j + 1; k < seq.getLength(); k++) {
					counts[j][k][seq.discreteVal(j)][seq.discreteVal(k)] += w;
					counts[k][j][seq.discreteVal(k)][seq.discreteVal(j)] += w;
				}
			}
		}
		return counts;
	}

	/**
	 * Computes the mutual information from <code>counts</code> counted on
	 * sequences with a total weight of <code>n</code>.
	 * 
	 * @param counts
	 *            the counts as defined in
	 *            {@link #getStatistics(DataSet, double[], int, double)}.
	 * @param n
	 *            the total weight
	 * 
	 * @return the mutual information
	 */
	protected static double[][] getMI(double[][][][] counts, double n) {
		double[][] mi = new double[counts.length][counts.length];

		for (int i = 0; i < counts.length; i++) {
			for (int j = 0; j < counts[i].length; j++) {
				if (i != j) {
					for (int a = 0; a < counts[i][j].length; a++) {
						for (int b = 0; b < counts[i][j][a].length; b++) {
							if(counts[i][j][a][b] > 0){
								mi[i][j] += (counts[i][j][a][b] / n)
										* Math
												.log(n
														* counts[i][j][a][b]
														/ (counts[i][i][a][a] * counts[j][j][b][b]));
							}
						}
					}
				} else {
					for (int a = 0; a < counts[i][i].length; a++) {
						if(counts[i][i][a][a] > 0){
							mi[i][i] -= (counts[i][i][a][a] / n)
									* Math.log((counts[i][i][a][a] / n));
						}
					}
				}
			}
		}
		return mi;
	}

	/**
	 * Computes the conditional mutual information from <code>fgStats</code> and
	 * <code>bgStats</code> counted on sequences with a total weight of
	 * <code>nFg</code> and <code>nBg</code>, respectively.
	 * 
	 * @param fgStats
	 *            the counts as defined in
	 *            {@link #getStatistics(DataSet, double[], int, double)} on the
	 *            foreground
	 * @param bgStats
	 *            the counts as defined in
	 *            {@link #getStatistics(DataSet, double[], int, double)} on the
	 *            background
	 * @param n
	 *            the total weight
	 * @param nFg
	 *            the total weight in the foreground
	 * @param nBg
	 *            the total weight in the background
	 * 
	 * @return the conditional mutual information
	 */
	protected static double[][] getCMI(double[][][][] fgStats,
			double[][][][] bgStats, double n, double nFg, double nBg) {
		double[][] cmi = new double[fgStats.length][fgStats.length];

		for (int i = 0; i < fgStats.length; i++) {
			for (int j = 0; j < fgStats[i].length; j++) {
				if (i != j) {
					for (int a = 0; a < fgStats[i][j].length; a++) {
						for (int b = 0; b < fgStats[i][j][a].length; b++) {
							cmi[i][j] += (fgStats[i][j][a][b] / n)
									* Math
											.log((fgStats[i][j][a][b] * (fgStats[j][j][b][b] + bgStats[j][j][b][b]))
													/ ((fgStats[i][j][a][b] + bgStats[i][j][a][b]) * (fgStats[j][j][b][b])));
							cmi[i][j] += (bgStats[i][j][a][b] / n)
									* Math
											.log((bgStats[i][j][a][b] * (fgStats[j][j][b][b] + bgStats[j][j][b][b]))
													/ ((fgStats[i][j][a][b] + bgStats[i][j][a][b]) * (bgStats[j][j][b][b])));
						}
					}
				} else {
					for (int a = 0; a < fgStats[i][i].length; a++) {
						cmi[i][i] += (fgStats[i][i][a][a] / n)
								* Math
										.log((fgStats[i][i][a][a] * n)
												/ ((fgStats[i][i][a][a] + bgStats[i][i][a][a]) * (nFg)));
						cmi[i][i] += (bgStats[i][i][a][a] / n)
								* Math
										.log((bgStats[i][i][a][a] * n)
												/ ((fgStats[i][i][a][a] + bgStats[i][i][a][a]) * (nBg)));
					}
				}
			}
		}

		return cmi;
	}

	/**
	 * Computes the explaining away residual from <code>fgStats</code> and
	 * <code>bgStats</code> counted on sequences with a total weight of
	 * <code>nFg</code> and <code>nBg</code>, respectively.
	 * 
	 * @param fgStats
	 *            the counts as defined in
	 *            {@link #getStatistics(DataSet, double[], int, double)} on the
	 *            foreground
	 * @param bgStats
	 *            the counts as defined in
	 *            {@link #getStatistics(DataSet, double[], int, double)} on the
	 *            background
	 * @param nFg
	 *            the total weight in the foreground
	 * @param nBg
	 *            the total weight in the background
	 * 
	 * @return the explaining away residual
	 */
	protected static double[][] getEAR(double[][][][] fgStats,
			double[][][][] bgStats, double nFg, double nBg) {
		double pcFg = nFg / (nFg + nBg);
		double pcBg = nBg / (nFg + nBg);

		double[][] ear = new double[fgStats.length][fgStats.length];

		for (int i = 0; i < fgStats.length; i++) {
			for (int j = 0; j < fgStats[i].length; j++) {
				if (i != j) {
					for (int a = 0; a < fgStats[i][j].length; a++) {
						for (int b = 0; b < fgStats[i][j][a].length; b++) {
							ear[i][j] += pcFg
									* ((fgStats[i][j][a][b] / nFg) * Math
											.log(nFg
													* fgStats[i][j][a][b]
													/ (fgStats[i][i][a][a] * fgStats[j][j][b][b])));
							ear[i][j] += pcBg
									* ((bgStats[i][j][a][b] / nBg) * Math
											.log(nBg
													* bgStats[i][j][a][b]
													/ (bgStats[i][i][a][a] * bgStats[j][j][b][b])));
							ear[i][j] -= ((fgStats[i][j][a][b] + bgStats[i][j][a][b]) / (nFg + nBg))
									* Math
											.log(((nFg + nBg) * (fgStats[i][j][a][b] + bgStats[i][j][a][b]))
													/ ((fgStats[i][i][a][a] + bgStats[i][i][a][a]) * (fgStats[j][j][b][b] + bgStats[j][j][b][b])));
						}
					}
				} else {
					ear[i][i] = 0;
				}
			}
		}
		return ear;
	}

	/**
	 * Computes the sum of all elements in the array <code>ar</code>.
	 * 
	 * @param ar
	 *            the array
	 * 
	 * @return the sum of the elements of the array
	 */
	protected static double sum(double[] ar) {
		double sum = 0;
		for (int i = 0; i < ar.length; i++) {
			sum += ar[i];
		}
		return sum;
	}

	/**
	 * Linearizes the arrays in the two-dimensional array <code>ar</code> to
	 * form a new, one-dimensional array.
	 * 
	 * @param ar
	 *            the two-dimensional array
	 * 
	 * @return the linearized one-dimensional array
	 */
	protected static double[] union(double[][] ar) {
		int num = 0;
		for (int i = 0; i < ar.length; i++) {
			num += ar[i].length;
		}
		double[] un = new double[num];
		num = 0;
		for (int i = 0; i < ar.length; i++) {
			System.arraycopy(ar[i], 0, un, num, ar[i].length);
			num += ar[i].length;
		}
		return un;
	}

	/**
	 * Indicates if {@link Measure} supports shifts.
	 * 
	 * @return if {@link Measure} supports shifts
	 */
	public boolean isShiftable() {
		return false;
	}
	
	/**
	 * This class is the super class of any {@link de.jstacs.parameters.ParameterSet} that can be used to instantiate a {@link Measure}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static abstract class MeasureParameterSet extends InstanceParameterSet<Measure> {

		/**
		 * Creates a new empty {@link MeasureParameterSet} for the given sub-class
		 * of {@link Measure},
		 * @param clazz the sub-class
		 */
		protected MeasureParameterSet( Class<? extends Measure> clazz ) {
			super( clazz );
		}
		
		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Recreates a {@link MeasureParameterSet} from its XML representation as
		 * returned by {@link #toXML()}.
		 * 
		 * @param xml  the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException if the XML code could not be parsed
		 */
		protected MeasureParameterSet( StringBuffer xml ) throws NonParsableException {
			super( xml );
		}
	}
}
