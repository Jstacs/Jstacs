/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.models.discrete.inhomogeneous;

import java.util.ArrayList;
import java.util.Arrays;

import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.graphs.DAG;
import de.jstacs.algorithms.graphs.MST;
import de.jstacs.algorithms.graphs.tensor.SymmetricTensor;
import de.jstacs.algorithms.graphs.tensor.Tensor;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.models.discrete.ConstraintManager;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This class can be used to learn the structure of any discrete model.
 * 
 * @author Jens Keilwagen
 */
public class StructureLearner {

	/**
	 * This <code>enum</code> defines the different types of learning that are
	 * possible with the {@link StructureLearner}.
	 * 
	 * @author Jens Keilwagen
	 */
	public enum LearningType {
		/**
		 * This value indicates that the model should learn its parameters using
		 * the <b>m</b>aximum <b>l</b>ikelihood approach (ML) (if ess = 0) or
		 * the <b>m</b>aximum <b>a</b> <b>p</b>osteriori approach (MAP)
		 * (otherwise).
		 */
		ML_OR_MAP,
		/**
		 * This value indicates that the model should learn its parameters using
		 * <b>B</b>ayesian <b>m</b>odel <b>a</b>veraging (BMA).
		 */
		BMA;
	}

	/**
	 * This <code>enum</code> defines the different types of models that can be
	 * learned with the {@link StructureLearner}.
	 * 
	 * @author Jens Keilwagen
	 */
	public enum ModelType {
		/**
		 * This value indicates that the model is an <b>i</b>nhomogeneous
		 * <b>M</b>arkov <b>m</b>odel (iMM).
		 */
		IMM,
		/**
		 * This value indicates that the model should compute a <b>p</b>ermuted
		 * <b>M</b>arkov <b>m</b>odel (pMM).
		 */
		PMM,
		/**
		 * This value indicates that the model should compute a <b>B</b>ayesian
		 * <b>n</b>etwork (BN).
		 */
		BN;
	}

	private AlphabetContainer con;

	private int length;

	private double ess;

	private int[] alphabetLength;

	/**
	 * Creates a new {@link StructureLearner} for a given
	 * {@link AlphabetContainer}, a given length and a given <b>e</b>quivalent
	 * <b>s</b>ample <b>s</b>ize (ess).
	 * 
	 * @param con
	 *            the alphabets this instance should use
	 * @param length
	 *            the length
	 * @param ess
	 *            the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize, has to
	 *            be non-negative)
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link AlphabetContainer} is not discrete, the length
	 *             is not matching with the {@link AlphabetContainer} or the ess
	 *             is below 0
	 */
	public StructureLearner( AlphabetContainer con, int length, double ess ) throws IllegalArgumentException {
		if( !con.isDiscrete() ) {
			throw new IllegalArgumentException( "The instance of AlphabetContainer has to be discrete." );
		}
		int i = con.getPossibleLength();
		if( i != 0 && i != length ) {
			throw new IllegalArgumentException( "The instance of AlphabetContainer and length are not matching." );
		}
		this.con = con;
		this.length = length;
		alphabetLength = new int[length];
		for( i = 0; i < length; i++ ) {
			alphabetLength[i] = (int)con.getAlphabetLengthAt( i );
		}
		setESS( ess );
	}

	/**
	 * Creates a {@link StructureLearner} with <b>e</b>quivalent <b>s</b>ample
	 * <b>s</b>ize (ess) = 0.
	 * 
	 * @param con
	 *            the alphabets this instance should use
	 * @param length
	 *            the length
	 * 
	 * @throws IllegalArgumentException
	 *             if the {@link AlphabetContainer} is not discrete or the
	 *             length is not matching with the {@link AlphabetContainer}
	 * 
	 * @see StructureLearner#StructureLearner(AlphabetContainer, int, double)
	 */
	public StructureLearner( AlphabetContainer con, int length ) throws IllegalArgumentException {
		this( con, length, 0 );
	}

	/**
	 * This method returns the {@link AlphabetContainer} of the
	 * {@link StructureLearner}.
	 * 
	 * @return the {@link AlphabetContainer} of the {@link StructureLearner}
	 */
	public AlphabetContainer getAlphabetContainer() {
		return con;
	}

	/**
	 * This method returns the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize)
	 * of the {@link StructureLearner}.
	 * 
	 * @return the ess of the {@link StructureLearner}
	 */
	public double getEss() {
		return ess;
	}

	/**
	 * This method sets the ess (<b>e</b>quivalent <b>s</b>ample <b>s</b>ize) of
	 * the {@link StructureLearner}.
	 * 
	 * @param ess
	 *            the ess of the {@link StructureLearner}
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>ess < 0</code>
	 */
	public void setESS( double ess ) throws IllegalArgumentException {
		if( ess < 0 ) {
			throw new IllegalArgumentException( "The value for ess has to be non-negative." );
		}
		this.ess = ess;
	}

	/**
	 * This method finds the optimal structure of a model by using a given
	 * learning method (in some sense).
	 * 
	 * @param data
	 *            the {@link DataSet}
	 * @param weights
	 *            the weights
	 * @param model
	 *            the kind of model
	 * @param order
	 *            the Markov order
	 * @param method
	 *            the learning method
	 * 
	 * @return the optimal structure of the specified model (in some sense)
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see StructureLearner#getTensor(DataSet, double[], byte, LearningType)
	 * @see StructureLearner#getStructure(Tensor, ModelType, byte)
	 */
	public int[][] getStructure( DataSet data, double[] weights, ModelType model, byte order, LearningType method ) throws Exception {
		int[][] dep;
		int counter1, counter2, counter3, idx;
		if( order == 0 ) {
			// the most simple case: PWM == IMM(0) == PMM(0) == BN(0)
			dep = new int[length][1];
			for( counter1 = 0; counter1 < length; counter1++ ) {
				dep[counter1][0] = counter1;
			}
		} else if( model == ModelType.IMM ) {
			// inhomogeneous markov models have a fixed structure
			dep = new int[length][];
			for( counter3 = 0, counter1 = 1; counter1 <= order; counter1++, counter3++ ) {
				dep[counter3] = new int[counter1];
				for( counter2 = 0; counter2 < counter1; counter2++ ) {
					dep[counter3][counter2] = counter2;
				}
			}
			for( counter1 = 0; counter3 < length; counter1++, counter3++ ) {
				dep[counter3] = new int[order + 1];
				for( counter2 = counter1, idx = 0; counter2 <= counter3; counter2++, idx++ ) {
					dep[counter3][idx] = counter2;
				}
			}
		} else {
			// learn structure of either BN or PMM
			Tensor t = getTensor( data, weights, order, method );
			dep = getStructure( t, model, order );
		}
		return dep;
	}

	/**
	 * This method can be used to determine the optimal structure of a model.
	 * 
	 * @param t
	 *            the tensor containing all relevant weights (includes the
	 *            learning method for the structure)
	 * @param model
	 *            the model type
	 * @param order
	 *            the model order
	 * 
	 * @return the optimal structure of the specified model
	 * 
	 * @throws Exception
	 *             if something in the algorithm went wrong
	 * 
	 * @see StructureLearner#getTensor(DataSet, double[], byte, LearningType)
	 */
	public static int[][] getStructure( Tensor t, ModelType model, byte order ) throws Exception {
		int length = t.getNumberOfNodes(), counter1, counter2, counter3;
		int[][] dep;
		if( model == ModelType.BN ) {
			if( order == 1 ) {
				// extract edge weights
				double[][] w = new double[length][];
				for( counter1 = 0; counter1 < w.length; counter1++ ) {
					w[counter1] = new double[length - 1 - counter1];
					counter2 = counter1;
					for( counter3 = 0, counter2++; counter3 < w[counter1].length; counter2++, counter3++ ) {
						w[counter1][counter3] = t.getValue( order, counter1, new int[]{ counter2 } );
						//System.out.println( counter1 + " -- " + counter2 + "\t" + w[counter1][counter3] );
					}
				}

				// compute MST
				int[][] dep2 = MST.kruskal( w );

				// re-structure
				dep = new int[length][];
				ArrayList<int[]> edges = new ArrayList<int[]>( dep2.length );

				for( counter3 = 0; counter3 < dep2.length; counter3++ ) {
					edges.add( dep2[counter3] );
					//System.out.println( Arrays.toString( dep2[counter3] ) );
				}

				int[] help;
				boolean[] used = new boolean[length];
				Arrays.fill( used, false );
				dep[0] = new int[]{ (int)0 };
				used[0] = true;

				do {
					for( counter3 = 0; counter3 < edges.size(); ) {
						help = edges.get( counter3 );
						if( used[help[0]] || used[help[1]] ) {
							if( used[help[1]] ) {
								counter2 = help[1];
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
			} else {
				// Bayesian networks with order > 1
				dep = DAG.computeMaximalKDAG( t );
			}
		} else {
			// PMM
			dep = DAG.getStructureFromPath( DAG.computeMaximalHP( t ), t );
		}

		// System.out.println( "score: " + DAG.getScore( t, dep ) );
		return dep;
	}

	private double[][] getSummands( DataSet data, double[] weights, byte order, LearningType method, double[] extra ) throws IllegalArgumentException,
			WrongAlphabetException {
		double help;
		if( method == LearningType.BMA && ess == 0 ) {
			throw new IllegalArgumentException( "The ESS has to be strict positive for BMA." );
		}
		InhCondProb[][] constr = new InhCondProb[order + 1][];
		ArrayList<InhCondProb> list = new ArrayList<InhCondProb>();
		byte counter1, counter3;
		int counter2;
		long l;

		// create all relavant constraints
		CombinationIterator com = new CombinationIterator( length, (byte)( order + 1 ) );
		for( counter1 = 1, counter3 = 0; counter3 <= order; counter1++, counter3++ ) {
			com.setCurrentLength( counter1 );

			l = com.getNumberOfCombinations( counter1 );
			if( l > Integer.MAX_VALUE ) {
				throw new IllegalArgumentException();
			}
			counter2 = (int)l;
			constr[counter3] = new InhCondProb[counter2];
			counter2--;
			while( counter2 >= 0 ) {
				constr[counter3][counter2] = new InhCondProb( com.getCombination(), alphabetLength, false );
				list.add( constr[counter3][counter2] );
				com.next();
				counter2--;
			}
		}
		// fill constraints
		double sum = ConstraintManager.countInhomogeneous( con, length, data, weights, true, list.toArray( new InhCondProb[0] ) ), all = sum + ess;

		// compute summand
		double[][] h = new double[order + 1][];
		for( counter1 = 0; counter1 <= order; counter1++ ) {
			h[counter1] = new double[constr[counter1].length];
			for( counter2 = 0; counter2 < h[counter1].length; counter2++ ) {
				if( method == LearningType.ML_OR_MAP ) {
					constr[counter1][counter2].estimateUnConditional( ess, sum );
					h[counter1][counter2] = all * ConstraintManager.getEntropy( constr[counter1][counter2] );
					if( ess > 0 ) {
						help = constr[counter1][counter2].getNumberOfSpecificConstraints();
						h[counter1][counter2] -= help * Gamma.logOfGamma( ess / help );
					}
				} else {
					h[counter1][counter2] = ConstraintManager.getLogGammaSum( constr[counter1][counter2], ess );
				}
			}
		}
		if( ess > 0 ) {
			extra[0] = Gamma.logOfGamma( ess );
		} else {
			extra[0] = 0;
		}
		if( method == LearningType.BMA ) {
			extra[0] -= Gamma.logOfGamma( all );
		}
		return h;
	}

	/**
	 * This method can be used to compute a {@link Tensor} that can be used to
	 * determine the optimal structure.
	 * 
	 * @param data
	 *            the data
	 * @param weights
	 *            the weights
	 * @param order
	 *            the Markov order
	 * @param method
	 *            the learning type
	 * 
	 * @return a tensor containing all values the are necessary to determine the
	 *         optimal structure (includes the learning method for the
	 *         structure)
	 * 
	 * @throws IllegalArgumentException
	 *             if something is wrong with the given arguments
	 * @throws WrongAlphabetException
	 *             if the {@link AlphabetContainer} of the data is not correct
	 * 
	 * @see StructureLearner#getStructure(Tensor, ModelType, byte)
	 */
	public SymmetricTensor getTensor( DataSet data, double[] weights, byte order, LearningType method ) throws IllegalArgumentException,
			WrongAlphabetException {
		double[] extra = new double[1];
		return fillTensor( getSummands( data, weights, order, method, extra ), order, extra[0] );
	}

	/**
	 * Fills the {@link Tensor} with the edge weights.
	 * 
	 * @param summands
	 *            the summands for the edge weights ( edge weight:
	 *            <code>w(x,\vector{y}) = summands(x) + summands(y) -
	 *            summands(x,y) - extra</code> )
	 * @param order
	 *            the maximal order of the {@link Tensor}
	 * @param extra
	 *            an extra value that is subtracted from all edges
	 * 
	 * @see SymmetricTensor
	 */
	private SymmetricTensor fillTensor( double[][] summands, byte order, double extra ) {
		SymmetricTensor t = new SymmetricTensor( length, order );
		int counter1, counter2, idx, swap, help, ind;
		byte counter3 = 0;
		CombinationIterator com = new CombinationIterator( length, (byte)( order + 1 ) );

		// compute edge weight:
		// w(x,\vector{y}) = summands(x) + summands(y) - summands(x,y) - extra

		int[] parents, comb, comb2;
		long l;
		boolean used[] = new boolean[length];
		for( counter3 = 1; counter3 <= order; counter3++ ) {
			com.setCurrentLength( counter3 );
			parents = new int[counter3];
			l = com.getNumberOfCombinations( counter3 );
			if( l > Integer.MAX_VALUE ) {
				throw new IllegalArgumentException();
			}
			counter2 = (int)l;
			counter2--;
			comb2 = new int[counter3 + 1];
			while( counter2 >= 0 ) {
				comb = com.getCombination();
				Arrays.fill( used, false );
				for( counter1 = 0; counter1 < counter3; counter1++ ) {
					parents[counter1] = comb[counter1];
					used[parents[counter1]] = true;
				}
				l = com.getIndex( comb );
				if( l > Integer.MAX_VALUE ) {
					throw new IllegalArgumentException();
				}
				idx = (int)l;
				System.arraycopy( comb, 0, comb2, 1, counter3 );
				swap = 0;
				for( counter1 = 0; counter1 < length; counter1++ ) {
					if( !used[counter1] ) {
						comb2[swap] = counter1;
						while( swap < counter3 && comb2[swap] > comb2[swap + 1] ) {
							help = comb2[swap];
							comb2[swap] = comb2[++swap];
							comb2[swap] = help;
						}
						//System.out.println( summands[0][length - 1 - counter1] + " + " + summands[counter3 - 1][idx] + " - " + summands[counter3][com.getIndex( comb2 )] + " - " + extra );
						l = com.getIndex( comb2 );
						if( l > Integer.MAX_VALUE ) {
							throw new IllegalArgumentException();
						}
						ind = (int)l;
						t.setValue( counter3, summands[0][length - 1 - counter1] + summands[counter3 - 1][idx]
												- summands[counter3][ind]
												- extra, counter1, parents );
					}
				}
				com.next();
				counter2--;
			}
		}
		return t;
	}
}
