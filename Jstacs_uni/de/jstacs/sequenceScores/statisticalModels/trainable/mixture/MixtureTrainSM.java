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

package de.jstacs.sequenceScores.statisticalModels.trainable.mixture;

import java.text.NumberFormat;
import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.utils.random.MRGParams;
import de.jstacs.utils.random.MultivariateRandomGenerator;

/**
 * The class for a mixture model of any {@link TrainableStatisticalModel}s.
 * 
 * <br>
 * <br>
 * 
 * If you use Gibbs sampling temporary files will be created in the Java temp
 * folder. These files will be deleted if no reference to the current instance
 * exists and the Garbage Collector is called. Therefore it is recommended to
 * call the Garbage Collector explicitly at the end of any application.
 * 
 * @author Jens Keilwagen, Berit Haldemann
 */
public class MixtureTrainSM extends AbstractMixtureTrainSM {

	/**
	 * Creates a new {@link MixtureTrainSM}. This constructor can be used for any
	 * algorithm since it takes all necessary values as parameters.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link MixtureTrainSM}, if the
	 *            model is trained using
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param estimateComponentProbs
	 *            the switch for estimating the component probabilities in the
	 *            algorithm or to hold them fixed; if the component parameters
	 *            are fixed, the values of <code>weights</code> will be used,
	 *            otherwise the <code>componentHyperParams</code> will be
	 *            incorporated in the adjustment
	 * @param componentHyperParams
	 *            the hyperparameters for the component assignment prior
	 *            <ul>
	 *            <li>will only be used if
	 *            <code>estimateComponentProbs == true</code>
	 *            <li>the array has to be <code>null</code> or has to have
	 *            length <code>models.length</code>
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == models.length</code>)
	 * @param algorithm
	 *            either {@link AbstractMixtureTrainSM.Algorithm#EM} or
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}
	 * @param alpha
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link AbstractMixtureTrainSM.Parameterization#THETA} or
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code> <li><code>dimension &lt; 1</code> <li>
	 *             <code>weights != null &amp;&amp; weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an <code>i
	 *             </code> where <code>weights[i] &lt; 0</code> <li><code>starts
	 *             &lt; 1</code> <li><code>componentHyperParams</code> are not
	 *             correct <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 */
	protected MixtureTrainSM( int length, TrainableStatisticalModel[] models, int starts, boolean estimateComponentProbs, double[] componentHyperParams,
							double[] weights, Algorithm algorithm, double alpha, TerminationCondition tc, Parameterization parametrization, //EM
							int initialIteration, int stationaryIteration, BurnInTest burnInTest ) //GIBBS_SAMPLING
																									throws IllegalArgumentException,
																									WrongAlphabetException,
																									CloneNotSupportedException {
		super( length,
				models,
				null,
				models.length,
				starts,
				estimateComponentProbs,
				componentHyperParams,
				weights,
				algorithm,
				alpha,
				tc,
				parametrization,
				initialIteration,
				stationaryIteration,
				burnInTest );
	}

	/**
	 * Creates an instance using EM and estimating the component probabilities.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link MixtureTrainSM}, if the
	 *            model is trained using
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param componentHyperParams
	 *            the hyperparameters for the component assignment prior
	 *            <ul>
	 *            <li>will only be used if
	 *            <code>estimateComponentProbs == true</code>
	 *            <li>the array has to be <code>null</code> or has to have
	 *            length <code>models.length</code>
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param alpha
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the type of the component probability parameterization
	 *            <ul>
	 *            <li>{@link AbstractMixtureTrainSM.Parameterization#THETA} or
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null &amp;&amp; weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see MixtureTrainSM#MixtureTrainSM(int, de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel[], int,
	 *      boolean, double[], double[],
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest)
	 * @see AbstractMixtureTrainSM.Algorithm#EM
	 */
	public MixtureTrainSM( int length, TrainableStatisticalModel[] models, int starts, double[] componentHyperParams, double alpha, TerminationCondition tc,
							Parameterization parametrization ) throws IllegalArgumentException, WrongAlphabetException,
																CloneNotSupportedException {
		this( length, models, starts, true, componentHyperParams, null, Algorithm.EM, alpha, tc, parametrization, 0, 0, null );
	}

	/**
	 * Creates an instance using EM and fixed component probabilities.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link MixtureTrainSM}, if the
	 *            model is trained using
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == models.length</code>)
	 * @param alpha
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link AbstractMixtureTrainSM.Parameterization#THETA} or
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null &amp;&amp; weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see MixtureTrainSM#MixtureTrainSM(int, de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel[], int,
	 *      boolean, double[], double[],
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest)
	 * @see AbstractMixtureTrainSM.Algorithm#EM
	 */
	public MixtureTrainSM( int length, TrainableStatisticalModel[] models, double[] weights, int starts, double alpha, TerminationCondition tc,
							Parameterization parametrization ) throws IllegalArgumentException, WrongAlphabetException,
																CloneNotSupportedException {
		this( length, models, starts, false, null, weights, Algorithm.EM, alpha, tc, parametrization, 0, 0, null );
	}

	/**
	 * Creates an instance using Gibbs Sampling and sampling the component
	 * probabilities.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link MixtureTrainSM}, if the
	 *            model is trained using
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param componentHyperParams
	 *            the hyperparameters for the component assignment prior
	 *            <ul>
	 *            <li>will only be used if
	 *            <code>estimateComponentProbs == true</code>
	 *            <li>the array has to be <code>null</code> or has to have
	 *            length <code>models.length</code>
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null &amp;&amp; weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see MixtureTrainSM#MixtureTrainSM(int, de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel[], int,
	 *      boolean, double[], double[],
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest)
	 * @see AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING
	 */
	public MixtureTrainSM( int length, TrainableStatisticalModel[] models, int starts, double[] componentHyperParams, int initialIteration,
							int stationaryIteration, BurnInTest burnInTest ) throws IllegalArgumentException, WrongAlphabetException,
																			CloneNotSupportedException {
		this( length,
				models,
				starts,
				true,
				componentHyperParams,
				null,
				Algorithm.GIBBS_SAMPLING,
				0,
				null,
				Parameterization.LAMBDA,
				initialIteration,
				stationaryIteration,
				burnInTest );
	}

	/**
	 * Creates an instance using Gibbs Sampling and fixed component
	 * probabilities.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link MixtureTrainSM}, if the
	 *            model is trained using
	 *            {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param weights
	 *            <code>null</code> or the weights for the components (than
	 *            <code>weights.length == models.length</code>)
	 * @param initialIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null &amp;&amp; weights.length != dimension</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see MixtureTrainSM#MixtureTrainSM(int, de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel[], int,
	 *      boolean, double[], double[],
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest)
	 * @see AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING
	 */
	public MixtureTrainSM( int length, TrainableStatisticalModel[] models, double[] weights, int starts, int initialIteration, int stationaryIteration,
							BurnInTest burnInTest ) throws IllegalArgumentException, WrongAlphabetException, CloneNotSupportedException {
		this( length,
				models,
				starts,
				false,
				null,
				weights,
				Algorithm.GIBBS_SAMPLING,
				0,
				null,
				Parameterization.LAMBDA,
				initialIteration,
				stationaryIteration,
				burnInTest );
	}

	/**
	 * The constructor for the interface {@link de.jstacs.Storable}. Creates a
	 * new {@link MixtureTrainSM} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation of the model as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} is not parsable
	 */
	public MixtureTrainSM( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#emitDataSetUsingCurrentParameterSet(int, int[])
	 */
	@Override
	protected Sequence[] emitDataSetUsingCurrentParameterSet( int n, int... lengths ) throws Exception {
		int[] numbers = new int[dimension];
		Arrays.fill( numbers, 0 );
		int counter = 0, no = 0, k = 0;
		// sample how many sequences each model should generate
		for( ; no < n; no++ ) {
			numbers[AbstractMixtureTrainSM.draw( weights, 0 )]++;
		}

		no = 0;
		DataSet help;
		Sequence[] seqs = new Sequence[n];
		if( length == 0 ) {
			// homogenous case
			for( ; counter < dimension; counter++ ) {
				if( numbers[counter] > 0 ) {
					if( lengths.length == 1 ) {
						help = model[counter].emitDataSet( n, lengths );
					} else {
						int[] array = new int[numbers[counter]];
						System.arraycopy( lengths, k, array, 0, numbers[counter] );
						help = model[counter].emitDataSet( n, array );
					}
					for( k = 0; k < help.getNumberOfElements(); k++ ) {
						seqs[no] = help.getElementAt( k );
					}
				}
			}
		} else {
			// inhomogenous case
			if( lengths == null || lengths.length == 0 ) {
				// System.out.println( Arrays.toString( weights ) );
				// System.out.println( Arrays.toString( numbers ) );

				// generate sequences
				for( ; counter < dimension; counter++ ) {
					if( numbers[counter] > 0 ) {
						if( model[counter].getLength() == 0 ) {
							help = model[counter].emitDataSet( numbers[counter], length );
						} else {
							help = model[counter].emitDataSet( numbers[counter], lengths );
						}
						for( no = 0; no < numbers[counter]; no++, k++ ) {
							seqs[k] = help.getElementAt( no );
						}
					}
				}
			} else {
				throw new Exception( "This is an inhomogeneous model. Please check parameter lengths." );
			}
		}
		return seqs;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#doFirstIteration(double[], de.jstacs.utils.random.MultivariateRandomGenerator, de.jstacs.utils.random.MRGParams[])
	 */
	@Override
	protected double[][] doFirstIteration( double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception {
		int counter1, counter2, d = sample[0].getNumberOfElements();
		double[][] seqweights = createSeqWeightsArray();
		double[] w = new double[dimension];
		initWithPrior( w );
		double[] help = new double[dimension];
		if( dataWeights == null ) {
			for( counter1 = 0; counter1 < d; counter1++ ) {
				help = m.generate( dimension, params[counter1] );
				for( counter2 = 0; counter2 < dimension; counter2++ ) {
					seqweights[counter2][counter1] = help[counter2];
					w[counter2] += help[counter2];
				}
			}
		} else {
			for( counter1 = 0; counter1 < d; counter1++ ) {
				help = m.generate( dimension, params[counter1] );
				for( counter2 = 0; counter2 < dimension; counter2++ ) {
					seqweights[counter2][counter1] = dataWeights[counter1] * help[counter2];
					w[counter2] += seqweights[counter2][counter1];
				}
			}
		}
		getNewParameters( 0, seqweights, w );
		return seqweights;
	}

	/**
	 * This method enables you to train a mixture model with a fixed start
	 * partitioning. This is useful to compare implementations or if one has a
	 * hypothesis how the components should look like.
	 * 
	 * @param data
	 *            the data set of sequences
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the data set
	 * @param partitioning
	 *            a kind of partitioning
	 *            <ol>
	 *            <li> <code>partitioning.length</code> has to be
	 *            <code>data.getNumberofElements()</code>
	 *            <li>for all i: <code>partitioning[i].length</code> has to be
	 *            <code>getNumberOfModels()</code>
	 *            <li>{@latex.inline $\\forall i:\\;\\sum_j partitioning[i][j] \\stackrel{!}{=}1$}
	 *            </ol>
	 * 
	 * @return the weighting array used to initialize, this array can be reused
	 *         in the following iterations
	 * 
	 * @throws Exception
	 *             if something went wrong or if the number of components is 1
	 */
	public double[][] doFirstIteration( DataSet data, double[] dataWeights, double[][] partitioning ) throws Exception {
		setTrainData( data );
		if( dimension > 1 ) {
			int counter1, counter2, d = data.getNumberOfElements();
			double[][] seqweights = createSeqWeightsArray();
			double[] w = new double[dimension];
			initWithPrior( w );
			double sum;
			for( counter1 = 0; counter1 < d; counter1++ ) {
				if( partitioning[counter1].length != dimension ) {
					throw new IllegalArgumentException( "The partitioning for sequence " + counter1 + " was wrong. (number of parts)" );
				}
				sum = 0;
				for( counter2 = 0; counter2 < dimension; counter2++ ) {
					if( partitioning[counter1][counter2] < 0 || partitioning[counter1][counter2] > 1 ) {
						throw new IllegalArgumentException( "The partitioning for sequence " + counter1
															+ " was wrong. (part "
															+ counter2
															+ "was incorrect)" );
					}
					seqweights[counter2][counter1] = ( ( dataWeights == null ) ? 1d : dataWeights[counter1] ) * partitioning[counter1][counter2];
					sum += partitioning[counter1][counter2];
					w[counter2] += seqweights[counter2][counter1];
				}
				if( sum != 1 ) {
					throw new IllegalArgumentException( "The partitioning for sequence " + counter1 + " was wrong. (sum of parts not 1)" );
				}
			}
			getNewParameters( 0, seqweights, w );
			return seqweights;
		} else {
			throw new OperationNotSupportedException();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#getLogProbUsingCurrentParameterSetFor(int, de.jstacs.data.Sequence, int, int)
	 */
	@Override
	protected double getLogProbUsingCurrentParameterSetFor( int component, Sequence s, int start, int end ) throws Exception {
		return logWeights[component] + model[component].getLogProbFor( s, start, end );
	}

	/* 
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString(NumberFormat nf) {
		StringBuffer sb = new StringBuffer( model.length * 100000 );
		sb.append( "Mixture model with parameter estimation by " + getNameOfAlgorithm() + ": \n" );
		sb.append( "number of starts:\t" + starts + "\n" );
		switch( algorithm ) {
			case EM:
				for( int i = 0; i < dimension; i++ ) {
					sb.append( nf.format( weights[i] ) + "\t" + model[i].getInstanceName() + "\n" + model[i].toString(nf) + "\n" );
				}
				break;
			case GIBBS_SAMPLING:
				sb.append( "burn in test              :\t" + burnInTest.getInstanceName() + "\n" );
				sb.append( "length of stationary phase:\t" + stationaryIteration + "\n" );

				sb.append( "Mixture model components:\n" );
				for( int i = 0; i < dimension; i++ ) {
					sb.append( ( i + 1 ) + ". component: " + model[i].getInstanceName() + "\n" );
				}
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		return sb.toString();
	}

	/**
	 * Computes sequence weights and returns the score.
	 */
	@Override
	protected double getNewWeights( double[] dataWeights, double[] w, double[][] seqweights ) throws Exception {
		double L = 0, currentWeight = 1;
		int counter1, counter2 = 0;
		Sequence seq;
		initWithPrior( w );
		double[] help = new double[dimension];

		for( counter1 = 0; counter1 < seqweights[0].length; counter1++ ) {
			seq = sample[0].getElementAt( counter1 );
			if( dataWeights != null ) {
				currentWeight = dataWeights[counter1];
			}
			for( counter2 = 0; counter2 < dimension; counter2++ ) {
				help[counter2] = model[counter2].getLogProbFor( seq ) + logWeights[counter2];
			}
			L += modifyWeights( help ) * currentWeight;
			for( counter2 = 0; counter2 < dimension; counter2++ ) {
				seqweights[counter2][counter1] = help[counter2] * currentWeight;
				w[counter2] += seqweights[counter2][counter1];
			}
		}
		return L;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#setTrainData(de.jstacs.data.DataSet)
	 */
	@Override
	protected void setTrainData( DataSet data ) {
		sample = new DataSet[]{ data };
	}
}
