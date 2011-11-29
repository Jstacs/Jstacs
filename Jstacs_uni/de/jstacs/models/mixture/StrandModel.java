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

package de.jstacs.models.mixture;

import java.util.Random;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.models.Model;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.utils.random.MRGParams;
import de.jstacs.utils.random.MultivariateRandomGenerator;

/**
 * This model handles sequences that can either lie on the forward strand or on
 * the reverse complementary strand. Therefore it is recommended to use this model only for
 * DNA, but it is not restricted to DNA.
 * 
 * <br>
 * <br>
 * 
 * If you use Gibbs Sampling temporary files will be created in the Java temp
 * folder. These files will be deleted if no reference to the current instance
 * exists and the Garbage Collector is called. Therefore it is recommended to
 * call the Garbage Collector explicitly at the end of any application.
 * 
 * @see de.jstacs.models.Model
 * 
 * @author Jens Keilwagen
 */
public class StrandModel extends AbstractMixtureModel {

	/**
	 * Creates a new {@link StrandModel}. This constructor can be used for any
	 * algorithm since it takes all necessary values as parameters.
	 * 
	 * @param model
	 *            the model building the basis of the {@link StrandModel}, if
	 *            the instance is trained using
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param estimateComponentProbs
	 *            the switch for estimating the component probabilities in the
	 *            algorithm or to hold them fixed; if the component parameters
	 *            are fixed, the value <code>forwardStrandProb</code> will be
	 *            used, otherwise the <code>componentHyperParams</code> will be
	 *            incorporated in the adjustment
	 * @param componentHyperParams
	 *            the hyperparameters for the component assignment prior
	 *            <ul>
	 *            <li>will only be used if
	 *            <code>estimateComponentProbs == true</code>
	 *            <li>the array has to be <code>null</code> or has to have
	 *            length 2
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param forwardStrandProb
	 *            the probability for the forward strand
	 * @param algorithm
	 *            either {@link AbstractMixtureModel.Algorithm#EM} or
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}
	 * @param alpha
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
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
	 *             <code>weights != null && weights.length != dimension</code>
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
	 */
	protected StrandModel( Model model, int starts, boolean estimateComponentProbs, double[] componentHyperParams,
							double forwardStrandProb, Algorithm algorithm, double alpha, TerminationCondition tc, Parameterization parametrization, //EM
							int initialIteration, int stationaryIteration, BurnInTest burnInTest ) //GIBBS_SAMPLING
																									throws CloneNotSupportedException,
																									IllegalArgumentException,
																									WrongAlphabetException {
		super( model.getLength(),
				new Model[]{ model },
				null,
				2,
				starts,
				estimateComponentProbs,
				componentHyperParams,
				new double[]{ forwardStrandProb, 1 - forwardStrandProb },
				algorithm,
				alpha,
				tc,
				parametrization,
				initialIteration,
				stationaryIteration,
				burnInTest );
		if( !alphabets.isReverseComplementable() ) {
			throw new WrongAlphabetException( "The given model uses an AlphabetContainer that can not be used for building a reverse complement." );
		}
	}

	/**
	 * Creates an instance using EM and estimating the component probabilities.
	 * 
	 * @param model
	 *            the model building the basis of the {@link StrandModel}, if
	 *            the instance is trained using
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
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
	 *            length 2
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param alpha
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the type of the component probability parameterization
	 *            <ul>
	 *            <li>{@link AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null && weights.length != dimension</code>
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
	 * @see StrandModel#StrandModel(de.jstacs.models.Model, int, boolean,
	 *      double[], double,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest )
	 * @see AbstractMixtureModel.Algorithm#EM
	 */
	public StrandModel( Model model, int starts, double[] componentHyperParams, double alpha, TerminationCondition tc, Parameterization parametrization )
																																			throws CloneNotSupportedException,
																																			IllegalArgumentException,
																																			WrongAlphabetException {
		this( model, starts, true, componentHyperParams, 0.5, Algorithm.EM, alpha, tc, parametrization, 0, 0, null );
	}

	/**
	 * Creates an instance using EM and fixed component probabilities.
	 * 
	 *@param model
	 *            the model building the basis of the {@link StrandModel}, if
	 *            the instance is trained using
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param forwardStrandProb
	 *            the probability for the forward strand
	 * @param alpha
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link AbstractMixtureModel.Algorithm#EM}<br>
	 *            the type of the component probability parameterization
	 *            <ul>
	 *            <li>{@link AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code>
	 *             <li><code>dimension &lt; 1</code>
	 *             <li>
	 *             <code>weights != null && weights.length != dimension</code>
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
	 * @see StrandModel#StrandModel(de.jstacs.models.Model, int, boolean,
	 *      double[], double,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest )
	 * @see AbstractMixtureModel.Algorithm#EM
	 */
	public StrandModel( Model model, int starts, double forwardStrandProb, double alpha, TerminationCondition tc, Parameterization parametrization )
																																		throws CloneNotSupportedException,
																																		IllegalArgumentException,
																																		WrongAlphabetException {
		this( model, starts, false, null, forwardStrandProb, Algorithm.EM, alpha, tc, parametrization, 0, 0, null );
	}

	/**
	 * Creates an instance using Gibbs Sampling and sampling the component
	 * probabilities.
	 * 
	 * @param model
	 *            the model building the basis of the {@link StrandModel}, if
	 *            the instance is trained using
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
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
	 *            length 2
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
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
	 *             <code>weights != null && weights.length != dimension</code>
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
	 * @see StrandModel#StrandModel(de.jstacs.models.Model, int, boolean,
	 *      double[], double,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest )
	 * @see AbstractMixtureModel.Algorithm#GIBBS_SAMPLING
	 */
	public StrandModel( Model model, int starts, double[] componentHyperParams, int initialIteration, int stationaryIteration,
						BurnInTest burnInTest ) throws CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException {
		this( model,
				starts,
				true,
				componentHyperParams,
				0.5,
				Algorithm.GIBBS_SAMPLING,
				0d,
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
	 *@param model
	 *            the model building the basis of the {@link StrandModel}, if
	 *            the instance is trained using
	 *            {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param forwardStrandProb
	 *            the probability for the forward strand
	 * @param initialIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
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
	 *             <code>weights != null && weights.length != dimension</code>
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
	 * @see StrandModel#StrandModel(de.jstacs.models.Model, int, boolean,
	 *      double[], double,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double,
	 *      TerminationCondition,
	 *      de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int,
	 *      int, de.jstacs.sampling.BurnInTest )
	 * @see AbstractMixtureModel.Algorithm#GIBBS_SAMPLING
	 */
	public StrandModel( Model model, int starts, double forwardStrandProb, int initialIteration, int stationaryIteration,
						BurnInTest burnInTest ) throws CloneNotSupportedException, IllegalArgumentException, WrongAlphabetException {
		this( model,
				starts,
				false,
				null,
				forwardStrandProb,
				Algorithm.GIBBS_SAMPLING,
				0d,
				null,
				Parameterization.LAMBDA,
				initialIteration,
				stationaryIteration,
				burnInTest );
	}

	/**
	 * The constructor for the interface {@link de.jstacs.Storable}. Creates a
	 * new {@link StrandModel} out of its XML representation.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} containing the XML representation of
	 *            the model
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} could not be parsed
	 */
	public StrandModel( StringBuffer stringBuff ) throws NonParsableException {
		super( stringBuff );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#setTrainData(de.jstacs.data.Sample)
	 */
	@Override
	public void setTrainData( Sample s ) throws Exception {
		int i = 0, n = s.getNumberOfElements();
		Sequence[] seq = new Sequence[2 * n];
		for( ; i < n; i++ ) {
			seq[2 * i] = s.getElementAt( i );
			seq[2 * i + 1] = seq[2 * i].reverseComplement();
		}
		sample = new Sample[]{ new Sample( "sample of both strands from " + s.getAnnotation(), seq ) };
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#doFirstIteration(double[], de.jstacs.utils.random.MultivariateRandomGenerator, de.jstacs.utils.random.MRGParams[])
	 */
	@Override
	protected double[][] doFirstIteration( double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception {
		int counter1, counter2, d = sample[0].getNumberOfElements();
		double[][] seqweights = createSeqWeightsArray();
		double[] w = new double[2];
		initWithPrior( w );
		d /= 2;
		if( dataWeights == null ) {
			for( counter1 = 0; counter1 < d; counter1++ ) {
				m.generate( seqweights[0], 2 * counter1, 2, params[counter1] );
				for( counter2 = 0; counter2 < 2; counter2++ ) {
					w[counter2] += seqweights[0][2 * counter1 + counter2];
				}
			}
		} else {
			double[] help = new double[2];
			for( counter1 = 0; counter1 < d; counter1++ ) {
				help = m.generate( 2, params[counter1] );
				for( counter2 = 0; counter2 < 2; counter2++ ) {
					seqweights[0][counter2 + 2 * counter1] = dataWeights[counter1] * help[counter2];
					w[counter2] += seqweights[0][counter2 + 2 * counter1];
				}
			}
		}
		getNewParameters( 0, seqweights, w );
		return seqweights;
	}

	/**
	 * Computes sequence weights and returns the score.
	 */
	@Override
	protected double getNewWeights( double[] dataWeights, double[] w, double[][] seqweights ) throws Exception {
		double L = 0, currentWeight = 1;
		int counter1 = 0, counter2 = 0;
		initWithPrior( w );
		double[] help = new double[2];

		while( counter1 < seqweights[0].length ) {
			if( dataWeights != null ) {
				currentWeight = dataWeights[counter2++];
			}
			help[0] = model[0].getLogProbFor( sample[0].getElementAt( counter1 ) ) + logWeights[0];
			help[1] = model[0].getLogProbFor( sample[0].getElementAt( counter1 + 1 ) ) + logWeights[1];

			L += modifyWeights( help ) * currentWeight;

			seqweights[0][counter1] = help[0] * currentWeight;
			w[0] += seqweights[0][counter1++];
			seqweights[0][counter1] = help[1] * currentWeight;
			w[1] += seqweights[0][counter1++];
		}
		return L;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer( model.length * 100000 );
		sb.append( "Strand model with parameter estimation by " + getNameOfAlgorithm() + ": \n" );
		sb.append( "number of starts:\t" + starts + "\n" );
		switch( algorithm ) {
			case EM:
				sb.append( weights[0] + "\tforward strand\n" );
				sb.append( weights[1] + "\tbackward strand\n\n" );
				sb.append( model[0].toString() );
				break;
			case GIBBS_SAMPLING:
				sb.append( "burn in test              :\t" + burnInTest.getInstanceName() + "\n" );
				sb.append( "length of stationary phase:\t" + stationaryIteration + "\n" );
				sb.append( "strand model component:" + model[0].getInstanceName() + "\n" );
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		return sb.toString();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#emitSampleUsingCurrentParameterSet(int, int[])
	 */
	@Override
	protected Sequence[] emitSampleUsingCurrentParameterSet( int n, int... lengths ) throws NotTrainedException, Exception {
		Sample nr = model[0].emitSample( n, lengths );
		Random r = new Random();
		Sequence[] seq = new Sequence[nr.getNumberOfElements()];
		for( int i = 0; i < seq.length; i++ ) {
			if( r.nextDouble() < weights[0] ) {
				seq[i] = nr.getElementAt( i );
			} else {
				seq[i] = nr.getElementAt( i ).reverseComplement();
			}
		}
		return seq;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#getLogProbUsingCurrentParameterSetFor(int, de.jstacs.data.Sequence, int, int)
	 */
	@Override
	protected double getLogProbUsingCurrentParameterSetFor( int component, Sequence s, int start, int end ) throws Exception {
		switch( component ) {
			case 0:
				return logWeights[0] + model[0].getLogProbFor( s, start, end );
			case 1:
				return logWeights[1] + model[0].getLogProbFor( s.reverseComplement(),
								( s.getLength() - end - 1 ),
								( s.getLength() - start - 1 ) );
			default:
				throw new IndexOutOfBoundsException( "component has to be in [0,1]; 0 = forward strand, 1 = backward strand" );
		}
	}
}
