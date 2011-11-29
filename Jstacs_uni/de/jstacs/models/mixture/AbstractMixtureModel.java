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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.Model;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sampling.GibbsSamplingModel;
import de.jstacs.sampling.SamplingComponent;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.RealTime;
import de.jstacs.utils.SafeOutputStream;
import de.jstacs.utils.Time;
import de.jstacs.utils.ToolBox;
import de.jstacs.utils.random.DirichletMRG;
import de.jstacs.utils.random.DirichletMRGParams;
import de.jstacs.utils.random.FastDirichletMRGParams;
import de.jstacs.utils.random.MRGParams;
import de.jstacs.utils.random.MultivariateRandomGenerator;
import de.jstacs.utils.random.SoftOneOfN;
import de.jtem.numericalMethods.calculus.specialFunctions.Gamma;

/**
 * This is the abstract class for all kinds of mixture models. It enables the
 * user to train the parameters using {@link Algorithm#EM} or
 * {@link Algorithm#GIBBS_SAMPLING}. If this instance is trained using
 * {@link Algorithm#GIBBS_SAMPLING} the internal models that will be adjusted
 * have to implement {@link SamplingComponent}. If you use Gibbs sampling
 * temporary files will be created in the Java temp folder. These files will be
 * deleted if no reference to the current instance exists and the Garbage
 * Collector is called. Therefore it is recommended to call the Garbage
 * Collector explicitly at the end of any application.
 * 
 * <br>
 * <br>
 * 
 * The model stores a reference to the last sample used in <code>train</code>.
 * This enables the user to estimate the parameters iteratively beginning with
 * the current set of parameters. Therefore you can use the method
 * {@link AbstractMixtureModel#continueIterations(double[], double[][], int, int)}
 * </code>.
 * 
 * <br>
 * <br>
 * 
 * The method {@link AbstractMixtureModel#setOutputStream(OutputStream)} enables
 * the user to get comments from the
 * {@link AbstractMixtureModel#train(Sample, double[])} method or to repress
 * them.
 * 
 * <br>
 * <br>
 * 
 * The method {@link AbstractMixtureModel#getScoreForBestRun()} enables the user
 * to optimize different instances of the same model (
 * {@link AbstractMixtureModel#clone()}) using the EM-algorithm on different
 * CPUs, to compare the results and to select the best trained model. This might
 * be useful to get the results faster (measured in real time).
 * 
 * <br>
 * <br>
 * 
 * <b>The reference to the internal sample is not stored if the model is stored
 * in a {@link StringBuffer}. So you can use these methods only after training
 * the parameters after (re)creating a model.</b>
 * 
 * @author Jens Keilwagen, Berit Haldemann
 * 
 * @see de.jstacs.sampling.SamplingComponent
 * @see System#gc()
 */
public abstract class AbstractMixtureModel extends AbstractModel {

	//private static final int NO_OF_UNSURE_DECIMAL_PLACES = 3;

	/**
	 * This <code>enum</code> defines the different types of parameterization
	 * for a probability that can be used in an {@link AbstractMixtureModel}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Parameterization {
		/**
		 * This value indicates that the component probabilities will be
		 * parameterized as {@latex.inline $\\theta_c = p(c)$}.
		 */
		THETA( -1 ),

		/**
		 * This value indicates that the component probabilities will be
		 * parameterized as {@latex.inline $\\lambda_c = \\log p(c)$}.
		 */
		LAMBDA( 0 );

		private double count;

		Parameterization( double count ) {
			this.count = count;
		}

		double getCount() {
			return count;
		}

	}

	/**
	 * This <code>enum</code> defines the different types of algorithms that can
	 * be used in an {@link AbstractMixtureModel}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static enum Algorithm {
		/**
		 * This value indicates that the model should be learned using the
		 * expectation maximization.
		 */
		EM,

		/**
		 * This value indicates that the model should be learned using Gibbs
		 * sampling.
		 */
		GIBBS_SAMPLING;
	}

	/**
	 * The probabilities for each component.
	 */
	protected double[] weights;

	/**
	 * The log probabilities for each component.
	 */
	protected double[] logWeights;

	/**
	 * The hyperparameters for estimating the probabilities of the components.
	 */
	protected double[] componentHyperParams;
	
	protected double ess;

	/**
	 * The model for the sequences.
	 */
	protected Model[] model;

	/**
	 * The alternative models for the EM.
	 */
	protected Model[] alternativeModel;

	/**
	 * The number of starts.
	 */
	protected int starts;

	/**
	 * The number of dimensions.
	 */
	protected int dimension;

	/**
	 * This field contains the value of objective function of the best start of the training. 
	 */
	protected double best;

	/**
	 * This is the stream for writing information while training.
	 */
	protected SafeOutputStream sostream;

	/**
	 * The sample that was used in the last training. Will not be stored in the
	 * {@link StringBuffer} when invoking {@link #toXML()}.
	 */
	protected Sample[] sample;

	/**
	 * The switch for estimating the component probabilities or not.
	 */
	protected boolean estimateComponentProbs;

	/**
	 * A switch for each model whether to optimize/adjust or not.
	 */
	protected boolean[] optimizeModel;

	/**
	 * The type of algorithm.
	 */
	protected Algorithm algorithm;

	/**
	 * A switch which indicates that the algorithm for determining the
	 * parameters has been run.
	 */
	protected boolean algorithmHasBeenRun;

	//EM
	/**
	 * The type of parameterization.
	 */
	private Parameterization parametrization;

	private double alpha = 1;
	private TerminationCondition tc;

	//GIBBS_SAMPLING

	/**
	 * The number of initial iterations.
	 */
	protected int initialIteration;

	/**
	 * The number of (stationary) iterations of the Gibbs Sampler.
	 */
	protected int stationaryIteration;

	/**
	 * The {@link BurnInTest} that is used to stop the sampling.
	 */
	protected BurnInTest burnInTest;

	/**
	 * Saving component probabilities in a file.
	 */
	protected BufferedWriter filewriter;

	/**
	 * Reading component probabilities from a file.
	 */
	protected BufferedReader filereader;

	/**
	 * The file in which the component probabilities are stored.
	 */
	protected File[] file;

	/**
	 * The current index of the parameter set while adjustment (optimization).
	 */
	protected int[] counter;

	/**
	 * The current index of the sampling.
	 */
	protected int samplingIndex;

	/**
	 * This array is used while training to avoid creating many new objects.
	 */
	protected double[] compProb;

	private double[][][] usedWeights;

	/**
	 * The weights of the (sub-)sequence used to train the components (internal models). The first dimension is used for the models, the second for the (sub-)sequences.
	 */
	protected double[][] seqWeights;

	/**
	 * Creates a new {@link AbstractMixtureModel}. This constructor can be used
	 * for any algorithm since it takes all necessary values as parameters.
	 * 
	 * @param length
	 *            the length used in this model
	 * @param models
	 *            the single models building the {@link AbstractMixtureModel},
	 *            if the model is trained using {@link Algorithm#GIBBS_SAMPLING}
	 *            the models that will be adjusted have to implement
	 *            {@link SamplingComponent}
	 * @param optimizeModel
	 *            an array of switches to determine whether a model should be
	 *            optimized or not
	 * @param dimension
	 *            the number of components
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
	 *            length <code>dimension</code>
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            then ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == dimension</code>)
	 * @param algorithm
	 *            either {@link Algorithm#EM} or
	 *            {@link Algorithm#GIBBS_SAMPLING}
	 * @param alpha
	 *            only for {@link Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link Parameterization#THETA} or
	 *            {@link Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use {@link Parameterization#LAMBDA}
	 *            <ul>
	 * @param initialIteration
	 *            only for {@link Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used for approximation
	 * @param burnInTest
	 *            only for {@link Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of length
	 *             <code>length</code> <li><code>dimension &lt; 1</code> <li>
	 *             <code>weights != null && weights.length != dimension</code>
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
	protected AbstractMixtureModel( int length, Model[] models, boolean[] optimizeModel, int dimension, int starts,
									boolean estimateComponentProbs, double[] componentHyperParams, double[] weights, Algorithm algorithm,
									double alpha, TerminationCondition tc, Parameterization parametrization, //EM parameters
									int initialIteration, int stationaryIteration, BurnInTest burnInTest )   //GIBBS_SAMPLING parameters
																											throws CloneNotSupportedException,
																											IllegalArgumentException,
																											WrongAlphabetException {
		super( models[0].getAlphabetContainer(), length );
		if( dimension < 1 ) {
			throw new IllegalArgumentException( "The dimension has to be at least 1." );
		}
		this.dimension = dimension;
		set( ArrayHandler.clone( models ),
				optimizeModel,
				starts,
				weights,
				estimateComponentProbs,
				componentHyperParams,
				algorithm,
				alpha,
				tc,
				parametrization,
				initialIteration,
				stationaryIteration,
				burnInTest != null ? burnInTest.clone() : null );
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
		algorithmHasBeenRun = false;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link AbstractMixtureModel} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation of the model as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	protected AbstractMixtureModel( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#clone()
	 */
	@Override
	public AbstractMixtureModel clone() throws CloneNotSupportedException {
		try {
			AbstractMixtureModel clone = (AbstractMixtureModel)super.clone();
			clone.weights = null;

			clone.set( ArrayHandler.clone( model ),
					optimizeModel,
					starts,
					weights,
					estimateComponentProbs,
					componentHyperParams,
					algorithm,
					alpha,
					tc,
					parametrization,
					initialIteration,
					stationaryIteration,
					burnInTest != null ? burnInTest.clone() : null );

			if( file != null ) {
				clone.counter = counter.clone();
				clone.file = new File[file.length];
				try {
					for( int i = 0; i < file.length; i++ ) {
						if( file[i] != null ) {
							clone.file[i] = File.createTempFile( "pi-", ".dat", null );
							FileManager.copy( file[i].getAbsolutePath(), clone.file[i].getAbsolutePath() );
						}
					}
				} catch ( IOException e ) {
					CloneNotSupportedException c = new CloneNotSupportedException( e.getMessage() );
					c.setStackTrace( e.getStackTrace() );
					throw c;
				}
			}

			// as long as parseParameterSet(int) and parseNextParameterSet() are protected methods
			clone.filereader = null;
			clone.filewriter = null;

			clone.setOutputStream( sostream.doesNothing() ? null : SafeOutputStream.DEFAULT_STREAM );
			clone.best = best;
			return clone;
		} catch ( IllegalArgumentException e ) {
			throw getCloneNotSupportedException( e );
		} catch ( WrongAlphabetException e ) {
			throw getCloneNotSupportedException( e );
		}
	}

	/**
	 * Only used in {@link #clone()} to throw a
	 * {@link CloneNotSupportedException}.
	 * 
	 * @param e
	 *            the original {@link Exception}
	 * 
	 * @return a {@link CloneNotSupportedException} containing all informations
	 */
	private static CloneNotSupportedException getCloneNotSupportedException( Exception e ) {
		CloneNotSupportedException ex = new CloneNotSupportedException( "impossible Exception in method clone in class AbstractMixtureModel: " + e.getMessage() );
		ex.setStackTrace( e.getStackTrace() );
		return ex;
	}

	/**
	 * This method creates the multivariate random generator that will be used
	 * during initialization.
	 * 
	 * @return a multivariate random generator
	 * 
	 * @see AbstractMixtureModel#getMRGParams()
	 */
	protected MultivariateRandomGenerator getMRG() {
		switch( algorithm ) {
			case EM:
				return DirichletMRG.DEFAULT_INSTANCE;
			case GIBBS_SAMPLING:
				return new SoftOneOfN();
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * This method creates the parameters used in a multivariate random
	 * generator while initialization.
	 * 
	 * @return the parameters for the multivariate random generator
	 * 
	 * @see AbstractMixtureModel#getMRG()
	 */
	protected MRGParams getMRGParams() {
		switch( algorithm ) {
			case EM:
				return new FastDirichletMRGParams( alpha );
			case GIBBS_SAMPLING:
				return null;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#train(de.jstacs.data.Sample, double[])
	 */
	public void train( Sample data, double[] dataWeights ) throws Exception {
		setNewAlphabetContainerInstance( data.getAlphabetContainer() );
		sample = null;
		System.gc();
		setTrainData( data );
		MultivariateRandomGenerator rg = getMRG();
		MRGParams[] params = new MRGParams[data.getNumberOfElements()];
		Arrays.fill( params, getMRGParams() );
		int i;
		switch( algorithm ) {
			case EM:
				double max = Double.NEGATIVE_INFINITY,
				current;
				double[] p = weights.clone();
				if( alternativeModel == null ) {
					alternativeModel = ArrayHandler.clone( model );
					for( i = 0; i < model.length; i++ ) {
						alternativeModel[i].setNewAlphabetContainerInstance( alphabets );
					}
				}
				for( i = 0; i < starts; i++ ) {
					current = iterate( i, dataWeights, rg, params );
					if( max < current ) {
						// swap models, ...
						swap();
						// swap weights
						p = weights.clone();
						// set new value
						max = current;
					}
				}
				swap();
				setWeights( p );
				p = null;
				best = max;
				sostream.writeln( "best = " + max );
				break;
			case GIBBS_SAMPLING:
				burnInTest.resetAllValues();
				initModelForSampling( starts );
				for( i = 0; i < starts; i++ ) {
					current = iterate( i, dataWeights, rg, params );
				}
				int anz,
				m,
				burnIn;
				boolean finished; //have all samplings passed the burn-in phase?
				do {
					burnIn = burnInTest.getLengthOfBurnIn();
					anz = m = 0;
					finished = true;
					for( i = 0; i < starts; i++ ) {
						m = counter[i] - burnIn;
						if( m > 0 ) {
							anz += m;
						} else {
							finished = false;
						}
					}
					//compute how many steps shall be done 
					anz = (int)Math.ceil( ( stationaryIteration - anz ) / (double)starts );

					if( anz > 0 ) {
						//go on with sampling
						for( i = 0; i < starts; i++ ) {
							sostream.writeln( "=== extend start: " + i + " ==========" );
							continueIterations( dataWeights, seqWeights, anz, i );
						}
					}
				} while( !finished );
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		rg = null;
		params = null;

		System.gc();
	}

	/**
	 * This method swaps the current component models with the alternative
	 * model.
	 * 
	 * <br>
	 * <br>
	 * 
	 * This method should <b>NOT</b> be made public and should <b>ONLY</b> be
	 * used in the <code>train</code>-method.
	 */
	protected void swap() {
		Model[] helpM = alternativeModel;
		alternativeModel = model;
		model = helpM;
	}

	/**
	 * This method is invoked by the <code>train</code>-method and sets for a
	 * given sample the sample that should be used for <code>train</code>.
	 * 
	 * @param data
	 *            the given sample of sequences
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected abstract void setTrainData( Sample data ) throws Exception;

	/**
	 * Creates an array that can be used for weighting sequences in the
	 * algorithm.
	 * 
	 * @return an array that can be used for weighting sequences in the
	 *         algorithm
	 */
	protected double[][] createSeqWeightsArray() {
		return new double[model.length][sample[0].getNumberOfElements()];
	}

	/**
	 * This method runs the train algorithm for the current model.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param dataWeights
	 *            the weights for each sequence or <code>null</code>
	 * @param m
	 *            the random generator for initiating the algorithm
	 * @param params
	 *            the parameters for the sequences
	 * 
	 * @return the score
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractMixtureModel#doFirstIteration(Sample, double[],
	 *      MultivariateRandomGenerator, MRGParams[])
	 * @see AbstractMixtureModel#continueIterations(double[], double[][])
	 * @see AbstractMixtureModel#continueIterations(double[], double[][], int,
	 *      int)
	 */
	public double iterate( Sample data, double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception {
		sample = null;
		System.gc();
		setTrainData( data );
		return iterate( 0, dataWeights, m, params );
	}

	/**
	 * This method runs the train algorithm for the current model and the
	 * internal data set.
	 * 
	 * @param start
	 *            the index of the training
	 * @param dataWeights
	 *            the weights for each sequence or <code>null</code>
	 * @param m
	 *            the random generator for initiating the algorithm
	 * @param params
	 *            the parameters for the sequences
	 * 
	 * @return the score
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see AbstractMixtureModel#doFirstIteration(Sample, double[],
	 *      MultivariateRandomGenerator, MRGParams[])
	 * @see AbstractMixtureModel#continueIterations(double[], double[][])
	 * @see AbstractMixtureModel#continueIterations(double[], double[][], int,
	 *      int)
	 */
	protected double iterate( int start, double[] dataWeights, MultivariateRandomGenerator m, MRGParams params[] ) throws Exception {
		sostream.writeln( "========== start: " + start + " ==========" );
		switch( algorithm ) {
			case EM:
				best = continueIterations( dataWeights, doFirstIteration( dataWeights, m, params ) );
				break;
			case GIBBS_SAMPLING:
				extendSampling( start );
				burnInTest.setCurrentSamplingIndex( start );
				seqWeights = doFirstIteration( dataWeights, m, params );
				samplingStopped();
				continueIterations( dataWeights, seqWeights, initialIteration, start );
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		algorithmHasBeenRun = true;
		return best;
	}

	/**
	 * This method will do the first step in the train algorithm for the current
	 * model. The initialization will be done by randomly setting the component
	 * membership. This is useful when nothing is known about the problem.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the sample
	 * 
	 * @return the weighting array used to initialize, this array can be reused
	 *         in the following iterations
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected double[][] doFirstIteration( Sample data, double[] dataWeights ) throws Exception {
		FastDirichletMRGParams[] params = new FastDirichletMRGParams[data.getNumberOfElements()];
		Arrays.fill( params, new FastDirichletMRGParams( alpha ) );
		return doFirstIteration( data, dataWeights, DirichletMRG.DEFAULT_INSTANCE, params );
	}

	/**
	 * This method will do the first step in the train algorithm for the current
	 * model. The initialization will be done by randomly setting the component
	 * membership. This is useful when nothing is known about the problem.
	 * 
	 * @param data
	 *            the sample of sequences
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the sample
	 * @param m
	 *            the multivariate random generator
	 * @param params
	 *            the parameters for the multivariate random generator
	 * 
	 * @return the weighting array used to initialize, this array can be reused
	 *         in the following iterations
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected double[][] doFirstIteration( Sample data, double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception {
		sample = null;
		System.gc();
		setTrainData( data );
		return doFirstIteration( dataWeights, m, params );
	}

	/**
	 * This method will do the first step in the train algorithm for the current
	 * model on the internal sample. The initialization will be done by randomly
	 * setting the component membership. This is useful when nothing is known
	 * about the problem.
	 * 
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the sample
	 * @param m
	 *            the multivariate random generator
	 * @param params
	 *            the parameters for the multivariate random generator
	 * 
	 * @return the weighting array used to initialize, this array can be reused
	 *         in the following iterations
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected abstract double[][] doFirstIteration( double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception;

	/**
	 * This method returns an instance of {@link Time} that is used for the {@link TerminationCondition} in the EM.
	 * 
	 * @return an instance of {@link Time} that is used during the EM
	 */
	protected Time getTime() {
		return new RealTime();
	}
	
	/**
	 * This method will run the train algorithm for the current model on the
	 * internal sample. The initialization will be done by using the models of
	 * the {@link AbstractMixtureModel}. So in this case the models have to be
	 * trained already. This method is useful for restarting the train algorithm
	 * at a certain point. The algorithm will stop if the difference between the
	 * optimized functions for two iterations is smaller than the specified
	 * threshold.
	 * 
	 * <br>
	 * <br>
	 * 
	 * If the difference becomes significant negative an exception is thrown.
	 * 
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the
	 *            internal sample (last sample the {@link AbstractMixtureModel}
	 *            was trained on)
	 * @param seqweights
	 *            <code>null</code> or an array for weighting the sequences, see
	 *            {@link #createSeqWeightsArray()}
	 * 
	 * @return a score for the model
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected double continueIterations( double[] dataWeights, double[][] seqweights ) throws Exception {
		if( sample == null ) {
			throw new OperationNotSupportedException( "There is no reference to an internal sample, so you can not go on with training." );
		}
		int i = 0;
		double[] w = new double[dimension];
		if( seqweights == null ) {
			seqweights = createSeqWeightsArray();
		}
		Time t = getTime();
		double pr = getLogPriorTerm(), L_old = Double.NEGATIVE_INFINITY, L_new = getNewWeights( dataWeights, w, seqweights );
		sostream.write( i + "\t" + L_new + "\t " + pr + "\t" );
		L_new += pr;
		sostream.writeln( L_new + "\t" + ( L_new - L_old ) );
		while( tc.doNextIteration(i, L_old, L_new, null, null, Double.NaN, t) ) {
			getNewParameters( ++i, seqweights, w );
			L_old = L_new;
			pr = getLogPriorTerm(); // in ML-case this should be 0
			L_new = getNewWeights( dataWeights, w, seqweights );
			sostream.write( i + "\t" + t.getElapsedTime() + "\t" + L_new + "\t " + pr + "\t" );
			L_new += pr;
			sostream.writeln( L_new + "\t" + ( L_new - L_old ) );
		}
		/*
		// if( sostream.doesNothing() ) System.out.println( "iteration: " + i + " \t" + L_new );
		if( L_new - L_old < 0 ) {
			// check negative abort
			// determine (approximately) the number of (possible) decimal places for the score (i) and
			// the smallest index of the decimal places of the difference that is not zero (d)
			String s = "" + L_new;
			int n = s.indexOf( "." );

			if( n >= 0 ) {
				n = 17 - n;
			} else {
				n = 17 - s.length();
			}

			L_old -= L_new;
			int d;
			if( L_old >= 1 ) {
				d = -(int)Math.ceil( Math.log( L_old ) / Math.log( 10 ) );
			} else {
				d = -(int)Math.floor( Math.log( L_old ) / Math.log( 10 ) );
			}

			if( n - NO_OF_UNSURE_DECIMAL_PLACES >= d ) {
				throw new EvaluationException( "The score decreases after " + i
												+ " steps! decimal places = "
												+ n
												+ ", delta = "
												+ ( -L_old ) );
			} else {
				s = "Negative abort after " + i + " steps caused by the limited precision of double.";
				if( sostream.doesNothing() ) {
					System.out.println( s );
				} else {
					sostream.writeln( s );
				}
			}
		}
		*/
		return L_new;
	}

	/**
	 * This method will run the train algorithm for the current model on the
	 * internal sample. The initialization will be done by using the models of
	 * the {@link AbstractMixtureModel}. So in this case the models have to be
	 * trained already. This method is useful for restarting the algorithm at a
	 * certain point. The algorithm will stop after the number of iterations.
	 * 
	 * @param dataWeights
	 *            <code>null</code> or the weights of each element of the
	 *            internal sample (last sample the {@link AbstractMixtureModel}
	 *            was trained on)
	 * @param seqweights
	 *            <code>null</code> or an array for weighting the sequences, see
	 *            {@link #createSeqWeightsArray()}
	 * @param iterations
	 *            the number of iterations that should be done
	 * @param start
	 *            the index of the run in a {@link Model#train(Sample)}-call
	 * 
	 * @return the current score (likelihood or posterior)
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected double continueIterations( double[] dataWeights, double[][] seqweights, int iterations, int start ) throws Exception {
		if( burnInTest != null ) {
			extendSampling( start );
			burnInTest.setCurrentSamplingIndex( start );
		}
		if( sample == null ) {
			throw new OperationNotSupportedException( "There is no reference to an internal sample, so you can not go on with training." );
		}
		double[] w = new double[dimension];
		if( seqweights == null ) {
			seqweights = createSeqWeightsArray();
		}
		double pr = getLogPriorTerm(), L_old = Double.NEGATIVE_INFINITY, L_new = getNewWeights( dataWeights, w, seqweights );
		for( int i = 0, j = burnInTest == null ? 0 : counter[samplingIndex]; i < iterations; i++, j++ ) {
			sostream.write( j + "\t" + L_new + "\t " + pr + "\t" );
			L_new += pr;
			sostream.writeln( L_new + "\t" + ( L_new - L_old ) );

			if( burnInTest != null ) {
				burnInTest.setValue( L_new );
			}

			getNewParameters( i, seqweights, w );
			L_old = L_new;
			pr = getLogPriorTerm(); // in ML-case this should be 0
			L_new = getNewWeights( dataWeights, w, seqweights );
		}
		if( burnInTest != null ) {
			samplingStopped();
		}
		return L_new+pr;
	}

	/**
	 * This method trains the internal models on the internal sample and the
	 * given weights.
	 * 
	 * @param iteration
	 *            the number of times this method has been invoked
	 * @param seqWeights
	 *            the weights for each model and sequence
	 * @param w
	 *            the weights for the components
	 * 
	 * @throws Exception
	 *             if the training of the internal models went wrong
	 */
	protected void getNewParameters( int iteration, double[][] seqWeights, double[] w ) throws Exception {
		for( int i = 0; i < seqWeights.length; i++ ) {
			getNewParametersForModel( i, iteration, 0, seqWeights[i] );
		}
		getNewComponentProbs( w );
	}

	/**
	 * This method trains the internal model with index <code>modelIndex</code>
	 * on the internal sample and the given weights.
	 * 
	 * @param modelIndex
	 *            the index of the model
	 * @param iteration
	 *            the number of times this method has been invoked for this
	 *            model
	 * @param sampleIndex
	 *            the index of the internal sample that should be used
	 * @param seqWeights
	 *            the weights for each sequence
	 * 
	 * @throws Exception
	 *             if the training of the internal model went wrong
	 */
	protected void getNewParametersForModel( int modelIndex, int iteration, int sampleIndex, double[] seqWeights ) throws Exception {
		if( optimizeModel[modelIndex] ) {
			switch( algorithm ) {
				case EM:
					if( model[modelIndex] instanceof AbstractMixtureModel ) {
						//runtime efficiency:
						//this allows us to use AbstractMixtureModels as internal models of an AbstractMixtureModel in an EM
						if( iteration == 0 ) {
							usedWeights[modelIndex] = ( (AbstractMixtureModel)model[modelIndex] ).doFirstIteration( sample[sampleIndex],
									seqWeights );
						} else {
							//TODO improve!?
							//maybe it is better to do more than one step
							( (AbstractMixtureModel)model[modelIndex] ).continueIterations( seqWeights, usedWeights[modelIndex], 1, 0 );
						}
					} else {
						model[modelIndex].train( sample[sampleIndex], seqWeights );
					}
					break;
				case GIBBS_SAMPLING:
					( (GibbsSamplingModel)model[modelIndex] ).drawParameters( sample[sampleIndex], seqWeights );
					( (GibbsSamplingModel)model[modelIndex] ).acceptParameters();
					break;
				default:
					throw new IllegalArgumentException( "The type of algorithm is unknown." );
			}
		}
	}

	/**
	 * Computes sequence weights and returns the score.
	 * 
	 * @param dataWeights
	 *            the weights for the internal sample (should not be changed)
	 * @param w
	 *            the array for the statistic of the component parameters (shall
	 *            be filled)
	 * @param seqweights
	 *            an array containing for each component the weights for each
	 *            sequence (shall be filled)
	 * 
	 * @return the score
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	protected abstract double getNewWeights( double[] dataWeights, double[] w, double[][] seqweights ) throws Exception;

	/**
	 * This method modifies the computed weights for one sequence and returns
	 * the score.
	 * 
	 * @param w
	 *            the weights
	 * 
	 * @return the score
	 */
	protected double modifyWeights( double[] w ) {
		switch( algorithm ) {
			case EM:
				return Normalisation.logSumNormalisation( w, 0, w.length, w, 0 );
			case GIBBS_SAMPLING:
				double l = Normalisation.logSumNormalisation( w, 0, w.length, w, 0 );
				int index = AbstractMixtureModel.draw( w, 0 );
				//double log_P_Of_X_and_U = l + Math.log( w[index] ); // likelihood of the sequence for the drawn component
				Arrays.fill( w, 0 );
				w[index] = 1;
				return l;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * This method sets the initial weights before counting the usage of each
	 * component. For ML the weights are set to 0 and for MAP they are set to
	 * the component hyperparameters.
	 * 
	 * @param w
	 *            the array of weights
	 */
	protected void initWithPrior( double[] w ) {
		System.arraycopy( componentHyperParams, 0, w, 0, dimension );
	}

	/**
	 * Returns the logarithmic probability for the sequence and the given
	 * component.
	 * 
	 * @param component
	 *            the index of the component
	 * @param s
	 *            the sequence
	 * 
	 * @return
	 *         <code>log P(s,component) = log P(s|component) + log P(component)</code>
	 * 
	 * @throws Exception
	 *             if the model was not trained yet or something else went wrong
	 * 
	 * @see AbstractMixtureModel#getNumberOfComponents()
	 */
	public double getLogProbFor( int component, Sequence s ) throws Exception {
		switch( algorithm ) {
			case EM:
				return getLogProbUsingCurrentParameterSetFor( component, s, 0, s.getLength() - 1 );
			case GIBBS_SAMPLING:
				//TODO label switching?
				throw new OperationNotSupportedException();
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * Returns the logarithmic probability for the sequence and the given
	 * component using the current parameter set.
	 * 
	 * @param component
	 *            the index of the component
	 * @param s
	 *            the sequence
	 * @param start
	 *            the start position in the sequence
	 * @param end
	 *            the end position in the sequence
	 * 
	 * @return
	 *         <code>log P(s,component) = log P(s|component) + log P(component)</code>
	 * 
	 * @throws Exception
	 *             if not trained yet or something else went wrong
	 * 
	 * @see AbstractMixtureModel#getNumberOfComponents()
	 */
	protected abstract double getLogProbUsingCurrentParameterSetFor( int component, Sequence s, int start, int end ) throws Exception;

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogProbFor(de.jstacs.data.Sequence, int, int)
	 */
	@Override
	public final double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws Exception {
		if( !isInitialized() ) {
			throw new NotTrainedException();
		}
		switch( algorithm ) {
			case EM:
				for( int i = 0; i < dimension; i++ ) {
					compProb[i] = getLogProbUsingCurrentParameterSetFor( i, sequence, startpos, endpos );
				}
				return Normalisation.getLogSum( compProb );
			case GIBBS_SAMPLING:
				int i,
				anz = 0,
				sampling = 0,
				burnIn = burnInTest.getLengthOfBurnIn();
				double res = Double.NEGATIVE_INFINITY;
				boolean b;
				for( ; sampling < starts; sampling++ ) {
					b = parseParameterSet( sampling, burnIn );
					while( b ) {
						for( i = 0; i < dimension; i++ ) {
							compProb[i] = getLogProbUsingCurrentParameterSetFor( i, sequence, startpos, endpos );
						}
						res = Normalisation.getLogSum( res, Normalisation.getLogSum( compProb ) );
						b = parseNextParameterSet();
						anz++;
					}
				}
				return res - Math.log( anz );
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getLogProbFor(de.jstacs.data.Sample)
	 */
	@Override
	public final double[] getLogScoreFor( Sample data ) throws Exception {
		if( !isInitialized() ) {
			throw new NotTrainedException();
		}
		switch( algorithm ) {
			case EM:
				return super.getLogScoreFor( data );
			case GIBBS_SAMPLING:
				int i,
				anz = 0,
				sampling = 0,
				k,
				burnIn = burnInTest.getLengthOfBurnIn();
				Sequence[] sequence = data.getAllElements();
				double[] res = new double[sequence.length];
				Arrays.fill( res, Double.NEGATIVE_INFINITY );
				boolean b;
				for( ; sampling < starts; sampling++ ) {
					b = parseParameterSet( sampling, burnIn );
					while( b ) {
						for( k = 0; k < sequence.length; k++ ) {
							for( i = 0; i < dimension; i++ ) {
								compProb[i] = getLogProbUsingCurrentParameterSetFor( i, sequence[k], 0, sequence[k].getLength() - 1 );
							}
							res[k] = Normalisation.getLogSum( res[k], Normalisation.getLogSum( compProb ) );
						}
						b = parseNextParameterSet();
						anz++;
					}
				}

				double d = Math.log( anz );
				for( k = 0; k < sequence.length; k++ ) {
					res[k] -= d;
				}
				return res;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getLogPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception {
		switch( algorithm ) {
			case GIBBS_SAMPLING:
				//TODO what should be returned???
				return 0;
			case EM:
				double erg = 0;
				for( int counter = 0; counter < model.length; counter++ ) {
					if( optimizeModel[counter] ) {
						erg += model[counter].getLogPriorTerm();
					}
				}
				return erg + getLogPriorTermForComponentProbs();
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * This method computes the part of the prior that comes from the component
	 * probabilities.
	 * 
	 * @return the part of the prior that comes from the component probabilities
	 */
	protected final double getLogPriorTermForComponentProbs() {
		double prior = 0, sum = 0;
		if( estimateComponentProbs && componentHyperParams[0] > 0 ) {
			for( int counter = 0; counter < dimension; counter++ ) {
				sum += componentHyperParams[counter];
				prior += ( componentHyperParams[counter] + parametrization.getCount() ) * logWeights[counter]
							- Gamma.logOfGamma( componentHyperParams[counter] );
			}
			prior += Gamma.logOfGamma( sum );
		}
		return prior;
	}

	/**
	 * Returns the value of the optimized function from the best run of the last
	 * training.
	 * 
	 * @return the value of the optimized function from the best run of the last
	 *         training
	 * 
	 * @throws NotTrainedException
	 *             if the training algorithm has not been run
	 * @throws OperationNotSupportedException
	 *             if this method is used for an instance that does not use the
	 *             EM
	 * 
	 * @see AbstractMixtureModel#train(Sample, double[])
	 * @see AbstractMixtureModel#algorithmHasBeenRun()
	 */
	public final double getScoreForBestRun() throws NotTrainedException, OperationNotSupportedException {
		if( algorithmHasBeenRun() ) {
			if( algorithm == Algorithm.EM ) {
				return best;
			} else {
				throw new OperationNotSupportedException();
			}
		} else {
			throw new NotTrainedException();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getInstanceName()
	 */
	public String getInstanceName() {
		StringBuffer erg = new StringBuffer( getClass().getSimpleName() + "(" );
		erg.append( model[0].getInstanceName() );
		for( int i = 1; i < model.length; i++ ) {
			erg.append( ", " );
			erg.append( model[i].getInstanceName() );
		}
		if( !estimateComponentProbs ) {
			erg.append( "; " + Arrays.toString( weights ) );
		}
		erg.append( ") " + getNameOfAlgorithm() );
		return erg.toString();
	}

	/**
	 * Returns the index <code>i</code> of the component with
	 * <code>P(i|s)<code> maximal. Therefore it computes 
	 * {@latex.ilb %preamble{\\usepackage{amsmath}}
	 * \\[i = \\operatorname{argmax}_j P(j|s)
	 * 	= \\operatorname{argmax}_j P(s|j) \\cdot P(j).\\]} 
	 * This method can be helpful for clustering.
	 * 
	 * @param s
	 *            the sequence
	 * 
	 * @return the index of the component
	 * 
	 * @throws Exception
	 *             if the model was not trained yet or something else went wrong
	 * 
	 * @see AbstractMixtureModel#getLogProbFor(int, Sequence)
	 */
	public int getIndexOfMaximalComponentFor( Sequence s ) throws Exception {
		switch( algorithm ) {
			case EM:
				double best = getLogProbFor( 0, s ),
				current;
				int index = 0,
				i = 1;
				while( i < dimension ) {
					current = getLogProbFor( i, s );
					if( current > best ) {
						best = current;
						index = i;
					}
					i++;
				}
				return index;
			case GIBBS_SAMPLING:
				//TODO see AbstractMixtureModel#getLogProbFor(int, Sequence)
				throw new OperationNotSupportedException();
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * Returns a deep copy of the models.
	 * 
	 * @return an array of {@link AbstractModel}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one model can not be cloned
	 * 
	 * @see AbstractMixtureModel#getModel(int)
	 */
	public final Model[] getModels() throws CloneNotSupportedException {
		return ArrayHandler.clone( model );
	}

	/**
	 * Returns a deep copy of the <code>i</code>-th model.
	 * 
	 * @param i
	 *            the index
	 * 
	 * @return a deep copy of the <code>i</code>-th model
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one model can not be cloned
	 * 
	 * @see AbstractMixtureModel#getModels()
	 */
	public final Model getModel( int i ) throws CloneNotSupportedException {
		return model[i].clone();
	}

	/**
	 * Returns the name of the used algorithm.
	 * 
	 * @return the name of the used algorithm
	 */
	public String getNameOfAlgorithm() {
		switch( algorithm ) {
			case EM:
				return "EM";
			case GIBBS_SAMPLING:
				return "Gibbs Sampling";
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * Returns the number of components the are modeled by this
	 * {@link AbstractMixtureModel}.
	 * 
	 * @return the number of components
	 */
	public final int getNumberOfComponents() {
		return dimension;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#getCharacteristics()
	 */
	@Override
	public ResultSet getCharacteristics() throws Exception {
		LinkedList<Result> infos = new LinkedList<Result>();
		ResultSet part;
		int i = 0, j;
		for( ; i < model.length; i++ ) {
			part = model[i].getCharacteristics();
			if( part != null && part.getNumberOfResults() > 0 ) {
				infos.add( new NumericalResult( "model number", "type of model " + model[i].getClass().getSimpleName(), new Integer( i ) ) );
				for( j = 0; j < part.getNumberOfResults(); j++ ) {
					infos.add( part.getResultAt( j ) );
				}
			}
		}
		infos.add( new StorableResult( "model", "the xml representation of the model", this ) );
		return new ResultSet( infos );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		LinkedList<NumericalResult> infos = new LinkedList<NumericalResult>();
		NumericalResultSet part;
		int i = 0, j;
		for( ; i < model.length; i++ ) {
			part = model[i].getNumericalCharacteristics();
			if( part != null && part.getNumberOfResults() > 0 ) {
				infos.add( new NumericalResult( "model number", "type of model " + model[i].getClass().getSimpleName(), new Integer( i ) ) );
				for( j = 0; j < part.getNumberOfResults(); j++ ) {
					infos.add( part.getResultAt( j ) );
				}
			}
		}
		return new NumericalResultSet( infos );
	}

	/**
	 * This method returns a deep copy of the weights for each component.
	 * 
	 * @return the weight for each component
	 */
	public final double[] getWeights() {
		return weights.clone();
	}

	/**
	 * This method indicates whether the parameters of the model has been
	 * determined by the internal algorithm.
	 * 
	 * @return <code>true</code> if the internal algorithm has been used to
	 *         determine the parameters of the model
	 */
	public boolean algorithmHasBeenRun() {
		return algorithmHasBeenRun;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.Model#isTrained()
	 */
	public boolean isInitialized() {
		switch( algorithm ) {
			case EM:
				int i = 0;
				while( i < model.length && model[i].isInitialized() ) {
					i++;
				}
				return i == model.length;
			case GIBBS_SAMPLING:
				return algorithmHasBeenRun;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/**
	 * Sets the parameter of the Dirichlet distribution which is used when you
	 * invoke <code>train</code> to init the gammas. It is recommended to use
	 * <code>alpha = 1</code> (uniform distribution on a simplex).
	 * 
	 * @param alpha
	 *            the parameter of the Dirichlet distribution with
	 *            <code>alpha &gt; 0</code>
	 * 
	 * @throws IllegalArgumentException
	 *             if <code>alpha  &lt;= 0</code>
	 */
	public final void setAlpha( double alpha ) throws IllegalArgumentException {
		if( alpha <= 0 ) {
			throw new IllegalArgumentException( "alpha has to be strict positive." );
		}
		this.alpha = alpha;
	}

	/**
	 * Sets the {@link OutputStream} that is used e.g. for writing information
	 * while training. It is possible to set <code>o=null</code>, than nothing
	 * will be written.
	 * 
	 * @param o
	 *            the {@link OutputStream}
	 */
	public final void setOutputStream( OutputStream o ) {
		sostream = SafeOutputStream.getSafeOutputStream( o );
	}

	/**
	 * Estimates the weights of each component.
	 * 
	 * @param weights
	 *            the array of weights, every element has to be non-negative and
	 *            the dimension has to be <code>dimension</code>
	 * 
	 * @throws Exception
	 *             a weight is less than 0
	 * 
	 * @see AbstractMixtureModel#getNumberOfComponents()
	 */
	protected void getNewComponentProbs( double[] weights ) throws Exception {
		if( estimateComponentProbs ) {
			double sum = 0;
			int i = 0;
			switch( algorithm ) {
				case EM:
					boolean map = componentHyperParams[0] != 0;
					for( ; i < dimension; i++ ) {
						if( map ) {
							weights[i] += parametrization.getCount();
						}
						sum += weights[i];
						if( weights[i] < 0 ) {
							throw new IllegalArgumentException( "Every weight has to be at least 0. Violate at position " + i + "." );
						}
					}
					for( i = 0; i < dimension; i++ ) {
						this.weights[i] = weights[i] / sum;
					}
					break;
				case GIBBS_SAMPLING:
					DirichletMRG.DEFAULT_INSTANCE.generate( this.weights, 0, dimension, new DirichletMRGParams( weights ) );

					filewriter.write( counter[samplingIndex] + "\t" );
					for( ; i < dimension; i++ ) {
						filewriter.write( this.weights[i] + "\t" );
					}
					filewriter.write( "\n" );
					filewriter.flush();
					break;
				default:
					throw new IllegalArgumentException( "The type of algorithm is unknown." );
			}
			for( i = 0; i < dimension; i++ ) {
				this.logWeights[i] = Math.log( this.weights[i] );
			}
		}

		if( algorithm == Algorithm.GIBBS_SAMPLING ) {
			counter[samplingIndex]++;
		}
	}

	/**
	 * Sets the weights of each component.
	 * 
	 * @param weights
	 *            every element has to be non-negative, the sum of all weights
	 *            has to be 1 and the dimension of <code>weights</code> has to
	 *            be <code>dimension</code>
	 * 
	 * @throws IllegalArgumentException
	 *             a weight is less than 0, the sum is not equal to 1 or the
	 *             dimension is incorrect
	 * 
	 * @see AbstractMixtureModel#getNumberOfComponents()
	 */
	protected void setWeights( double... weights ) throws IllegalArgumentException {
		if( weights.length != dimension ) {
			throw new IllegalArgumentException( "The number of weights is incorrect" );
		}
		double sum = 0;
		int i = 0;
		for( ; i < dimension; i++ ) {
			sum += weights[i];
			if( weights[i] < 0 ) {
				throw new IllegalArgumentException( "Every weight has to be at least 0. Violate at position " + i + "." );
			}
		}
		if( Math.abs( 1d - sum ) > 1E-9 ) {
			throw new IllegalArgumentException( "The weights do not sum to 1." );
		} else {
			if( this.weights == null ) {
				this.weights = new double[dimension];
				this.logWeights = new double[dimension];
			}
			for( i = 0; i < dimension; i++ ) {
				this.weights[i] = weights[i];
				this.logWeights[i] = Math.log( this.weights[i] );
			}
		}
	}

	public double getESS(){
		return ess;
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( xml, length, "length" );
		XMLParser.appendObjectWithTags( xml, dimension, "dimension" );
		XMLParser.appendObjectWithTags( xml, starts, "starts" );
		XMLParser.appendObjectWithTags( xml, estimateComponentProbs, "estimateComponentProbs" );
		XMLParser.appendObjectWithTags( xml, componentHyperParams, "componentHyperParams" );
		XMLParser.appendObjectWithTags( xml, model, "models" );
		XMLParser.appendObjectWithTags( xml, optimizeModel, "optimizeModel" );
		XMLParser.appendObjectWithTags( xml, algorithmHasBeenRun, "algorithmHasBeenRun" );
		XMLParser.appendObjectWithTags( xml, weights, "weights" );

		//algorithm specific values
		XMLParser.appendObjectWithTags( xml, algorithm, "algorithm" );
		switch( algorithm ) {
			case EM:
				XMLParser.appendObjectWithTags( xml, alpha, "alpha" );
				XMLParser.appendObjectWithTags( xml, tc, "terminationCondition" );
				XMLParser.appendObjectWithTags( xml, parametrization, "parametrization" );
				break;
			case GIBBS_SAMPLING:
				XMLParser.appendObjectWithTags( xml, initialIteration, "initialIteration" );
				XMLParser.appendObjectWithTags( xml, stationaryIteration, "stationaryIteration" );
				XMLParser.appendObjectWithTags( xml, burnInTest, "burnInTest" );
				XMLParser.appendObjectWithTags( xml, file != null, "hasParameterFiles" );
				if( file != null ) {
					XMLParser.appendObjectWithTags( xml, counter, "counter" );
					try {
						String content;
						for( int i = 0; i < counter.length; i++ ) {
							if( file[i] != null ) {
								content = FileManager.readFile( file[i] ).toString();
							} else {
								content = "";
							}
							XMLParser.appendObjectWithTagsAndAttributes( xml, content, "fileContent", "pos=\"" + i + "\"" );
						}
					} catch ( IOException e ) {
						RuntimeException r = new RuntimeException( e.getMessage() );
						r.setStackTrace( e.getStackTrace() );
						throw r;
					}
				}
				break;
		}

		XMLParser.appendObjectWithTags( xml, best, "best" );
		xml.append( getFurtherInformation() );
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
	}

	/**
	 * This method is used in the subclasses to append further information to
	 * the XML representation.
	 * 
	 * @return a part of the XML representation
	 * 
	 * @see AbstractMixtureModel#extractFurtherInformation(StringBuffer)
	 */
	protected StringBuffer getFurtherInformation() {
		return new StringBuffer( 1 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML( StringBuffer representation ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( representation, getClass().getSimpleName() );
		length = XMLParser.extractObjectForTags( xml, "length", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		dimension = XMLParser.extractObjectForTags( xml, "dimension", int.class );
		starts = XMLParser.extractObjectForTags( xml, "starts", int.class );
		estimateComponentProbs = XMLParser.extractObjectForTags( xml, "estimateComponentProbs", boolean.class );
		componentHyperParams = XMLParser.extractObjectForTags( xml, "componentHyperParams", double[].class );
		model = XMLParser.extractObjectForTags( xml, "models", Model[].class );
		optimizeModel = XMLParser.extractObjectForTags( xml, "optimizeModel", boolean[].class );
		algorithmHasBeenRun = XMLParser.extractObjectForTags( xml, "algorithmHasBeenRun", boolean.class );
		double[] w = XMLParser.extractObjectForTags( xml, "weights", double[].class );

		algorithm = XMLParser.extractObjectForTags( xml, "algorithm", Algorithm.class );
		try {
			switch( algorithm ) {
				case EM:
					parametrization = XMLParser.extractObjectForTags( xml, "parametrization", Parameterization.class );
					if( XMLParser.hasTag(xml, "epsilon", null, null) ) {
						tc = new SmallDifferenceOfFunctionEvaluationsCondition( XMLParser.extractObjectForTags( xml, "epsilon", double.class ) );
					} else {
						tc = XMLParser.extractObjectForTags( xml, "terminationCondition", TerminationCondition.class );
					}
					set( model, optimizeModel, starts, w, estimateComponentProbs, componentHyperParams, algorithm,
					//EM
							XMLParser.extractObjectForTags( xml, "alpha", double.class ),
							tc,
							parametrization,
							//Gibbs Sampling
							0,
							0,
							null );
					break;
				case GIBBS_SAMPLING:
					set( model, optimizeModel, starts, w, estimateComponentProbs, componentHyperParams, algorithm,
					//EM
							0d,
							null,
							Parameterization.LAMBDA,
							//Gibbs Sampling
							XMLParser.extractObjectForTags( xml, "initialIteration", int.class ),// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
							XMLParser.extractObjectForTags( xml, "stationaryIteration", int.class ),
							XMLParser.extractObjectForTags( xml, "burnInTest", BurnInTest.class ) );
					if( XMLParser.extractObjectForTags( xml, "hasParameterFiles", boolean.class ) ) {
						counter = XMLParser.extractObjectForTags( xml, "counter", int[].class );
						file = new File[counter.length];
						try {
							String content;
							Map<String,String> filter = new TreeMap<String, String>();
							for( int i = 0; i < counter.length; i++ ) {
								filter.clear();
								filter.put( "pos", ""+i );
								content = XMLParser.extractObjectAndAttributesForTags( xml, "fileContent", null, filter, String.class );
								if( !content.equalsIgnoreCase( "" ) ) {
									file[i] = File.createTempFile( "pi-", ".dat", null );
									FileManager.writeFile( file[i], new StringBuffer( content ) );
								}
							}
						} catch ( IOException e ) {
							NonParsableException r = new NonParsableException( e.getMessage() );
							r.setStackTrace( e.getStackTrace() );
							throw r;
						}
					} else {
						file = null;
					}
					break;
				default:
					throw new IllegalArgumentException( "The type of algorithm is unknown." );
			}
		} catch ( Exception e ) {
			NonParsableException n = new NonParsableException( e.getMessage() );
			n.setStackTrace( e.getStackTrace() );
			throw n;
		}

		best = XMLParser.extractObjectForTags( xml, "best", double.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		alphabets = model[0].getAlphabetContainer();
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
		extractFurtherInformation( xml );
	}

	/**
	 * This method is used in the subclasses to extract further information from
	 * the XML representation and to set these as values of the instance.
	 * 
	 * @param xml
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML representation is not parsable
	 * 
	 * @see AbstractMixtureModel#getFurtherInformation()
	 */
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#set(de.jstacs.data.AlphabetContainer)
	 */
	@Override
	protected void set( AlphabetContainer abc ) {
		for( int i = 0; i < model.length; i++ ) {
			model[i].setNewAlphabetContainerInstance( abc );
		}
	}

	/**
	 * This method is used in the constructor and in the methods
	 * {@link #clone()} and {@link #fromXML(StringBuffer)} to set all necessary
	 * values.
	 */
	private void set( Model[] model, boolean[] optimizeModel, int starts, double[] weights, boolean estimateComponentProbs,
			double[] componentHyperParams, Algorithm algorithm, double alpha, TerminationCondition tc, Parameterization parametrization,
			int initialIteration, int stationaryIteration, BurnInTest burnInTest ) throws IllegalArgumentException, WrongAlphabetException {
		int i = 0;
		if( starts < 1 ) {
			throw new IllegalArgumentException( "The number of iterations has to be at least 1." );
		}
		this.starts = starts;
		AlphabetContainer abc = model[0].getAlphabetContainer();
		for( ; i < model.length; i++ ) {
			if( i != 0 && !model[i].setNewAlphabetContainerInstance( abc ) ) {
				throw new WrongAlphabetException( "The models have to have the same alphabet like the AbstractMixtureModel. Violated at position " + i
													+ "." );
			}
			if( model[i] instanceof AbstractMixtureModel ) {
				( (AbstractMixtureModel)model[i] ).setOutputStream( null );
			}
			checkLength( i, model[i].getLength() );
		}
		if( optimizeModel == null ) {
			this.optimizeModel = new boolean[model.length];
			Arrays.fill( this.optimizeModel, true );
		} else {
			if( optimizeModel.length != model.length ) {
				throw new IllegalArgumentException( "The dimension of the switch whether the individual models should be optimized/adjusted has wrong dimension." );
			} else {
				this.optimizeModel = new boolean[model.length];
				System.arraycopy( optimizeModel, 0, this.optimizeModel, 0, optimizeModel.length );
			}
		}

		if( weights == null ) {
			weights = new double[dimension];
			Arrays.fill( weights, 1d / dimension );
		}
		setWeights( weights );

		this.model = model;
		this.alternativeModel = null;
		this.estimateComponentProbs = estimateComponentProbs;

		if(componentHyperParams == null){
			this.ess = 0;
		}else{
			this.ess = ToolBox.sum( componentHyperParams );
		}
		
		boolean minValueOfUsedHyperParamIsZero;
		if( !estimateComponentProbs || componentHyperParams == null ) {
			this.componentHyperParams = new double[dimension];
			if( !estimateComponentProbs ) {
				minValueOfUsedHyperParamIsZero = false;
			} else {
				minValueOfUsedHyperParamIsZero = true;
			}
		} else {
			if( componentHyperParams.length != dimension ) {
				throw new IllegalArgumentException( "The dimension of the component assignment hyperparameter is not correct." );
			}
			this.componentHyperParams = new double[dimension];
			minValueOfUsedHyperParamIsZero = componentHyperParams[0] == 0;
			for( i = 0; i < dimension; i++ ) {
				if( componentHyperParams[i] < 0 || ( minValueOfUsedHyperParamIsZero && componentHyperParams[i] > 0 )
					|| ( !minValueOfUsedHyperParamIsZero && componentHyperParams[i] == 0 ) ) {
					throw new IllegalArgumentException( "The " + i + "-th component assignment hyperparameter is not correct." );
				}
				this.componentHyperParams[i] = componentHyperParams[i];
			}
		}
		best = Double.NEGATIVE_INFINITY;
		compProb = new double[dimension];
		usedWeights = new double[model.length][][];

		switch( algorithm ) {
			case EM:
				if( parametrization == Parameterization.THETA || parametrization == Parameterization.LAMBDA ) {
					this.parametrization = parametrization;
				} else {
					throw new IllegalArgumentException( "The type of parametrization is unknown." );
				}
				setAlpha( alpha );
				if( tc == null ) {
					throw new NullPointerException();
				}
				if( !tc.isSimple() ) {
					throw new IllegalArgumentException( "The TerminationCondition has to be simple." );
				}
				this.tc = tc;
				break;
			case GIBBS_SAMPLING:
				if( minValueOfUsedHyperParamIsZero ) {
					throw new IllegalArgumentException( "The component hyper parameters have to be set to positive values." );
				}
				if( initialIteration <= 0 ) {
					throw new IllegalArgumentException( "The given number of intial iterations has to be at least 1." );
				}
				if( initialIteration * starts > stationaryIteration ) {
					throw new IllegalArgumentException( "The given number of intial iterations has to be most (stationaryIteration/starts)." );
				}
				if( stationaryIteration <= 0 ) {
					throw new IllegalArgumentException( "The given number of iterations has to be at least 1." );
				}
				if( burnInTest == null ) {
					throw new IllegalArgumentException( "You have to specify a burn in test." );
				}
				this.initialIteration = initialIteration;
				this.stationaryIteration = stationaryIteration;
				this.burnInTest = burnInTest;

				checkModelsForGibbsSampling();
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		this.algorithm = algorithm;
	}

	/**
	 * This method can be used to check whether the necessary models have
	 * implemented the {@link SamplingComponent}.
	 */
	protected void checkModelsForGibbsSampling() {
		for( int i = 0; i < model.length; i++ ) {
			if( optimizeModel[i] && !( model[i] instanceof GibbsSamplingModel ) ) {
				throw new IllegalArgumentException( "The model for component " + i
													+ " doesn't implement the interface GibbsSamplingComponent!" );
			}
		}
	}

	/**
	 * This method checks if the length <code>l</code> of the model with index
	 * <code>index</code> is capable for the current instance. Otherwise an
	 * {@link IllegalArgumentException} is thrown.
	 * 
	 * @param index
	 *            the index of the model
	 * @param l
	 *            the length of the model
	 * 
	 * @throws IllegalArgumentException
	 *             if the model instance can not be used
	 */
	protected void checkLength( int index, int l ) {
		if( l != 0 && length != l ) {
			throw new IllegalArgumentException( "The models have to use the same length like the AbstractMixtureModel. Violated at position " + index
												+ "." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.AbstractModel#emitSample(int, int[])
	 */
	@Override
	public Sample emitSample( int n, int... lengths ) throws Exception {
		if( !isInitialized() ) {
			throw new NotTrainedException();
		}
		Sequence[] seqs;
		switch( algorithm ) {
			case EM:
				seqs = emitSampleUsingCurrentParameterSet( n, lengths );
				break;
			case GIBBS_SAMPLING:
				int[] anz = new int[starts];
				int i = 0,
				all = 0,
				j,
				k,
				burnIn = burnInTest.getLengthOfBurnIn();
				for( ; i < starts; i++ ) {
					all = anz[i] = all + Math.max( 0, counter[i] - burnIn );
				}
				int[] no = new int[all],
				len;
				for( i = 0; i < n; i++ ) {
					no[r.nextInt( all )]++;
				}

				seqs = new Sequence[n];
				Sequence[] help;
				all = 0;
				for( i = 0, j = 0; i < starts; i++ ) {
					parseParameterSet( i, burnIn );
					for( ; all < anz[i]; all++ ) {
						if( no[all] > 0 ) {
							if( lengths == null || lengths.length <= 1 ) {
								len = lengths;
							} else {
								len = new int[no[all]];
								System.arraycopy( lengths, j, len, 0, no[all] );
							}
							help = emitSampleUsingCurrentParameterSet( no[all], len );
							for( k = 0; k < no[all]; k++, j++ ) {
								seqs[j] = help[k];
							}
						}
						parseNextParameterSet();
					}
				}
				break;
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
		return new Sample( "sampled from " + getInstanceName(), seqs );
	}

	/**
	 * The method returns an array of sequences using the current parameter set.
	 * 
	 * @param n
	 *            the number of sequences to be sampled
	 * @param lengths
	 *            the corresponding lengths
	 * 
	 * @return an array of sequences
	 * 
	 * @throws Exception
	 *             if it was impossible to sample the sequences
	 * 
	 * @see AbstractModel#emitSample(int, int...)
	 */
	protected abstract Sequence[] emitSampleUsingCurrentParameterSet( int n, int... lengths ) throws Exception;

	/**
	 * This method allows the user to parse the set of parameters with index
	 * <code>burnInIteration</code> of a specific <code>sampling</code> (from a
	 * file).
	 * 
	 * @param sampling
	 *            the index of the sampling
	 * @param burnInIteration
	 *            the number of iterations that should be skipped
	 * 
	 * @return <code>true</code> if the parameter set could be parsed
	 * 
	 * @throws Exception
	 *             if something went wrong while reading or parsing the
	 *             parameter set
	 */
	protected boolean parseParameterSet( int sampling, int burnInIteration ) throws Exception {
		boolean parsed = true;
		for( int i = 0; i < model.length; i++ ) {
			if( optimizeModel[i] ) {
				parsed &= ( (SamplingComponent)model[i] ).parseParameterSet( sampling, burnInIteration );
			}
		}
		parsed &= parseComponentParameterSet( sampling, burnInIteration );
		return parsed;
	}

	private boolean parseComponentParameterSet( int sampling, int burnInIteration ) throws IOException {
		if( filereader != null ) {
			filereader.close();
		}
		String str;
		filereader = new BufferedReader( new FileReader( file[sampling] ) );

		while( ( str = filereader.readLine() ) != null ) {
			if( Integer.parseInt( str.substring( 0, str.indexOf( "\t" ) ) ) == burnInIteration ) {
				parse( str );
				return true;
			}
		}
		return false;
	}

	private void parse( String str ) {
		String[] strarray = str.split( "\t" );
		for( int l = 1, i = 0; i < model.length; i++ ) {
			this.weights[i] = Double.parseDouble( strarray[l++] );
			this.logWeights[i] = Math.log( this.weights[i] );
		}
	}

	/**
	 * This method allows the user to parse the next set of parameters (from a
	 * file).
	 * 
	 * @return <code>true</code> if the parameter set could be parsed
	 * 
	 * @throws Exception
	 *             if something went wrong while reading or parsing the
	 *             parameter set
	 */
	protected boolean parseNextParameterSet() throws Exception {
		String str = filereader.readLine();
		if( str == null ) {
			return false;
		}
		parse( str );

		boolean parsed = true;
		for( int i = 0; i < model.length && parsed; i++ ) {
			if( optimizeModel[i] ) {
				parsed &= ( (SamplingComponent)model[i] ).parseNextParameterSet();
			}
		}
		return parsed;
	}

	/**
	 * This method initializes the model for the sampling. For instance this
	 * method can be used to create new files where all parameter sets will be
	 * stored.
	 * 
	 * @param starts
	 *            the number of sampling starts
	 * 
	 * @throws IOException
	 *             if the files could not be handled properly
	 */
	protected void initModelForSampling( int starts ) throws IOException {
		if( file != null && file.length == starts ) {
			FileOutputStream o;
			for( int i = 0; i < starts; i++ ) {
				if( file[i] != null ) {
					o = new FileOutputStream( file[i] );
					o.close();
				}
				counter[i] = 0;
			}
		} else {
			deleteParameterFiles();
			file = new File[starts];
			counter = new int[starts];
		}
		for( int i = 0; i < model.length; i++ ) {
			if( optimizeModel[i] ) {
				( (SamplingComponent)model[i] ).initForSampling( starts );
			}
		}
	}

	/**
	 * This method prepares the model to extend an existing sampling.
	 * 
	 * @param sampling
	 *            the index of the sampling
	 * 
	 * @throws Exception
	 *             if the internal files could not be handled properly
	 */
	protected void extendSampling( int sampling ) throws Exception {
		if( file[sampling] == null ) {
			file[sampling] = File.createTempFile( "pi-", ".dat", null );
		} else {
			parseComponentParameterSet( sampling, counter[sampling] - 1 );
			filereader.close();
			filereader = null;
		}
		filewriter = new BufferedWriter( new FileWriter( file[sampling], true ) );

		for( int i = 0; i < model.length; i++ ) {
			if( optimizeModel[i] ) {
				( (SamplingComponent)model[i] ).extendSampling( sampling, true );
			}
		}
		samplingIndex = sampling;
	}

	/**
	 * This method is the opposite of the method
	 * {@link #initModelForSampling(int)}. It can be used for closing any
	 * streams of writer, ...
	 * 
	 * @throws IOException
	 *             if the {@link FileWriter} could not be closed properly
	 */
	protected void samplingStopped() throws IOException {
		for( int i = 0; i < model.length; i++ ) {
			if( optimizeModel[i] ) {
				( (SamplingComponent)model[i] ).samplingStopped();
			}
		}
		filewriter.close();
		filewriter = null;
	}

	/**
	 * This method returns <code>true</code> if the object is currently used in
	 * a sampling, otherwise <code>false</code>.
	 * 
	 * @return <code>true</code> if the object is currently used in a sampling
	 */
	protected boolean isInSamplingMode() {
		int i = 0;
		while( i < model.length ) {
			if( !optimizeModel[i] || ( (SamplingComponent)model[i] ).isInSamplingMode() ) {
				i++;
			} else {
				break;
			}
		}
		return ( i == model.length ) && ( filewriter != null );
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		model = alternativeModel = null;
		weights = logWeights = componentHyperParams = compProb = null;
		sample = null;
		counter = null;
		optimizeModel = null;
		usedWeights = null;

		if( filereader != null ) {
			filereader.close();
		}
		if( filewriter != null ) {
			filewriter.close();
		}
		deleteParameterFiles();
		super.finalize();
	}

	private void deleteParameterFiles() {
		if( file != null ) {
			for( int i = 0; i < file.length; i++ ) {
				if( file[i] != null ) {
					file[i].delete();
				}
			}
		}
	}

	private final static Random r = new Random();

	/**
	 * This method draws an index of an array corresponding to the probabilities
	 * encoded in the entries of the array.
	 * 
	 * @param w
	 *            an array containing probabilities starting at position
	 *            <code>start</code>
	 * @param start
	 *            the start index
	 * 
	 * @return the drawn index
	 */
	public final static int draw( double[] w, int start ) {
		double p = r.nextDouble();
		int i = start;
		while( i < w.length && p > w[i] ) {
			p -= w[i++];
		}
		if( i == w.length ) {
			i--;
		}
		return i;
	}

	/**
	 * This method returns the index of a maximal entry in the array
	 * <code>w</code> between index <code>start</code> and <code>end</code>.
	 * 
	 * @param w
	 *            an array
	 * @param start
	 *            the start index (inclusive)
	 * @param end
	 *            the end index (exclusive)
	 * 
	 * @return the index of the maximal entry
	 */
	public final static int max( double[] w, int start, int end ) {
		int i = start + 1, max = start;
		for( ; i < end; i++ ) {
			if( w[i] > w[max] ) {
				max = i;
			}
		}
		return max;
	}
}
