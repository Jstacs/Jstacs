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

package de.jstacs.models.mixture.motif;

import java.util.Arrays;
import java.util.LinkedList;

import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.models.Model;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.positionprior.PositionPrior;
import de.jstacs.models.mixture.motif.positionprior.UniformPositionPrior;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.MRGParams;
import de.jstacs.utils.random.MultivariateRandomGenerator;

/**
 * This class enables the user to search for a single motif in a sequence. The
 * user is enabled to train the model either &quot;one occurrence per
 * sequence&quot; (=OOPS) or &quot;zero or one occurrence per sequence&quot;
 * (=ZOOPS).
 * 
 * <br>
 * <br>
 * 
 * If EM is used for training the parameters are trained in a MEME-like manner.
 * 
 * <br>
 * <br>
 * 
 * <b>Currently only EM is implemented.</b>
 * 
 * @author Jens Keilwagen
 */
public class SingleHiddenMotifMixture extends HiddenMotifMixture {

	private int[] refBgSample;
	private boolean trainOnlyMotifModel;
	private boolean correctPhaseShift;
	
	/**
	 * The order of the background model.
	 */
	protected byte bgMaxMarkovOrder;

	/**
	 * Creates a new {@link SingleHiddenMotifMixture}. This constructor can be
	 * used for any algorithm since it takes all necessary values as parameters.
	 * 
	 * @param motif
	 *            the motif model, if the model is trained using
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING} the
	 *            model has to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}.
	 * @param bg
	 *            the background model for the flanking sequences and for those
	 *            sequences that do not contain a binding site, if
	 *            <code>trainOnlyMotifModel == false</code> and
	 *            <code>algorithm == {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}</code>
	 *            the model has to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}.
	 *            The model has to be able to score sequences of arbitrary
	 *            length.
	 * @param trainOnlyMotifModel
	 *            a switch whether to train only the motif model
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param componentHyperParams
	 *            the hyperparameters for the component assignment prior
	 *            <ul>
	 *            <li>will only be used if
	 *            <code>estimateComponentProbs == true</code>
	 *            <li>the array has to be <code>null</code> or has to have
	 *            length <code>dimension</code>
	 *            <li><code>null</code> or an array with all values zero (0)
	 *            than ML
	 *            <li>otherwise (all values positive) a prior is used (MAP, MP,
	 *            ...)
	 *            <li>depends on the <code>parameterization</code>
	 *            </ul>
	 * @param weights
	 *            <code>null</code> or the weights for the components (then
	 *            <code>weights.length == dimension</code>)
	 * @param posPrior
	 *            this object determines the positional distribution that shall
	 *            be used
	 * @param algorithm
	 *            either {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM} or
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}
	 * @param alpha
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used in approximation
	 * @param burnInTest
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws CloneNotSupportedException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of the
	 *             corresponding length
	 *             <li><code>weights != null && weights.length != 2</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws IllegalArgumentException
	 *             if not all <code>models</code> work on the same simple
	 *             alphabet
	 * @throws WrongAlphabetException
	 *             if the <code>models</code> can not be cloned
	 */
	protected SingleHiddenMotifMixture( Model motif, Model bg, boolean trainOnlyMotifModel, int starts, double[] componentHyperParams,
										double[] weights, PositionPrior posPrior, Algorithm algorithm, double alpha, TerminationCondition tc,
										Parameterization parametrization, int initialIteration, int stationaryIteration,
										BurnInTest burnInTest ) throws CloneNotSupportedException, IllegalArgumentException,
																WrongAlphabetException {
		super( new Model[]{ motif, bg },
				new boolean[]{true,!trainOnlyMotifModel},
				2,
				starts,
				weights == null,
				componentHyperParams,
				weights,
				posPrior == null ? new UniformPositionPrior() : posPrior,
				algorithm,
				alpha,
				tc,
				parametrization,
				initialIteration,
				stationaryIteration,
				burnInTest );
		bgMaxMarkovOrder = model[1].getMaximalMarkovOrder();
		this.trainOnlyMotifModel = trainOnlyMotifModel;
		this.correctPhaseShift = true;
	}

	/**
	 * Creates a new {@link SingleHiddenMotifMixture} using EM and estimating
	 * the probability for finding a motif.
	 * 
	 * @param motif
	 *            the motif model
	 * @param bg
	 *            the background model for the flanking sequences and for those
	 *            sequences that do not contain a binding site. The model has to
	 *            be able to score sequences of arbitrary length.
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
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
	 * @param posPrior
	 *            this object determines the positional distribution that shall
	 *            be used
	 * @param trainOnlyMotifModel
	 *            a switch whether to train only the motif model
	 * @param alpha
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            the type of the component probability parameterization
	 *            <ul>
	 *            <li>{@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of the
	 *             corresponding length
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same simple
	 *             alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see SingleHiddenMotifMixture#SingleHiddenMotifMixture(de.jstacs.models.Model, de.jstacs.models.Model, boolean, int, double[], double[], de.jstacs.models.mixture.motif.positionprior.PositionPrior, de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double, de.jstacs.algorithms.optimization.termination.TerminationCondition,  de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int, int, de.jstacs.sampling.BurnInTest)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM
	 */
	public SingleHiddenMotifMixture( Model motif, Model bg, boolean trainOnlyMotifModel, int starts, double[] componentHyperParams,
										PositionPrior posPrior, double alpha, TerminationCondition tc, Parameterization parametrization )
																															throws CloneNotSupportedException,
																															IllegalArgumentException,
																															WrongAlphabetException {
		this( motif,
				bg,
				trainOnlyMotifModel,
				starts,
				componentHyperParams,
				null,
				posPrior,
				Algorithm.EM,
				alpha,
				tc,
				parametrization,
				0,
				0,
				null );
	}

	/**
	 * Creates a new {@link SingleHiddenMotifMixture} using EM and fixed
	 * probability for finding a motif.
	 * 
	 * @param motif
	 *            the motif model
	 * @param bg
	 *            the background model for the flanking sequences and for those
	 *            sequences that do not contain a binding site. The model has to
	 *            be able to score sequences of arbitrary length.
	 * @param starts
	 *            the number of times the algorithm will be started in the
	 *            <code>train</code>-method, at least 1
	 * @param motifProb
	 *            the probability of finding a motif in a sequence (in [0,1])
	 * @param posPrior
	 *            this object determines the positional distribution that shall
	 *            be used
	 * @param trainOnlyMotifModel
	 *            a switch whether to train only the motif model
	 * @param alpha
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            the type of the component probability parameterization
	 *            <ul>
	 *            <li>{@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#THETA} or
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            <li>the parameterization of in a component is determined by
	 *            the component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link de.jstacs.models.mixture.AbstractMixtureModel.Parameterization#LAMBDA}
	 *            </ul>
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of the
	 *             corresponding length
	 *             <li><code>motifProb &lt; 0</code> or
	 *             <code>motifProb &gt; 1</code>
	 *             <li><code>starts &lt; 1</code>
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same simple
	 *             alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 * 
	 * @see SingleHiddenMotifMixture#SingleHiddenMotifMixture(de.jstacs.models.Model, de.jstacs.models.Model, boolean, int, double[], double[], de.jstacs.models.mixture.motif.positionprior.PositionPrior, de.jstacs.models.mixture.AbstractMixtureModel.Algorithm, double, de.jstacs.algorithms.optimization.termination.TerminationCondition,  de.jstacs.models.mixture.AbstractMixtureModel.Parameterization, int, int, de.jstacs.sampling.BurnInTest)     
	 * @see de.jstacs.models.mixture.AbstractMixtureModel.Algorithm#EM
	 */
	public SingleHiddenMotifMixture( Model motif, Model bg, boolean trainOnlyMotifModel, int starts, double motifProb,
										PositionPrior posPrior, double alpha, TerminationCondition tc, Parameterization parametrization )
																															throws CloneNotSupportedException,
																															IllegalArgumentException,
																															WrongAlphabetException {
		this( motif,
				bg,
				trainOnlyMotifModel,
				starts,
				null,
				new double[]{ motifProb, 1d - motifProb },
				posPrior,
				Algorithm.EM,
				alpha,
				tc,
				parametrization,
				0,
				0,
				null );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SingleHiddenMotifMixture} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation of the model as a {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	public SingleHiddenMotifMixture( StringBuffer xml ) throws NonParsableException {
		super( xml );
		bgMaxMarkovOrder = model[1].getMaximalMarkovOrder();
		int i = 1;
		while( i < model.length && !optimizeModel[i] ) {
			i++;
		}
		trainOnlyMotifModel = i == model.length;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#setTrainData(de.jstacs.data.Sample)
	 */
	@Override
	protected void setTrainData( DataSet data ) throws Exception {
		LinkedList<Sequence> fg = new LinkedList<Sequence>();
		LinkedList<Sequence> bg = new LinkedList<Sequence>();

		int i = 0, l, start, motifLength = getMotifLength( 0 );

		boolean rev = model[0] instanceof StrandModel;
		refBgSample = new int[data.getNumberOfElements() + 1];
		Sequence s;
		while( i < data.getNumberOfElements() ) {
			s = data.getElementAt( i );
			if( rev ) {
				s.reverseComplement();
			}

			if( !trainOnlyMotifModel ) {
				bg.add( s );
			}
			l = s.getLength() - motifLength;
			for( start = 0; start <= l; start++ ) {
				fg.add( s.getSubSequence( start, motifLength ) );
				if( !trainOnlyMotifModel ) {
					bg.add( s.getSubSequence( 0, start ) );
					bg.add( s.getSubSequence( start + motifLength ) );
				}
			}
			refBgSample[++i] = trainOnlyMotifModel ? i : bg.size();
		}
		Sequence[] empty = new Sequence[0];
		sample = new DataSet[]{ new DataSet( "possible motifs", fg.toArray( empty ) ), data };
		if( bg.size() != 0 ) {
			sample[1] = new DataSet( "possible background", bg.toArray( empty ) );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#createSeqWeightsArray()
	 */
	@Override
	protected double[][] createSeqWeightsArray() {
		if( trainOnlyMotifModel ) {
			return new double[][]{ new double[sample[0].getNumberOfElements()] };
		} else {
			return new double[][]{ new double[sample[0].getNumberOfElements()], new double[sample[1].getNumberOfElements()] };
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#doFirstIteration(double[], de.jstacs.utils.random.MultivariateRandomGenerator, de.jstacs.utils.random.MRGParams[])
	 */
	@Override
	protected double[][] doFirstIteration( double[] dataWeights, MultivariateRandomGenerator m, MRGParams[] params ) throws Exception {
		int j, i = 0, fgStart = 0, bgStart = 0, l = refBgSample.length - 1, ml = getMotifLength( 0 ) - 1, len;
		double d;
		double[][] seqweights = createSeqWeightsArray();
		double[] helpArray, w = new double[2];
		initWithPrior( w );
		if( !estimateComponentProbs && weights[0] == 1 ) {
			// there has to be a motif in each sequence
			while( i < l ) {
				len = sample[1].getElementAt( refBgSample[i] ).getLength() - ml;
				m.generate( seqweights[0], fgStart, len, params[i] );
				d = ( dataWeights == null ) ? 1d : dataWeights[i];
				if( trainOnlyMotifModel ) {
					len = fgStart + len;
					while( fgStart < len ) {
						seqweights[0][fgStart++] *= d;
					}
				} else {
					seqweights[1][bgStart++] = 0;
					len = fgStart + len;
					for( ; fgStart < len; fgStart++ ) {
						seqweights[0][fgStart] *= d;
						seqweights[1][bgStart++] = seqweights[1][bgStart++] = seqweights[0][fgStart];
					}
				}
				i++;
			}
		} else {
			ml -= 1;
			while( i < l ) {
				len = sample[1].getElementAt( refBgSample[i] ).getLength() - ml;
				helpArray = m.generate( len, params[i] );
				d = ( dataWeights == null ) ? 1d : dataWeights[i];
				w[0] += ( 1d - helpArray[0] ) * d;
				w[1] += helpArray[0] * d;
				if( trainOnlyMotifModel ) {
					for( j = 1; j < helpArray.length; j++, fgStart++ ) {
						seqweights[0][fgStart] = helpArray[j] * d;
					}
				} else {
					seqweights[1][bgStart++] = helpArray[0] * d;
					for( j = 1; j < helpArray.length; j++, fgStart++ ) {
						seqweights[0][fgStart] = helpArray[j] * d;
						seqweights[1][bgStart++] = seqweights[1][bgStart++] = seqweights[0][fgStart];
					}
				}
				i++;
			}
		}
		getNewParameters( 0, seqweights, w );
		return seqweights;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#getNewWeights(double[], double[], double[][])
	 */
	@Override
	protected double getNewWeights( double[] dataWeights, double[] w, double[][] seqweights ) throws Exception {
		double ll = 0, logPSeq, currentWeight = 1;
		Sequence seq;
		int seqIndex = 0, motifLength = getMotifLength( 0 ), end, j, w0Index = 0, w1Index = 0, start;

		// simple solution

		// weights[0][w0Index] = logWeights[0] // prior for motif
		// + prior.getLogPriorForPositions( j ) // position prior
		// + model[1].getLogProbFor( seq, 0, j - 1 ) // front bg sequence
		// + model[0].getLogProbFor( seq, j ) // motif
		// + model[1].getLogProbFor( seq, j + motifLength ) // end bg sequence
		//

		// complex solution (only possible if for the homogeneous model the probability distribution is rewritable
		// as a product of position dependent terms.)
		// -> not possible for mixture models (e.g. a real mixture model, a cyclic model, ...)

		int b1, b2, l, s, e;
		double[] help = new double[2];
		initWithPrior( w );
		while( w0Index < sample[0].getNumberOfElements() ) {
			seq = sample[1].getElementAt( refBgSample[seqIndex] );
			if( dataWeights != null ) {
				currentWeight = dataWeights[seqIndex];
			}
			l = seq.getLength();
			logPSeq = model[1].getLogProbFor( seq, 0, l - 1 );

			end = l - motifLength;
			b1 = -bgMaxMarkovOrder;
			b2 = motifLength + bgMaxMarkovOrder - 1;
			start = w0Index;
			for( j = 0; j <= end; j++, b1++, b2++, w0Index++ ) {
				s = Math.max( b1, 0 );
				e = Math.min( b2, l - 1 );
				seqweights[0][w0Index] = posPrior.getLogPriorForPositions( l, j ) // position prior
											+ model[0].getLogProbFor( sample[0].getElementAt( w0Index ), 0, motifLength - 1 )// motif
											- model[1].getLogProbFor( seq, s, e ) // bg sequence around the motif
											+ model[1].getLogProbFor( seq, s, j - 1 ) // part of the bg sequence left of the motif
											+ model[1].getLogProbFor( seq, j + motifLength, e ); // part of the bg sequence right of the
				// motif
			}
			ll += currentWeight * ( logPSeq + modify( help, seqweights[0], start, w0Index ) );

			help[0] *= currentWeight;

			w[0] += help[0];
			w[1] += currentWeight * help[1];

			if( trainOnlyMotifModel ) {
				for( ; start < w0Index; start++ ) {
					seqweights[0][start] *= help[0];
				}
			} else {
				seqweights[1][w1Index++] = currentWeight * help[1];
				for( ; start < w0Index; start++ ) {
					seqweights[0][start] *= help[0];
					seqweights[1][w1Index++] = seqweights[1][w1Index++] = seqweights[0][start];
				}
			}
			seqIndex++;
		}
		return ll;
	}

	/**
	 * This method modifies the computed weights for one sequence and returns
	 * the score.
	 * 
	 * @param containsMotif
	 *            an array to return the weights for containing a motif (index
	 *            0) or containing no motif (index 1)
	 * @param startpos
	 *            the array containing the scores for each start position
	 *            (including no motif in the sequence)
	 * @param start
	 *            the start index
	 * @param end
	 *            the end index
	 * 
	 * @return the score
	 */
	protected double modify( double[] containsMotif, double[] startpos, int start, int end ) {
		switch( algorithm ) {
			case EM:
				containsMotif[0] = logWeights[0] + Normalisation.logSumNormalisation( startpos, start, end, startpos, start );
				containsMotif[1] = logWeights[1];
				return Normalisation.logSumNormalisation( containsMotif, 0, 2, containsMotif, 0 );
			case GIBBS_SAMPLING:
				throw new IllegalArgumentException( "Gibbs Sampling currently not implemented." );
			default:
				throw new IllegalArgumentException( "The type of algorithm is unknown." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.AbstractMixtureModel#getLogProbUsingCurrentParameterSetFor(int, de.jstacs.data.Sequence, int, int)
	 */
	@Override
	protected double getLogProbUsingCurrentParameterSetFor( int component, Sequence seq, int start, int end ) throws Exception {
		switch( component ) {
			case 0:
				int current = 0,
				l = end - start + 1,
				motifLength = getMotifLength( 0 ),
				m = l - motifLength;
				int s,
				e,
				b1 = -bgMaxMarkovOrder,
				b2 = motifLength + bgMaxMarkovOrder - 1;
				double all = model[1].getLogProbFor( seq, start, end ),
				res = Double.NEGATIVE_INFINITY;
				for( ; current <= m; current++, b1++, b2++ ) {
					s = Math.max( b1, 0 );
					e = Math.min( b2, l - 1 );
					res = Normalisation.getLogSum( res, // current value 
							posPrior.getLogPriorForPositions( l, current ) //position
									+ model[0].getLogProbFor( seq, start + current ) // motif
									- model[1].getLogProbFor( seq, start + s, start + e ) // bg sequence around the motif
									+ model[1].getLogProbFor( seq, start + s, start + current - 1 ) // part of the bg sequence left of the motif
									+ model[1].getLogProbFor( seq, start + current + motifLength, start + e ) ); // part of the bg sequence right of the motif
				}
				return all + logWeights[0] + res;
			case 1:
				return logWeights[1] + model[1].getLogProbFor( seq, start, end );
			default:
				throw new IndexOutOfBoundsException( "This model has only two components (0=motif, 1=no motif)." );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getProfileOfScoresFor(int, int, de.jstacs.data.Sequence, int, de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile)
	 */
	public double[] getProfileOfScoresFor( int component, int motif, Sequence sequence, int startpos, KindOfProfile kind ) throws Exception {
		if( component == 0 && motif == 0 ) {
			int motifLength = getMotifLength( motif ), l = sequence.getLength() - startpos, len = l - motifLength + 1;
			switch( algorithm ) {
				case EM:
					double d;
					if( kind == KindOfProfile.UNNORMALIZED_JOINT ) {
						d = logWeights[component];
					} else {
						d = 0;
					}
					double[] weights = new double[len];
					int s,
					e,
					b1 = -bgMaxMarkovOrder,
					b2 = motifLength + bgMaxMarkovOrder - 1;
					for( int current = 0; current < len; current++, b1++, b2++ ) {
						s = Math.max( b1, 0 );
						e = Math.min( b2, l - 1 );
						weights[current] = d + posPrior.getLogPriorForPositions( l, current + startpos ) //position
											+ model[0].getLogProbFor( sequence, current + startpos ) // motif
											- model[1].getLogProbFor( sequence, s, e ) // bg sequence around the motif
											+ model[1].getLogProbFor( sequence, s, current + startpos - 1 ) // part of the bg sequence left of the motif
											+ model[1].getLogProbFor( sequence, current + startpos + motifLength, e ); // part of the bg sequence right of the motif
					}
					if( kind == KindOfProfile.NORMALIZED_CONDITIONAL ) {
						d = Normalisation.getLogSum( 0, weights.length, weights );
						for( int current = 0; current < len; current++, b1++, b2++ ) {
							weights[current] -= d;
						}
					}
					return weights;
				case GIBBS_SAMPLING:
					throw new IllegalArgumentException( "Gibbs Sampling currently not implemented." );
				default:
					throw new IllegalArgumentException( "The type of algorithm is unknown." );
			}
		} else {
			//return null;
			throw new IndexOutOfBoundsException();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.HiddenMotifMixture#getMinimalSequenceLength()
	 */
	@Override
	public int getMinimalSequenceLength() {
		if( estimateComponentProbs || weights[1] != 0 ) {
			return 0;
		} else {
			return model[0].getLength();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getMotifLength(int)
	 */
	public int getMotifLength( int motif ) {
		if( motif == 0 ) {
			return model[0].getLength();
		} else {
			throw new IndexOutOfBoundsException();
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifs()
	 */
	public int getNumberOfMotifs() {
		return 1;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getNumberOfMotifsInComponent(int)
	 */
	public int getNumberOfMotifsInComponent( int component ) {
		switch( component ) {
			case 0:
				return 1; //component with motif
			case 1:
				return 0; //component without component
			default:
				throw new IndexOutOfBoundsException();
		}
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getStrandFor(int, int, de.jstacs.data.Sequence, int)
	 */
    public double[] getStrandProbabilitiesFor( int component, int motif, Sequence sequence, int startpos ) throws Exception {
        if( component == 0 && motif == 0 ) {
            if( model[0] instanceof StrandModel ) {
                Sequence help = sequence.getSubSequence(startpos, model[0].getLength());
                double[] logProbs = {
                     ( (StrandModel)model[0] ).getLogProbFor( 0, help ),
                     ( (StrandModel)model[0] ).getLogProbFor( 1, help )
                };
                Normalisation.logSumNormalisation(logProbs);
                return logProbs;
            } else {
                return new double[]{1.0,0.0};
            }
        } else {
            throw new IndexOutOfBoundsException();
        }
    }

	/* (non-Javadoc)
	 * @see de.jstacs.motifDiscovery.MotifDiscoverer#getGlobalIndexOfMotifInComponent(int, int)
	 */
	public int getGlobalIndexOfMotifInComponent( int component, int motif ) {
		if( component == 0 && motif == 0 ) {
			return 0;
		} else {
			throw new IndexOutOfBoundsException();
		}
	}
	
	/* (non-Javadoc)
	 * @see de.jstacs.models.mixture.motif.HiddenMotifMixture#trainBgModel(de.jstacs.data.Sample, double[])
	 */
	public void trainBgModel( DataSet data, double[] weights ) throws Exception {
		model[1].train( data, weights );
	}

/*
    public int getBestShift( Sample data, double[] dataWeights ) throws Exception {
//TODO remove
    	setTrainData( data );
    	System.out.println(this);
    	
    	double[][] seqWeights = createSeqWeightsArray();
    	
    	for( int modelIndex = 0; modelIndex < model.length; modelIndex++ ) {
    		if( model[modelIndex] instanceof AbstractMixtureModel ) {
    			( (AbstractMixtureModel)model[modelIndex] ).setTrainData( sample[modelIndex] );
    			usedWeights[modelIndex] = ( (AbstractMixtureModel)model[modelIndex] ).createSeqWeightsArray();
			}
    	}
    	if( algorithmHasBeenRun() ) {
    		sostream.writeln( "best:\t" + getScoreForBestRun() );
    	}
    	double ll = continueIterations(dataWeights, seqWeights, 0, 0);
    	sostream.writeln( "re-compute:\t" + ll );
//TODO end
    	
    	double[][] help = createSeqWeightsArray();
    	Model[] backup = ArrayHandler.clone(model);
    	int m = model[0].getLength()/2, idx=0;
    	double best = Double.NEGATIVE_INFINITY, current;
    	for( int i = -m; i <= m; i++ ) {
    		//shift weights & estimate
    		estimateShiftedParameters( i, seqWeights, help );
    		//compute
    		current = continueIterations(dataWeights, help, 0, 0);
    		//test
    		if( best < current ) {
    			idx = i;
    			best = current;
    		}
    		sostream.writeln( i + "\t" + current + "\t" + idx );
    		//reset to initial parameters
    		model = ArrayHandler.clone(backup);
    	}
    	return idx;
    }    
*/    
    private void estimateShiftedParameters( int shift, double[][] originalWeights, double[][] newWeights ) throws Exception {
    	//reset
    	for( int i = 0; i < newWeights.length; i++ ) {
    		Arrays.fill(newWeights[i], 0);
    	}
    	//set
		int motifLength = model[0].getLength(), end, i = 0, seqIndex = 0, start, bgIndex=0;
		double[] strandProbs;
    	while( i < sample[0].getNumberOfElements() ) {
    		start = i;
			end = start + sample[1].getElementAt( refBgSample[seqIndex++] ).getLength() - motifLength;
			for( int j = 0; i <= end; j++, i++ ) {
				strandProbs = getStrandProbabilitiesFor(0, 0, sample[0].getElementAt( i ), 0);
				
				//shift circular
								
				//motif
				newWeights[0][getIndexForCircularShift( start, end, i-shift )] += strandProbs[0]*originalWeights[0][i];
				newWeights[0][getIndexForCircularShift( start, end, i+shift )] += strandProbs[1]*originalWeights[0][i];
			}
			
			//bg
	    	if( !trainOnlyMotifModel ) {
	    		//TODO test
	    		newWeights[1][bgIndex] = originalWeights[1][bgIndex];
	    		bgIndex++;
				for( ; start < end; start++ ) {
					newWeights[1][bgIndex++] = newWeights[1][bgIndex++] = newWeights[0][start];
				}
	    	}			
    	}
    	//estimate
    	for( i = 0; i < newWeights.length; i++ ) {
			getNewParametersForModel( i, 1, i, newWeights[i] );
		}
    }
    
    private int getIndexForCircularShift( int start, int end, int proposal ) {
    	if( proposal < start ) {
    		proposal = end - (start-proposal)-1;
		}
		if( proposal > end ) {
			proposal = start + (proposal-end)-1;
		}
		return proposal;
    }
    
	@Override
	protected double iterate( int start, double[] dataWeights, MultivariateRandomGenerator m, MRGParams params[] ) throws Exception {
		sostream.writeln( "========== start: " + start + " ==========" );
		switch( algorithm ) {
			case EM:
				seqWeights = doFirstIteration( dataWeights, m, params );
				double[][] help = null;
				Model[] backup = null;
				double best, current;
				int ml = model[0].getLength()/2, shift;
				do {
					//optimize
					this.best = continueIterations( dataWeights, seqWeights );
					shift = 0;
					if(correctPhaseShift) {
						//test shift
						if( help == null ) {
							help = createSeqWeightsArray();
						}
						backup = ArrayHandler.clone(model);
						best = Double.NEGATIVE_INFINITY;
				    	for( int i = -ml; i <= ml; i++ ) {
				    		//shift weights & estimate
				    		estimateShiftedParameters( i, seqWeights, help );
				    		//compute
				    		current = continueIterations( dataWeights, help, 0, 0 );
				    		//test
				    		if( best < current ) {
				    			shift = i;
				    			best = current;
				    		}
				    		sostream.writeln( i + "\t" + current + "\t" + shift );
				    		//reset to initial parameters
				    		model = ArrayHandler.clone(backup);
				    	}
						
						//found relevant shift
						if( shift != 0 ) {
							estimateShiftedParameters( shift, seqWeights, help );
						}
					}
				}while( shift != 0 );
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
	 * Enables or disables the phase shift correction. By default, shift correction is enabled.
	 * @param correct switch that determines whether to correct shifts or not
	 */
	public void setShiftCorrection(boolean correct) {
		correctPhaseShift = correct;
	}
}
