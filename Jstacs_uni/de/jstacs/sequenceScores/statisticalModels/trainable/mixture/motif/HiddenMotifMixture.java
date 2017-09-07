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

package de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif;

import java.text.NumberFormat;
import java.util.Arrays;

import javax.naming.OperationNotSupportedException;

import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.motifDiscovery.MotifDiscoverer;
import de.jstacs.sampling.BurnInTest;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.mixture.motif.positionprior.PositionPrior;

/**
 * This is the main class that every generative motif discoverer should
 * implement.
 * 
 * @author Jens Keilwagen
 */
public abstract class HiddenMotifMixture extends AbstractMixtureTrainSM implements MotifDiscoverer {

	/**
	 * The prior for the positions.
	 */
	protected PositionPrior posPrior;

	/**
	 * Creates a new {@link HiddenMotifMixture}. This constructor can be used
	 * for any algorithm since it takes all necessary values as parameters.
	 * 
	 * @param models
	 *            the single models building the {@link HiddenMotifMixture}, if
	 *            the model is trained using
	 *            {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING} the
	 *            models that will be adjusted have to implement
	 *            {@link de.jstacs.sampling.SamplingComponent}.
	 *            The models that are used for the flanking sequences have to
	 *            be able to score sequences of arbitrary length.
	 * @param optimzeArray
	 *            a array of switches whether to train or not the corresponding model
	 * @param components
	 * 			  the number of components (e.g. for ZOOPS this is 2)
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
	 * @param posPrior
	 *            this object determines the positional distribution that shall
	 *            be used
	 * @param algorithm
	 *            either {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#EM} or
	 *            {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}
	 * @param alpha
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the positive parameter for the Dirichlet distribution which is
	 *            used when you invoke <code>train</code> to initialize the
	 *            gammas. It is recommended to use <code>alpha = 1</code>
	 *            (uniform distribution on a simplex).
	 * @param tc
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the {@link TerminationCondition} for stopping the EM-algorithm,
	 *            <code>tc</code> has to return <code>true</code> from {@link TerminationCondition#isSimple()}
	 * @param parametrization
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#EM}<br>
	 *            the type of the component probability parameterization;
	 *            <ul>
	 *            <li>{@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization#THETA} or
	 *            {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            <li>the parameterization of a component is determined by the
	 *            component model
	 *            <li>it is recommended to use the same parameterization for the
	 *            components and the component assignment probabilities
	 *            <li>it is recommended to use
	 *            {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Parameterization#LAMBDA}
	 *            </ul>
	 * @param initialIteration
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the initial sampling phase (at least 1,
	 *            at most <code>stationaryIteration/starts</code>)
	 * @param stationaryIteration
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the positive length of the stationary phase (at least 1)
	 *            (summed over all starts), i.e. the number of parameter sets
	 *            that is used in approximation
	 * @param burnInTest
	 *            only for {@link de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM.Algorithm#GIBBS_SAMPLING}<br>
	 *            the test that will be used to determine the length of the
	 *            burn-in phase
	 * 
	 * @throws IllegalArgumentException
	 *             if
	 *             <ul>
	 *             <li>the models are not able to score the sequence of the
	 *             corresponding length
	 *             <li><code>weights != null &amp;&amp; weights.length != 2</code>
	 *             <li><code>weights != null</code> and it exists an
	 *             <code>i</code> where <code>weights[i] &lt; 0</code>
	 *             <li><code>starts &lt; 1</code>
	 *             <li><code>componentHyperParams</code> are not correct
	 *             <li>the algorithm specific parameters are not correct
	 *             </ul>
	 * @throws WrongAlphabetException
	 *             if not all <code>models</code> work on the same simple
	 *             alphabet
	 * @throws CloneNotSupportedException
	 *             if the <code>models</code> can not be cloned
	 */
	protected HiddenMotifMixture( TrainableStatisticalModel[] models, boolean[] optimzeArray, int components, int starts, boolean estimateComponentProbs, double[] componentHyperParams,
									double[] weights, PositionPrior posPrior, Algorithm algorithm,
									double alpha, TerminationCondition tc, Parameterization parametrization, //EM parameters
									int initialIteration, int stationaryIteration, BurnInTest burnInTest ) //GIBBS_SAMPLING parameters
																											throws CloneNotSupportedException,
																											IllegalArgumentException,
																											WrongAlphabetException {
		super( posPrior.getLength(),
				models,
				optimzeArray,
				components,
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
		if( !alphabets.isSimple() ) {
			throw new WrongAlphabetException( "The AlphabetContainer has to be simple." );
		}
		this.posPrior = posPrior.clone();
		this.posPrior.setMotifLength( getMotifLength( 0 ) );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link HiddenMotifMixture} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation of the model as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} can not be parsed
	 */
	protected HiddenMotifMixture( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#clone()
	 */
	@Override
	public HiddenMotifMixture clone() throws CloneNotSupportedException {
		HiddenMotifMixture clone = (HiddenMotifMixture)super.clone();
		clone.posPrior = posPrior.clone();
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#getFurtherInformation()
	 */
	@Override
	protected StringBuffer getFurtherInformation() {
		StringBuffer erg = new StringBuffer( 1000 );
		XMLParser.appendObjectWithTags( erg, posPrior, "posPrior" );
		return erg;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#extractFurtherInformation(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		posPrior = XMLParser.extractObjectForTags( xml, "posPrior", PositionPrior.class );
		posPrior.setMotifLength( getMotifLength( 0 ) );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#train(de.jstacs.data.DataSet, double[])
	 */
	@Override
	public void train( DataSet data, double[] weights ) throws Exception {
		if( data.getMinimalElementLength() < getMinimalSequenceLength() ) {
			throw new IllegalArgumentException( "The data set contains sequence that are not allowed in this MotifDiscoverer. The minimal length is " + getMinimalSequenceLength()
												+ "." );
		}
		super.train( data, weights );
	}

	protected void getNewParameters( int iteration, double[][] seqWeights, double[] w ) throws Exception {
		for( int i = 0; i < seqWeights.length; i++ ) {
			getNewParametersForModel( i, iteration, i, seqWeights[i] );
		}
		getNewComponentProbs( w );
	}
	
	/**
	 * This method trains the background model. This can be useful if the
	 * background model is not trained during the EM-algorithm.
	 * 
	 * @param data
	 *            the data set
	 * @param weights
	 *            the weights
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public abstract void trainBgModel( DataSet data, double[] weights ) throws Exception;

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#checkLength(int, int)
	 */
	@Override
	protected void checkLength( int index, int l ) {
		if( index == 0 ) {
			if( length != 0 && length < l ) {
				throw new IllegalArgumentException( "The motif length is bigger than the length of the sequences the should be modeled." );
			}
		} else {
			if( l != 0 ) {
				throw new IllegalArgumentException( "All models accept the motif model have to be homogeneous. Violated at position " + index
													+ "." );
			}
		}
	}

	/**
	 * Returns the minimal length a sequence respectively a data set has to have.
	 * 
	 * @return the minimal length a sequence respectively a data set has to have
	 */
	public abstract int getMinimalSequenceLength();

	/* (non-Javadoc)
	 * @see de.jstacs.sequenceScores.statisticalModels.trainable.mixture.AbstractMixtureTrainSM#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		StringBuffer erg = new StringBuffer( getClass().getSimpleName() + "(" );
		erg.append( model[0].getInstanceName() );
		for( int i = 1; i < model.length; i++ ) {
			erg.append( ", " );
			erg.append( model[i].getInstanceName() );
		}
		erg.append( "; " + posPrior.getInstanceName() );
		if( !estimateComponentProbs ) {
			erg.append( "; " + Arrays.toString( weights ) );
		}
		erg.append( ")" );
		return erg.toString();
	}

	/* 
	 * (non-Javadoc)
	 * @see de.jstacs.sequenceScores.SequenceScore#toString(java.text.NumberFormat)
	 */
	@Override
	public String toString( NumberFormat nf ) {
		StringBuffer sb = new StringBuffer( 100000 );
		sb.append( nf.format(weights[0]) + "\tmotif\n" );
		sb.append( nf.format(weights[1]) + "\tno motif\n\n" );
		sb.append( "position prior: " + posPrior.getInstanceName() + "\n\n" );
		for( int i = 0; i < dimension; i++ ) {
			sb.append( model[i].getInstanceName() + "\n" + model[i].toString(nf) + "\n" );
		}
		return sb.toString();
	}

	/**
	 * Standard implementation throwing an
	 * {@link OperationNotSupportedException}.
	 */
	@Override
	protected Sequence[] emitDataSetUsingCurrentParameterSet( int n, int... lengths ) throws Exception {
		throw new OperationNotSupportedException();
	}
}
