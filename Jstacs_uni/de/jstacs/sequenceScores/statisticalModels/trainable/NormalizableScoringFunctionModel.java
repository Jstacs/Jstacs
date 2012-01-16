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

package de.jstacs.sequenceScores.statisticalModels.trainable;

import java.io.OutputStream;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.algorithms.optimization.StartDistanceForecaster;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.LogGenDisMixFunction;
import de.jstacs.classifier.differentiableSequenceScoreBased.logPrior.CompositeLogPrior;
import de.jstacs.classifier.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.DataSet.WeightedDataSetFactory;
import de.jstacs.data.DataSet.WeightedDataSetFactory.SortOperation;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;
import de.jstacs.sequenceScores.statisticalModels.differentiable.IndependentProductDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.UniformDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.homogeneous.UniformHomogeneousDiffSM;
import de.jstacs.utils.SafeOutputStream;

/**
 * This model can be used to use a DifferentiableStatisticalModel as model.
 * It enables the user to train the DifferentiableStatisticalModel in a generative way.
 * 
 * @author Jens Keilwagen
 * 
 * @see DifferentiableStatisticalModel
 * @see LogGenDisMixFunction
 */
public class NormalizableScoringFunctionModel extends AbstractTrainSM
{
	private SafeOutputStream out;
	
	/**
	 * The internally used {@link DifferentiableStatisticalModel}.
	 */
	protected DifferentiableStatisticalModel nsf;
	private double logNorm, lineps, startD;
	private AbstractTerminationCondition tc;
	private byte algo;
	private int threads;

	/**
	 * The main constructor that creates an instance with the user given parameters.
	 * 
	 * @param nsf the {@link DifferentiableStatisticalModel} that should be used
	 * @param threads the number of threads that should be used for optimization
	 * @param algo the algorithm that should be used for the optimization
	 * @param tc the {@link AbstractTerminationCondition} for stopping the optimization
	 * @param lineps the line epsilon for stopping the line search in the optimization
	 * @param startD the start distance that should be used initially
	 * 
	 * @throws CloneNotSupportedException if <code>nsf</code> can not be cloned
	 */
	public NormalizableScoringFunctionModel( DifferentiableStatisticalModel nsf, int threads, byte algo, AbstractTerminationCondition tc, double lineps, double startD ) throws CloneNotSupportedException
	{
		super( nsf.getAlphabetContainer(), nsf.getLength() );
		if( threads < 1 )
		{
			throw new IllegalArgumentException( "The number of threads has to be positive." );
		}
		this.threads = threads;
		this.tc = tc.clone();
		if( lineps < 0 )
		{
			throw new IllegalArgumentException( "The value of lineps has to be non-negative." );
		}
		this.lineps = lineps;
		if( startD <= 0 )
		{
			throw new IllegalArgumentException( "The value of startD has to be positive." );
		}
		this.startD = startD;
		this.algo = algo;
		this.nsf = (DifferentiableStatisticalModel) nsf.clone();
		if( isInitialized() )
		{
			logNorm = nsf.getLogNormalizationConstant();
		}
		else
		{
			logNorm = Double.NEGATIVE_INFINITY;
		}
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link NormalizableScoringFunctionModel} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             is thrown if the {@link StringBuffer} could not be parsed
	 */
	public NormalizableScoringFunctionModel( StringBuffer stringBuff ) throws NonParsableException
	{
		super( stringBuff );
	}

	public NormalizableScoringFunctionModel clone() throws CloneNotSupportedException
	{
		NormalizableScoringFunctionModel clone = (NormalizableScoringFunctionModel) super.clone();
		clone.nsf = (DifferentiableStatisticalModel) nsf.clone();
		clone.tc = tc.clone();
		clone.setOutputStream( out.doesNothing() ? null : SafeOutputStream.DEFAULT_STREAM );
		return clone;
	}
	
	public void train( DataSet data, double[] weights ) throws Exception
	{
		if( !data.getAlphabetContainer().checkConsistency( alphabets ) )
		{
			throw new WrongAlphabetException( "The AlphabetConatainer of the sample and the model do not match." );
		}
		if( length != 0 && length != data.getElementLength() )
		{
			throw new WrongLengthException( "The length of the elements of the sample is not suitable for the model." );
		}
		
		if( nsf instanceof IndependentProductDiffSM ) {
			IndependentProductDiffSM ipsf = (IndependentProductDiffSM) nsf;
			DifferentiableStatisticalModel[] nsfs = ArrayHandler.cast( DifferentiableStatisticalModel.class, ipsf.getFunctions() );
			DataSet[] part = new DataSet[1], packedData = { data };
			double[][] partWeights, packedWeights = { weights };
			for( int a, i = 0; i < nsfs.length; i++ ) {
				a = ipsf.extractSequenceParts( i, packedData, part );
				partWeights = ipsf.extractWeights( a, packedWeights );
				nsfs[i] = train( part[0], partWeights[0], nsfs[i] );
			}
			nsf = new IndependentProductDiffSM( ipsf.getESS(), true, nsfs, ipsf.getIndices(), ipsf.getPartialLengths(), ipsf.getReverseSwitches() );
		} else {
			nsf = train( data, weights, nsf );
		}
	}
	
	private DifferentiableStatisticalModel train( DataSet data, double[] weights, DifferentiableStatisticalModel nsf ) throws Exception {
		if( !(nsf instanceof UniformDiffSM || nsf instanceof UniformHomogeneousDiffSM ) ) {
			WeightedDataSetFactory wsf = new WeightedDataSetFactory( SortOperation.NO_SORT, data, weights );
			DataSet small = wsf.getDataSet();
			double[] smallWeights = wsf.getWeights();  
			
			double[] params;
			DifferentiableStatisticalModel best = null;
			double current, max = Double.NEGATIVE_INFINITY, fac = data.getNumberOfElements(), ess = nsf.getESS();
			fac = fac / (ess+ fac) * (ess == 0 ? 1d : 2d);
			
			DifferentiableStatisticalModel[] score = { (DifferentiableStatisticalModel) nsf.clone() };
			LogPrior prior = new CompositeLogPrior();
			double[] beta = LearningPrinciple.getBeta( ess == 0 ? LearningPrinciple.ML : LearningPrinciple.MAP );
			LogGenDisMixFunction f = new LogGenDisMixFunction( threads, score, new DataSet[]{small}, new double[][]{smallWeights}, prior, beta, true, false );
			NegativeDifferentiableFunction minusF = new NegativeDifferentiableFunction( f );
			StartDistanceForecaster sd =
				//new ConstantStartDistance( startD*fac );
				new LimitedMedianStartDistance( 5, startD*fac );
			for( int i = 0; i < nsf.getNumberOfRecommendedStarts(); i++ )
			{
				out.writeln( "start: " + i );
				//TODO freeParams???
				score[0].initializeFunction( 0, false, new DataSet[]{small}, new double[][]{smallWeights} );
				f.reset( score );
				params = f.getParameters( KindOfParameter.PLUGIN );
				sd.reset();
				Optimizer.optimize( algo, minusF, params, tc, lineps*fac, sd, out );
				current = f.evaluateFunction( params );
				if( current > max )
				{
					best = score[0];
					max = current;
				}
				score[0] = (DifferentiableStatisticalModel) nsf.clone();
			}
			out.writeln( "best: " + max );
			nsf = best;
			logNorm = nsf.getLogNormalizationConstant();
			f.stopThreads();
			System.gc();
		}
		return nsf;
	}

	public double getLogProbFor( Sequence sequence, int startpos, int endpos ) throws NotTrainedException, Exception
	{
		if( !isInitialized() )
		{
			throw new NotTrainedException();
		}
		if( !sequence.getAlphabetContainer().checkConsistency(alphabets) )
		{
			throw new WrongAlphabetException( "The AlphabetContainer of the sequence and the model do not match." );
		}
		if( startpos < 0 )
		{
			throw new IllegalArgumentException( "Check start position." );
		}
		if( endpos+1 < startpos || endpos >= sequence.getLength() )
		{
			throw new IllegalArgumentException( "Check end position." );
		}
		if( length != 0 && length != endpos-startpos+1 )
		{
			throw new WrongLengthException( "Check length of the sequence." );
		}
		return nsf.getLogScoreFor( sequence, startpos ) - logNorm;
	}
	
	public double getLogPriorTerm() throws Exception
	{
		return nsf.getLogPriorTerm() - nsf.getESS()*logNorm;
	}

	public String getInstanceName()
	{
		return "model using " + nsf.getInstanceName();
	}

	public boolean isInitialized()
	{
		return nsf.isInitialized();
	}

	public NumericalResultSet getNumericalCharacteristics() throws Exception
	{
		return null;
	}

	public String toString()
	{
		return nsf.toString();
	}

	private static final String XML_TAG = "NormalizableScoringFunctionModel";
	
	protected void fromXML( StringBuffer xml ) throws NonParsableException
	{
		StringBuffer rep = XMLParser.extractForTag( xml, XML_TAG );
		nsf = XMLParser.extractObjectForTags( rep, "DifferentiableStatisticalModel", DifferentiableStatisticalModel.class );
		threads = XMLParser.extractObjectForTags( rep, "threads", int.class );
		algo = XMLParser.extractObjectForTags( rep, "algorithm", byte.class );
		if( XMLParser.hasTag( rep, "tc", null, null ) ) {
			tc = (AbstractTerminationCondition) XMLParser.extractObjectForTags( rep, "tc" );
		} else {
			try {
				tc = new SmallDifferenceOfFunctionEvaluationsCondition( XMLParser.extractObjectForTags( rep, "eps", double.class ) );
			} catch (Exception e) {
				NonParsableException n = new NonParsableException( e.getMessage() );
				throw n;
			}
		}
		lineps = XMLParser.extractObjectForTags( rep, "lineps", double.class );
		startD = XMLParser.extractObjectForTags( rep, "startDistance", double.class );
		if( isInitialized() )
		{
			logNorm = nsf.getLogNormalizationConstant();
		}
		else
		{
			logNorm = Double.NEGATIVE_INFINITY;
		}
		alphabets = nsf.getAlphabetContainer();
		length = nsf.getLength();
		setOutputStream( SafeOutputStream.DEFAULT_STREAM );
	}
	
	public StringBuffer toXML()
	{
		StringBuffer xml = new StringBuffer( 100000 );
		XMLParser.appendObjectWithTags( xml, nsf, "DifferentiableStatisticalModel" );
		XMLParser.appendObjectWithTags( xml, threads, "threads" );
		XMLParser.appendObjectWithTags( xml, algo, "algorithm" );
		XMLParser.appendObjectWithTags( xml, tc, "tc" );
		XMLParser.appendObjectWithTags( xml, lineps, "lineps" );
		XMLParser.appendObjectWithTags( xml, startD, "startDistance" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
	
	/**
	 * Sets the OutputStream that is used e.g. for writing information while training. It is possible to set
	 * <code>o=null</code>, than nothing will be written.
	 * 
	 * @param o
	 *            the OutputStream
	 */
	public final void setOutputStream( OutputStream o )
	{
		out = SafeOutputStream.getSafeOutputStream( o );
	}
	
	/**
	 * Returns a copy of the internally used {@link DifferentiableStatisticalModel}.
	 * @return a copy of the internally used {@link DifferentiableStatisticalModel}
	 * @throws CloneNotSupportedException if the internal instance could not be cloned
	 */
	public DifferentiableStatisticalModel getFunction() throws CloneNotSupportedException
	{
		return (DifferentiableStatisticalModel) nsf.clone();
	}
}
