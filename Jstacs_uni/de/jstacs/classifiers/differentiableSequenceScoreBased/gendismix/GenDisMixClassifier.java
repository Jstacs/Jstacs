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

package de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix;

import java.util.Arrays;

import de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.data.DataSet;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;

/**
 * This class implements a classifier the optimizes the following function
 * 
 * {@latex.ilb \\[f(\\underline{\\lambda}|C,D,\\underline{\\alpha},\\underline{\\beta})
 * := \\beta_0 \\log p(C|D,\\underline{\\lambda}) + \\beta_1 \\log p(C,D|\\underline{\\lambda}) + \\beta_2 \\log p(\\underline{\\lambda}|\\underline{\\alpha})
 * .\\]}
 * 
 * The weights {@latex.inline $\\beta_i$} have to sum to 1. For special weights the optimization turns out to be
 * well known
 * <ul>
 * <li> if the weights are (0,1,0), one obtains maximum likelihood,
 * <li> if the weights are (0,0.5,0.5), one obtains maximum a posteriori,
 * <li> if the weights are (1,0,0), one obtains maximum conditional likelihood,
 * <li> if the weights are (0.5,0,0.5), one obtains maximum supervised posterior,
 * <li> if the {@latex.inline $\\beta_2$}=0, one obtains the generative-discriminative trade-off,
 * <li> if the {@latex.inline $\\beta_2$}=0.5, one obtains the penalized generative-discriminative trade-off.
 * </ul>
 * 
 * Of course, there are also some very interesting cases with other weights.
 * 
 * @author Jens Keilwagen
 * 
 * @see LearningPrinciple
 * @see LogGenDisMixFunction
 */
public class GenDisMixClassifier extends ScoreClassifier
{
	/**
	 * The prior that is used in this classifier.
	 */
	protected LogPrior prior;
	
	/**
	 * The function that is optimized in this classifier. 
	 */
	protected LogGenDisMixFunction function;
	
	/**
	 * The weights that determine a point of the simplex that corresponds to the function that is optimized by this classifier.
	 */
	protected double[] beta;

	/**
	 * This constructor creates an instance and sets the value of the last (external) optimization. This constructor should only be used
	 * if it is clear that now generative part is used in the objective function.
	 * 
	 * @param params
	 *            the common parameters
	 * @param prior
	 *            the prior for the parameters of the classifier
	 * @param lastScore
	 *            the score of the last optimization
	 * @param beta
	 *            the weights for log likelihood, log conditional likelihood and log prior
	 * @param score
	 *            the functions for the classes
	 * 
	 * @throws CloneNotSupportedException if an instance could not be clone
	 */	
	protected GenDisMixClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, double lastScore, double[] beta,
			DifferentiableSequenceScore... score ) throws CloneNotSupportedException
	{
		super( params, lastScore, score );
		setWeights( beta );
		setPrior( prior );
	}
	
	/**
	 * This constructor creates an instance and sets the value of the last (external) optimization.
	 * 
	 * @param params
	 *            the common parameters
	 * @param prior
	 *            the prior for the parameters of the classifier
	 * @param lastScore
	 *            the score of the last optimization
	 * @param beta
	 *            the weights for log likelihood, log conditional likelihood and log prior
	 * @param score
	 *            the functions for the classes
	 * 
	 * @throws CloneNotSupportedException if an instance could not be clone
	 */	
	public GenDisMixClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, double lastScore, double[] beta,
			DifferentiableStatisticalModel... score ) throws CloneNotSupportedException
	{
		this( params, prior, lastScore, beta, (DifferentiableSequenceScore[]) score );
	}
	
	/**
	 * The main constructor.
	 * 
	 * @param params
	 *            the common parameters
	 * @param prior
	 *            the prior for the parameters of the classifier
	 * @param beta
	 *            the weights for log likelihood, log conditional likelihood and log prior
	 * @param score
	 *            the functions for the classes
	 * 
	 * @throws CloneNotSupportedException if an instance could not be clone
	 */	
	public GenDisMixClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, double[] beta,
			DifferentiableStatisticalModel... score ) throws CloneNotSupportedException
	{
		this( params, prior, NOT_TRAINED_VALUE, beta, score );
	}

	/**
	 * This convenience constructor agglomerates the <code>genBeta, disBeta,</code> and <code>priorBeta</code> into an array and calls the main constructor.
	 * 
	 * @param params
	 *            the common parameters
	 * @param prior
	 *            the prior for the parameters of the classifier
	 * @param genBeta
	 *            the weight for the log likelihood
	 * @param disBeta
	 *            the weight for the log conditional likelihood
	 * @param priorBeta
	 *            the weight for the log prior
	 * @param score
	 *            the functions for the classes
	 * 
	 * @throws CloneNotSupportedException if an instance could not be clone
	 */	
	public GenDisMixClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, double genBeta, double disBeta,
			double priorBeta, DifferentiableStatisticalModel... score ) throws CloneNotSupportedException
	{
		this( params, prior, NOT_TRAINED_VALUE, new double[]{ disBeta, genBeta, priorBeta }, score );
	}

	/**
	 * This convenience constructor  creates an array of weights for an elementary learning principle and calls the main constructor.
	 * 
	 * @param params
	 *            the common parameters
	 * @param prior
	 *            the prior for the parameters of the classifier
	 * @param key
	 * 			  the key for an elementary {@link LearningPrinciple}
	 * @param score
	 *            the functions for the classes
	 * 
	 * @throws CloneNotSupportedException if an instance could not be clone
	 * 
	 * @see LearningPrinciple#getBeta(LearningPrinciple)
	 */	
	public GenDisMixClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, LearningPrinciple key, DifferentiableStatisticalModel... score )
			throws CloneNotSupportedException
	{
		this( params, prior, NOT_TRAINED_VALUE, LearningPrinciple.getBeta( key ), score );
	}
	

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public GenDisMixClassifier( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}
	
	public GenDisMixClassifier clone() throws CloneNotSupportedException
	{
		GenDisMixClassifier clone = (GenDisMixClassifier) super.clone();
		clone.prior = prior.getNewInstance();
		clone.beta = beta.clone();
		return clone;
	}

	protected LogGenDisMixFunction getFunction( DataSet[] data, double[][] weights ) throws Exception
	{
		GenDisMixClassifierParameterSet p = (GenDisMixClassifierParameterSet) params;
		if( data.length > 1 ) {
			return new LogGenDisMixFunction( p.getNumberOfThreads(), score,
				data, weights, prior, beta, p.shouldBeNormalized(), p.useOnlyFreeParameter() );
		} else {
			return new OneDataSetLogGenDisMixFunction( p.getNumberOfThreads(), score,
				data[0], weights, prior, beta, p.shouldBeNormalized(), p.useOnlyFreeParameter() );
		}
	}

	/**
	 * This method set a new prior that should be used for optimization. Since it could not be ensured that the
	 * classifier is optimal now <code>isTrained</code> will return <code>false</code> after invoking this method.
	 * 
	 * @param prior
	 *            the new prior
	 */
	public void setPrior( LogPrior prior )
	{
		if( prior != null )
		{
			this.prior = prior;
		}
		else
		{
			this.prior = DoesNothingLogPrior.defaultInstance;
		}
		hasBeenOptimized = false;
	}
	
	/**
	 * This method set the weights for the summand of the function. Since it could not be ensured that the classifier is optimal
	 * now <code>isTrained</code> will return <code>false</code> after invoking this method.
	 * 
	 * @param beta
	 *            the weights
	 *            
	 * @throws IllegalArgumentException if the weights array is not correct
	 * 
	 * @see LearningPrinciple#getBeta(LearningPrinciple)
	 */
	public void setWeights( double... beta ) throws IllegalArgumentException
	{
		this.beta = LearningPrinciple.checkWeights( beta );
		hasBeenOptimized = false;
	}

	private static final String XML_TAG = "gendismix-classifier";

	protected String getXMLTag()
	{
		return XML_TAG;
	}

	protected StringBuffer getFurtherClassifierInfos()
	{
		StringBuffer xml = super.getFurtherClassifierInfos();
		XMLParser.appendObjectWithTags( xml, beta, "beta" );
		if( !(prior instanceof DoesNothingLogPrior) )
		{
			StringBuffer pr = new StringBuffer( 1000 );
			pr.append( "<prior>\n" );
			XMLParser.appendObjectWithTags( pr, prior.getClass(), "class" );
			pr.append( prior.toXML() );
			pr.append( "\t</prior>\n" );
			xml.append( pr );
		}
		return xml;
	}

	protected void extractFurtherClassifierInfosFromXML( StringBuffer xml ) throws NonParsableException
	{
		super.extractFurtherClassifierInfosFromXML( xml );
	
		beta = LearningPrinciple.checkWeights( XMLParser.extractObjectForTags( xml, "beta", double[].class ) );
		
		StringBuffer pr = XMLParser.extractForTag( xml, "prior" );
		if( pr != null )
		{
			Class clazz = XMLParser.extractObjectForTags( pr, "class", Class.class );
			try
			{
				prior = (LogPrior) clazz.getConstructor( new Class[]{ StringBuffer.class } ).newInstance( pr );
			}
			catch( NoSuchMethodException e )
			{
				NonParsableException n = new NonParsableException( "You must provide a constructor " + clazz.getSimpleName() + "(StringBuffer)." );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
			catch( Exception e )
			{
				NonParsableException n = new NonParsableException( "problem at " + clazz.getSimpleName() + ": " + e.getMessage() );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
		}
		else
		{
			prior = DoesNothingLogPrior.defaultInstance;
		}
		if( beta[LearningPrinciple.PRIOR_INDEX] > 0 ) {
			try
			{
				prior.set( ((GenDisMixClassifierParameterSet) params).useOnlyFreeParameter(), score );
			}
			catch( Exception e )
			{
				NonParsableException n = new NonParsableException( "problem when setting the kind of parameter: " + e.getMessage() );
				n.setStackTrace( e.getStackTrace() );
				throw n;
			}
		}
	}

	/**
	 * This method creates an array of GenDisMixClassifiers by using the cross-product of the given
	 * {@link DifferentiableStatisticalModel}s.
	 * 
	 * @param params
	 *            the parameters that will be used in all classifiers
	 * @param prior
	 *            the prior that will be used in all classifiers
	 * @param weights
	 *            the weights that will be used in all classifiers
	 * @param functions
	 *            the {@link DifferentiableStatisticalModel}s
	 *            <ol>
	 *            <li> functions[i] are the {@link DifferentiableStatisticalModel}s that can be used for class i
	 *            <li> functions.length has to be at least 2
	 *            </ol>
	 * 
	 * @return an array of GenDisMixClassifiers
	 * 
	 * @throws CloneNotSupportedException
	 *             if the some item could not be cloned
	 * 
	 * @see GenDisMixClassifier#GenDisMixClassifier(GenDisMixClassifierParameterSet, LogPrior, double[],
	 *      DifferentiableStatisticalModel[])
	 */
	public static GenDisMixClassifier[] create( GenDisMixClassifierParameterSet params, LogPrior prior, double[] weights,
			DifferentiableStatisticalModel[]... functions ) throws CloneNotSupportedException
	{
		int anz = 1, counter2;
		int[] current = new int[functions.length], max = new int[functions.length];
		DifferentiableStatisticalModel[] sf = new DifferentiableStatisticalModel[functions.length];
		for( int counter1 = 0; counter1 < functions.length; counter1++ )
		{
			anz *= functions[counter1].length;
			max[counter1] = functions[counter1].length - 1;
		}

		GenDisMixClassifier[] erg = new GenDisMixClassifier[anz];
		anz = sf.length - 1;

		for( int counter1 = 0; counter1 < erg.length; counter1++ )
		{
			for( counter2 = 0; counter2 < sf.length; counter2++ )
			{
				sf[counter2] = functions[counter2][current[counter2]];
			}

			erg[counter1] = new GenDisMixClassifier( params, prior, weights, sf );

			counter2 = 0;
			while( counter2 < anz && current[counter2] == max[counter2] )
			{
				current[counter2++] = 0;
			}
			current[counter2]++;
		}

		return erg;
	}

	public String getInstanceName()
	{
		return super.getInstanceName() + " with weights=" + Arrays.toString( beta )
			+ (prior == null || prior == DoesNothingLogPrior.defaultInstance ? "" : (" and with " + prior.getInstanceName()));
	}
	
	/**
	* This method returns the number of used threads while optimization.
	* 
	* @return the number of used threads while optimization
	*/
	public int getNumberOfThreads() {
		return ((GenDisMixClassifierParameterSet) params).getNumberOfThreads();
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer( score.length * 5000 );
		String heading = "function ";
		for( int i =0; i < score.length; i++ ) {
			sb.append( heading + i );
			sb.append( "\n" + score[i].toString() + "\n" );
		}
		sb.append( "class weights: " );
		for(int i=0;i<getNumberOfClasses();i++){
			sb.append( getClassWeight( i )+" " );
		}
		sb.append( "\n" );
		return sb.toString();
	}
	
	/**
	 * This method allows to set the number of threads used while optimization.
	 * 
	 * @param threads the number of threads
	 * 
	 * @throws IllegalValueException if the number of threads is not possible
	 */
	public void setNumberOfThreads( int threads ) throws IllegalValueException
	{
		((GenDisMixClassifierParameterSet) params).setNumberOfThreads( threads );
	}
}
