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

package de.jstacs.classifiers.differentiableSequenceScoreBased.msp;

import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifierParameterSet;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.LearningPrinciple;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior;
import de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior;
import de.jstacs.io.NonParsableException;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;

/**
 * This class implements a classifier that allows the training via MCL or MSP principle.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class MSPClassifier extends GenDisMixClassifier {

	/**
	 * This convenience constructor creates an {@link MSPClassifier} that used MCL principle for training.
	 * 
	 * @param params
	 *            the parameter set for the classifier
	 * @param score
	 *            the {@link DifferentiableSequenceScore}s for the classes
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one {@link DifferentiableSequenceScore} could not be cloned
	 * 
	 * @see DoesNothingLogPrior
	 * @see MSPClassifier#MSPClassifier(GenDisMixClassifierParameterSet, LogPrior, DifferentiableSequenceScore...)
	 */
	public MSPClassifier( GenDisMixClassifierParameterSet params, DifferentiableSequenceScore... score ) throws CloneNotSupportedException {
		this( params, null, score );
	}

	/**
	 * The default constructor that creates a new {@link MSPClassifier} from a
	 * given parameter set, a prior and {@link DifferentiableSequenceScore}s for the
	 * classes.
	 * 
	 * @param params
	 *            the parameter set for the classifier
	 * @param prior
	 *            the prior that shall be used; the learning principle depends on the prior: if the prior is <code>null</code> or {@link de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior} then MCL is used otherwise MSP
	 * @param score
	 *            the {@link DifferentiableSequenceScore}s for the classes
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one {@link DifferentiableSequenceScore} could not be cloned
	 * 
	 * @see MSPClassifier#MSPClassifier(GenDisMixClassifierParameterSet, LogPrior, double, DifferentiableSequenceScore...)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier#NOT_TRAINED_VALUE
	 * @see MSPClassifier#setPrior(LogPrior)
	 */
	public MSPClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, DifferentiableSequenceScore... score ) throws CloneNotSupportedException {
		this( params, prior, NOT_TRAINED_VALUE, score );
	}
	
	/**
	 * This constructor that creates a new {@link MSPClassifier} from a
	 * given parameter set, a prior and {@link DifferentiableSequenceScore}s for the
	 * classes. Additionally, the value <code>lastScore</code> can be set to determine the return value of
	 * {@link de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier#getLastScore()}. This might be useful if the parameters of the classifier are
	 * determined by an external procedure. 
	 * 
	 * @param params
	 *            the parameter set for the classifier
	 * @param prior
	 *            the prior that shall be used; the learning principle depends on the prior: if the prior is <code>null</code> or {@link de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.DoesNothingLogPrior} then MCL is used otherwise MSP
	 * @param lastScore
	 * 			  the score of the last optimization  
	 * @param score
	 *            the {@link DifferentiableSequenceScore}s for the classes
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one {@link DifferentiableSequenceScore} could not be cloned
	 * 
	 * @see GenDisMixClassifier#GenDisMixClassifier(GenDisMixClassifierParameterSet, LogPrior, LearningPrinciple, de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel...)
	 * @see LearningPrinciple#MSP
	 * @see MSPClassifier#setPrior(LogPrior)
	 */
	public MSPClassifier( GenDisMixClassifierParameterSet params, LogPrior prior, double lastScore, DifferentiableSequenceScore... score ) throws CloneNotSupportedException {
		super( params, prior, lastScore, LearningPrinciple.getBeta( LearningPrinciple.MSP ), score );
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}.
	 * 
	 * @param xml the xml representation
	 * 
	 * @throws NonParsableException if the representation could not be parsed.
	 */
	public MSPClassifier( StringBuffer xml ) throws NonParsableException
	{
		super( xml );
	}
	
	private static final String XML_TAG = "cll-classifier";

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier#getXMLTag()
	 */
	@Override
	protected String getXMLTag() {
		return XML_TAG;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.ScoreClassifier#getInstanceName()
	 */
	@Override
	public String getInstanceName() {
		return getClass().getSimpleName() + ( prior == null || prior == DoesNothingLogPrior.defaultInstance ? ""
																										: " with " + prior.getInstanceName() );
	}
}
