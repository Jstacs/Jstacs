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

package de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore;
import de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableStatisticalModel;

/**
 * Abstract class for priors that penalize each parameter value independently
 * and have some variances (and possible means) as hyperparameters. Such priors
 * are e.g. the {@link SeparateGaussianLogPrior} or the
 * {@link SeparateLaplaceLogPrior}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class SeparateLogPrior extends LogPrior implements Cloneable {

	/**
	 * The {@link DifferentiableSequenceScore}s using the parameters that shall be
	 * penalized.
	 */
	protected DifferentiableStatisticalModel[] funs;

	/**
	 * The base variances for the parameters of the {@link DifferentiableSequenceScore}s of
	 * each class, the means of the non-class parameters should be 0.
	 */
	protected double[] vars;

	/**
	 * The variances for the class parameters, as specified by the user.
	 */
	protected double[] classVars;

	/**
	 * The means for the class parameters, as specified by the user.
	 */
	protected double[] classMus;

	/**
	 * Indicates, if only free parameters shall be used and hence penalized.
	 */
	protected boolean freeParameters;

	/**
	 * Creates a new {@link SeparateLogPrior} using the class-specific base
	 * variances <code>vars</code>, the variances <code>classVars</code> and the
	 * means <code>classMus</code> for the class parameters.
	 * 
	 * @param vars
	 *            the base variances for each class
	 * @param classVars
	 *            the variances for the class parameters
	 * @param classMus
	 *            the means for the class parameters
	 */
	public SeparateLogPrior( double[] vars, double[] classVars, double[] classMus ) {
		this.vars = vars;
		this.classMus = classMus;
		this.classVars = classVars;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link SeparateLogPrior} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SeparateLogPrior} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public SeparateLogPrior( StringBuffer xml ) throws NonParsableException {
		StringBuffer content = XMLParser.extractForTag( xml, this.getClass().getSimpleName() );
		vars = XMLParser.extractObjectForTags( content, "variances", double[].class );
		classVars = XMLParser.extractObjectForTags( content, "classVariances", double[].class );
		classMus = XMLParser.extractObjectForTags( content, "classMeans", double[].class );
		unset();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#set(boolean, de.jstacs.sequenceScores.statisticalModels.differentiable.DifferentiableSequenceScore[])
	 */
	@Override
	public void set( boolean freeParameters, DifferentiableSequenceScore... funs ) throws Exception {
		this.funs = new DifferentiableStatisticalModel[funs.length];
		for( int i = 0; i < funs.length; i++ ) {
			if( !( funs[i] instanceof DifferentiableStatisticalModel ) ) {
				throw new Exception( "Only DirectedGraphicalModels allowed." );
			} else {
				this.funs[i] = (DifferentiableStatisticalModel)funs[i];
			}
		}
		this.freeParameters = freeParameters;
		this.unset();
	}

	/**
	 * Resets all internally pre-computed values, e.g. the hyperparameters for
	 * each parameter.
	 */
	protected abstract void unset();

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#getNewInstance()
	 */
	@Override
	public SeparateLogPrior getNewInstance() throws CloneNotSupportedException {
		SeparateLogPrior clone = (SeparateLogPrior)super.clone();
		clone.vars = vars.clone();
		clone.unset();
		clone.funs = null;
		clone.classVars = classVars.clone();
		clone.classMus = classMus.clone();
		return clone;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifiers.differentiableSequenceScoreBased.logPrior.LogPrior#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 20 * ( vars.length + classVars.length + classMus.length ) + 1000 );
		XMLParser.appendObjectWithTags( xml, vars, "variances" );
		XMLParser.appendObjectWithTags( xml, classVars, "classVariances" );
		XMLParser.appendObjectWithTags( xml, classMus, "classMeans" );
		XMLParser.addTags( xml, this.getClass().getSimpleName() );
		return xml;
	}
}
