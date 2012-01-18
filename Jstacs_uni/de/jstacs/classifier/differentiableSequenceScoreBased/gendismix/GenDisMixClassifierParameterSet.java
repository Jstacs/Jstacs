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

package de.jstacs.classifier.differentiableSequenceScoreBased.gendismix;

import de.jstacs.DataType;
import de.jstacs.classifier.differentiableSequenceScoreBased.ScoreClassifier;
import de.jstacs.classifier.differentiableSequenceScoreBased.ScoreClassifierParameterSet;
import de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * This class contains the parameters for the {@link GenDisMixClassifier}.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class GenDisMixClassifierParameterSet extends ScoreClassifierParameterSet {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link GenDisMixClassifierParameterSet} out of its XML
	 * representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link GenDisMixClassifierParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} could not be parsed)
	 * 
	 * @see ScoreClassifierParameterSet#ScoreClassifierParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public GenDisMixClassifierParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
		if( parameters != null && parameters.size() == 7 ) {//TODO remove later
			parameters.add( getThreadsParameter() );
		}
	}

	/**
	 * The default constructor that constructs a new
	 * {@link GenDisMixClassifierParameterSet}.
	 * 
	 * @param alphabet
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the sequences
	 * @param algo
	 *            the algorithm that shall be used for optimization
	 * @param eps
	 *            the epsilon for stopping the optimization
	 * @param lineps
	 *            the epsilon for stopping the line search
	 * @param startD
	 *            the start distance for the line search
	 * @param free
	 *            the switch for using only the free or all parameters in a
	 *            {@link de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore}
	 * @param kind
	 *            indicates the kind of class parameter initialization
	 * @param norm
	 *            the switch for using a normalization while optimization
	 * @param threads
	 * 			  the {@link NumberFormatException} of threads used during an optimization
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see KindOfParameter
	 * @see GenDisMixClassifierParameterSet#GenDisMixClassifierParameterSet(Class,
	 *      de.jstacs.data.AlphabetContainer, int, byte, double, double, double, boolean,
	 *      de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter, boolean,int)
	 */
	public GenDisMixClassifierParameterSet( AlphabetContainer alphabet, int length, byte algo, double eps, double lineps, double startD,
										boolean free, KindOfParameter kind, boolean norm, int threads ) throws Exception {
		this( GenDisMixClassifier.class, alphabet, length, algo, eps, lineps, startD, free, kind, norm, threads );
	}

	/**
	 * The default constructor that constructs a new
	 * {@link GenDisMixClassifierParameterSet}.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param alphabet
	 *            the {@link AlphabetContainer}
	 * @param length
	 *            the length of the sequences
	 * @param algo
	 *            the algorithm that shall be used for optimization
	 * @param eps
	 *            the epsilon for stopping the optimization
	 * @param lineps
	 *            the epsilon for stopping the line search
	 * @param startD
	 *            the start distance for the line search
	 * @param free
	 *            the switch for using only the free or all parameters in a
	 *            {@link de.jstacs.sequenceScores.differentiable.DifferentiableSequenceScore}
	 * @param kind
	 *            indicates the kind of class parameter initialization
	 * @param norm
	 *            the switch for using a normalization while optimization
	 * @param threads
	 * 			  the {@link NumberFormatException} of threads used during an optimization
	 * 
	 * @throws Exception
	 *             if something went wrong
	 * 
	 * @see KindOfParameter
	 * @see ScoreClassifierParameterSet#ScoreClassifierParameterSet(Class,
	 *      de.jstacs.data.AlphabetContainer, int, byte, double, double, double, boolean,
	 *      de.jstacs.classifier.differentiableSequenceScoreBased.OptimizableFunction.KindOfParameter)
	 */
	protected GenDisMixClassifierParameterSet( Class<? extends ScoreClassifier> instanceClass, AlphabetContainer alphabet, int length, byte algo,
											double eps, double lineps, double startD, boolean free, KindOfParameter kind, boolean norm, int threads )
																																		throws Exception {
		super( instanceClass, alphabet, length, algo, eps, lineps, startD, free, kind );
		parameters.add( new SimpleParameter( DataType.BOOLEAN,
				"Normalize",
				"If true the conditional likelihood will be normalized to the number of samples.",
				true,
				new Boolean( true ) ) );
		parameters.add( getThreadsParameter() );
		getParameterForName( "Normalize" ).setValue( norm );
		setNumberOfThreads( threads );
	}
	
	private static Parameter getThreadsParameter() {
		try {
			return new SimpleParameter( DataType.INT,
					"Threads",
					"The number of threads that is used during an optimization.",
					true,
					new NumberValidator<Integer>(1,128),
					1 );
		} catch( Exception e ) {
			throw new RuntimeException( e.getMessage() );
		}
	}

	/**
	 * This method indicates if a normalization shall be used while
	 * optimization. The normalization is done through division by the number of
	 * sequences.
	 * 
	 * @return <code>true</code> if a normalization shall be used while
	 *         optimization, <code>false</code> otherwise
	 */
	public boolean shouldBeNormalized() {
		return (Boolean)getParameterForName( "Normalize" ).getValue();
	}
	
	/**
	 * This method returns the number of threads that should be used during optimization.
	 * 
	 * @return the number of threads that should be used during optimization
	 */
	public int getNumberOfThreads() {
		return (Integer)getParameterForName( "Threads" ).getValue();
	}
	
	/***
	 * This method set the number of threads used during optimization.
	 * 
	 * @param threads the number of threads used during optimization
	 * 
	 * @throws IllegalValueException if the value could not be set
	 */
	public void setNumberOfThreads( int threads ) throws IllegalValueException {
		getParameterForName( "Threads" ).setValue( threads );
	}
}
