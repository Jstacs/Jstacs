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

package de.jstacs.classifiers.assessment;

import java.util.ArrayList;
import java.util.Collection;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.results.MeanResultSet;
import de.jstacs.results.Result;

/**
 * This class is the superclass used by all
 * {@link ClassifierAssessmentAssessParameterSet}s. It is a container of user
 * specified parameters that are necessary to define a run of an
 * <code>assess( ... )</code>-method in a {@link ClassifierAssessment} instance.
 * 
 * <br>
 * <br>
 * 
 * It is recommended to extend this class for each subclass of
 * {@link ClassifierAssessment}.
 * 
 * @see ClassifierAssessment#assess(de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet, ClassifierAssessmentAssessParameterSet, de.jstacs.data.DataSet...)
 * @see ClassifierAssessment#assess(de.jstacs.classifiers.performanceMeasures.NumericalPerformanceMeasureParameterSet, ClassifierAssessmentAssessParameterSet, de.jstacs.utils.ProgressUpdater, de.jstacs.data.DataSet...)
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com)
 * 
 */
public class ClassifierAssessmentAssessParameterSet extends ParameterSet {

	private void extend() {
		try {
			parameters.add(0, new SimpleParameter( DataType.BOOLEAN,
					"storeAll",
					"store all performance measures in each iteration seperately",
					true, false )
			);
		} catch ( Exception e ) {
			throw new RuntimeException( e );
		}
	}
	
	/**
	 * Constructs a new {@link ClassifierAssessmentAssessParameterSet} with
	 * empty parameter values.
	 * 
	 * @throws UnsupportedOperationException
	 *             if the {@link ClassifierAssessmentAssessParameterSet} could
	 *             not be constructed or the parameters could not be loaded
	 * 
	 * @see ParameterSet#ParameterSet()
	 */
	public ClassifierAssessmentAssessParameterSet() throws UnsupportedOperationException {
		super();
		try{
			extend();
			//0-subSequenceLength
			this.parameters.add( new SimpleParameter( DataType.INT,
					"elementLength",
					"Defines the lengths of overlapping windows " + "of data, that should be used to train and test " + "classifiers/models.",
					true ) );

			//1-excpetionIfMPNotComputable
			this.parameters.add( new SimpleParameter( DataType.BOOLEAN,
					"exceptionIfMeasureParamaterNotComputable",
					"True causes ClassiefierAssessment to throw " + "an error if measure-parameters should be computed "
					+ "that are not computable for the given classifiers.",
					true ) );
		}catch(DatatypeNotValidException doesnothappen){
			throw new RuntimeException( doesnothappen );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ClassifierAssessmentAssessParameterSet} out of its
	 * XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * @throws NonParsableException
	 *             if the {@link ClassifierAssessmentAssessParameterSet} could
	 *             not be reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see ParameterSet#ParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public ClassifierAssessmentAssessParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
		if( parameters.get("storeAll")==null ) {
			extend();
		}
	}

	/**
	 * Constructs a new {@link ClassifierAssessmentAssessParameterSet} with
	 * given parameter values.
	 * 
	 * @param elementLength
	 *            defines the length of elements (sequences) the classifiers to
	 *            be assessed are able to classify
	 * @param exceptionIfMPNotComputable
	 *            A {@link ClassifierAssessmentAssessParameterSet} is used in
	 *            combination with an
	 *            {@link de.jstacs.classifiers.performanceMeasures.AbstractPerformanceMeasure}-object to call
	 *            <code>assess( ... )</code>-methods of
	 *            {@link ClassifierAssessment}s. If
	 *            <code>exceptionIfMPNotComputable==true</code> an
	 *            {@link Exception} is thrown if a user selected measure
	 *            parameters that could not be computed.
	 * @throws IllegalValueException
	 *             is thrown in case of out-of-range or invalid given parameters
	 * @throws UnsupportedOperationException  if the {@link ClassifierAssessmentAssessParameterSet} could
	 *             not be constructed or the parameters could not be loaded
	 * 
	 * @see ParameterSet#ParameterSet()
	 */
	public ClassifierAssessmentAssessParameterSet( int elementLength, boolean exceptionIfMPNotComputable ) throws IllegalValueException, UnsupportedOperationException {
		this();
		this.parameters.get( "elementLength" ).setValue( new Integer( elementLength ) );
		this.parameters.get( "exceptionIfMeasureParamaterNotComputable" ).setValue( new Boolean( exceptionIfMPNotComputable ) );
	}

	//	**********************
	//	member methods
	//	**********************


	/**
	 * Returns the length of elements (sequences) defined by this
	 * {@link ClassifierAssessmentAssessParameterSet}.
	 * 
	 * @return the element length defined by this
	 *         {@link ClassifierAssessmentAssessParameterSet}
	 */
	public int getElementLength() {
		return ( (Integer)( this.getParameterForName( "elementLength" ).getValue() ) ).intValue();
	}

	/**
	 * Returns the flag defined by this
	 * {@link ClassifierAssessmentAssessParameterSet}.
	 * 
	 * @return the flag defined by this
	 *         {@link ClassifierAssessmentAssessParameterSet} (<code>true</code>
	 *         : an exception is thrown if a user selected measure parameters
	 *         that could not be computed, <code>false</code>: no exception is
	 *         thrown)
	 */
	public boolean getExceptionIfMPNotComputable() {
		return ( (Boolean)( getParameterForName( "exceptionIfMeasureParamaterNotComputable" ).getValue() ) ).booleanValue();
	}
	
	/**
	 * Returns the flag for storing all performance measures in each iteration defined by this
	 * {@link ClassifierAssessmentAssessParameterSet}.
	 * 
	 * @return the flag defined by this
	 *         {@link ClassifierAssessmentAssessParameterSet} 
	 *         (<code>true</code>: all performance measures in each iteration are stored separately,
	 *         <code>false</code>: only values need for mean and standard error are stored)
	 *         
	 * @see MeanResultSet#MeanResultSet(boolean, de.jstacs.results.SimpleResult...)
	 * @see MeanResultSet#getAllValues(int)
	 */
	public boolean getStoreAll() {
		return ( (Boolean)( getParameterForName( "storeAll" ).getValue() ) ).booleanValue();
	}

	/**
	 * Returns a {@link Collection} of parameters containing informations about
	 * this {@link ClassifierAssessmentAssessParameterSet}.
	 * 
	 * @return a {@link Collection} of parameters containing informations about
	 *         this {@link ClassifierAssessmentAssessParameterSet}
	 */
	public Collection<Result> getAnnotation() {
		return new ArrayList<Result>( 0 );
	}

	/**
	 * This method allows to set the switch for storing all individual performance measure values of each iteration of the {@link ClassifierAssessment}.
	 * @param b the value to be set
	 * @throws IllegalValueException if the parameter could not be set
	 */
	public void setStoreAll( boolean b ) throws IllegalValueException {
		 getParameterForName( "storeAll" ).setValue(b);		
	}
}
