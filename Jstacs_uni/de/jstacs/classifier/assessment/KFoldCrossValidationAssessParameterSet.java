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

package de.jstacs.classifier.assessment;

import java.util.ArrayList;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.DataSet.PartitionMethod;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class implements a {@link ClassifierAssessmentAssessParameterSet} that
 * must be used to call method <code>assess( ... )</code> of a
 * {@link KFoldCrossValidation}.<br>
 * It contains user specific parameters necessary for a run of a
 * {@link KFoldCrossValidation}.
 * 
 * @author Andre Gohr (bioinf (nospam:.) ag (nospam:@) googlemail (nospam:.)
 *         com)
 * 
 */
public class KFoldCrossValidationAssessParameterSet extends ClassifierAssessmentAssessParameterSet {

	/**
	 * Constructs a new {@link KFoldCrossValidationAssessParameterSet} with empty parameter
	 * values. This constructor should only be used to create &quot;filled&quot;
	 * {@link KFoldCrossValidationAssessParameterSet}s, i.e. to create
	 * {@link KFoldCrossValidationAssessParameterSet}s from a set of values and not to fill
	 * it from the platform user interface.
	 * 
	 * @throws UnsupportedOperationException
	 *             if the {@link KFoldCrossValidationAssessParameterSet} could not be
	 *             constructed or the parameters could not be loaded
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet()
	 */
	public KFoldCrossValidationAssessParameterSet() throws UnsupportedOperationException {
		super();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link KFoldCrossValidationAssessParameterSet} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * @throws NonParsableException
	 *             if the {@link KFoldCrossValidationAssessParameterSet} could not be
	 *             reconstructed out of the XML representation (the
	 *             {@link StringBuffer} <code>representation</code> could not be
	 *             parsed)
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public KFoldCrossValidationAssessParameterSet( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Constructs a new {@link KFoldCrossValidationAssessParameterSet} with given parameter
	 * values.
	 * 
	 * @param dataSplitMethod
	 *            defines the method used to split user supplied data into
	 *            <code>k</code> mutually exclusive random-splits (available
	 *            options are:
	 *            {@link PartitionMethod#PARTITION_BY_NUMBER_OF_ELEMENTS} and
	 *            {@link PartitionMethod#PARTITION_BY_NUMBER_OF_SYMBOLS})
	 * @param elementLength
	 *            defines the length of elements (sequences) the classifiers to
	 *            be assessed are able to classify
	 * @param exceptionIfMPNotComputable
	 *            a {@link KFoldCrossValidationAssessParameterSet} is used in combination
	 *            with a {@link de.jstacs.classifier.performanceMeasures.NumericalPerformanceMeasureParameters}-object
	 *            to call <code>assess( ... )</code>-methods of
	 *            {@link KFoldCrossValidation}s, if
	 *            <code>exceptionIfMPNotComputable==true</code> an
	 *            {@link Exception} is thrown in case of a user selected measure
	 *            parameters that could not be computed.
	 * @param k
	 *            defines the number of mutually exclusive random-splits of user
	 *            supplied data. Each part is used once as a test dataset and
	 *            the union of the remaining k-1 parts is once used as train
	 *            dataset. Thus <code>k</code> also defines how many (
	 *            <code>k</code>) repeated classifier trainings and classifier
	 *            evaluations (tests) are performed.
	 * 
	 * @throws IllegalValueException
	 *             in case of out-of-range or invalid given parameters
	 * 
	 * @see ClassifierAssessmentAssessParameterSet#ClassifierAssessmentAssessParameterSet(int,
	 *      boolean)
	 * @see de.jstacs.data.DataSet.PartitionMethod
	 */
	public KFoldCrossValidationAssessParameterSet( PartitionMethod dataSplitMethod, int elementLength, boolean exceptionIfMPNotComputable, int k )
																																		throws IllegalValueException {
		super( elementLength, exceptionIfMPNotComputable );

		( this.parameters.get( 2 ) ).setValue( new Integer( k ) );

		( this.parameters.get( 3 ) ).setValue( dataSplitMethod );
	}

	//	**********************
	//	member methods
	//	**********************

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#initializeMyParametersArrayList()
	 */
	@Override
	protected void initializeMyParametersArrayList() {
		initParameterList( 4 );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#loadParameters()
	 */
	@Override
	protected void loadParameters() throws Exception {
		super.loadParameters();

		//2-k
		this.parameters.add( new SimpleParameter( DataType.INT,
				"k",
				"Defines the folds of a KFoldCrossValidation.",
				true,
				new NumberValidator<Integer>( 2, Integer.MAX_VALUE ) ) );

		//3-dataSplitMethod
		this.parameters.add( new EnumParameter( PartitionMethod.class, "The method used to compute the percentages of the partitions", true ) );
	}

	/**
	 * Returns the number of mutually exclusive random-splits of user supplied
	 * data defined by this {@link KFoldCrossValidationAssessParameterSet}.
	 * 
	 * @return the number of mutually exclusive random-splits of user-supplied
	 *         data defined by this {@link KFoldCrossValidationAssessParameterSet}
	 */
	public int getK() {
		return ( (Integer)( this.getParameterAt( 2 ).getValue() ) ).intValue();
	}

	/**
	 * Returns the {@link PartitionMethod} defining how the mutually exclusive
	 * random-splits of user supplied data are generated.
	 * 
	 * @return the {@link PartitionMethod} defining how the mutually exclusive
	 *         random-splits of user supplied data are generated
	 * 
	 * @see de.jstacs.data.DataSet.PartitionMethod
	 */
	public PartitionMethod getDataSplitMethod() {
		return (PartitionMethod)( (EnumParameter)getParameterAt( 3 ) ).getValue();
	}

	/* (non-Javadoc)
	 * @see de.jstacs.classifier.assessment.ClassifierAssessmentAssessParameterSet#getAnnotation()
	 */
	@Override
	public ArrayList<Result> getAnnotation() {
		ArrayList<Result> l = new ArrayList<Result>( 2 );
		l.add( new NumericalResult( "k", "The folds of a KFoldCrossValidation.", getK() ) );
		l.add( new CategoricalResult( "dataSplitMethod",
				"Describes how data should be splitted in ClassifierAssessment.evaluateClassifier())",
				getDataSplitMethod().name() ) );
		return l;
	}

}
