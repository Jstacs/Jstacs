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

package de.jstacs.classifier;

import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.MultiSelectionCollectionParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.RangeParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.RangeParameter.Scale;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.CategoricalResult;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.Result;

/**
 * This class holds the parameters for the <code>evaluate</code>-methods of a
 * classifier.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see AbstractClassifier#evaluate(MeasureParameters, boolean,
 *      de.jstacs.data.Sample...)
 * @see AbstractClassifier#evaluateAll(MeasureParameters, boolean,
 *      de.jstacs.data.Sample...)
 */
public final class MeasureParameters extends MultiSelectionCollectionParameter {

	/**
	 * This <code>enum</code> defines all measures that are currently
	 * implemented in Jstacs.
	 * 
	 * @author Jan Grau
	 */
	public enum Measure {

		/**
		 * The ratio of correctly classified sequences and all sequences.
		 */
		ClassificationRate( "Classification rate", "The ratio of correctly classified sequences and all sequences", true ),
		/**
		 * The sensitivity TP / (TP + FN) for a fixed, previously chosen
		 * specificity TN / (FP + TN).
		 */
		Sensitivity( "Sensitivity for fixed specificity",
				"The sensitivity TP / (TP + FN) for a fixed, previously chosen specificity TN / (FP + TN)", true ),
		/**
		 * The false positive rate FP / (FP + TN) for a fixed, previously chosen
		 * sensitivity TP / (TP + FN).
		 */
		FalsePositiveRate( "False positive rate for fixed sensitivity",
				"The false positive rate FP / (FP + TN) for a fixed, previously chosen sensitivity TP / (TP + FN)", true ),
		/**
		 * The positive predictive value TP / (TP+FP) for a fixed, previously
		 * chosen sensitivity TP / (TP + FN).
		 */
		PositivePredictiveValue( "Positive predictive value for fixed sensitivity",
				"The positive predictive value TP / (TP+FP) for a fixed, previously chosen sensitivity TP / (TP + FN)", true ),
		/**
		 * The area under the receiver operating characteristic (ROC) curve.
		 */
		AreaUnderROCCurve( "Area under ROC curve", "The area under the receiver operating characteristic (ROC) curve", true ),
		/**
		 * The area under the precision recall (PR) curve.
		 */
		AreaUnderPrecisionRecallCurve( "Area under PR curve", "The area under the precision recall (PR) curve", true),
		/**
		 * The maximum of all possible correlation coefficients.
		 */
		MaximumCorrelationCoefficient( "Maximum correlation coefficient", "The maximum of all possible correlation coefficients", true ),
		/**
		 * A part of the ROC curve on a grid between minimum and maximum specificity.
		 */
		PartialROCCurve( "Partial ROC curve", "A part of the ROC curve between minimum and maximum specificity", true ),
		/**
		 * The points of the receiver operating characteristic (ROC) curve.
		 */
		ReceiverOperatingCharacteristicCurve( "Receiver operating characteristic curve",
				"The points of the receiver operating characteristic (ROC) curve", false ),
		/**
		 * The points of the precision recall (PR) curve.
		 */
		PrecisionRecallCurve( "Precision recall curve", "The points of the precision recall (PR) curve", false );

		private String name, comment;
		private boolean isNumeric;

		Measure( String name, String comment, boolean isNumeric ) {
			this.name = name;
			this.comment = comment;
			this.isNumeric = isNumeric;
		}
		
		/**
		 * This method returns <code>true</code> if the measure returns a scalar that can be used in a {@link NumericalResult}.
		 * @return <code>true</code> if the measure returns a scalar that can be used in a {@link NumericalResult}
		 * 
		 * @see MeasureParameters#isNumeric()
		 */
		public boolean isNumeric() {
			return isNumeric;
		}

		/**
		 * Returns the name of the {@link Measure}.
		 * 
		 * @return the name
		 */
		public String getNameString() {
			return name;
		}

		/**
		 * Returns a comment on the {@link Measure}.
		 * 
		 * @return the comment
		 */
		public String getCommentString() {
			return comment;
		}

		private static String[] getMeasureNames() {
			Measure[] vals = Measure.values();
			String[] names = new String[vals.length];
			for( int i = 0; i < names.length; i++ ) {
				names[i] = vals[i].getNameString();
			}
			return names;
		}

		private static String[] getMeasureComments() {
			Measure[] vals = Measure.values();
			String[] comments = new String[vals.length];
			for( int i = 0; i < comments.length; i++ ) {
				comments[i] = vals[i].getCommentString();
			}
			return comments;
		}

	}

	private static final String[] getStringArray( String[] template, boolean all ) {
		if( all ) {
			return template;
		} else {
			String[] copy = new String[template.length - 2];
			System.arraycopy( template, 0, copy, 0, copy.length );
			return copy;
		}
	}

	private static final Object[] getValuesArray( boolean all, Double sp, Double snForFPR, Double snForPPV, Double minSpec, Double maxSpec,
			Integer steps, Scale scale ) throws ParameterException {
		Object[] values = new Object[8 + ( all ? 2 : 0 )];

		NumberValidator<Double> validator = new NumberValidator<Double>( new Double( 0 ), new Double( 1 ) );
		values[0] = values[4] = values[5] = values[6] = new SimpleParameterSet( new Parameter[]{} );
		if( sp == null ) {
			values[1] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Specificity",
					"The fixed specificity",
					true,
					validator ) } );
		} else {
			values[1] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Specificity",
					"The fixed specificity",
					true,
					validator,
					sp ) } );
		}
		if( snForFPR == null ) {
			values[2] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Sensitivity",
					"The fixed sensitivity",
					true,
					validator ) } );
		} else {
			values[2] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Sensitivity",
					"The fixed sensitivity",
					true,
					validator,
					snForFPR ) } );
		}
		if( snForPPV == null ) {
			values[3] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Sensitivity",
					"The fixed sensitivity",
					true,
					validator ) } );
		} else {
			values[3] = new SimpleParameterSet( new Parameter[]{ new SimpleParameter( DataType.DOUBLE,
					"Sensitivity",
					"The fixed sensitivity",
					true,
					validator,
					snForPPV ) } );
		}
		try {
			RangeParameter rp = new RangeParameter( new SimpleParameter( DataType.DOUBLE, "Specificity", "The specificity", true, validator ) );
			values[7] = new SimpleParameterSet( new Parameter[]{ rp } );

			if( minSpec != null && maxSpec != null && steps != null && scale != null ) {
				rp.setValues( minSpec, steps, maxSpec, scale );
			}
		} catch ( Exception e ) {
			e.printStackTrace();
			throw new ParameterException( e.getMessage() );
		}
		if( all ) {
			values[values.length - 2] = values[values.length - 1] = values[0];
		}

		return values;
	}

	/**
	 * Creates a new empty instance of {@link MeasureParameters}.
	 * 
	 * @param evaluateAll
	 *            indicates if the instance contains the parameters for the
	 *            method
	 *            {@link AbstractClassifier#evaluateAll(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 *            or only the parameters for the method
	 *            {@link AbstractClassifier#evaluate(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 * 
	 * @throws ParameterException
	 *             if something went wrong while constructing the parameters
	 * 
	 * @see Measure
	 * @see MultiSelectionCollectionParameter#MultiSelectionCollectionParameter(DataType,
	 *      String[], Object[], String[], String, String, boolean)
	 */
	public MeasureParameters( boolean evaluateAll ) throws ParameterException {
		super( DataType.PARAMETERSET,
				getStringArray( Measure.getMeasureNames(), evaluateAll ),
				getValuesArray( evaluateAll, null, null, null, null, null, null, null ),
				getStringArray( Measure.getMeasureComments(), evaluateAll ),
				"Performance measures",
				"Choose the desired performance measures",
				true );
	}

	/**
	 * Creates a new instance of {@link MeasureParameters}. All measures are
	 * switched on and the parameters of specificity and sensitivity are set to
	 * the given values.
	 * 
	 * @param evaluateAll
	 *            indicates if the instance contains the parameters for the
	 *            method
	 *            {@link AbstractClassifier#evaluateAll(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 *            or only the parameters for the method
	 *            {@link AbstractClassifier#evaluate(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 * @param sp
	 *            the (fixed) specificity
	 * @param snForFPR
	 *            the (fixed) sensitivity for the FPR
	 * @param snForPPV
	 *            the (fixed) sensitivity for the PPV
	 * 
	 * @throws ParameterException
	 *             if something went wrong while constructing the parameters
	 * 
	 * @see Measure
	 * @see MultiSelectionCollectionParameter#MultiSelectionCollectionParameter(DataType,
	 *      String[], Object[], String[], String, String, boolean)
	 */
	public MeasureParameters( boolean evaluateAll, double sp, double snForFPR, double snForPPV ) throws ParameterException {
		super( DataType.PARAMETERSET,
				getStringArray( Measure.getMeasureNames(), evaluateAll ),
				getValuesArray( evaluateAll, sp, snForFPR, snForPPV, null, null, null, null ),
				getStringArray( Measure.getMeasureComments(), evaluateAll ),
				"Performance measures",
				"Choose the desired performance measures",
				true );
		setValue( getStringArray( Measure.getMeasureNames(), evaluateAll ) );
		// may be changed in the future
		setSelected( Measure.PartialROCCurve, false );
	}

	/**
	 * Creates a new instance of {@link MeasureParameters}. All measures are
	 * switched on and the parameters of specificity and sensitivity are set to
	 * the given values.
	 * 
	 * @param evaluateAll
	 *            indicates if the instance contains the parameters for the
	 *            method
	 *            {@link AbstractClassifier#evaluateAll(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 *            or only the parameters for the method
	 *            {@link AbstractClassifier#evaluate(MeasureParameters, boolean, de.jstacs.data.Sample...)}
	 * @param sp
	 *            the (fixed) specificity
	 * @param snForFPR
	 *            the (fixed) sensitivity for the FPR
	 * @param snForPPV
	 *            the (fixed) sensitivity for the PPV
	 * @param minSpec
	 *            only used for {@link Measure#PartialROCCurve}: the minimal specificity
	 * @param maxSpec
	 *            only used for {@link Measure#PartialROCCurve}: the maximal specificity
	 * @param steps
	 *            only used for {@link Measure#PartialROCCurve}: the number of steps
	 * @param scale
	 *            only used for {@link Measure#PartialROCCurve}: see {@link RangeParameter}
	 * 
	 * @throws ParameterException
	 *             if something went wrong while constructing the parameters
	 * 
	 * @see Scale
	 * @see Measure
	 * @see RangeParameter
	 * @see MultiSelectionCollectionParameter#MultiSelectionCollectionParameter(DataType,
	 *      String[], Object[], String[], String, String, boolean)
	 */
	public MeasureParameters( boolean evaluateAll, double sp, double snForFPR, double snForPPV, double minSpec, double maxSpec, int steps,
								Scale scale ) throws ParameterException {
		super( DataType.PARAMETERSET,
				getStringArray( Measure.getMeasureNames(), evaluateAll ),
				getValuesArray( evaluateAll, sp, snForFPR, snForPPV, minSpec, maxSpec, steps, scale ),
				getStringArray( Measure.getMeasureNames(), evaluateAll ),
				"Performance measures",
				"Choose the desired performance measures",
				true );
		setValue( getStringArray( Measure.getMeasureNames(), evaluateAll ) );
	}
	
	/**
	 * This method returns <code>true</code> if all selected measures returns scalars. Hence, the result can be returned a {@link de.jstacs.results.NumericalResultSet}.
	 * @return <code>true</code> if all selected measure returns scalars
	 * 
	 * @see Measure#isNumeric()
	 * @see de.jstacs.results.NumericalResultSet
	 * @see de.jstacs.results.ResultSet
	 */
	public boolean isNumeric() {
		Measure[] m = Measure.values();
		for( int i = 0; i < m.length; i++ ) {
			if( !m[i].isNumeric() && isSelected( m[i] ) ) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Selects or deselects an option <code>sel</code> depending on the new
	 * selection <code>b</code>.
	 * 
	 * @param sel
	 *            the option
	 * @param b
	 *            the new selection
	 * 
	 * @see MeasureParameters#setSelected(String, boolean)
	 */
	public void setSelected( Measure sel, boolean b ) {
		this.setSelected( sel.getNameString(), b );
	}

	/**
	 * Indicates if the option <code>sel</code> is selected.
	 * 
	 * @param sel
	 *            the option
	 * 
	 * @return <code>true</code> if the option <code>sel</code> is selected,
	 *         <code>false</code> otherwise
	 * 
	 * @see MeasureParameters#isSelected(String)
	 */
	public boolean isSelected( Measure sel ) {
		return isSelected( sel.getNameString() );
	}

	/**
	 * Deselects all measures, i.e. calls {@link #setSelected(Measure, boolean)}
	 * for all {@link Measure}s.
	 * 
	 * @see MeasureParameters#setSelected(Measure, boolean)
	 */
	public void deselectAll() {
		Measure[] vals = Measure.values();
		for( int i = 0; i < vals.length; i++ ) {
			this.setSelected( vals[i], false );
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new instance of {@link MeasureParameters} out of its XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link MeasureParameters} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see MultiSelectionCollectionParameter#MultiSelectionCollectionParameter(StringBuffer)
	 * @see de.jstacs.Storable
	 */
	public MeasureParameters( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}

	/**
	 * Returns the selected parameters and their values as a {@link LinkedList}
	 * of {@link Result}s.
	 * 
	 * @return a {@link LinkedList} of all selected parameters and their values
	 */
	public LinkedList<Result> getAnnotation() {
		LinkedList<Result> list = new LinkedList<Result>();
		SimpleParameterSet sps;
		Parameter p;
		int i = 0, j, l = parameters.getNumberOfParameters();
		for( ; i < l; i++ ) {
			if( isSelected( i ) ) {
				p = parameters.getParameterAt( i );
				sps = (SimpleParameterSet)p.getValue();
				list.add( new CategoricalResult( p.getName(), p.getComment(), true ) );
				for( j = 0; j < sps.getNumberOfParameters(); j++ ) {
					p = sps.getParameterAt( j );
					if( p instanceof RangeParameter ) {
						list.add( new CategoricalResult( p.getName(), p.getComment(), ( (RangeParameter)p ).valuesToString() ) );
					} else {
						list.add( new NumericalResult( p.getName(), p.getComment(), ( (Number)p.getValue() ).doubleValue() ) );
					}
				}
			}
		}
		return list;
	}
}
