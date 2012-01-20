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

package de.jstacs.data.alphabets;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * Class for a continuous alphabet.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class ContinuousAlphabet extends Alphabet {

	private static final String XML_TAG = "ContinuousAlphabet";

	/**
	 * The minimal and the maximal value.
	 */
	private double min, max;
	
	/**
	 * Allow NaN values?
	 */
	private boolean allowNaN;

	private ContinuousAlphabetParameterSet parameters;

	/**
	 * The constructor for the {@link de.jstacs.InstantiableFromParameterSet}
	 * interface. Creates a new {@link ContinuousAlphabet} from a given set of
	 * parameters.
	 * 
	 * @param parameters
	 *            the parameter set for the {@link ContinuousAlphabet}
	 * 
	 * @see de.jstacs.InstantiableFromParameterSet
	 */
	public ContinuousAlphabet( ContinuousAlphabetParameterSet parameters ) {
		this( (Double)parameters.getParameterAt( 0 ).getValue(), (Double)parameters.getParameterAt( 1 ).getValue(), (Boolean)parameters.getParameterAt( 2 ).getValue() );
		try {
			this.parameters = (ContinuousAlphabetParameterSet)parameters.clone();
		} catch ( Exception e ) {
			this.parameters = null;
		}
	}

	/**
	 * Creates a new {@link ContinuousAlphabet} from a minimal and a maximal
	 * value. NaN values are not allowed.
	 * 
	 * @param min
	 *            the minimal value
	 * @param max
	 *            the maximal value
	 * 
	 * @throws IllegalArgumentException
	 *             if the minimum or the maximum could not be set
	 */
	public ContinuousAlphabet( double min, double max ) throws IllegalArgumentException {
		this(min,max,false);
	}
	
	/**
	 * Creates a new {@link ContinuousAlphabet} from a minimal and a maximal
	 * value.
	 * 
	 * @param min
	 *            the minimal value
	 * @param max
	 *            the maximal value
	 * @param allowNaN
	 * 			  true if NaN values are allowed
	 * 
	 * @throws IllegalArgumentException
	 *             if the minimum or the maximum could not be set
	 */
	public ContinuousAlphabet( double min, double max, boolean allowNaN ) throws IllegalArgumentException {
		if( Double.isInfinite( min ) || Double.isNaN( min ) || Double.isNaN( max ) || Double.isInfinite( max ) ) {
			throw new IllegalArgumentException( "min and max have to be numbers (not infinity, NaN, ...)" );
		}
		if( min >= max ) {
			throw new IllegalArgumentException( "constraint violated: min < max" );
		}
		this.min = min;
		this.max = max;
		this.allowNaN = allowNaN;
	}

	/**
	 * Creates a new {@link ContinuousAlphabet} with minimum and maximum value
	 * being -{@link Double#MAX_VALUE} and {@link Double#MAX_VALUE},
	 * respectively. NaN values are not allowed.
	 * 
	 * @throws IllegalArgumentException
	 *             if the minimum or the maximum could not be set
	 * 
	 * @see ContinuousAlphabet#ContinuousAlphabet(double, double)
	 */
	public ContinuousAlphabet() {
		this(false);
	}
	
	/**
	 * Creates a new {@link ContinuousAlphabet} with minimum and maximum value
	 * being -{@link Double#MAX_VALUE} and {@link Double#MAX_VALUE},
	 * respectively.
	 * 
	 * @param allowNaN
	 * 			  true if NaN values are allowed
	 * 
	 * @throws IllegalArgumentException
	 *             if the minimum or the maximum could not be set
	 * 
	 * @see ContinuousAlphabet#ContinuousAlphabet(double, double)
	 */
	public ContinuousAlphabet(boolean allowNaN) {
		this( -Double.MAX_VALUE, Double.MAX_VALUE, allowNaN );
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ContinuousAlphabet} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ContinuousAlphabet} could not be reconstructed
	 *             out of the XML representation (the {@link StringBuffer} could
	 *             not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public ContinuousAlphabet( StringBuffer xml ) throws NonParsableException {
		StringBuffer help = XMLParser.extractForTag( xml, XML_TAG );
		min = XMLParser.extractObjectForTags( help, "MIN", double.class );
		max = XMLParser.extractObjectForTags( help, "MAX", double.class );
		Object o = XMLParser.extractObjectForTags( help, "allowNaN", boolean.class );
		if(o != null){
			allowNaN = (Boolean)o;
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public ContinuousAlphabetParameterSet getCurrentParameterSet() throws Exception {
		if( parameters != null ) {
			return (ContinuousAlphabetParameterSet)parameters.clone();
		} else {
			return new ContinuousAlphabetParameterSet( min, max, allowNaN );
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 100 );
		xml.append( "\t" );
		XMLParser.appendObjectWithTags( xml, min, "MIN" );
		xml.append( "\t" );
		XMLParser.appendObjectWithTags( xml, max, "MAX" );
		xml.append( "\t" );
		XMLParser.appendObjectWithTags( xml, allowNaN, "allowNaN" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo( Alphabet a ) {
		if( !getClass().equals( a.getClass() ) ) {
			return getClass().getName().compareTo( a.getClass().getName() );
		} else {
			ContinuousAlphabet b = (ContinuousAlphabet)a;
			int s = (int)Math.signum( min - b.min );
			if( s == 0 ) {
				return (int)Math.signum( max - b.max );
			} else {
				return s;
			}
		}
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.Alphabet#getMin()
	 */
	@Override
	public double getMin() {
		return min;
	}

	/**
	 * Returns the maximal value of this alphabet.
	 * 
	 * @return the maximal value of this alphabet
	 */
	public double getMax() {
		return max;
	}

	/**
	 * Indicates if <code>candidat</code> is an element of the internal
	 * interval.
	 * 
	 * @param candidat
	 *            the value to be tested
	 * 
	 * @return <code>true</code> if <code>candidat</code> is an element of the
	 *         internal interval, <code>false</code> otherwise
	 */
	public final boolean isEncodedSymbol( double candidat ) {
		return (Double.isNaN( candidat ) && allowNaN) || ( ( min <= candidat ) && ( candidat <= max ) );
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.Alphabet#length()
	 */
	@Override
	public double length() {
		return max - min;
	}

	/* (non-Javadoc)
	 * @see de.jstacs.data.alphabets.Alphabet#toString()
	 */
	@Override
	public String toString() {
		return "[" + min + ", " + max + "]";
	}

	/**
	 * Class for the {@link de.jstacs.parameters.ParameterSet} of a
	 * {@link ContinuousAlphabet}.
	 * 
	 * @author Jan Grau
	 */
	public static class ContinuousAlphabetParameterSet extends AlphabetParameterSet {

		/**
		 * Creates a new {@link ContinuousAlphabetParameterSet} with empty
		 * values.
		 * @see de.jstacs.data.alphabets.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(Class) Alphabet.AlphabetParameterSet#AlphabetParameterSet(Class)
		 */
		public ContinuousAlphabetParameterSet() {
			super( ContinuousAlphabet.class );
			try {
				parameters.add( new SimpleParameter( DataType.DOUBLE, "Minimum", "The minimal value of the alphabet.", true ) );
				parameters.add( new SimpleParameter( DataType.DOUBLE, "Maximum", "The maximum value of the alphabet", true ) );
				parameters.add( new SimpleParameter( DataType.BOOLEAN, "Allow NaN", "Allow NaN values", true, false ) );
			} catch ( DatatypeNotValidException doesnothappen ) {
				throw new RuntimeException( doesnothappen );
			}catch ( IllegalValueException doesnothappen ) {
				throw new RuntimeException( doesnothappen );
			}
		}

		/**
		 * Creates a new {@link ContinuousAlphabetParameterSet} from a minimum
		 * and a maximum value. NaN values are not allowed.
		 * 
		 * @param min
		 *            the minimum
		 * @param max
		 *            the maximum
		 * 
		 * @throws Exception
		 *             if minimum or maximum could not be set
		 * 
		 * @see de.jstacs.data.alphabets.ContinuousAlphabet.ContinuousAlphabetParameterSet#ContinuousAlphabet.ContinuousAlphabetParameterSet() ContinuousAlphabet.ContinuousAlphabetParameterSet#ContinuousAlphabetParameterSet()
		 */
		public ContinuousAlphabetParameterSet( double min, double max ) throws Exception {
			this(min,max,false);
		}
		
		/**
		 * Creates a new {@link ContinuousAlphabetParameterSet} from a minimum
		 * and a maximum value.
		 * 
		 * @param min
		 *            the minimum
		 * @param max
		 *            the maximum
		 * @param allowNaN
		 * 			  true if NaN values are allowed
		 * 
		 * @throws Exception
		 *             if minimum or maximum could not be set
		 * 
		 * @see de.jstacs.data.alphabets.ContinuousAlphabet.ContinuousAlphabetParameterSet#ContinuousAlphabet.ContinuousAlphabetParameterSet() ContinuousAlphabet.ContinuousAlphabetParameterSet#ContinuousAlphabetParameterSet()
		 */
		public ContinuousAlphabetParameterSet( double min, double max, boolean allowNaN ) throws Exception {
			this();
			parameters.get( 0 ).setValue( new Double( min ) );
			parameters.get( 1 ).setValue( new Double( max ) );
			parameters.get( 2 ).setValue( allowNaN );
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Creates a new {@link ContinuousAlphabetParameterSet} out of its XML
		 * representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link ContinuousAlphabetParameterSet} could not
		 *             be reconstructed out of the XML representation (the
		 *             {@link StringBuffer} <code>representation</code> could
		 *             not be parsed)
		 * 
		 * @see de.jstacs.data.alphabets.Alphabet.AlphabetParameterSet#Alphabet.AlphabetParameterSet(StringBuffer) Alphabet.AlphabetParameterSet#AlphabetParameterSet(StringBuffer)
		 * @see de.jstacs.Storable
		 */
		public ContinuousAlphabetParameterSet( StringBuffer representation ) throws NonParsableException {
			super( representation );
		}

		/* (non-Javadoc)
		 * @see de.jstacs.parameters.InstanceParameterSet#getInstanceComment()
		 */
		@Override
		public String getInstanceComment() {
			return "An alphabet that consists of real values between a minimum and a maximum value.";
		}
	}
}