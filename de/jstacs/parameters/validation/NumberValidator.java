/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.parameters.validation;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;

/**
 * Class that validates all subclasses of {@link Number} that implement
 * {@link Comparable} (e.g. <code>Double, Long, Float</code>) for compliance
 * with a specified lower and upper bound.
 * 
 * @author Jan Grau
 * 
 * @param <E>
 *            the subclass of {@link Number} to be validated
 */
public class NumberValidator<E extends Comparable<? extends Number>> implements
		ParameterValidator, GalaxyConvertible {

	/**
	 * The lower bound to check against
	 */
	private E lowerBound;
	/**
	 * The upper bound to check against
	 */
	private E upperBound;
	/**
	 * The class of <code>E</code>
	 */
	private Class clazz;

	/**
	 * The error message, <code>null</code> if no error occurred
	 */
	private String errorMessage;

	/**
	 * Constructs a {@link NumberValidator} for a given upper and lower bound.
	 * 
	 * @param lowerBound
	 *            the lower bound
	 * @param upperBound
	 *            the upper bound
	 */
	public NumberValidator(E lowerBound, E upperBound) {
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.clazz = lowerBound.getClass();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link NumberValidator} out of a XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the representation could not be parsed
	 */
	public NumberValidator(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	@SuppressWarnings("unchecked")
	public NumberValidator clone() {
		NumberValidator clone = new NumberValidator(lowerBound, upperBound);
		clone.errorMessage = errorMessage;
		return clone;

	}

	/**
	 * Returns the lower bound of the {@link NumberValidator}.
	 * 
	 * @return the lower bound of the {@link NumberValidator}
	 */
	public E getLowerBound() {
		return lowerBound;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.ParameterValidator#getErrorMessage()
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.ParameterValidator#checkValue
	 * (java.lang.Object)
	 */
	@SuppressWarnings("unchecked")
	public boolean checkValue(Object value) {
		if (value == null) {
			errorMessage = "Value is null.";
			return false;
		} else if (!clazz.isInstance(value)
				|| !(value instanceof Comparable && value instanceof Number)) {
			errorMessage = "Value is not of the correct numeric datatype but of "
					+ value.getClass() + ".";
			return false;
		} else {
			Comparable c = (Comparable) value;
			if (c.compareTo(lowerBound) < 0) {
				errorMessage = "Value " + c.toString()
						+ " is less than the lower bound " + lowerBound + ".";
				return false;
			} else if (c.compareTo(upperBound) > 0) {
				errorMessage = "Value " + c.toString()
						+ " is greater than the upper bound " + upperBound
						+ ".";
				return false;
			} else {
				errorMessage = null;
				return true;
			}
		}
	}

	/**
	 * Parses a {@link NumberValidator} from the XML representation as returned
	 * by {@link NumberValidator#toXML()}.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	@SuppressWarnings("unchecked")
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,"NumberValidator");
		String className = XMLParser.extractObjectForTags(representation, "className", String.class );
		try {
			clazz = Class.forName(className);

			String lower = XMLParser.extractObjectForTags(representation, "lowerBound", String.class );
			String upper = XMLParser.extractObjectForTags(representation, "upperBound", String.class );

			if (clazz.equals(Double.class)) {
				lowerBound = (E) new Double(lower);
				upperBound = (E) new Double(upper);
			} else if (clazz.equals(Float.class)) {
				lowerBound = (E) new Float(lower);
				upperBound = (E) new Float(upper);
			} else if (clazz.equals(Byte.class)) {
				lowerBound = (E) new Byte(lower);
				upperBound = (E) new Byte(upper);
			} else if (clazz.equals(Short.class)) {
				lowerBound = (E) new Short(lower);
				upperBound = (E) new Short(upper);
			} else if (clazz.equals(Integer.class)) {
				lowerBound = (E) new Integer(lower);
				upperBound = (E) new Integer(upper);
			} else if (clazz.equals(Long.class)) {
				lowerBound = (E) new Long(lower);
				upperBound = (E) new Long(upper);
			} else {
				throw new NonParsableException();
			}

		} catch (Exception e) {

		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();

		XMLParser.appendObjectWithTags(buf, clazz.getName(), "className");
		XMLParser.appendObjectWithTags(buf, lowerBound.toString(), "lowerBound");
		XMLParser.appendObjectWithTags(buf, upperBound.toString(), "upperBound");
		XMLParser.addTags(buf, "NumberValidator");

		return buf;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return "valid range = [" + lowerBound + ", " + upperBound + "]";
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception {
		XMLParser.addIndentation(descBuffer, indentation);
		descBuffer.append( "<validator type=\"in_range\" min=\""+lowerBound+"\" max=\""+upperBound+"\" message=\"Value is not in the specified range ["+lowerBound+", "+upperBound+"].\"/>\n" );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
				
	}

	@Override
	public void toGalaxyTest(String namePrefix, int depth, StringBuffer testBuffer, int indentation) throws Exception {}
}
