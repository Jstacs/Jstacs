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

package de.jstacs.results;

import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * Abstract class for a {@link Result} with a value of a primitive data type or
 * {@link String}.
 * 
 * @see Result
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class SimpleResult extends Result implements
		Comparable<SimpleResult> {
	/**
	 * The result.
	 */
	protected Comparable result;

	/**
	 * The main constructor which takes the main information of a result.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            the comment for the result
	 * @param datatype
	 *            the data type of the result
	 */
	protected SimpleResult(String name, String comment, DataType datatype) {
		super(name, comment, datatype);
	}

	/**
	 * This is the constructor for {@link de.jstacs.Storable}. Creates a new
	 * {@link SimpleResult} out of the XML representation as returned by
	 * {@link #toXML()}.
	 * 
	 * @param rep
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the representation could not be parsed.
	 */
	protected SimpleResult(StringBuffer rep) throws NonParsableException {
		super(rep);
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		XMLParser.appendObjectWithTags(buf, result.toString(), "result");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		switch (datatype) {
			case STRING:
				result = XMLParser.extractObjectForTags(representation, "result", String.class );
				break;
			case BOOLEAN:
				result = XMLParser.extractObjectForTags(representation, "result", boolean.class );
				break;
			case INT:
				result = XMLParser.extractObjectForTags(representation, "result", int.class );
				break;
			case LONG:
				result = XMLParser.extractObjectForTags(representation, "result", long.class );
				break;
			case DOUBLE:
				result = XMLParser.extractObjectForTags(representation, "result", double.class );
				break;
			default:
				throw new NonParsableException("Result not of expected datatype");
		}
	}

	
	/**
	 * Sets the value of the result to <code>newValue</code>.
	 * 
	 * @param newValue
	 *            the new value of this result
	 * 
	 * @throws IllegalValueException
	 *             if <code>newValue</code> is not of the expected data type or
	 *             out of range
	 * 
	 * @deprecated
	 */
	@Deprecated
	void setResult(Comparable newValue) throws IllegalValueException {
		set(newValue);
	}

	/**
	 * Sets the value of this {@link SimpleResult} to <code>result</code>.
	 * 
	 * @param result
	 *            the result
	 * 
	 * @throws IllegalValueException
	 *             if <code>result</code> is an illegal value for this result
	 */
	private void set(Comparable result) throws IllegalValueException {
		if (datatype == DataType.DOUBLE
				&& (result instanceof Double || result instanceof Float)) {
			if (result instanceof Float) {
				this.result = new Double(((Float) result).doubleValue());
			} else {
				this.result = result;
			}
		} else if (datatype == DataType.INT
				&& (result instanceof Integer || result instanceof Byte || result instanceof Short)) {
			if (result instanceof Byte) {
				this.result = new Integer(((Byte) result).intValue());
			} else if (result instanceof Short) {
				this.result = new Integer(((Short) result).intValue());
			} else {
				this.result = result;
			}
		} else if (datatype == DataType.LONG && result instanceof Long) {
			this.result = result;
		} else if (datatype == DataType.BOOLEAN && result instanceof Boolean) {
			this.result = result;
		} else if (datatype == DataType.STRING && result instanceof String) {
			this.result = result;
		} else if (result instanceof String) {
			if (datatype == DataType.DOUBLE) {
				this.result = new Double((String) result);
			} else if (datatype == DataType.INT) {
				this.result = new Integer((String) result);
			} else if (datatype == DataType.BOOLEAN) {
				this.result = new Boolean((String) result);
			} else if (datatype == DataType.LONG) {
				this.result = new Long((String) result);
			} else {
				throw new IllegalValueException(
						"Value not of the expected datatype!");
			}
		} else {
			throw new IllegalValueException(
					"Value not of the expected datatype!");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#getResult()
	 */
	@Override
	public Comparable getValue() {
		return result;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return name + ": " + result.toString();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(SimpleResult r) {
		if (datatype != r.datatype) {
			return datatype.ordinal() - r.datatype.ordinal();
		}
		int c = name.compareTo(r.name);
		if (c != 0) {
			return c;
		}
		c = comment.compareTo(r.comment);
		if (c != 0) {
			return c;
		}
		return result.compareTo(r.result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object arg) {
		if (arg instanceof SimpleResult) {
			return compareTo((SimpleResult) arg) == 0;
		} else {
			return false;
		}
	}
}
