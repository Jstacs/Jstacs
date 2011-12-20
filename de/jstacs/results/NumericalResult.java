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
import de.jstacs.NonParsableException;

/**
 * Class for numerical {@link Result} values. Possible data types are the
 * numerical data types of {@link DataType}.
 * 
 * @author Jan Grau
 */
public class NumericalResult extends SimpleResult {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link NumericalResult} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public NumericalResult(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/**
	 * Creates a {@link NumericalResult} of a primitive numerical data type.
	 * 
	 * @param datatype
	 *            the data type
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 */
	protected NumericalResult(DataType datatype, String name, String comment, Comparable result) {
		super(name, comment, datatype);
		try {
			setResult(result);
		} catch (Exception e) {
			// impossible
			IllegalArgumentException n = new IllegalArgumentException(e
					.getMessage());
			n.setStackTrace(e.getStackTrace());
			throw n;
		}
	}

	/**
	 * The simplified constructor for the primitive type <code>double</code>.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 * 
	 * @see #NumericalResult(DataType, String, String, Comparable)
	 */
	public NumericalResult(String name, String comment, double result) {
		this(DataType.DOUBLE, name, comment, new Double(result));
	}

	/**
	 * The simplified constructor for the primitive type <code>int</code>.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 * 
	 * @see #NumericalResult(DataType, String, String, Comparable)
	 */
	public NumericalResult(String name, String comment, int result) {
		this(DataType.INT, name, comment, new Integer(result));
	}

	/**
	 * The simplified constructor for the type {@link Integer}.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 * 
	 * @see #NumericalResult(DataType, String, String, Comparable)
	 */
	public NumericalResult(String name, String comment, Integer result) {
		this(DataType.INT, name, comment, result);
	}

	/**
	 * The simplified constructor for the primitive type <code>long</code>.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 * 
	 * @see #NumericalResult(DataType, String, String, Comparable)
	 */
	public NumericalResult(String name, String comment, long result) {
		this(DataType.LONG, name, comment, new Long(result));
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.SimpleResult#toString()
	 */
	@Override
	public String toString() {
		return getResultString() + " \t= " + name + " \t(" + comment + ")";
	}

	private String getResultString() {
		String erg = "" + result;
		int l = 20 - erg.length();
		for (int i = 0; i < l; i++) {
			erg += " ";
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "numericalResult";
	}
}
