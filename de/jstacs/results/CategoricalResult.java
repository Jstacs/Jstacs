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
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * A class for categorical results (i.e. non-numerical results) for primitives
 * and {@link String}s.
 * 
 * @author Jan Grau
 */
public class CategoricalResult extends SimpleResult {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link CategoricalResult} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public CategoricalResult(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/**
	 * Creates a result of a primitive categorical data type or a {@link String}
	 * .
	 * 
	 * @param datatype
	 *            the primitive data type
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 * 
	 * @throws IllegalValueException
	 *             if the result value is not of the expected data type
	 */
	public CategoricalResult(DataType datatype, String name, String comment,
			Comparable result) throws IllegalValueException {
		super(name, comment, datatype);
		if (!(datatype == DataType.STRING || datatype == DataType.BOOLEAN)) {
			throw new IllegalValueException(name, "Datatype must be categorical!");
		}
		setResult(result);
	}

	/**
	 * Creates a result of a {@link String}.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 */
	public CategoricalResult(String name, String comment, String result) {
		super(name, comment, DataType.STRING);
		this.result = result;
	}

	/**
	 * Creates a result of a <code>boolean</code>.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            a comment on the result
	 * @param result
	 *            the result itself
	 */
	public CategoricalResult(String name, String comment, boolean result) {
		super(name, comment, DataType.BOOLEAN);
		this.result = result;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "categoricalResult";
	}
}
