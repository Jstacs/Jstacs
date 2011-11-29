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

import java.util.LinkedList;

import de.jstacs.NonParsableException;

/**
 * Class for a set of numerical result values, which are all of the type
 * {@link NumericalResult}.
 * 
 * @author Jan Grau
 */
public class NumericalResultSet extends ResultSet {

	/**
	 * Constructs a {@link NumericalResultSet} containing one
	 * {@link NumericalResult}.
	 * 
	 * @param result
	 *            the {@link NumericalResult} to be contained
	 */
	public NumericalResultSet(NumericalResult result) {
		super(result);
	}

	/**
	 * Constructs a {@link NumericalResultSet} from some arrays of
	 * {@link NumericalResult}s.
	 * 
	 * @param results
	 *            the results
	 */
	public NumericalResultSet(NumericalResult[]... results) {
		super(results);
	}

	/**
	 * Constructs a {@link NumericalResultSet} from a {@link LinkedList} of
	 * {@link NumericalResult}s.
	 * 
	 * @param results
	 *            the {@link LinkedList} of {@link NumericalResult}s
	 */
	public NumericalResultSet(LinkedList<? extends NumericalResult> results) {
		super(results);
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link NumericalResultSet} from the corresponding XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed.
	 */
	public NumericalResultSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.ResultSet#getResultAt(int)
	 */
	@Override
	public NumericalResult getResultAt(int index) {
		return (NumericalResult) results[index];
	}
}
