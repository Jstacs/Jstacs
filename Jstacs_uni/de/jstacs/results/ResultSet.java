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

import java.util.Collection;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.XMLParser;

/**
 * Class for a set of {@link Result}s which provides methods to access single
 * {@link Result}s in the set, to retrieve the number of {@link Result}s in the
 * set, to get a {@link String} representation or an XML representation of all
 * the {@link Result}s in the set.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class ResultSet implements Storable {
	/**
	 * The set of results as an array of {@link Result}s.
	 */
	protected Result[] results;

	/**
	 * Constructs a new {@link ResultSet} containing one {@link Result}.
	 * 
	 * @param result
	 *            the {@link Result} to be contained
	 */
	public ResultSet(Result result) {
		this.results = new Result[] { result };
	}

	/**
	 * Constructs a new {@link ResultSet} from a two-dimensional array of
	 * {@link Result}s.
	 * 
	 * @param results
	 *            the two-dimensional array of {@link Result}s
	 */
	public ResultSet(Result[]... results) {
		if (results == null) {
			this.results = new Result[0];
		} else {
			int c = 0, i;
			for( i = 0; i < results.length; c += (results[i] == null) ? 0 : results[i].length, i++ );

			this.results = new Result[c];
			c = 0;
			for (i = 0; i < results.length; i++) {
				if (results[i] != null) {
					System.arraycopy(results[i], 0, this.results, c, results[i].length);
					c += results[i].length;
				}
			}
		}
	}

	/**
	 * Constructs a new {@link ResultSet} from a {@link Collection} of type
	 * {@link Result}.
	 * 
	 * @param results
	 *            a {@link Collection} of {@link Result}s
	 */
	public ResultSet(Collection<? extends Result> results) {
		this.results = results.toArray(new Result[0]);
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ResultSet} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public ResultSet(StringBuffer representation) throws NonParsableException {
		fromXML(representation);
	}

	/**
	 * Returns {@link Result} number <code>index</code> in this
	 * {@link ResultSet}.
	 * 
	 * @param index
	 *            the index of the {@link Result}
	 * 
	 * @return the {@link Result} at <code>index</code>
	 */
	public Result getResultAt(int index) {
		return results[index];
	}

	/**
	 * Returns all internal results as an array of {@link Result}s.
	 * 
	 * @return all internal results as an array of {@link Result}s
	 */
	public Result[] getResults() {
		Result[] res = new Result[results.length];
		System.arraycopy(results, 0, res, 0, res.length);
		return res;
	}

	/**
	 * Returns the number of {@link Result}s in this {@link ResultSet}
	 * 
	 * @return the number of {@link Result}s
	 */
	public int getNumberOfResults() {
		return results.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags( buf, results, "resStrings" );
		XMLParser.addTags(buf, "resultSet");
		return buf;
	}

	/**
	 * Parses the contents of a {@link ResultSet} from its XML representation as
	 * returned by {@link #toXML()}.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	protected void fromXML(StringBuffer representation) throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "resultSet" );
		results = XMLParser.extractObjectForTags( representation, "resStrings", Result[].class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuffer help = new StringBuffer(100 * results.length);
		for (int i = 0; i < results.length; i++) {
			help.append(results[i].toString() + "\n");
		}
		return help.toString();
	}

	/**
	 * This method enables you to search for a column. It returns the index of
	 * the column, if it could be found, otherwise -1.
	 * 
	 * @param columnName
	 *            the name of the column
	 * 
	 * @return the index of the column, if it could be found, otherwise -1.
	 */
	public int findColumn(String columnName) {
		for (int i = 0; i < results.length; i++) {
			if (results[i].getName().equals(columnName)) {
				return i;
			}
		}
		return -1;
	}
}
