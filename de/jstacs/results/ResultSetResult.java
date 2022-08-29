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

package de.jstacs.results;

import de.jstacs.io.NonParsableException;

/**
 * Class for a {@link Result} that contains a single {@link ResultSet}. Useful for building hierarchies of {@link Result}s and
 * {@link ResultSet}s.
 * 
 * @author Jan Grau
 *
 */
public class ResultSetResult extends ListResult {

	/**
	 * Creates a new {@link ResultSetResult} with given name, comment, annotation, and content.
	 * @param name the name of the result
	 * @param comment a comment on the result
	 * @param annotation annotation (may be <code>null</code>)
	 * @param result the result
	 */
	public ResultSetResult( String name, String comment, ResultSet annotation, ResultSet result ) {
		super( name, comment, annotation, result );
	}

	/**
	 * Creates a {@link ResultSetResult} from its XML representation.
	 * @param representation the XML representation
	 * @throws NonParsableException if XML could not be parsed
	 */
	public ResultSetResult( StringBuffer representation ) throws NonParsableException {
		super( representation );
	}
	
	public ResultSet getResultSet() {
		return list[0];
	}

}
