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

package de.jstacs;

/**
 * This is the root interface for all immutable objects that must be stored in
 * e.g. a file or a database. Classes that implement this interface <b>must</b>
 * provide a constructor with a single parameter of type {@link StringBuffer}.
 * 
 * <br>
 * <br>
 * 
 * The recommended way to store the objects is an XML representation, because it
 * is human-readable and flexible enough to describe all objects properly. For
 * convenience {@link de.jstacs.io.XMLParser} provides methods to store all
 * primitive data types, <code>String</code>s, {@link Storable}s and their array
 * types.
 * 
 * <br>
 * <br>
 * 
 * For writing or reading a {@link StringBuffer} to or from a file you can use
 * {@link de.jstacs.io.FileManager}.
 * 
 * @see de.jstacs.io.XMLParser
 * @see de.jstacs.io.FileManager
 * 
 * @author Jan Grau, Jens Keilwagen
 */

public interface Storable {

	/**
	 * This method returns an XML representation as {@link StringBuffer} of an
	 * instance of the implementing class.
	 * 
	 * @return the XML representation
	 */
	public StringBuffer toXML();
}
