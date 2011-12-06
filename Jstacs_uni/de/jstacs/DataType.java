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
 * This <code>enum</code> defines a number of data types that can be used for
 * {@link de.jstacs.parameters.Parameter}s and {@link de.jstacs.results.Result}
 * s.
 * 
 * @author Jens Keilwagen, Jan Grau
 * 
 * @see de.jstacs.parameters.Parameter
 * @see de.jstacs.results.Result
 */
public enum DataType {

	/**
	 * This value indicates the data type <code>boolean</code>.
	 */
	BOOLEAN,

	/**
	 * This value indicates the data type <code>char</code>.
	 */
	CHAR,

	/**
	 * This value indicates the data type <code>byte</code>.
	 */
	BYTE,

	/**
	 * This value indicates the data type <code>short</code>.
	 */
	SHORT,

	/**
	 * This value indicates the data type <code>int</code>.
	 */
	INT,

	/**
	 * This value indicates the data type <code>long</code>.
	 */
	LONG,

	/**
	 * This value indicates the data type <code>float</code>.
	 */
	FLOAT,

	/**
	 * This value indicates the data type <code>double</code>.
	 */
	DOUBLE,

	/**
	 * This value indicates the data type {@link String}.
	 */
	STRING,

	/**
	 * This value indicates the data type HTML.
	 */
	HTML,

	/**
	 * This value indicates the data type png.
	 */
	PNG,

	/**
	 * This value indicates the data type {@link Storable}.
	 */
	STORABLE,

	/**
	 * This value indicates the data type {@link de.jstacs.data.DataSet}.
	 */
	DATASET,

	/**
	 * This value indicates the data type {@link de.jstacs.results.ListResult}.
	 */
	LIST,

	/**
	 * This value indicates the data type
	 * {@link de.jstacs.parameters.ParameterSet}.
	 */
	PARAMETERSET,

	/**
	 * This value indicates the data type
	 * {@link de.jstacs.parameters.FileParameter.FileRepresentation}.
	 */
	FILE;

	/**
	 * Checks if the {@link DataType} <code>from</code> can be casted to the
	 * {@link DataType} <code>to</code> without losing information. The
	 * following casts are allowed:
	 * <ul>
	 * <li> {@link DataType#BYTE} -> {@link DataType#BYTE},
	 * {@link DataType#SHORT}, {@link DataType#INT}, {@link DataType#LONG},
	 * {@link DataType#DOUBLE}</li>
	 * <li> {@link DataType#SHORT} -> {@link DataType#SHORT},
	 * {@link DataType#INT}, {@link DataType#LONG}, {@link DataType#DOUBLE}</li>
	 * <li> {@link DataType#INT} -> {@link DataType#INT}, {@link DataType#LONG},
	 * {@link DataType#DOUBLE}</li>
	 * <li> {@link DataType#LONG} -> {@link DataType#LONG}</li>
	 * <li> {@link DataType#DOUBLE} -> {@link DataType#DOUBLE}</li>
	 * </ul>
	 * 
	 * @param from
	 *            the {@link DataType} to cast from
	 * @param to
	 *            the {@link DataType} to cast to
	 * 
	 * @return <code>true</code> if the cast is possible, <code>false</code>
	 *         otherwise
	 */
	public static boolean canBeCastedFromTo( DataType from, DataType to ) {
		switch( from ) {
			case BYTE:
				if( to == BYTE ) {
					return true;
				}
			case SHORT:
				if( to == SHORT ) {
					return true;
				}
			case INT:
				if( to == INT ) {
					return true;
				}
			case DOUBLE:
				if( to == DOUBLE ) {
					return true;
				} else if( from == DOUBLE && to == LONG ) {
					return false;
				}
			case LONG:
				if( to == LONG ) {
					return true;
				} else if( from == LONG && to == DOUBLE ) {
					return false;
				}
			default:
				return false;
		}
	}
}
