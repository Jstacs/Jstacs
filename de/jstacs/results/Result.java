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

import java.awt.image.BufferedImage;

import de.jstacs.AnnotatedEntity;
import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * The abstract class for any result. Each result should be immutable. In
 * analogy to the {@link de.jstacs.parameters.Parameter} classes, the {@link Result} classes provide
 * the possibility to return the results of a computation together with some
 * annotation in a standardized way.
 * 
 * @author Jan Grau
 */
public abstract class Result extends AnnotatedEntity {

	private String originalName;
	
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
	protected Result(String name, String comment, DataType datatype) {
		super( name, comment, datatype );
	}

	/**
	 * The standard constructor for the interface {@link Storable}. Creates a
	 * new {@link Result} out of its XML representation.
	 * 
	 * @param rep
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation is not parsable
	 * 
	 * @see Storable
	 */
	protected Result(StringBuffer rep) throws NonParsableException {
		super( rep );
	}

	/**
	 * Returns <code>true</code> if the {@link Result} <code>test</code> and the
	 * current object have the same data type, name and comment for the result.
	 * 
	 * <br>
	 * 
	 * <b>The method does NOT answer whether the method {@link #getValue()}
	 * returns an instance of {@link Comparable}.</b>
	 * 
	 * @param test
	 *            the {@link Result} to be tested
	 * 
	 * @return <code>true</code> if the {@link Result} <code>test</code> and the
	 *         current object have the same data type, name and comment for the
	 *         result
	 */
	public boolean isComparableResult(Result test) {
		return datatype == test.datatype && name.equals(test.name)
				&& comment.equals(test.comment);
	}

	/**
	 * Returns <code>true</code> if the data type of the {@link Result}
	 * <code>test</code> can be casted to that of this instance and both have
	 * the same name and comment for the {@link Result}.
	 * 
	 * <br>
	 * 
	 * <b>The method does NOT answer whether the method {@link #getValue()}
	 * returns an instance of {@link Comparable}.</b>
	 * 
	 * @param test
	 *            the {@link Result} to be tested
	 * 
	 * @return <code>true</code> if the {@link Result} <code>test</code> and the
	 *         current object have a castable data type, name and comment for
	 *         the result
	 */
	public boolean isCastableResult(Result test) {
		return DataType.canBeCastedFromTo(test.datatype, datatype)
				&& name.equals(test.name) && comment.equals(test.comment);
	}

	/**
	 * This method provides the possibility to check the compliance of some
	 * result <code>value</code> with one of the pre-defined {@link DataType}s
	 * before creating a new {@link Result} and possibly running into an
	 * {@link Exception}.
	 * 
	 * @param datatype
	 *            the data type the <code>value</code> should comply to
	 * @param value
	 *            the value
	 * 
	 * @return if <code>value</code> can be stored in a {@link Result} of
	 *         {@link DataType} <code>datatype</code>
	 */
	public static boolean checkDatatype(DataType datatype, Object value) {
		if (value instanceof String) {
			String val = (String) value;
			try {
				switch (datatype) {
				case DOUBLE:
					Double.parseDouble(val);
					break;
				case LONG:
					Long.parseLong(val);
					break;
				case INT:
					Integer.parseInt(val);
					break;
				case STRING:
					break;
				case BOOLEAN:
					Boolean.parseBoolean(val);
					break;
				default:
					return false;
				}
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		} else {
			switch (datatype) {
			case DOUBLE:
				return value instanceof Double || value instanceof Float;
			case LONG:
				return value instanceof Long;
			case INT:
				return value instanceof Integer || value instanceof Byte
						|| value instanceof Short;
			case BOOLEAN:
				return value instanceof Boolean;
			case PNG:
				return value instanceof BufferedImage;
			case STORABLE:
				return value instanceof Storable;
			case HTML:
			case LIST:
			default:
				return false;
			}
		}
	}

	/**
	 * Factory method to create a new {@link Result}.
	 * 
	 * @param name
	 *            the name of the {@link Result} as chosen by the user
	 * @param comment
	 *            the comment on the {@link Result} as chosen by the user
	 * @param datatype
	 *            the data type of <code>value</code>
	 * @param value
	 *            the value of the {@link Result}
	 * 
	 * @return the new {@link Result} instance
	 * 
	 * @throws IllegalValueException
	 *             if <code>datatype</code> is {@link DataType#HTML},
	 *             {@link DataType#LIST} or another {@link DataType} not
	 *             implemented, yet
	 */
	public static Result createResult(String name, String comment,
			DataType datatype, Object value) throws IllegalValueException {
		switch (datatype) {
		case DOUBLE:
		case LONG:
		case INT:
			return new NumericalResult(datatype, name, comment,
					(Comparable) value);
		case STRING:
		case BOOLEAN:
			return new CategoricalResult(datatype, name, comment,
					(Comparable) value);
		case PNG:
			return new ImageResult(name, comment, (BufferedImage) value);
		case STORABLE:
			return new StorableResult(name, comment, (Storable) value);
		case HTML:
		case LIST:
		default:
			throw new IllegalValueException(name,"wrong datatype");
		}
	}
	
	/**
	 * Renames this Result. The original name will be stored and is available from {@link #getOriginalName()}. 
	 * @param newName the new name
	 */
	public void rename(String newName){
		if(this.originalName == null){
			this.originalName = this.name;
		}
		this.name = newName;
	}
	
	/**
	 * Returns the original name (i.e., the name upon object creation) of this {@link Result}, 
	 * which may be just the name if {@link #rename(String)} has not been called on this object, yet.
	 * @return the original name
	 */
	public String getOriginalName(){
		if(this.originalName == null){
			return this.name;
		}else{
			return this.originalName;
		}
	}

	@Override
	protected void appendFurtherInfos(StringBuffer buf) {
		XMLParser.appendObjectWithTags(buf, originalName, "originalName");		
	}

	@Override
	protected void extractFurtherInfos(StringBuffer buf) throws NonParsableException {
		originalName = (String) XMLParser.extractObjectForTags(buf, "originalName");
	}
	
	

}
