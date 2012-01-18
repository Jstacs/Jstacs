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

package de.jstacs.parameters;

import de.jstacs.AnnotatedEntity;
import de.jstacs.DataType;
import de.jstacs.io.NonParsableException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * Abstract class for a parameter that shall be used as the parameter of some
 * method, constructor, etc. This class provides annotation of the parameter,
 * e.g. a name and a comment, the data type, and additional properties for
 * non-atomic parameters. The intention of the {@link Parameter} classes is to
 * link this annotation directly to the parameter values, which could not be
 * achieved by primitive data types.<br/> {@link Parameter}s should always be
 * embedded into a {@link ParameterSet} and then passed to the respective
 * method.
 * 
 * @see de.jstacs.parameters.ParameterSet
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 */
public abstract class Parameter extends AnnotatedEntity implements Cloneable {

	/**
	 * If this {@link Parameter} is enclosed in a {@link ParameterSet}, this
	 * variable holds a reference to that {@link ParameterSet}.
	 */
	protected ParameterSet parent;

	/**
	 * The main constructor which takes the main information of a {@link Parameter}.
	 * 
	 * @param name
	 *            the name of the result
	 * @param comment
	 *            the comment for the result
	 * @param type
	 *            the data type of the result
	 */
	public Parameter( String name, String comment, DataType type ) {
		super( name, comment, type );
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}. Creates a
	 * new {@link Parameter} out of its XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML representation is not parsable
	 * 
	 * @see de.jstacs.Storable
	 * @see #extractFurtherInfos(StringBuffer)
	 */
	public Parameter( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}


	/**
	 * Returns <code>true</code> if the {@link Parameter} is required,
	 * <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if the {@link Parameter} is required,
	 *         <code>false</code> otherwise
	 */
	public abstract boolean isRequired();

	/**
	 * Checks the value for correctness, e.g. for numerical parameters this
	 * might be checking if the value is within specified bounds.
	 * 
	 * @param value
	 *            the value to be checked
	 * 
	 * @return <code>true</code> if the value is valid, <code>false</code>
	 *         otherwise
	 */
	public abstract boolean checkValue(Object value);

	/**
	 * Sets the value of this {@link Parameter} to <code>value</code>.
	 * 
	 * @param value
	 *            the new value of the {@link Parameter}
	 * 
	 * @throws IllegalValueException
	 *             if the specified value is not valid for this
	 *             {@link Parameter}
	 */
	public abstract void setValue(Object value) throws IllegalValueException;

	/**
	 * Returns <code>true</code> if the parameter either has a default value or
	 * the value was set by the user, <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if value has a default value or was set,
	 *         <code>false</code> otherwise
	 */
	public abstract boolean hasDefaultOrIsSet();

	/**
	 * Returns <code>true</code> if the parameter was set by the user,
	 * <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if the parameter was set, <code>false</code>
	 *         otherwise
	 */
	public abstract boolean isSet();

	/**
	 * Returns <code>true</code> if the parameter is of an atomic data type,
	 * <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if the parameter is atomic, <code>false</code>
	 *         otherwise
	 */
	public abstract boolean isAtomic();

	/**
	 * If a value could not be set successfully this method returns the
	 * corresponding error message.
	 * 
	 * @return an error message if a value could not be set successfully
	 */
	public abstract String getErrorMessage();

	/**
	 * Resets the parameter and its contents to the default values or
	 * <code>null</code> if no defaults are given.
	 */
	public abstract void reset();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public Parameter clone() throws CloneNotSupportedException {
		Parameter clone = (Parameter) super.clone();
		clone.parent = null;
		return clone;
	}

	/**
	 * Sets the default value of the {@link Parameter} to
	 * <code>defaultValue</code>. This method also sets the current
	 * value of this {@link Parameter} to the default.
	 * 
	 * @param defaultValue
	 *            the default value
	 * 
	 * @throws Exception
	 *             if the default value is not an appropriate value for this
	 *             {@link Parameter}
	 */
	public abstract void setDefault(Object defaultValue) throws Exception;

	/**
	 * Sets the reference of the enclosing {@link ParameterSet} of this
	 * {@link Parameter} to <code>parent</code>.
	 * 
	 * @param parent
	 *            the new enclosing {@link ParameterSet}
	 */
	public void setParent(ParameterSet parent) {
		this.parent = parent;
	}

	/**
	 * Returns a reference to the {@link ParameterSet} enclosing this
	 * {@link Parameter}.
	 * 
	 * @return the reference to the enclosing {@link ParameterSet}
	 */
	public ParameterSet getParent() {
		return parent;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos(StringBuffer source) throws NonParsableException {
		
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		
	}
	
	/**
	 * This method checks whether the given {@link Parameter} is comparable to the current instance, i.e. whether
	 * the {@link Class}, the {@link DataType}, the name and the comment are identical. If necessary, all these
	 * characteristics are checked recursively.  
	 * 
	 * In other words, the method returns <code>true</code> if the parameters only differ in their specific raw values.
	 * 
	 * @param p the {@link Parameter} for the comparison
	 * 
	 * @return <code>true</code> if the parameters only differ in their values, otherwise <code>false</code>
	 * 
	 * @see Object#getClass()
	 * @see #getDatatype()
	 * @see #getName()
	 * @see #getComment()
	 * @see DataType#PARAMETERSET
	 * @see ParameterSet#isComparable(ParameterSet)
	 */
	public boolean isComparable( Parameter p ) {
		boolean res = getClass().equals(p.getClass()) && getDatatype() == p.getDatatype() && getName().equals(p.getName()) && getComment().equals( p.getComment() );
		if( res && getDatatype() == DataType.PARAMETERSET ) {
			return ((ParameterSet)getValue()).isComparable( (ParameterSet) p.getValue() );
		} else {
			return res;
		}
	}
}