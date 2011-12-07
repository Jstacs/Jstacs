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
import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;
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
 * @author Jan Grau
 * 
 */
public abstract class Parameter extends AnnotatedEntity implements Cloneable {

	/**
	 * In cases, when the validity of some {@link ParameterSet} depends on the
	 * value of this {@link Parameter}, this variable holds a reference to that
	 * {@link ParameterSet}.
	 */
	protected ParameterSet neededReference;

	/**
	 * In addition to the reference, the id of the {@link ParameterSet} is
	 * saved, which assists reconstruction from XML.
	 * 
	 * @see Parameter#neededReference
	 */
	protected Long neededReferenceId;

	/**
	 * If this {@link Parameter} is enclosed in a {@link ParameterSet}, this
	 * variable holds a reference to that {@link ParameterSet}.
	 */
	protected ParameterSet parent;

	private long id;

	/**
	 * Creates a new {@link Parameter} and generates the internal id.
	 */
	public Parameter() {
		id = System.currentTimeMillis() + this.hashCode();
	}

	/**
	 * Returns the id of this {@link Parameter}.
	 * 
	 * @return the id of this {@link Parameter}
	 */
	public long getId() {
		return id;
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
	 * Returns a reference to a {@link ParameterSet} whose
	 * {@link ParameterSet#hasDefaultOrIsSet()}-method depends on the value of
	 * this {@link Parameter}. If no such {@link ParameterSet} exists,
	 * <code>null</code> is returned.
	 * 
	 * @return the reference to the {@link ParameterSet}
	 * 
	 * @see ParameterSet#hasDefaultOrIsSet()
	 */
	public ParameterSet getNeededReference() {
		return neededReference;
	}

	/**
	 * Returns the id of the {@link ParameterSet} that would be returned by
	 * {@link Parameter#getNeededReference()}.
	 * 
	 * @return the id of the {@link ParameterSet}
	 * 
	 * @see Parameter#getNeededReference()
	 */
	public Long getNeededReferenceId() {
		return neededReferenceId;
	}

	/**
	 * Sets an internal reference to a {@link ParameterSet} whose validity
	 * depends on the value of this {@link Parameter}.
	 * 
	 * @param reference
	 *            to the {@link ParameterSet} depending on the value of this
	 *            {@link Parameter}
	 */
	public void setNeededReference(ParameterSet reference) {
		this.neededReference = reference;
		this.neededReferenceId = reference.getId();
	}

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
	 * Simplifies the {@link Parameter} and its contents to the relevant
	 * information. This could be e.g. to reset the contents of those values of
	 * a {@link CollectionParameter} that are not selected.
	 */
	public abstract void simplify();

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
		clone.neededReference = null;
		clone.parent = null;
		return clone;
	}

	/**
	 * Sets the default value of the {@link Parameter} to
	 * <code>defaultValue</code>.
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

	/**
	 * Parses a {@link Parameter} from a XML representation as returned by
	 * {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 * 
	 * @see Parameter#toXML()
	 */
	protected void fromXML(StringBuffer source) throws NonParsableException {
		source = XMLParser.extractForTag(source, "parameter");
		this.id = XMLParser.extractObjectForTags(source, "id", long.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		try {
			this.neededReferenceId = XMLParser.extractObjectForTags(source, "neededReferenceId", long.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		} catch (NonParsableException e) {
			this.neededReferenceId = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, id, "id");
		XMLParser.addTags(buf, "parameter");
		if (neededReference != null) {
			XMLParser.appendObjectWithTags(buf, neededReferenceId,
					"neededReferenceId");
		}
		return buf;
	}

}