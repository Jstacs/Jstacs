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

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.NonParsableException;
import de.jstacs.io.ParameterSetParser;
import de.jstacs.io.XMLParser;
import de.jstacs.io.ParameterSetParser.NotInstantiableException;

/**
 * Container class for a set of {@link Parameter}s that can be used to
 * instantiate another class.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 */
public abstract class InstanceParameterSet extends ParameterSet {
	/**
	 * The class that can be instantiated using this {@link ParameterSet}.
	 */
	private Class instanceClass;

	/**
	 * Constructs an {@link InstanceParameterSet} from the class that can be
	 * instantiated using this {@link InstanceParameterSet}. 
	 * 
	 * @param instanceClass
	 *            the class to be instantiated
	 */
	public InstanceParameterSet(Class instanceClass) {
		super();
		if (instanceClass == null) {
			throw new IllegalArgumentException(	"The instanceClass can not be null." );
		}
		this.instanceClass = instanceClass;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link InstanceParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link InstanceParameterSet} could not be
	 *             reconstructed out of the {@link StringBuffer}
	 *             <code>representation</code>
	 */
	public InstanceParameterSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/**
	 * Returns the class of the instances that can be constructed using this
	 * set.
	 * 
	 * @return the class of the instances
	 */
	public Class getInstanceClass() {
		return this.instanceClass;
	}

	/**
	 * Returns a new instance of the class of {@link #getInstanceClass()} that
	 * was created using this {@link ParameterSet}.
	 * 
	 * @return the new instance of the class
	 * 
	 * @throws NotInstantiableException
	 *             if the class could not be instantiated
	 * 
	 * @see de.jstacs.io.ParameterSetParser
	 */
	public InstantiableFromParameterSet getInstance()
			throws NotInstantiableException {
		return ParameterSetParser.getInstanceFromParameterSet(this);
	}

	/**
	 * Returns a comment (a textual description) of the class that can be
	 * constructed using this {@link ParameterSet}.
	 * 
	 * @return the comment of the class
	 */
	public abstract String getInstanceComment();

	/**
	 * Returns the name of an instance of the class that can be constructed
	 * using this {@link ParameterSet}.
	 * 
	 * @return the name of the class
	 */
	public abstract String getInstanceName();

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superParameterSet");
		XMLParser.appendObjectWithTags(buf, instanceClass,"instanceClass");
		XMLParser.addTags(buf, "instanceParameterSet");
		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "instanceParameterSet" );
		super.fromXML( XMLParser.extractForTag( representation, "superParameterSet" ) );
		instanceClass = XMLParser.extractObjectForTags( representation, "instanceClass", Class.class );
	}
}
