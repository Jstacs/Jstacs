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
package de.jstacs.algorithms.optimization.termination;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.InstanceParameterSet;

/**
 * This class is the abstract super class of many {@link TerminationCondition}s. It implements several methods for
 * loading, saving, cloning, ... allowing for a fast implementation of a new {@link TerminationCondition}.
 * 
 * @author Jens Keilwagen
 */
public abstract class AbstractTerminationCondition implements TerminationCondition {

	/**
	 * The internally used parameters.
	 */
	protected AbstractTerminationConditionParameterSet parameter;
	
	/**
	 * This is the main constructor creating an instance from a given parameter set. 
	 * 
	 * @param parameter the set of parameters
	 * 
	 * @throws CloneNotSupportedException if <code>parameter</code> can not be cloned properly.
	 */
	protected AbstractTerminationCondition( AbstractTerminationConditionParameterSet parameter ) throws CloneNotSupportedException {
		this.parameter = parameter.clone();
		set();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link AbstractTerminationCondition} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link AbstractTerminationCondition} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	@SuppressWarnings("unchecked")
	protected AbstractTerminationCondition( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getXmlTag() );
		parameter = (AbstractTerminationConditionParameterSet) XMLParser.extractObjectForTags( xml, "parameter" );
		if( !parameter.getInstanceClass().isAssignableFrom(this.getClass()) ) {
			throw new NonParsableException( "The loaded ParameterSet can not be used for this class" );
		}
		set();
		
	}
	
	public AbstractTerminationCondition clone() throws CloneNotSupportedException {
		AbstractTerminationCondition clone = (AbstractTerminationCondition) super.clone();
		clone.parameter = parameter.clone();
		clone.set();
		return clone;
	}
	
	/**
	 * This method returns the xml tag that is used in the method {@link AbstractTerminationCondition#toXML()} and
	 * in the constructor {@link AbstractTerminationCondition#AbstractTerminationCondition(StringBuffer)}. 
	 * 
	 * @return the xml tag of the instance
	 */
	protected abstract String getXmlTag();

	/**
	 * This method sets internal member variables from {@link #parameter}.
	 * It is used in the constructors.
	 */
	protected abstract void set();
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public final StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, parameter, "parameter" );
		XMLParser.addTags( xml, getXmlTag() );
		return xml;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.InstantiableFromParameterSet#getCurrentParameterSet()
	 */
	public AbstractTerminationConditionParameterSet getCurrentParameterSet() throws Exception {
		return (AbstractTerminationConditionParameterSet) parameter.clone();
	}
	
	/**
	 * This class implements the super class of all parameter sets of instances from {@link AbstractTerminationCondition}.
	 * 
	 * @author Jens Keilwagen
	 */
	public static abstract class AbstractTerminationConditionParameterSet extends InstanceParameterSet<AbstractTerminationCondition> {
		
		/**
		 * Constructs an {@link AbstractTerminationConditionParameterSet} from the class that can be
		 * instantiated using this {@link AbstractTerminationConditionParameterSet}.
		 * 
		 * @param instanceClass
		 *            the class to be instantiated
		 * 
		 * @throws IllegalArgumentException
		 *             if <code>instanceClass</code> is null
		 */
		public AbstractTerminationConditionParameterSet(Class<? extends AbstractTerminationCondition> instanceClass)
				throws IllegalArgumentException {
			super(instanceClass);
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}.
		 * Constructs an {@link AbstractTerminationConditionParameterSet} out of an XML representation.
		 * 
		 * @param representation
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link AbstractTerminationConditionParameterSet} could not be
		 *             reconstructed out of the {@link StringBuffer}
		 *             <code>representation</code>
		 */
		public AbstractTerminationConditionParameterSet(StringBuffer representation)
				throws NonParsableException {
			super(representation);
		}
		
		public AbstractTerminationConditionParameterSet clone() throws CloneNotSupportedException {
			return (AbstractTerminationConditionParameterSet) super.clone();
		}
	}
}
