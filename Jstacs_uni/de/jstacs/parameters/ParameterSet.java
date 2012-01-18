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

import java.util.ArrayList;
import java.util.Iterator;

import de.jstacs.AnnotatedEntityList;
import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * (Container) class for a set of {@link Parameter}s. This is the base class for
 * other {@link ParameterSet}s that provide more specialized methods, e.g. for
 * expanding the set of {@link Parameter}s at runtime (
 * {@link ExpandableParameterSet} ) or defining an array of {@link ParameterSet}
 * s from a common template ( {@link ArrayParameterSet} ).
 * 
 * @author Jan Grau
 */
public abstract class ParameterSet implements Storable, Cloneable, GalaxyConvertible {
	
	//XXX
	public static String getName( Class<? extends ParameterSet> c ) {
		return c.getSimpleName();
	}
	
	public static String getComment( Class<? extends ParameterSet> c ) {
		return "";
	}
	
	public static String getName( ParameterSet p ) {
		return getName( p.getClass() );
	}
	
	public static String getComment( ParameterSet p ) {
		return getComment( p.getClass() );
	}
	
	/**
	 * This method tries to find the correct name ({@link String}) for your
	 * choice. This method is useful if you handle {@link AbstractSelectionParameter}s.
	 * 
	 * @param names
	 *            the names
	 * @param values
	 *            the values that can be set
	 * @param current
	 *            the value to be set
	 * @param hasAlternative
	 *            indicates whether the last entry of names is an alternative
	 *            parameter
	 * 
	 * @return the {@link StringIndexOutOfBoundsException} that has to be used
	 * 
	 * @throws IllegalArgumentException
	 *             if no match could be found
	 * 
	 * @see AbstractSelectionParameter
	 */
	public static int getIndex(String[] names, Object[] values,
			Comparable current, boolean hasAlternative)
			throws IllegalArgumentException {
		int i = 0, l = values.length;
		if (hasAlternative) {
			l--;
		}
		while (i < l && !current.equals(values[i])) {
			i++;
		}
		if (i == values.length) {
			throw new IllegalArgumentException(
					"Could not find a matching constant.");
		} else {
			return i;
		}
	}
	
	/**
	 * The set of parameters
	 */
	protected ParameterList parameters;

	/**
	 * The error message of the last error or <code>null</code>
	 */
	protected String errorMessage;

	/**
	 * If this {@link ParameterSet} is contained in a
	 * {@link ParameterSetContainer}, this variable holds a reference to that
	 * {@link ParameterSetContainer}.
	 */
	protected ParameterSetContainer parent;

	/**
	 * Constructs a new {@link ParameterSet} with empty parameter values.
	 */
	public ParameterSet() {
		initParameterList();
	}

	/**
	 * Creates a full clone (deep copy) of this {@link ParameterSet}. As a
	 * convenience-method the user can use
	 * <code>fillWithStandardFieldsForClone(ParameterSet)</code> on a newly
	 * created instance of a subclass of {@link ParameterSet} to obtain a
	 * clone/copy of all standard member variables (those already defined in
	 * {@link ParameterSet}) in the passed {@link ParameterSet}. Using this
	 * method, the cloning process becomes merely three-step:<br>
	 * <ul>
	 * <li>Create a new instance of your subclass of {@link ParameterSet}, most
	 * likely with an empty constructor or the one taking just the instance
	 * class.
	 * <li>Call <code>this.fillWithStandardFieldsForClone</code> on this
	 * instance.
	 * <li>Return the instance.
	 * </ul>
	 * This method fulfills the conventions of {@link Object}'s method
	 * {@link Object#clone()}.
	 * 
	 * @return a deep clone/copy of this {@link ParameterSet}
	 */
	@Override
	public ParameterSet clone() throws CloneNotSupportedException {

		ParameterSet ret = (ParameterSet) super.clone();

		if (this.parameters != null) {
			ret.initParameterList(this.parameters.size());
			for (int i = 0; i < this.parameters.size(); i++) {
				ret.parameters.add(i, this.parameters.get(i).clone());
			}
		}
		return ret;
	}

	/**
	 * Initializes the internal set of {@link Parameter}s, which is a
	 * {@link ParameterList}.
	 */
	protected final void initParameterList() {
		this.parameters = new ParameterList();
	}

	/**
	 * Initializes the internal set of {@link Parameter}s, which is a
	 * {@link ParameterList}, with an initial number of {@link Parameter}s of
	 * <code>initCapacity</code>.
	 * 
	 * @param initCapacity
	 *            the initial number of {@link Parameter}s
	 */
	protected final void initParameterList(int initCapacity) {
		this.parameters = new ParameterList(initCapacity);
	}
	
	/**
	 * Returns <code>true</code> if the parameters of this {@link ParameterSet}
	 * have been loaded.
	 * 
	 * @return the state of the parameters (loaded or not)
	 */
	//TODO remove?
	public boolean parametersLoaded() {
		return (parameters != null);
	}
	
	/**
	 * Returns the index of the first parameter in the set of required
	 * parameters that has not been set. If all required parameters have been
	 * set, it returns -1.
	 * 
	 * @return the index of the required parameter that has not been set
	 */
	/*
	 * public int getIndexOfRequiredParameterNotSet() { if( requiredParameters
	 * == null ) { return -1; } else { for( int i = 0; i <
	 * requiredParameters.length; i++ ) { if( !requiredParameters[i].isSet() ) {
	 * return i; } } return -1; } }
	 */

	/**
	 * Returns <code>true</code> if this {@link ParameterSet} contains only
	 * atomic parameters, i.e. the parameters do not contain
	 * {@link ParameterSet}s themselves.
	 * 
	 * @return <code>true</code> if the {@link ParameterSet} contains only
	 *         atomic parameters, <code>false</code> otherwise
	 */
	public boolean isAtomic() {
		for (int i = 0; i < parameters.size(); i++) {
			if (!parameters.get(i).isAtomic()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns <code>true</code> if all parameters in this {@link ParameterSet}
	 * are either set by the user or have default values. If any additional
	 * constraints are required on your parameters you should either specify
	 * some {@link de.jstacs.parameters.validation.ParameterValidator} on these parameters or implement these
	 * constraints by overriding this method in your implementation of
	 * {@link ParameterSet}. It is recommended to specify a useful remark which
	 * constraint failed in the member-variable <code>errorMessage</code>, which
	 * will be displayed to the user. In the overriding method
	 * <code>super.</code>{@link #hasDefaultOrIsSet()} should be called prior to
	 * checking specific constraints.
	 * 
	 * @return <code>true</code> if all parameters have some allowed value set,
	 *         <code>false</code> otherwise
	 */
	public boolean hasDefaultOrIsSet() {
		for (int i = 0; i < parameters.size(); i++) {
			if (parameters.get(i).isRequired()
					&& (!parameters.get(i).hasDefaultOrIsSet())) {
				return false;
			}
		}
		errorMessage = null;
		return true;
	}

	/**
	 * Returns the message of the last error that occurred. If no error occurred
	 * this method returns <code>null</code>.
	 * 
	 * @return the error message of the last error or <code>null</code>
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/**
	 * Constructs a {@link ParameterSet} out of an array of {@link Parameter}s.
	 * The {@link Parameter}s are not cloned, but passed by reference.
	 * 
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	protected ParameterSet(Parameter[] parameters) {
		initParameterList(parameters.length);
		for (int i = 0; i < parameters.length; i++) {
			this.parameters.add(parameters[i]);
		}
	}

	/**
	 * Constructs a {@link ParameterSet} out of an {@link ArrayList} of
	 * {@link Parameter}s. The {@link Parameter}s are not cloned, but passed by
	 * reference.
	 * 
	 * @param parameters
	 *            the {@link Parameter}s
	 */
	protected ParameterSet(ArrayList<Parameter> parameters) {
		initParameterList( parameters.size() );
		Iterator<Parameter> parIt = parameters.iterator();
		while (parIt.hasNext()) {
			this.parameters.add(parIt.next());
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link ParameterSet} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>representation</code>
	 */
	public ParameterSet(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/**
	 * Returns the number of parameters in the {@link ParameterSet}.
	 * 
	 * @return the number of parameters in the {@link ParameterSet}
	 */
	public int getNumberOfParameters() {
		if (parameters != null) {
			return parameters.size();
		} else {
			return 0;
		}
	}

	/**
	 * Returns the names of all {@link Parameter}s in this {@link ParameterSet}.
	 * 
	 * @return an array containing the names of all {@link Parameter}s in this {@link ParameterSet}
	 */
	public String[] getAllParameterNames(){
		return parameters.getNames();
	}
	
	/**
	 * Returns the {@link Parameter} with name <code>name</code>.
	 * 
	 * @param name
	 *            the name of the {@link Parameter}
	 * 
	 * @return the {@link Parameter} with name <code>name</code>
	 */
	public Parameter getParameterForName(String name){
		return parameters.get(name);
	}
	
	/**
	 * Returns the {@link Parameter} at position <code>i</code>.
	 * 
	 * @param i
	 *            the position in the {@link ParameterSet}
	 * 
	 * @return the {@link Parameter} at position <code>i</code>
	 */
	public Parameter getParameterAt(int i) {
		return parameters.get(i);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		if (parameters != null) {
			StringBuffer buf2 = new StringBuffer();
			XMLParser.appendObjectWithTags(buf2, parameters.size(), "numberOfParameters");
			for (int i = 0; i < parameters.size(); i++) {
				XMLParser.appendObjectWithTags(buf2, parameters.get( i ), "parameter");
				/*
				if (parameters.get(i) == null) {
					XMLParser.appendObjectWithTags(buf2, "null", "parameter");
				} else {
					StringBuffer buf3 = new StringBuffer();
					XMLParser.appendObjectWithTags(buf3, parameters.get(i).getClass().getName(), "className");
					buf3.append(parameters.get(i).toXML());
					XMLParser.addTags(buf3, "parameter");
					buf2.append(buf3);
				}
				*/
			}
			XMLParser.addTags(buf2, "set");
			buf.append(buf2);
		} else {
			XMLParser.appendObjectWithTags(buf, null, "set");
		}
		XMLParser.addTags(buf, "parameterSet");

		return buf;
	}

	/**
	 * Resets all {@link Parameter}s in this {@link ParameterSet} to their
	 * default values or <code>null</code> if not default value was provided.
	 * 
	 * @see Parameter#reset()
	 */
	public void reset() {
		parameters = null;
	}

	/**
	 * Parses the instance fields of a {@link ParameterSet} from the XML
	 * representation as returned by {@link #toXML()}.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 * 
	 * @see ParameterSet#toXML()
	 */
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation, "parameterSet");	
		StringBuffer buf = XMLParser.extractForTag(representation, "set");
		if (!buf.toString().equalsIgnoreCase("null")) {
			int numPars = XMLParser.extractObjectForTags(buf, "numberOfParameters", int.class );
			parameters = new ParameterList(numPars);
			for (int i = 0; i < numPars; i++) {
				parameters.add( (Parameter) XMLParser.extractObjectForTags( buf, "parameter" ) );
			}
		}
	}

	/**
	 * Returns the enclosing {@link ParameterSetContainer} of this
	 * {@link ParameterSet} or <code>null</code> if none exists.
	 * 
	 * @return the enclosing {@link ParameterSetContainer} (<code>parent</code>)
	 */
	public ParameterSetContainer getParent() {
		return parent;
	}

	/**
	 * Sets the enclosing {@link ParameterSetContainer} of this
	 * {@link ParameterSet} to <code>parent</code>.
	 * 
	 * @param parent
	 *            the new parent
	 */
	public void setParent(ParameterSetContainer parent) {
		this.parent = parent;
	}

	
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		for(int i=0;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).toGalaxy( namePrefix+"_ps", configPrefix, depth+1, descBuffer, configBuffer );
			descBuffer.append( "\n" );
			configBuffer.append( "\n" );
		}	
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		for(int i=0;i<getNumberOfParameters();i++){
			((GalaxyConvertible)getParameterAt( i )).fromGalaxy( namePrefix+"_ps", command );
		}
	}
	
	/**
	 * This method checks whether the given {@link ParameterSet} is comparable to the current instance, i.e. whether
	 * the {@link Class} and the number of parameters are identical, and the individual {@link Parameter}s are comparable.
	 * 
	 * In other words, the method returns <code>true</code> if the {@link ParameterSet}s only differ in their specific raw values.
	 * 
	 * @param p the {@link ParameterSet} for the comparison
	 * 
	 * @return <code>true</code> if the {@link ParameterSet}s only differ in their values, otherwise <code>false</code>
	 * 
	 * @see Object#getClass()
	 * @see #getNumberOfParameters()
	 * @see Parameter#isComparable(Parameter)
	 */
	public boolean isComparable( ParameterSet p ){
		boolean res = getClass().equals(p.getClass());
		if( res ) {
			int n = getNumberOfParameters(), i = 0;
			res = n == p.getNumberOfParameters();
			while( res && i < n && getParameterAt(i).isComparable( p.getParameterAt(i) ) ) {
				i++;
			}
			return res && i < n;
		}
		return res;
	}
	
	/**
	 * Class for a {@link AnnotatedEntityList} that automatically sets
	 * the {@link Parameter#parent} field to the enclosing {@link ParameterSet}.
	 * @author Jan Grau
	 *
	 */
	protected class ParameterList extends AnnotatedEntityList<Parameter>{

		/**
		 * 
		 */
		public ParameterList() {
			super();
		}

		/**
		 * @param initialSize
		 */
		public ParameterList( int initialSize ) {
			super( initialSize );
		}

		@Override
		public void set( int idx, Parameter p ) {
			p.setParent( ParameterSet.this );
			super.set( idx, p );
		}

		@Override
		public void add( int idx, Parameter p ) {
			p.setParent( ParameterSet.this );
			super.add( idx, p );
		}

		@Override
		public void add( Parameter... p ) {
			for(int i=0;i<p.length;i++){
				p[i].setParent( ParameterSet.this );
			}
			super.add( p );
		}
	}
}
