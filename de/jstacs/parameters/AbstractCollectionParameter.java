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

import java.util.HashSet;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a collection parameter, i.e. a parameter that provides some
 * collection of possible values the user can choose from. Instances of this
 * class can be used as parameters in a {@link ParameterSet}.
 * 
 * @see de.jstacs.parameters.ParameterSet
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public abstract class AbstractCollectionParameter extends Parameter implements Rangeable, GalaxyConvertible {

	/**
	 * The internal {@link ParameterSet} that holds the possible values
	 */
	protected ParameterSet parameters;

	/**
	 * <code>true</code> if the user has selected an item
	 */
	protected boolean userSelected;

	/**
	 * <code>true</code> if the current parameter is required, false otherwise
	 */
	private boolean required;

	/**
	 * If a value was illegal for the collection parameter, this field holds the
	 * error message.
	 */
	protected String errorMessage;

	/**
	 * Indicates if this {@link AbstractCollectionParameter} shall be rangeable
	 */
	private boolean rangeable;

	// default constructor
	private AbstractCollectionParameter( DataType datatype, String name, String comment, boolean required ) {
		super( name, comment, datatype );
		this.required = required;
		this.userSelected = false;
		this.rangeable = true;
	}
	
	protected abstract void init();

	/**
	 * Constructor for a {@link AbstractCollectionParameter} of {@link SimpleParameter}s.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if this {@link AbstractCollectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and the
	 *             <code>values</code> is different or the collection is
	 *             inconsistent for some other reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see #AbstractCollectionParameter(DataType, String[], Object[], String[], String, String, boolean)
	 * @see SimpleParameter
	 */
	public AbstractCollectionParameter(DataType datatype, String[] keys, Object[] values, String name, String comment, boolean required)
			throws InconsistentCollectionException, IllegalValueException,
			DatatypeNotValidException {
		this( datatype, keys, values, null, name, comment, required );
	}

	/**
	 * Constructor for a {@link AbstractCollectionParameter}.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param comments
	 *            the comments on the possible values
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the {@link AbstractCollectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and <code>values</code> is
	 *             different or the collection is inconsistent for some other
	 *             reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see #createParameterSet(Object[], String[], String[])
	 */
	public AbstractCollectionParameter(DataType datatype, String[] keys, Object[] values, String[] comments, String name, String comment,
			boolean required) throws InconsistentCollectionException, IllegalValueException, DatatypeNotValidException {
		this(datatype, name, comment, required);

		if ( !(values instanceof Parameter[]) && (keys == null || keys.length != values.length || (comments != null && keys.length != comments.length)) ) {
			throw new InconsistentCollectionException( "You have to define the same number of keys and values for a AbstractCollectionParameter!");
		}

		createParameterSet(values, keys, comments);
	}

	/**
	 * Constructor for a {@link AbstractCollectionParameter} from an array of
	 * {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link AbstractCollectionParameter} that lets the user select from a list of
	 * possible options that all require an own set of {@link Parameter}s to be
	 * instantiated.
	 * 
	 * @param values
	 *            the array of {@link ParameterSet}s
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param comments
	 *            the comments on the possible values
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @see #createParameterSet(Object[], String[], String[])
	 * @see ParameterSet#getName(ParameterSet)
	 * @see ParameterSet#getComment(ParameterSet)
	 */
	public AbstractCollectionParameter( String name, String comment, boolean required, ParameterSet... values) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		this(DataType.PARAMETERSET, name, comment, required);
		//XXX try catch?
		createParameterSet(values, null, null);
	}

	public AbstractCollectionParameter( String name, String comment, boolean required, Class<? extends ParameterSet>... values) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		this(DataType.PARAMETERSET, name, comment, required);
		//XXX try catch?
		createParameterSet(values, null, null);
	}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores an instance of {@link AbstractCollectionParameter} from a XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public AbstractCollectionParameter(StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/**
	 * Creates a new {@link ParameterSet} from an array of values, an array of
	 * names and an array of comments.
	 * 
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param comments
	 *            the comments on the possible values
	 * 
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the data type is not allowed, i.e. is not one of the
	 *             primitive data types, {@link String} or
	 *             {@link DataType#PARAMETERSET}
	 */
	@SuppressWarnings("unchecked")
	protected void createParameterSet(Object[] values, String[] keys, String[] comments) throws DatatypeNotValidException, IllegalValueException, InconsistentCollectionException {
		Parameter[] pars = new Parameter[values.length];
		HashSet<String> hash = new HashSet<String>();
		for (int i = 0; i < pars.length; i++) {
			if( values[i] instanceof Parameter ) {
				pars[i] = (Parameter) values[i];
			} else if( values[i] instanceof ParameterSet ) {
				pars[i] = new ParameterSetContainer( (ParameterSet) values[i] );
			} else if( values[i] instanceof Class && ParameterSet.class.isAssignableFrom( (Class) values[i] )  ) {
				pars[i] = new ParameterSetContainer( (Class<? extends ParameterSet>) values[i] );
			} else {
				if (keys == null || keys[i] == null) {
					throw new IllegalArgumentException( "You have to state the key for entity " + i);
				}
				String c;
				if (comments != null) {
					c = comments[i];
				} else {
					c = "";
				}
				pars[i] = new SimpleParameter(datatype, keys[i], c, false);
				pars[i].setValue(values[i]);
			}
			if( !hash.contains( pars[i].getName() ) ) {
				hash.add( pars[i].getName() );
			} else {
				throw new InconsistentCollectionException( "The key \"" + pars[i].getName() +"\" is used multiple times." );
			}
		}

		parameters = new SimpleParameterSet(pars);
		init();
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public AbstractCollectionParameter clone() throws CloneNotSupportedException {
		AbstractCollectionParameter clone = (AbstractCollectionParameter) super.clone();
		clone.parameters = parameters == null ? null : parameters.clone();
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		// return !(getValue() instanceof ParameterSet);
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getValue() instanceof ParameterSet) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns <code>true</code>, if this {@link AbstractCollectionParameter} has a
	 * default value.
	 * 
	 * @return if this {@link AbstractCollectionParameter} has a default value
	 * 
	 * @see #setDefault(Object)
	 */
	public abstract boolean hasDefault();

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return isUserSelected();
	}



	/**
	 * Returns the possible values in this collection.
	 * 
	 * @return the values/parameters in this collection
	 */
	public ParameterSet getParametersInCollection() {
		return parameters;
	}

	/**
	 * Sets the value returned by {@link #isRangeable()} to
	 * <code>rangeable</code>.
	 * 
	 * @param rangeable
	 *            the new value
	 */
	public void setRangeable(boolean rangeable) {
		this.rangeable = rangeable;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#isRangeable()
	 */
	public boolean isRangeable() {
		return rangeable;
	}
	
	@SuppressWarnings("unchecked")
	protected int check(Object value) {
		if (value == null) {
			return -1;
		}
		
		Object val2 = value;
		if( value instanceof ParameterSet ) {
			val2 = ParameterSet.getName( (ParameterSet)value );
		} else if( value instanceof Class && ParameterSet.class.isAssignableFrom( (Class) value )  ) {
			val2 = ParameterSet.getName( (Class<? extends ParameterSet>)value );
		} else if( value instanceof Parameter ) {
			val2 = ((Parameter)value).getName();
		}
		
		if ( val2 instanceof String ) {
			for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
				if (parameters.getParameterAt(i).getName().equals(val2)) {
					errorMessage = null;
					return i;
				}
			}
		}
		errorMessage = "The value is not in the set of defined values: " + value + ".";
		return -1;
	}

	/**
	 * Returns <code>true</code> if the key specified by <code>value</code> is
	 * in the set of keys of this {@link AbstractCollectionParameter}.
	 * 
	 * @return <code>true</code> if the key specified by <code>value</code> is
	 *         in the set of keys of this {@link AbstractCollectionParameter},
	 *         <code>false</code> otherwise
	 */
	@Override
	public boolean checkValue(Object value) {
		return check(value) >= 0 ;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		return errorMessage;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, userSelected, "userSelected");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, rangeable, "rangeable");
		XMLParser.appendObjectWithTags(buf, parameters, "collection");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		
		required = XMLParser.extractObjectForTags(representation, "required", boolean.class );
		userSelected = XMLParser.extractObjectForTags(representation, "userSelected", boolean.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		StringBuffer help = XMLParser.extractForTag(representation, "rangeable");
		if (help == null) {
			rangeable = false;
		} else {
			rangeable = Boolean.parseBoolean(help.toString());
		}
		parameters = new SimpleParameterSet(XMLParser.extractForTag(representation,"collection"));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return required;
	}

	/**
	 * Returns <code>true</code> if the option at position <code>idx</code> is
	 * selected.
	 * 
	 * @param idx
	 *            the position
	 * 
	 * @return <code>true</code> if the option at position <code>idx</code> is
	 *         selected, <code>false</code> otherwise
	 */
	public abstract boolean isSelected(int idx);

	/**
	 * Returns <code>true</code> if the value was selected by the user.
	 * 
	 * @return <code>true</code> if the value was selected by the user.
	 */
	public boolean isUserSelected() {
		return userSelected;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o2) {
		if (o2 instanceof AbstractCollectionParameter) {
			ParameterSet parSet2 = ((AbstractCollectionParameter) o2)
					.getParametersInCollection();
			if (parSet2.getNumberOfParameters() != parameters
					.getNumberOfParameters()) {
				return false;
			} else {
				if (!(((AbstractCollectionParameter) o2).getName().equals(name) && ((AbstractCollectionParameter) o2)
						.getComment().equals(comment))) {
					return false;
				}

				for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
					if (!parameters.getParameterAt(i).equals(
							parSet2.getParameterAt(i))) {
						return false;
					}
				}
				return true;
			}
		} else {
			return false;
		}
	}

	/**
	 * This exception is thrown if the {@link AbstractCollectionParameter} is
	 * inconsistent for some reason.
	 * 
	 * @author Jan Grau
	 */
	public class InconsistentCollectionException extends ParameterException {

		private static final long serialVersionUID = -2703514434545861722L;

		/**
		 * Constructs a new {@link InconsistentCollectionException} with message
		 * <code>message</code>.
		 * 
		 * @param message
		 *            the error message
		 */
		public InconsistentCollectionException(String message) {
			super(message);
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#getRangedInstance()
	 */
	@Override
	public MultiSelectionParameter getRangedInstance() throws Exception {
		//boolean[] selected = new boolean[parameters.getNumberOfParameters()];
		//selected[0] = true;
		Parameter[] p = new Parameter[parameters.getNumberOfParameters()];
		for( int i = 0; i < p.length; i++ ) {
			p[i] = parameters.getParameterAt(i);
		}
		MultiSelectionParameter par = 
			new MultiSelectionParameter( datatype, null, p, null, getName(), getComment(), required );
		// new MultiSelectionParameter( this.parameters.clone(), selected, new boolean[parameters.getNumberOfParameters()], false, getName(), getComment(), required, datatype, errorMessage, 0 );
		return par;
	}
	

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		StringBuffer buf = new StringBuffer();
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			if(isSelected( i )){
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\" selected=\"true\"", false );
			}else{
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\"", false );
			}
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"select\" format=\"text\" name=\""+namePrefix+"\" label=\""+getName()+"\" optional=\""+(!isRequired())+"\" help=\""+getComment()+"\"" );
	
		StringBuffer buf3 = new StringBuffer("${"+configPrefix+(isAtomic() ? "" : namePrefix+"_cond.")+namePrefix+"}");
		XMLParser.addTags( buf3, namePrefix );
		configBuffer.append( buf3 );
		if(!isAtomic()){
			for(int i=0;i<parameters.getNumberOfParameters();i++){
				StringBuffer temp = new StringBuffer();
				StringBuffer temp2 = new StringBuffer();
				((GalaxyConvertible)parameters.getParameterAt( i )).toGalaxy( namePrefix+"_opt"+i, configPrefix+namePrefix+"_cond.", depth+1, temp, temp2 );
				XMLParser.addTagsAndAttributes( temp, "when", "value=\""+parameters.getParameterAt( i ).getName()+"\"" );
				buf.append( temp );
				configBuffer.append( (i== 0 ? "#if " : "#elif ")+"$"+configPrefix+namePrefix+"_cond."+namePrefix+" == \""+parameters.getParameterAt( i ).getName()+"\"\n" );
				configBuffer.append( temp2 );
				configBuffer.append( "\n" );
			}
			configBuffer.append( "#end if\n" );
		}
		
		
		if(!isAtomic()){
			XMLParser.addTagsAndAttributes( buf, "conditional", "name=\""+namePrefix+"_cond\"" );
		}
		
		descBuffer.append( buf );
	}
}